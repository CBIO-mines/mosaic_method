#!/usr/bin/env python3

"""
global_fitting.py – Block coordinate descent fitting of kappa, {theta_i, log10_xi_i}

The single global parameter is kappa,
the ratio of the rate of the least and most conserved loci of a comparison.
All per-pair parameters (theta_i, log10_xi_i) are conditionally independent
given kappa and are fitted as independent 2-D problems.

Structure
---------
1. GRID SCAN over kappa:
   For each grid point, solve all per-pair 2-D problems independently.

2. BLOCK COORDINATE DESCENT:
   Alternate between:
     - optimizing kappa with pairwise parameters fixed
     - optimizing all pairwise parameters with kappa fixed
   until convergence.

3. Final joint L-BFGS-B polish over all parameters simultaneously,
   warm-started from the best refined grid point(s).
"""

import argparse
import itertools
import os
from functools import partial, reduce
import matplotlib.pyplot as plt
import seaborn as sns
import jax
import jax.numpy as jnp
from jax import value_and_grad
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import pearsonr, spearmanr

import yaml

from mosaic_method.parsing import get_genome_comp, sum_mlds, bin_mld

jax.config.update("jax_enable_x64", True)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
SMAL_DIF = 0.1
DELTA = 0.82
MUS = 3.64e-9
BOUNDS_KAPPA = (1.1, 30)  # kappa in (1, inf)
BOUNDS_THETA = (0.01, 3.0)
BOUNDS_LOG10_XI = (-5.0, -0.1)
N_PAIR_STARTS = (
    6  # number of multi-start points for fit_one_pair; tune for cost/robustness
)


def make_diverse_starts(
    n_starts, bounds_theta=BOUNDS_THETA, bounds_log10_xi=BOUNDS_LOG10_XI
):
    """
    n_starts deterministic (theta, log10_xi) starting points spanning
    bounds_theta (log-spaced) x bounds_log10_xi (linearly spaced).  The two
    axes are paired with a fixed offset (rather than the diagonal) so the
    points spread across the box instead of collapsing to a line.
    """
    log10_theta = np.linspace(
        np.log10(bounds_theta[0]), np.log10(bounds_theta[1]), n_starts
    )
    log10_xi = np.linspace(*bounds_log10_xi, n_starts)
    log10_xi = np.roll(log10_xi, max(1, n_starts // 3))
    return [(10.0**lt, lx) for lt, lx in zip(log10_theta, log10_xi)]


DIVERSE_STARTS = make_diverse_starts(N_PAIR_STARTS)


# ---------------------------------------------------------------------------
# Physics kernel
# ---------------------------------------------------------------------------


def theoretical_mld_vectorized(
    thetas, xis_raw, smal_dif, ml_safe, kappa_raw, delta, L0s
):
    """
    Vectorised theoretical MLD, fully JAX-differentiable.

    theta              : (n,)  values
    xis_raw            : (n,)  log10 values
    ml_safe            : (n, max_len)
    L0s                : (n,)
    Returns mt         : (n, max_len)
    """
    xis = (10.0**xis_raw)[:, None]
    kappa = 1 / kappa_raw
    L0 = L0s[:, None]

    ml_low = ml_safe - smal_dif
    ml_hi = ml_safe + smal_dif
    thetas_2d = thetas[:, None]  # (n,1) mus·tau, uncapped
    thetas_capped = jnp.minimum(delta, thetas_2d)  # (n,1) mua·tau = min(δ, mus·tau)

    def _mc_term(x):
        return (
            (1.0 + x * kappa * thetas_2d) * jnp.exp(-x * kappa * thetas_2d)
            - (1.0 + x * thetas_capped) * jnp.exp(-x * thetas_capped)
        ) / (x**2)

    mc = (
        (2.0 * L0)
        / ((1 - kappa**2) * thetas_2d**2)
        * (_mc_term(ml_low) + _mc_term(ml_hi) - 2.0 * _mc_term(ml_safe))
        / smal_dif**2
    )

    def _mh_low_branch(x):
        return (
            (1 - kappa) * thetas_2d * x
            + jnp.exp(-x * thetas_2d)
            - jnp.exp(-x * kappa * thetas_2d)
        ) / (x**2)

    def _mh_high_branch(x):
        return (
            (1 - kappa) * thetas_2d * x
            + (1.0 + x * (delta - thetas_2d)) * jnp.exp(-x * delta)
            - jnp.exp(-x * kappa * thetas_2d)
        ) / (x**2)

    use_low = (thetas < delta)[:, None]

    def _mh(x):
        return jnp.where(use_low, _mh_low_branch(x), _mh_high_branch(x))

    mh = (
        (L0 * xis)
        / (thetas_2d * (1 - kappa**2))
        * (_mh(ml_low) + _mh(ml_hi) - 2.0 * _mh(ml_safe))
        / smal_dif**2
    )

    return mh + mc


# ---------------------------------------------------------------------------
# Loss functions
# ---------------------------------------------------------------------------


def global_loss(opt_pars, precomputed, smal_dif=SMAL_DIF, delta=DELTA):
    """
    Full joint loss.
    opt_pars: [kappa, theta_0, log10_xi_0, ...]
    """
    _, ml_pad, emp_pad, mask, L0s = precomputed
    ml_safe = jnp.where(mask, ml_pad, 1.0)

    kappa = opt_pars[0]
    thetas = opt_pars[1::2]
    log10_xis = opt_pars[2::2]

    mt = theoretical_mld_vectorized(
        thetas, log10_xis, smal_dif, ml_safe, kappa, delta, L0s
    )
    denom = mt + emp_pad
    ratio = jnp.where(denom > 0.0, ((mt - emp_pad) / denom) ** 2, 0.0)
    ratio = jnp.where(mask, ratio, 0.0)
    return (ratio.sum(axis=1) / mask.sum(axis=1)).sum()


def pair_loss(p, kappa, ml_row, emp_row, mask_row, L0s_1, smal_dif, delta):
    """
    Loss for a single pair.  All data is passed explicitly so that
    jax.jit sees one function signature for every pair and compiles once.
    """
    theta, log10_xi = p[0], p[1]
    ml_safe = jnp.where(mask_row, ml_row, 1.0)
    mt = theoretical_mld_vectorized(
        jnp.array([theta]),
        jnp.array([log10_xi]),
        smal_dif,
        ml_safe[None, :],
        kappa,
        delta,
        L0s_1,
    )[0]
    denom = mt + emp_row
    ratio = jnp.where(denom > 0.0, ((mt - emp_row) / denom) ** 2, 0.0)
    ratio = jnp.where(mask_row, ratio, 0.0)
    return ratio.sum() / mask_row.sum()


pair_lag = jax.jit(value_and_grad(pair_loss, argnums=0))


# ---------------------------------------------------------------------------
# Precomputation
# ---------------------------------------------------------------------------


def precompute_pairs(binned_mlds, L0_df):
    pairs = list(binned_mlds.keys())
    n = len(pairs)
    max_len = max(len(v) for v in binned_mlds.values())

    ml_pad = np.zeros((n, max_len))
    emp_pad = np.zeros((n, max_len))
    mask = np.zeros((n, max_len), dtype=bool)
    L0s = np.zeros(n)

    L0_lookup = {(row.bac1, row.bac2): row.L0 for row in L0_df.itertuples(index=False)}

    for i, pair in enumerate(pairs):
        df = binned_mlds[pair]
        ml, emp, L = df["match_length"].values, df["freq"].values, len(df)
        ml_pad[i, :L] = ml
        emp_pad[i, :L] = emp
        mask[i, :L] = True
        L0s[i] = L0_lookup[pair]

    return (
        pairs,
        jnp.array(ml_pad),
        jnp.array(emp_pad),
        jnp.array(mask),
        jnp.array(L0s),
    )


# ---------------------------------------------------------------------------
# Inner solver
# ---------------------------------------------------------------------------
def fit_one_pair(
    i, precomputed, kappa, smal_dif=SMAL_DIF, delta=DELTA, warm_start=None
):
    """
    Optimise (theta, log10_xi) for pair i at fixed kappa.

    DIVERSE_STARTS (N_PAIR_STARTS fixed points spanning BOUNDS_THETA x
    BOUNDS_LOG10_XI) are always tried, plus warm_start (previous grid-point
    / BCD-iteration solution) if provided — exploits smoothness across
    kappa/iterations without depending on any single starting point to land
    in a good basin. The result with the lowest loss is returned.
    """
    _, ml_pad, emp_pad, mask, L0s = precomputed
    ml_row = ml_pad[i]
    emp_row = emp_pad[i]
    mask_row = mask[i]
    L0s_1 = L0s[i : i + 1]

    def _obj(p):
        loss, grad = pair_lag(
            jnp.array(p), kappa, ml_row, emp_row, mask_row, L0s_1, smal_dif, delta
        )
        loss = float(loss)
        grad = np.array(grad, dtype=np.float64)
        if not np.isfinite(loss):
            print(f"[pair {i}] NONFINITE LOSS")
        if not np.all(np.isfinite(grad)):
            print(f"[pair {i}] NONFINITE GRADIENT")
            print(np.where(~np.isfinite(grad))[0])
            # Replace NaN/Inf gradient entries with 0 so L-BFGS-B internal
            # state is not corrupted; the zero tells the optimizer not to move
            # in those directions.
            grad = np.where(np.isfinite(grad), grad, 0.0)
        return loss, grad

    _bounds = [BOUNDS_THETA, BOUNDS_LOG10_XI]
    _opts = {"ftol": 1e-9, "gtol": 1e-8, "maxls": 50}

    inits = [np.array(p) for p in DIVERSE_STARTS]
    if warm_start is not None:
        inits.append(np.asarray(warm_start))

    best = None
    for init in inits:
        res = minimize(
            _obj, init, method="L-BFGS-B", jac=True, bounds=_bounds, options=_opts
        )
        if best is None or res.fun < best.fun:
            best = res

    return float(best.x[0]), float(best.x[1]), float(best.fun)


def fit_pair_block(
    precomputed, kappa, smal_dif=SMAL_DIF, delta=DELTA, warm_starts=None
):
    """
    warm_starts : (n, 2) array of per-pair [thetas, log10_xis], or None.
                  Passed as a second starting point alongside the kappa-proportional init.
    """
    n = len(precomputed[0])
    thetas = np.empty(n)
    log10_xis = np.empty(n)
    total_loss = 0.0

    for i in range(n):
        ws = None if warm_starts is None else warm_starts[i]
        t, lx, lv = fit_one_pair(i, precomputed, kappa, smal_dif, delta, warm_start=ws)
        thetas[i] = t
        log10_xis[i] = lx
        total_loss += lv

    return thetas, log10_xis, total_loss


def kappa_only_fit(
    kappa,
    thetas,
    log10_xis,
    precomputed,
    bounds_kappa,
    smal_dif=SMAL_DIF,
    delta=DELTA,
    verbose=True,
):
    """
    Optimize kappa only, keeping all pairwise parameters fixed.
    """
    fixed_pairs = np.array([thetas, log10_xis], dtype=np.float64).T.ravel()

    _lag = jax.jit(
        value_and_grad(
            lambda kk: global_loss(
                jnp.concatenate([jnp.atleast_1d(kk), jnp.asarray(fixed_pairs)]),
                precomputed,
                smal_dif,
                delta,
            )
        )
    )

    def _obj(p):
        l, g = _lag(jnp.asarray(p))
        return np.array(l, dtype=np.float64), np.array(g, dtype=np.float64)

    res = minimize(
        _obj,
        np.array([kappa], dtype=np.float64),
        method="L-BFGS-B",
        jac=True,
        bounds=[bounds_kappa],
        options={"maxfun": int(1e7), "maxcor": 30},
    )

    if verbose:
        print("\n=== kappa-only fit ===")
        print(f"success : {res.success}")
        print(f"status  : {res.status}")
        print(f"message : {repr(res.message)}")
        print(f"nit     : {res.nit}")
        print(f"final f : {res.fun:.12f}")
        print(f"start kappa : {kappa:.6f}")
        print(f"end kappa   : {float(res.x[0]):.6f}")
        print("===============================\n")

    return res


# ---------------------------------------------------------------------------
# Grid scan over kappa
# ---------------------------------------------------------------------------


def grid_scan(precomputed, kappa_grid, smal_dif=SMAL_DIF, delta=DELTA, verbose=True):
    """
    Evaluate min_{pairs} L(kappa) at every point of the kappa grid.
    Returns a list of result dicts sorted by loss (best first).
    """
    results = []
    warm_starts = None  # populated from previous grid point after first solve

    for k, kappa in enumerate(kappa_grid):
        thetas, log10_xis, loss = fit_pair_block(
            precomputed, kappa, smal_dif, delta, warm_starts=warm_starts
        )
        results.append(
            {
                "kappa": kappa,
                "thetas": thetas,
                "log10_xis": log10_xis,
                "loss": loss,
            }
        )
        if verbose:
            print(
                f"  grid {k+1}/{len(kappa_grid)}  "
                f"kappa={kappa:.2e}  loss={loss:.6f}"
            )
        warm_starts = np.stack([thetas, log10_xis], axis=1)

    results.sort(key=lambda r: r["loss"])
    return results


# ---------------------------------------------------------------------------
# Block descent polish
# ---------------------------------------------------------------------------


def block_coordinate_descent_refine(
    kappa,
    thetas,
    log10_xis,
    precomputed,
    bounds_kappa=BOUNDS_KAPPA,
    n_iter=10,
    tol=1e-10,
    smal_dif=SMAL_DIF,
    delta=DELTA,
    verbose=True,
):
    """
    Block coordinate descent:
      1) optimize kappa with pairwise parameters fixed
      2) optimize pairwise parameters with kappa fixed
    Repeat until convergence, then return the refined state.
    """
    cur_kappa = float(kappa)
    cur_thetas = np.asarray(thetas, dtype=np.float64)
    cur_log10_xis = np.asarray(log10_xis, dtype=np.float64)

    cur_x = np.concatenate(
        [[cur_kappa], np.ravel(np.column_stack([cur_thetas, cur_log10_xis]))]
    )
    prev_loss = float(global_loss(cur_x, precomputed, smal_dif, delta))

    best_loss = prev_loss
    best_kappa = cur_kappa
    best_thetas = cur_thetas.copy()
    best_log10_xis = cur_log10_xis.copy()

    for it in range(n_iter):
        kres = kappa_only_fit(
            cur_kappa,
            cur_thetas,
            cur_log10_xis,
            precomputed,
            bounds_kappa=bounds_kappa,
            smal_dif=smal_dif,
            delta=delta,
            verbose=verbose,
        )
        cur_kappa = float(kres.x[0])

        cur_thetas, cur_log10_xis, cur_loss = fit_pair_block(
            precomputed,
            cur_kappa,
            smal_dif=smal_dif,
            delta=delta,
            warm_starts=np.stack([cur_thetas, cur_log10_xis], axis=1),
        )

        if verbose:
            print(f"[BCD] iter={it+1}  kappa={cur_kappa:.6f}  loss={cur_loss:.12f}")

        if cur_loss > prev_loss + tol * (1.0 + abs(prev_loss)):
            print(
                f"[BCD] WARNING: loss increased at iter={it+1} "
                f"({prev_loss:.12f} -> {cur_loss:.12f})"
            )

        if cur_loss < best_loss:
            best_loss = cur_loss
            best_kappa = cur_kappa
            best_thetas = cur_thetas.copy()
            best_log10_xis = cur_log10_xis.copy()

        if abs(prev_loss - cur_loss) <= tol * (1.0 + abs(prev_loss)):
            break
        prev_loss = cur_loss

    return best_kappa, best_thetas, best_log10_xis, best_loss


# ---------------------------------------------------------------------------
# Joint polish
# ---------------------------------------------------------------------------


def joint_polish(
    kappa,
    thetas,
    log10_xis,
    precomputed,
    bounds_kappa,
    smal_dif=SMAL_DIF,
    delta=DELTA,
    verbose=True,
):
    """
    Full joint L-BFGS-B over all parameters simultaneously.
    """
    n = len(thetas)
    pair_flat = np.array([thetas, log10_xis]).T.ravel()
    warm = np.concatenate([[kappa], pair_flat])
    bounds = [bounds_kappa] + [BOUNDS_THETA, BOUNDS_LOG10_XI] * n

    _lag = jax.jit(
        value_and_grad(lambda p: global_loss(p, precomputed, smal_dif, delta))
    )

    def _obj(p):
        l, g = _lag(jnp.array(p))
        return np.array(l, dtype=np.float64), np.array(g, dtype=np.float64)

    res = minimize(
        _obj,
        warm,
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        options={"maxfun": int(1e7), "maxcor": 30},
    )
    if verbose:
        print("\n=== Joint polish diagnostics ===")
        print(f"success : {res.success}")
        print(f"status  : {res.status}")
        print(f"message : {repr(res.message)}")
        print(f"nit     : {res.nit}")
        print(f"nfev    : {res.nfev}")
        print(f"njev    : {res.njev}")
        print(f"final f : {res.fun:.12f}")

        if hasattr(res, "jac"):
            gnorm = np.linalg.norm(res.jac)
            gmax = np.max(np.abs(res.jac))
            print(f"|g|     : {gnorm:.3e}")
            print(f"max|g|  : {gmax:.3e}")

        print("===============================\n")
    return res


# ---------------------------------------------------------------------------
# Extract results
# ---------------------------------------------------------------------------


def extract_fitted_params(res, precomputed):
    pairs, _, _, _, L0s = precomputed
    x = res.x
    kappa = x[0]
    thetas = x[1::2]
    log10_xis = x[2::2]
    return [
        {
            "species_1": pair[0],
            "species_2": pair[1],
            "kappa": kappa,
            "theta": thetas[i],
            "xi": 10.0 ** log10_xis[i],
            "L0": float(L0s[i]),
            "minimum": res.fun,
        }
        for i, pair in enumerate(pairs)
    ]


# ---------------------------------------------------------------------------
# Pair subsetting
# ---------------------------------------------------------------------------


def subset_precomputed(precomputed, indices):
    """Return a precomputed tuple restricted to the given pair indices."""
    pairs, ml_pad, emp_pad, mask, L0s = precomputed
    idx = jnp.array(indices)
    return (
        [pairs[i] for i in indices],
        ml_pad[idx],
        emp_pad[idx],
        mask[idx],
        L0s[idx],
    )


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------


def plot_kappa_tranche(
    res, precomputed, smal_dif=SMAL_DIF, delta=DELTA, out_path="kappa_tranche.png"
):
    """Plots the marginal objective function at the minimum for kappa values."""

    x = res.x
    kappa = float(x[0])
    pair_pars = x[1:]
    kappa_ranges = np.linspace(BOUNDS_KAPPA[0], BOUNDS_KAPPA[1], 1000)
    obj_fun = np.zeros_like(kappa_ranges)
    for i, kp in enumerate(kappa_ranges):
        obj_fun[i] = global_loss(
            np.concatenate([np.array([kp]), pair_pars]), precomputed, smal_dif, delta
        )
    fig = plt.Figure(dpi=300)
    ax = fig.subplots()
    ax.plot(kappa_ranges, obj_fun)
    ax.axvline(kappa, ls="--", label="Optimum")
    ax.set_title("κ Slice")
    fig.savefig(out_path)


def hessian_diagnostics(res, precomputed, smal_dif=SMAL_DIF, delta=DELTA):
    """
    Analyse the objective landscape at the optimum using the bordered
    block-diagonal structure of the Hessian.  Never builds the full
    (2N+1)×(2N+1) matrix.

    For each pair i:
      H_i = d²L_i / d(theta_i, log10_xi_i)²              (2×2)
      c_i = d²L_i / d(theta_i, log10_xi_i) d(kappa)  (2,)
      a_i = d²L_i / d(kappa)²                              (scalar)

    The Schur complement  a_total - Σ cᵢᵀ Hᵢ⁻¹ cᵢ  is the effective
    curvature in the kappa direction after marginalising the pair params.

    Parameters
    ----------
    res         : OptimizeResult from joint_polish
    precomputed : output of precompute_pairs

    Returns
    -------
    dict with keys H_all, eigenvalues, cond_numbers, c_all, a_all, schur_kappa
    """
    pairs, ml_pad, emp_pad, mask, L0s = precomputed
    n = len(pairs)

    x = res.x
    kappa = float(x[0])
    thetas = x[1::2]
    log10_xis = x[2::2]

    p_all = jnp.array(np.stack([thetas, log10_xis], axis=1))  # (N, 2)
    lm = jnp.array(kappa)  # scalar 0-d
    L0s_2d = L0s[:, None]  # (N, 1)

    # pair_loss already takes kappa — no wrapper needed.
    # Freeze smal_dif and delta so the vmapped functions share one compilation.
    _pair_loss_fixed = partial(pair_loss, smal_dif=smal_dif, delta=delta)

    _in_axes = (0, None, 0, 0, 0, 0)

    H_fn = jax.jit(jax.vmap(jax.hessian(_pair_loss_fixed, argnums=0), in_axes=_in_axes))
    c_fn = jax.jit(
        jax.vmap(
            jax.jacobian(jax.grad(_pair_loss_fixed, argnums=0), argnums=1),
            in_axes=_in_axes,
        )
    )
    a_fn = jax.jit(jax.vmap(jax.hessian(_pair_loss_fixed, argnums=1), in_axes=_in_axes))

    args = (p_all, lm, ml_pad, emp_pad, mask, L0s_2d)

    print("Computing per-pair Hessians …", flush=True)
    H_all = np.array(H_fn(*args))  # (N, 2, 2)
    print("Computing cross-derivatives …", flush=True)
    c_all = np.array(c_fn(*args))  # (N, 2)
    print("Computing kappa curvature …", flush=True)
    a_all = np.array(a_fn(*args))  # (N,)

    # Per-pair eigenvalues (ascending)
    eigenvalues = np.linalg.eigvalsh(H_all)  # (N, 2)
    cond_numbers = eigenvalues[:, 1] / np.maximum(eigenvalues[:, 0], 1e-30)

    # Schur complement for kappa
    H_inv_c = np.array(
        [np.linalg.lstsq(H_all[i], c_all[i], rcond=None)[0] for i in range(n)]
    )  # (N, 2)
    schur_kappa = float(a_all.sum()) - float((c_all * H_inv_c).sum())

    # ------------------------------------------------------------------ report
    print("\n=== Hessian diagnostics at optimum ===")
    print(f"n_pairs={n}   kappa={kappa:.4f}")

    lam_min, lam_max = eigenvalues[:, 0], eigenvalues[:, 1]
    print("\n-- Per-pair (theta, log10_xi) blocks --")
    print(
        f"  min eigenvalue : min={lam_min.min():.2e}  "
        f"p5={np.percentile(lam_min, 5):.2e}  "
        f"median={np.median(lam_min):.2e}"
    )
    print(
        f"  max eigenvalue : median={np.median(lam_max):.2e}  "
        f"p95={np.percentile(lam_max, 95):.2e}  "
        f"max={lam_max.max():.2e}"
    )
    print(
        f"  condition number: median={np.median(cond_numbers):.1e}  "
        f"p95={np.percentile(cond_numbers, 95):.1e}  "
        f"max={cond_numbers.max():.1e}"
    )
    n_flat = int((lam_min < 1e-6).sum())
    print(f"  pairs with min eigenvalue < 1e-6 : {n_flat}/{n}")

    worst = np.argsort(cond_numbers)[-5:][::-1]
    print(f"\n  5 worst-conditioned pairs:")
    for idx in worst:
        s1, s2 = pairs[idx]
        print(
            f"    {s1[:24]:24s}  {s2[:24]:24s}  "
            f"cond={cond_numbers[idx]:.2e}  "
            f"eig=({lam_min[idx]:.2e}, {lam_max[idx]:.2e})"
        )

    print("\n-- kappa curvature --")
    print(f"  raw  Σ d²L_i/d(kappa)²  : {a_all.sum():.3e}")
    print(f"  Schur complement             : {schur_kappa:.3e}")
    if schur_kappa < 0:
        print("  WARNING: non-positive Schur complement — kappa is poorly identified")

    return {
        "H_all": H_all,
        "eigenvalues": eigenvalues,
        "cond_numbers": cond_numbers,
        "c_all": c_all,
        "a_all": a_all,
        "schur_kappa": schur_kappa,
    }


def scale_mus(results_df, mus):
    """
    Add columns of tau, rho and muc and mus based on a value of mus
    """
    results_df["tau"] = results_df["theta"] / mus
    results_df["log10tau"] = np.log10(results_df["tau"])
    results_df["muc"] = 1 / results_df["kappa"] * mus
    results_df["mus"] = mus
    results_df["rho"] = results_df["xi"] * mus
    results_df["log10rho"] = np.log10(results_df["rho"])


def compare_fits(results_dic, columns=("theta", "xi", "kappa")):
    """
    Merge any number of fit CSVs on (species_1, species_2).

    Parameters
    ----------
    results_dic : dict[str, str]
        Mapping from run name to CSV path.
    columns : iterable of str
        Columns to carry from each CSV, suffixed with the run name.

    Returns
    -------
    pd.DataFrame with one row per species pair and one column per
    (column, run-name) combination, suitable for further analysis or
    passing to plot_tau_comparison.
    """
    frames = []
    for name, csv_path in results_dic.items():
        df = pd.read_csv(csv_path)  # , index_col=0)
        keep = ["species_1", "species_2"] + [c for c in columns if c in df.columns]
        df = df[keep].rename(
            columns={c: f"{c}_{name}" for c in columns if c in df.columns}
        )
        frames.append(df)

    return reduce(
        lambda left, right: pd.merge(
            left, right, on=["species_1", "species_2"], how="inner"
        ),
        frames,
    )


def plot_distance_comparison(merged_df, run_names, out_path=None, distance="tau"):
    """
    Pairwise scatter plots of fitted theta for all combinations of run_names.

    Parameters
    ----------
    merged_df : pd.DataFrame
        Output of compare_fits (must contain theta_{name} columns).
    run_names : list[str]
        Run names to compare; all N*(N-1)/2 pairs are plotted.
    out_path : str or None
        If given, save the figure there; otherwise show interactively.
    """
    pairs = list(itertools.combinations(run_names, 2))
    n = len(pairs)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4), squeeze=False)

    for ax, (name_x, name_y) in zip(axes[0], pairs):
        col_x = f"{distance}_{name_x}"
        col_y = f"{distance}_{name_y}"
        x = merged_df[col_x].values
        y = merged_df[col_y].values

        pear = pearsonr(x, y).statistic
        spear = spearmanr(x, y).statistic

        ax.scatter(x, y, alpha=0.3, s=10, rasterized=True)
        lims = [min(x.min(), y.min()), max(x.max(), y.max())]
        ax.plot(lims, lims, "k--", linewidth=0.8)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(f"{distance}  [{name_x}]")
        ax.set_ylabel(f"{distance}  [{name_y}]")
        ax.set_title(f"Pearson r={pear:.3f}  Spearman ρ={spear:.3f}")

    fig.tight_layout()
    if out_path:
        fig.savefig(out_path, dpi=300)
    else:
        plt.show()
    return fig


def load_mlds_from_config(cfg, override_censor=0):
    """
    Load and bin all pairwise MLDs described by a YAML config file.

    Parameters
    ----------
    cfg : dict
        dict with keys: results_dir, database_name,
        taxon_csv, cluster_name, censor.

    Returns
    -------
    binned_mlds : dict[tuple, pd.DataFrame]
    L0_df       : pd.DataFrame
    ncomps      : list[str, str, int]
    """
    if override_censor > 0:
        censor = override_censor
    else:
        censor = float(cfg["censor"])
    L0_df = pd.read_csv(os.path.join(cfg["results_dir"], "L0.csv"))
    taxon_df = pd.read_csv(cfg["taxon_csv"], index_col=0)
    levels = list(
        itertools.combinations(sorted(taxon_df[cfg["cluster_name"]].unique()), 2)
    )

    binned_mlds = {}
    ncomps = {}
    for level in levels:
        genome_comps = get_genome_comp(level, taxon_df, cfg["cluster_name"])
        summed_mld = sum_mlds(
            genome_comps,
            os.path.join(cfg["results_dir"], cfg["database_name"]),
        )
        binned_mlds[level] = bin_mld(
            summed_mld,
            linear_bin_width=3,
            limit_size=30.5,
            power_increment=0.1,
            ncomp=len(genome_comps),
            censor=censor,
        )
        ncomps[level] = len(genome_comps)
    ncomps_list = [[k1, k2, nc] for (k1, k2), nc in ncomps.items()]

    return binned_mlds, L0_df, ncomps_list


def convert_taurho(results_df, muc_def=6e-11, mus_def=5e-9):
    """
    Converts results computed with time to dimensionless units.
    """
    try:
        muc = results_df["muc"]
    except KeyError:
        muc = muc_def
    try:
        mus = results_df["mus"]
    except KeyError:
        mus = mus_def
    if not "tau" in results_df.columns:
        results_df["tau"] = 10 ** results_df["log10tau"]
        results_df["rho"] = 10 ** results_df["log10rho"]

    results_df["xi"] = results_df["rho"] / mus
    results_df["theta"] = results_df["tau"] * mus
    results_df["kappa"] = mus / muc
    return results_df


# ---------------------------------------------------------------------------
# High-level fitting pipeline
# ---------------------------------------------------------------------------


def run_global_fitting(
    precomputed,
    kappa_grid_n=10,
    top_k=1,
    cond_threshold=np.inf,
    bounds_kappa=BOUNDS_KAPPA,
    smal_dif=SMAL_DIF,
    delta=DELTA,
    verbose=True,
    out_path=".",
):
    """
    Full fitting pipeline: grid scan → block coordinate descent → joint polish → optional bad-pair exclusion.

    Parameters
    ----------
    precomputed       : output of precompute_pairs
    kappa_grid_n      : number of grid points along kappa
    top_k             : number of best grid points to joint-polish
    cond_threshold    : exclude pairs whose (theta, log10_xi) Hessian condition
                        number exceeds this value and re-polish the rest
    bounds_kappa: (lo, hi) bounds for kappa
    smal_dif, delta   : physics constants

    Returns
    -------
    list[dict]  — one dict per pair, with keys:
        species_1, species_2, kappa, theta, xi, L0, minimum, excluded
    """
    kappa_grid = np.linspace(*bounds_kappa, kappa_grid_n)

    if verbose:
        print(f"n_pairs={len(precomputed[0])}  kappa bounds={bounds_kappa}")
        print("\n=== Grid scan over kappa ===")
    grid_results = grid_scan(precomputed, kappa_grid, smal_dif, delta, verbose=verbose)

    k = min(top_k, len(grid_results))
    if verbose:
        print(f"\nRefining top {k} grid point(s) with block coordinate descent:")
        for rank, g in enumerate(grid_results[:k]):
            print(f"  rank {rank+1}: kappa={g['kappa']:.2e}  grid loss={g['loss']:.6f}")

    final_res = None
    for rank, g in enumerate(grid_results[:k]):
        if verbose:
            print(
                f"\n=== Block coordinate descent {rank+1}/{k} "
                f"(kappa={g['kappa']:.2e}) ==="
            )
        kappa_bcd, thetas_bcd, log10_xis_bcd, _ = block_coordinate_descent_refine(
            g["kappa"],
            g["thetas"],
            g["log10_xis"],
            precomputed,
            bounds_kappa=bounds_kappa,
            smal_dif=smal_dif,
            delta=delta,
            verbose=verbose,
        )
        if verbose:
            print(f"\n=== Joint polish {rank+1}/{k} ===")
        res = joint_polish(
            kappa_bcd,
            thetas_bcd,
            log10_xis_bcd,
            precomputed,
            bounds_kappa=bounds_kappa,
            smal_dif=smal_dif,
            delta=delta,
            verbose=verbose,
        )
        if final_res is None or res.fun < final_res.fun:
            final_res = res
        print("==================================================")
        print("==================================================")
        print(f"Kappa grid: {g['kappa']:.3f}")
        print(f"Kappa BCD: {kappa_bcd:.3f}")
        print(f"Kappa JP: {res.x[0]:.3f}")
        print("==================================================")
        print("==================================================")

    if verbose:
        print(f"\nBest joint-polish loss: {final_res.fun:.6f}")

    diag = hessian_diagnostics(final_res, precomputed, smal_dif, delta)

    rows = extract_fitted_params(final_res, precomputed)
    plot_kappa_tranche(final_res, precomputed, out_path=out_path)
    for row in rows:
        row["excluded"] = False

    if np.isfinite(cond_threshold):
        bad_idx = np.where(diag["cond_numbers"] > cond_threshold)[0]
        good_idx = np.where(diag["cond_numbers"] <= cond_threshold)[0]
        if len(bad_idx) > 0:
            pairs = precomputed[0]
            if verbose:
                print(
                    f"\nExcluding {len(bad_idx)} pair(s) with condition > {cond_threshold:.1e}:"
                )
                for idx in bad_idx:
                    print(
                        f"  {pairs[idx][0]} vs {pairs[idx][1]}"
                        f"  cond={diag['cond_numbers'][idx]:.2e}"
                    )

            x = final_res.x
            precomputed_good = subset_precomputed(precomputed, good_idx)
            if verbose:
                print("\n=== Final polish (excluding poorly-conditioned pairs) ===")
            res_good = joint_polish(
                float(x[0]),
                x[1::2][good_idx],
                x[2::2][good_idx],
                precomputed_good,
                bounds_kappa=bounds_kappa,
                smal_dif=smal_dif,
                delta=delta,
                verbose=verbose,
            )

            good_rows = extract_fitted_params(res_good, precomputed_good)
            for row in good_rows:
                row["excluded"] = False
            for idx in bad_idx:
                rows[idx]["excluded"] = True
            for j, i in enumerate(good_idx):
                rows[i] = good_rows[j]
        elif verbose:
            print(
                f"\nNo pairs exceed cond_threshold={cond_threshold:.1e}, "
                f"skipping exclusion step."
            )

    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Global MLD fitting: grid scan over kappa + block coordinate descent + joint polish."
    )
    parser.add_argument("config", type=str, help="Path to YAML config.")
    parser.add_argument("csv_out", type=str, help="Output CSV filename.")
    parser.add_argument(
        "--kappa_grid_n",
        type=int,
        default=20,
        help="Grid points along kappa (default 10).",
    )
    parser.add_argument(
        "--top_k",
        type=int,
        default=1,
        help="Joint-polish the k best grid points and keep the lowest loss (default 1).",
    )
    parser.add_argument(
        "--cond_threshold",
        type=float,
        default=np.inf,
        help="Exclude pairs whose (theta, log10_xi) Hessian condition number "
        "exceeds this value and run a final joint polish on the remaining pairs.",
    )
    parser.add_argument(
        "--censor_grid_n",
        type=int,
        default=0,
        help="test how many points should be censored from the mld (default 0).",
    )

    args = parser.parse_args()

    with open(args.config) as cf:
        cfg = yaml.safe_load(cf)

    censor_list = [0]
    if args.censor_grid_n != 0:
        linear_bin_width = 3
        limit_size = np.round(30).astype(int)
        linear_ind = pd.IntervalIndex.from_arrays(
            np.arange(
                1, (limit_size // linear_bin_width) * linear_bin_width, linear_bin_width
            ),
            np.arange(
                linear_bin_width,
                ((limit_size // linear_bin_width) + 1) * linear_bin_width,
                linear_bin_width,
            ),
            closed="both",
        )
        bins_center = [np.sqrt(x.right * x.left) for x in linear_ind]
        censor_list += bins_center[: args.censor_grid_n]

    with open(args.config, "r") as f:
        cfg = yaml.safe_load(f)
    for censor in censor_list:
        binned_mlds, L0_df, ncomps = load_mlds_from_config(cfg, censor)
        precomputed = precompute_pairs(binned_mlds, L0_df)

        rows = run_global_fitting(
            precomputed,
            kappa_grid_n=args.kappa_grid_n,
            top_k=args.top_k,
            cond_threshold=args.cond_threshold,
        )

        out_path = os.path.join(
            cfg["results_dir"], f"censor_{censor:.1f}_{args.csv_out}"
        )
        res_df = pd.DataFrame(rows)
        scale_mus(res_df, MUS)
        res_df["censor"] = censor
        res_df.to_csv(out_path)
        print(f"\nResults written to {out_path}")
