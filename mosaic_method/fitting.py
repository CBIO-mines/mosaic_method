#!/usr/bin/env python3

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.optimize import dual_annealing


def theoretical_mld(
    opt_pars, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit=False
):
    """
    Computes the theoretical match length distribution according to the paper.
    """
    ml_low = match_lengths - smal_dif
    ml_hi = match_lengths + smal_dif
    tau = np.double(10.0 ** opt_pars[0])
    rho = np.double(10.0 ** opt_pars[1])
    if L0_fit:
        L0 = np.double(10.0 ** opt_pars[2])
    mua = min(delta / tau, mus)

    def _mc_term(x):
        return (
            2
            * (
                (1 + x * mua * tau) * np.exp(-x * mua * tau)
                - (1 + x * muc * tau) * np.exp(-x * muc * tau)
            )
            / (x**2 * (muc**2 - mus**2) * tau**2)
        )

    mc = _mc_term(match_lengths)
    mc_low = _mc_term(ml_low)
    mc_hi = _mc_term(ml_hi)

    mc_diff = L0 * (mc_low + mc_hi - 2 * mc) / smal_dif**2
    np.nan_to_num(mc, copy=False)

    def _mh_low_branch(x):
        # used when tau < delta/mus
        return (
            x * (mus - muc) * tau + np.exp(-x * mus * tau) - np.exp(-x * muc * tau)
        ) / (x**2 * (mus**2 - muc**2) * tau)

    def _mh_high_branch(x):
        # used when tau >= delta/mus
        return (
            x * (mus - muc) * tau
            + (1 + x * (delta - mus * tau)) * np.exp(-x * delta)
            - np.exp(-x * muc * tau)
        ) / (x**2 * (mus**2 - muc**2) * tau)

    if tau < delta / mus:
        mh = _mh_low_branch(match_lengths)
        mh_low = _mh_low_branch(ml_low)
        mh_hi = _mh_low_branch(ml_hi)
    else:
        mh = _mh_high_branch(match_lengths)
        mh_low = _mh_high_branch(ml_low)
        mh_hi = _mh_high_branch(ml_hi)

    mh_diff = L0 * rho * (mh_low + mh_hi - 2 * mh) / smal_dif**2

    return mh_diff, mc_diff


def Lllocal(
    opt_pars, empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit=False
):
    """
    The squared relative difference to minimize.
    """
    mh_calc, mc_calc = theoretical_mld(
        opt_pars, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit
    )
    mt_calc = mh_calc + mc_calc
    return np.mean(((mt_calc - empirical_mld) / (mt_calc + empirical_mld)) ** 2.0)


def minus3Lllocal(
    opt_pars, empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit=False
):
    """
    The squared relative difference to minimize, only mh.
    """
    mh_calc, _ = theoretical_mld(
        opt_pars, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit
    )
    return np.mean(((mh_calc - empirical_mld) / (mh_calc + empirical_mld)) ** 2.0)


def minus4Lllocal(
    opt_pars, empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit=False
):
    """
    The squared relative difference to minimize, only mc.
    """
    _, mc_calc = theoretical_mld(
        opt_pars, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit
    )
    return np.mean(((mc_calc - empirical_mld) / (mc_calc + empirical_mld)) ** 2.0)


def _optimize(objective, opt_method, init_pars, args):
    if opt_method == "dual-annealing":
        return dual_annealing(objective, bounds=[(4, 10), (-12, -9)], args=args)
    else:
        return minimize(
            objective,
            init_pars,
            method=opt_method,
            args=args,
            tol=1e-8,
            options={"disp": False},
        )


def fit_params(
    opt_method,
    init_pars,
    empirical_mld,
    smal_dif,
    match_lengths,
    mus,
    muc,
    delta,
    L0,
    only_minus4=False,
):
    """
    Interface to minimize from scipy
    """
    if opt_method not in [
        "Nelder-Mead",
        "BFGS",
        "L-BFGS-B",
        "Powell",
        "COBYLA",
        "dual-annealing",
    ]:
        sys.exit("Unexistent/unimplemented optimization method requested")

    L0_fit = not L0
    args = (empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0, L0_fit)

    res_opt_minus4 = _optimize(minus4Lllocal, opt_method, init_pars, args)
    if only_minus4:
        return None, None, res_opt_minus4

    res_opt_full = _optimize(Lllocal, opt_method, init_pars, args)
    res_opt_minus3 = _optimize(minus3Lllocal, opt_method, init_pars, args)
    return res_opt_full, res_opt_minus3, res_opt_minus4


def write_results(res_opt, out_pars, L0, res_minus3_opt=None):
    """
    Writes the results of the fit to specified files.
    """
    opted_pars = res_opt.x
    obj_fun = res_opt.fun
    with open(out_pars, "w") as outfile:
        outfile.write("log10tau,log10rho,L0,minimum")
        if res_minus3_opt is not None:
            outfile.write(",minimum_minus3")
        outfile.write("\n")
        outfile.write(f"{opted_pars[0]}")
        for par in opted_pars[1:]:
            outfile.write(f",{par}")
        if L0:
            outfile.write(f",{L0}")
        outfile.write(f",{obj_fun}")
        if res_minus3_opt is not None:
            outfile.write(f",{res_minus3_opt.fun}")
        outfile.write("\n")


def plot_surface(
    min_logtau,
    max_logtau,
    min_logrho,
    max_logrho,
    num_points,
    output_file,
    empirical_mld,
    smal_dif,
    match_lengths,
    mus,
    muc,
    delta,
    L0,
    fitted_params=None,
):
    """Plots the Lllocal surface in a given region of the parameters to optimize."""

    x_range = np.linspace(min_logtau, max_logtau, num_points)
    y_range = np.linspace(min_logrho, max_logrho, num_points)
    x_vals, y_vals = np.meshgrid(x_range, y_range)

    # Calculate the corresponding Z values using Lllocal
    z_vals = np.zeros_like(x_vals)
    for i in range(len(x_range)):
        for j in range(len(y_range)):
            opt_pars = (x_vals[i, j], y_vals[i, j])
            z_vals[i, j] = Lllocal(
                opt_pars, empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0
            )

    # Create a 3D surface plot
    fig = plt.figure(figsize=(10, 8))

    for i, azim in enumerate(range(0, 280, 90)):
        ax = fig.add_subplot(221 + i, projection="3d", computed_zorder=False)
        if fitted_params is not None:
            ax.scatter(
                fitted_params[0],
                fitted_params[1],
                Lllocal(
                    fitted_params,
                    empirical_mld,
                    smal_dif,
                    match_lengths,
                    mus,
                    muc,
                    delta,
                    L0,
                ),
                color="red",
                s=100,
                label="Fitted Parameters",
                zorder=10,
            )
            ax.legend()
        ax.plot_surface(x_vals, y_vals, z_vals, cmap="viridis", zorder=1)
        ax.view_init(azim=azim, elev=60)
        ax.set_title(f"View with Azimuth = {azim}")
        # Set labels for the axes
        ax.set_xlabel("logtau")
        ax.set_ylabel("logrho")
        ax.set_zlabel("Lllocal")

    fig.tight_layout()

    # Save the plot to the specified output file
    fig.savefig(output_file, dpi=300)


def plot_residuals(
    opt_pars, empirical_mld, smal_dif, match_lengths, mus, muc, delta, L0, output_file
):
    """
    Plots the residuals of the fit, using Tommaso's normalization
    """
    mh_calc, mc_calc = theoretical_mld(
        opt_pars, smal_dif, match_lengths, mus, muc, delta, L0
    )
    mt_calc = mh_calc + mc_calc
    normalized_residuals = (mt_calc - empirical_mld) / np.sqrt(empirical_mld)
    fig = plt.Figure()
    ax = fig.subplots()
    ax.plot(match_lengths, normalized_residuals)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Match Length")
    ax.set_ylabel("Normalized Residuals")
    fig.savefig(output_file, dpi=300)
