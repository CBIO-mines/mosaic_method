#!/usr/bin/env python3
"Plots the fit of the mosaic model to the MLDs. Translation of plot_mld_fit.R."

import os

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import string

from mosaic_method.fitting import theoretical_mld


def plot_mld_fit(
    binned_mld,
    muc,
    mus,
    delta,
    fitted_params,
    level,
    L0,
    outfile=None,
    ax=None | plt.Axes,
    panel=None | str,
):
    """Plots the fit of the mosaic model to the MLDs.

    If `ax` is provided, it is cleared and reused instead of creating a new
    Figure/Axes. Callers looping over many pairs should pass a shared `ax` to
    avoid piling up unclosed Figure/Axes/Legend reference cycles.
    """
    max_x = binned_mld["match_length"].max()
    min_y = binned_mld[binned_mld["freq"] > 0]["freq"].min()
    if binned_mld["freq"].sum() == 0:
        # mld is empty, do not plot
        return
    r = np.logspace(0, np.log10(max_x), 1000)
    mh, mc = theoretical_mld(fitted_params, 0.1, r, mus, muc, delta, L0, False)

    max_y = max([binned_mld["freq"].max(), mh.max(), mc.max()])
    owns_fig = ax is None
    if owns_fig:
        fig = plt.Figure()
        ax = fig.subplots()
    else:
        ax.clear()
        fig = ax.figure
    ax.plot(
        binned_mld["match_length"],
        binned_mld["freq"],
        "o",
        label="Observed",
        color="black",
    )
    ax.plot(r, mh, label="Transferred part", color="red")
    ax.plot(r, mc, label="Conserved part", color="blue")
    ax.plot(r, mc + mh, label="Complete fit", color="black")
    ax.set_ylim(min_y / 10, max_y * 10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_title(
        f"MLD fit for {level[0]} vs {level[1]}",
        fontsize=10,
    )
    ax.set_xlabel("Match length")
    ax.set_ylabel("Frequency")
    if panel is not None:
        ax.text(
            -0.15,
            1,
            panel,
            transform=ax.transAxes,
            fontweight="bold",
            va="top",
            ha="left",
            fontsize=10,
        )

    if outfile:
        fig.savefig(outfile, dpi=300)
        if owns_fig:
            plt.close(fig)
    else:
        plt.show()


def plot_some_mlds(
    comps: list[tuple[str]],
    res_df: pd.DataFrame,
    binned_mlds_path: str,
    outpath: str,
    muc=6e-11,
    mus=5e-9,
    delta=0.25,
    wrap=3,
    L0_df=None | pd.DataFrame,
):

    plot_df = pd.concat(
        [
            res_df[(res_df["species_1"] == s1) & (res_df["species_2"] == s2)]
            for s1, s2 in comps
        ]
    )
    nplots = len(comps)
    nrow = nplots // wrap
    fig, axes = plt.subplots(
        nrow, wrap, layout="constrained", figsize=(wrap * 4, nrow * 4)
    )
    try:
        muc = plot_df["muc"].values[0]
    except KeyError:
        pass
    try:
        mus = plot_df["mus"].values[0]
    except KeyError:
        pass

    for panel, ax, (pair, pair_df) in zip(
        list(string.ascii_uppercase)[:nplots],
        axes.flatten(),
        plot_df.groupby(["species_1", "species_2"]),
    ):
        binned_mld = pd.read_csv(
            os.path.join(
                binned_mlds_path,
                f"{pair[0]}_{pair[1]}.csv",
            )
        )
        fitted_params = pair_df[["log10tau", "log10rho"]].values[0]
        try:
            L0 = plot_df["L0"].values[0]
        except KeyError:
            L0 = L0_df.loc[
                (L0_df["bac1"] == pair[0]) & (L0_df["bac2"] == pair[1]), "L0"
            ].values[0]
        plot_mld_fit(
            binned_mld,
            muc,
            mus,
            delta,
            fitted_params,
            pair,
            L0,
            ax=ax,
            panel=panel,
        )
    fig.savefig(outpath, dpi=300)


def plot_global_pairwise(row, ax, binned_mlds_path, delta, muc, mus):
    binned_mld = pd.read_csv(
        os.path.join(binned_mlds_path, f"{row['species_1']}_{row['species_2']}.csv")
    )
    max_x = binned_mld["match_length"].max()
    min_y = binned_mld[binned_mld["freq"] > 0]["freq"].min()
    r = np.logspace(0, np.log10(max_x), 1000)

    global_muc = row["muc"]
    try:
        global_mus = row["mus"]
    except KeyError:
        global_mus = 3.64e-9
        global_log10tau = np.log10(row["tau"])
        global_log10rho = np.log10(row["rho"])
        pair_log10tau = row["log10tau"]
        pair_log10rho = row["log10rho"]

    mh_pair, mc_pair = theoretical_mld(
        [pair_log10tau, pair_log10rho], 0.1, r, mus, muc, delta, row["L0"], False
    )
    mt_pair = mh_pair + mc_pair
    mh_global, mc_global = theoretical_mld(
        [global_log10tau, global_log10rho],
        0.1,
        r,
        global_mus,
        global_muc,
        delta,
        row["L0"],
        False,
    )
    mt_global = mh_global + mc_global
    max_y = max([binned_mld["freq"].max(), mt_pair.max(), mt_global.max()])

    ax.plot(
        binned_mld["match_length"],
        binned_mld["freq"],
        "o",
        label="Observed",
        color="black",
    )
    ax.plot(r, mt_pair, label="Theta ratio from litterature", color="tab:orange")
    ax.plot(r, mt_global, label="fitted Theta ratio", color="tab:blue")
    ax.set_ylim(min_y / 10, max_y * 10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()
    ax.set_title(f"MLD fit for {row['species_1']} vs {row['species_2']}")
    ax.set_xlabel("Match length")
    ax.set_ylabel("Frequency")


def plot_global_fit(
    global_fit_csv, results_csv, muc, mus, binned_mlds_path, delta, outdir
):
    """plots the global fits and data, alongside the pairwise fit ?"""
    os.makedirs(outdir, exist_ok=True)
    global_fit_df = pd.read_csv(global_fit_csv).drop(["L0"], axis=1)
    results_df = pd.read_csv(results_csv)
    global_pair_df = pd.merge(
        global_fit_df,
        results_df,
        left_on=["species_1", "species_2"],
        right_on=["species_1", "species_2"],
    )
    fig, ax = plt.subplots()
    for _, row in global_pair_df.iterrows():
        ax.clear()
        plot_global_pairwise(row, ax, delta, binned_mlds_path, muc, mus)
        fig.savefig(
            os.path.join(outdir, f"{row['species_1']}_{row['species_2']}.png"), dpi=300
        )
