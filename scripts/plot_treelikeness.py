#!/usr/bin/env python3

"""
Loading ape::delta.plot result in python, then computing and plotting delta plots.
"""

import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch
import seaborn as sns
import seaborn.objects as so
import yaml
import jax.numpy as jnp
from mosaic_method.fitting import theoretical_mld
from mosaic_method.global_fitting_reparam_kappa_inv import (
    load_mlds_from_config,
    convert_taurho,
    precompute_pairs,
    pair_lag,
    SMAL_DIF,
    DELTA,
)

COLORS = {"mosaic-pair": "tab:orange", "mosaic-tree": "tab:blue"}
PAIR_MUC = 6e-11
PAIR_MUS = 5e-9
GLOBAL_MUS = 3.64e-9


def compute_residuals_clade(results_list, smal_dif=SMAL_DIF, delta=DELTA):
    """
    Compares per-pair fit residuals across clades as a boxplot.

    Parameters
    ----------
    results_list : list[tuple]
        list of (clade name, config_path, csv_path, fitting_procedure).
        ``config_path`` is passed to :func:`load_mlds_from_config`;
        ``csv_path`` is the output CSV written by the fitting pipeline
        (must contain columns kappa/muc&mus, theta/tau, xi/rho, species_1, species_2).
        example:
                results_list = [
                    ("Enterobacteriales", "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz_1_kappa/main_config_entero.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz_1_kappa/results.csv", "mosaic-tree"),
                    ("Enterobacteriales", "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz/main_config_entero.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz/results.csv", "mosaic-pair"),
                    ("Bacillales", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_1_kappa/main_config_bacillales_species_w_staph_global.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_1_kappa/results.csv", "mosaic-tree"),
                    ("Bacillales", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_w_staph_w_outgroup/main_config_bacillales_species_w_staph.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_w_staph_w_outgroup/results.csv", "mosaic-pair"),
                    ("Methanobacteriota", "/home/paulimer/Documents/results_bacteria_mlds/results_methano_1_kappa/main_config_methanobacterium_MLD.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_methano_1_kappa/results.csv", "mosaic-tree"),
                    ("Methanobacteriota", "/home/paulimer/Documents/results_bacteria_mlds/results_methano/main_config_methanobacterium_MLD.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_methano/results.csv", "mosaic-pair")
                ]

    Returns
    -------
    matplotlib.figure.Figure
    """
    clade_losses = []

    for clade_name, config_path, csv_path, fitting_procedure in results_list:
        with open(config_path, "r") as f:
            cfg = yaml.safe_load(f)
        cfg["censor"] = 0
        binned_mlds, L0_df, _ = load_mlds_from_config(cfg)
        precomputed = precompute_pairs(binned_mlds, L0_df)
        pairs, ml_pad, emp_pad, mask, L0s = precomputed

        fit_df = pd.read_csv(csv_path)
        if not "kappa" in fit_df.columns:
            fit_df = convert_taurho(fit_df)
        kappa = float(fit_df["kappa"].iloc[0])
        pair_index = {pair: i for i, pair in enumerate(pairs)}

        for row in fit_df.itertuples(index=False):
            pair_resid = {}
            key = (row.species_1, row.species_2)
            pair_resid["pair"] = str(key)
            pair_resid["clade"] = clade_name
            pair_resid["fitting procedure"] = fitting_procedure
            if key not in pair_index:
                print(f"key {key} not in pair_index")
                continue
            i = pair_index[key]
            p = jnp.array([row.theta, np.log10(row.xi)])
            loss_val, _ = pair_lag(
                p,
                kappa,
                ml_pad[i],
                emp_pad[i],
                mask[i],
                L0s[i : i + 1],
                smal_dif,
                delta,
            )
            pair_resid["residuals"] = float(loss_val)
            clade_losses.append(pair_resid)

    clade_losses_df = pd.DataFrame(clade_losses)
    return clade_losses_df


def plot_residuals_clades(
    sf: mpl.figure.SubFigure, residuals_df: pd.DataFrame, panel="X"
):
    """
    Plots the boxplots of pairwise summed residuals for each clades, for different fitting procedures.

    Parameters
    ----------
    sf: mpl.figure.SubFigure
    The subfigure to plot on.

    residuals_df: pd.DataFrame
    Dataframe with columns: ["residuals", "comparison", "clade", "fitting procedure"]
    """
    residuals_df.rename(
        {"fitting procedure": "Fitting procedure"}, inplace=True, axis=1
    )
    ax = sf.subplots()
    sns.boxplot(
        residuals_df,
        x="clade",
        y="residuals",
        hue="Fitting procedure",
        hue_order=["mosaic-pair", "mosaic-tree"],
        ax=ax,
        palette=COLORS,
        saturation=1,
    )

    target_rgba = to_rgba(COLORS["mosaic-pair"])
    pair_boxes = [
        p for p in ax.patches if np.allclose(p.get_facecolor(), target_rgba, atol=1e-6)
    ]

    for box in pair_boxes:
        box.set_hatch("//")
        box.set_edgecolor(
            "white"
        )  # hatch is drawn in the edge color, so make it visible
        box.set_linewidth(0.8)

    ax.set_ylabel("Per-pair summed residual")
    ax.set_xlabel("Clade")
    ax.set_title("Fit quality across clades")
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
    legend_handles = [
        Patch(
            facecolor=COLORS["mosaic-pair"],
            hatch="//",
            edgecolor="white",
            label="mosaic-pair",
        ),
        Patch(facecolor=COLORS["mosaic-tree"], label="mosaic-tree"),
    ]
    ax.legend(handles=legend_handles, title="Fitting procedure")
    # plt.xticks(rotation=30, ha="right")


def delta_plot_hist(
    sf: mpl.figure.SubFigure, res_dir: str, title: str = "title", panel="X"
):
    delta_hist_df = pd.read_csv(os.path.join(res_dir, "pair_global_hist_delta.csv"))
    delta_hist_df.rename(
        {"fitting_procedure": "Fitting procedure"}, inplace=True, axis=1
    )
    delta_hist_df["Fitting procedure"] = delta_hist_df["Fitting procedure"].astype(
        "category"
    )
    delta_hist_df["Fitting procedure"] = delta_hist_df[
        "Fitting procedure"
    ].cat.rename_categories({"pairwise": "mosaic-pair", "global": "mosaic-tree"})

    ax = sf.subplots()
    sns.lineplot(
        delta_hist_df,
        x="delta",
        y="counts",
        hue="Fitting procedure",
        hue_order=["mosaic-pair", "mosaic-tree"],
        style="Fitting procedure",
        style_order=["mosaic-tree", "mosaic-pair"],
        palette=COLORS,
        ax=ax,
    )
    ax.set_ylabel("Quartet counts")
    ax.set_xlabel("\u03b4")
    ax.set_title(title)
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
    # handles, labels = ax.get_legend_handles_labels()
    # order = ["mosaic-pair, mosaic-tree"]
    # ax.legend([handles[idx] for idx in order], [labels[idx] for idx in order])


def overall_figure(
    results_list: list[tuple], out_path: str, residuals_df: pd.DataFrame | None = None
):
    """
    Plots the figure with residuals and treelikeness.

    Parameters
    ----------
    results_list : list[tuple]
        list of (clade name, config_path, csv_path, fitting_procedure).
        ``config_path`` is passed to :func:`load_mlds_from_config`;
        ``csv_path`` is the output CSV written by the fitting pipeline
        (must contain columns kappa/muc&mus, theta/tau, xi/rho, species_1, species_2).
        example:
                results_list = [
                    ("Enterobacteriales", "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz_1_kappa/main_config_entero.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz_1_kappa/results.csv", "mosaic-tree"),
                    ("Enterobacteriales", "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz/main_config_entero.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz/results.csv", "mosaic-pair"),
                    ("Bacillales", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_1_kappa/main_config_bacillales_species_w_staph_global.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_1_kappa/results.csv", "mosaic-tree"),
                    ("Bacillales", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_w_staph_w_outgroup/main_config_bacillales_species_w_staph.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_w_staph_w_outgroup/results.csv", "mosaic-pair"),
                    ("Methanobacteriota", "/home/paulimer/Documents/results_bacteria_mlds/results_methano_1_kappa/main_config_methanobacterium_MLD.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_methano_1_kappa/results.csv", "mosaic-tree"),
                    ("Methanobacteriota", "/home/paulimer/Documents/results_bacteria_mlds/results_methano/main_config_methanobacterium_MLD.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_methano/results.csv", "mosaic-pair")
                ]


    """
    if residuals_df is None:
        residuals_df = compute_residuals_clade(results_list)
    with plt.style.context("seaborn-v0_8-paper"):
        sns.set_style("ticks", {"axes.grid": True})
        fig = plt.figure(layout="constrained", figsize=(8, 8), dpi=300)
        sf1, sf2, sf3, sf4 = fig.subfigures(2, 2, wspace=0.07).flatten()
        plot_residuals_clades(sf1, residuals_df, panel="A")
        delta_plot_hist(
            sf2,
            "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_1_kappa/",
            "Bacillales treelikeness",
            panel="B",
        )
        delta_plot_hist(
            sf3,
            "/home/paulimer/Documents/results_bacteria_mlds/results_methano_1_kappa/",
            "Methanobacteriota treelikeness",
            panel="C",
        )
        delta_plot_hist(
            sf4,
            "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz_1_kappa/",
            "Enterobacteriales treelikeness",
            panel="D",
        )

        fig.savefig(out_path)
        plt.close(fig)


if __name__ == "__main__":
    res_dir = "/home/paulimer/Documents/results_bacteria_mlds/results_methano"
    fig = plt.figure(layout="constrained", figsize=(10, 4))
    sf1, sf2 = fig.subfigures(1, 2, wspace=0.07)
    residuals_df = pd.read_csv(
        "/home/paulimer/Documents/results_bacteria_mlds/residuals_clades_fit.csv"
    )
    residuals_df["fitting procedure"] = residuals_df["fitting procedure"].astype(
        "category"
    )
    residuals_df["fitting procedure"] = residuals_df[
        "fitting procedure"
    ].cat.rename_categories({"pairwise": "mosaic-pair", "global": "mosaic-tree"})
    plot_residuals_clades(sf2, residuals_df)
    fig.savefig("/home/paulimer/Downloads/test.png")
