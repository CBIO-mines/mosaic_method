#!/usr/bin/env python3
"""Both for fig 3 and in general"""

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
from plot_treelikeness import compute_residuals_clade
from mosaic_method.plot_mld_fit import plot_mld_fit

COLORS = {"mosaic-pair": "tab:orange", "mosaic-tree": "tab:blue"}
PAIR_MUC = 6e-11
PAIR_MUS = 5e-9
GLOBAL_MUS = 3.64e-9


def choose_mlds(residuals_df: pd.DataFrame, by_diff: bool = False):
    if by_diff:
        try:
            residuals_diff = residuals_df.drop("Unnamed: 0", axis=1).pivot(
                index=["pair", "clade"], columns=["fitting procedure"]
            )
        except KeyError:
            residuals_diff = residuals_df.pivot(
                index=["pair", "clade"], columns=["fitting procedure"]
            )
            res_dif = (
                residuals_diff["residuals"]["mosaic-tree"]
                - residuals_diff["residuals"]["mosaic-pair"]
            )
        res_dif_df = pd.DataFrame(res_dif.reset_index())
        res_dif_df = res_dif_df.rename({0: "diff"}, axis=1)
        res_dif_df["pct_rank"] = res_dif_df["diff"].rank(pct=True)
        plot_pairs = {}
        for q in [0.25, 0.5, 0.75]:
            plot_pairs[q] = eval(
                res_dif_df.iloc[(res_dif_df["pct_rank"] - q).abs().argsort()[:1]][
                    "pair"
                ].values[0]
            )
        return plot_pairs

    else:
        global_residuals = residuals_df.loc[
            residuals_df["fitting procedure"] == "mosaic-tree"
        ].copy()
        global_residuals["pct_rank"] = global_residuals["residuals"].rank(pct=True)
        plot_pairs = {}
        for q in [0.1, 0.5, 0.9]:
            plot_pairs[q] = eval(
                global_residuals[global_residuals["pct_rank"] == q]["pair"].values[0]
            )
        return plot_pairs


def plot_example_fits(
    row_pair: dict,
    row_global: dict,
    ax: plt.Axes,
    binned_mlds_path: str,
    delta: float,
    muc: float,
    mus: float,
    quantile: float | None = None,
    panel: str = "X",
):
    binned_mld = pd.read_csv(
        os.path.join(
            binned_mlds_path,
            "binned_mlds",
            f"{row_pair['species_1']}_{row_pair['species_2']}.csv",
        )
    )
    max_x = binned_mld["match_length"].max()
    min_y = binned_mld[binned_mld["freq"] > 0]["freq"].min()
    match_length = np.logspace(0, np.log10(max_x), 1000)

    global_muc = row_global["muc"]
    global_mus = row_global["mus"]
    global_log10tau = row_global["log10tau"]
    global_log10rho = row_global["log10rho"]
    pair_log10tau = row_pair["log10tau"]
    pair_log10rho = row_pair["log10rho"]

    mh_pair, mc_pair = theoretical_mld(
        [pair_log10tau, pair_log10rho],
        0.1,
        match_length,
        mus,
        muc,
        delta,
        row_global["L0"],
        False,
    )
    mt_pair = mh_pair + mc_pair
    mh_global, mc_global = theoretical_mld(
        [global_log10tau, global_log10rho],
        0.1,
        match_length,
        global_mus,
        global_muc,
        delta,
        row_global["L0"],
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
    ax.plot(match_length, mt_pair, label="mosaic-pair fit", color="tab:orange", ls="--")
    ax.plot(match_length, mt_global, label="mosaic-tree fit", color="tab:blue", ls="-")
    ax.set_ylim(min_y / 10, max_y * 10)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.text(
        0.1,
        0.1,
        f"{row_pair['species_1']} vs \n{row_pair['species_2']}",
        style="italic",
        transform=ax.transAxes,
    )
    if quantile is not None:
        in_words = {25: "Most", 50: "Average", 75: "Least"}
        ax.set_title(f"{in_words[int(quantile)]} improved fits")
    ax.set_xlabel("Match length")
    ax.set_ylabel("Normalized counts")
    if panel != "X":
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
    ax.legend()


def plot_cases(
    results_list: list[tuple], out_path: str, residuals_df: pd.DataFrame | None = None
):
    """
    Plots the figure with global and not global fit.

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
    plot_pairs = choose_mlds(residuals_df, True)
    res_df_list = []
    res_dir_list = []
    for clade_name, config_path, csv_path, fitting_procedure in results_list:
        res_df = pd.read_csv(csv_path)
        res_dir = os.path.split(config_path)[0]
        res_dir_list.append(
            {
                "fitting procedure": fitting_procedure,
                "clade name": clade_name,
                "res dir": res_dir,
            }
        )
        res_df["fitting procedure"] = fitting_procedure
        res_df["clade name"] = clade_name
        res_df_list.append(res_df)
    big_res_df = pd.concat(res_df_list).reset_index(drop=True)
    res_dir_df = pd.DataFrame(res_dir_list)
    plot_pair_indexes = {}
    for q, pair in plot_pairs.items():
        plot_pair_indexes[q] = np.where(
            (big_res_df["species_1"] == pair[0]) & (big_res_df["species_2"] == pair[1])
        )

    fig, axes = plt.subplots(1, len(plot_pairs), layout="constrained", figsize=(12, 4))
    for ax, (q, plot_pair_ind), panel in zip(
        axes, plot_pair_indexes.items(), ["A", "B", "C"]
    ):
        plot_df = big_res_df.iloc[plot_pair_ind].copy()[
            [
                "species_1",
                "species_2",
                "log10tau",
                "log10rho",
                "fitting procedure",
                "muc",
                "mus",
                "L0",
                "clade name",
            ]
        ]
        row_global = (
            plot_df[plot_df["fitting procedure"] == "mosaic-tree"].iloc[0].to_dict()
        )
        row_pair = (
            plot_df[plot_df["fitting procedure"] == "mosaic-pair"].iloc[0].to_dict()
        )
        binned_mld_path = res_dir_df[
            (res_dir_df["fitting procedure"] == "mosaic-pair")
            & (res_dir_df["clade name"] == row_global["clade name"])
        ]["res dir"].values[0]

        plot_example_fits(
            row_pair,
            row_global,
            ax,
            binned_mld_path,
            DELTA,
            PAIR_MUC,
            PAIR_MUS,
            q * 100,
            panel,
        )
    fig.savefig(out_path, dpi=300)


def plot_all_of_one(
    results_list: list, which_species: str, out_path: str, wrap: int = 10
):
    """plots all mlds involving one species, as an example

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
                    ("Enterobacteriales", "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz/main_config_entero.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz/results.csv", "mosaic-pair")
                ]
        or:
                results_list = [
                    ("Bacillales", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_1_kappa/main_config_bacillales_species_w_staph_global.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_1_kappa/results.csv", "mosaic-tree"),
                    ("Bacillales", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_w_staph_w_outgroup/main_config_bacillales_species_w_staph.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_w_staph_w_outgroup/results.csv", "mosaic-pair"),
                ]
        or:
                results_list = [
                    ("Methanobacteriota", "/home/paulimer/Documents/results_bacteria_mlds/results_methano_1_kappa/main_config_methanobacterium_MLD.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_methano_1_kappa/results.csv", "mosaic-tree"),
                    ("Methanobacteriota", "/home/paulimer/Documents/results_bacteria_mlds/results_methano/main_config_methanobacterium_MLD.yaml", "/home/paulimer/Documents/results_bacteria_mlds/results_methano/results.csv", "mosaic-pair")
                ]
    """
    res_df_list = []
    res_dir_list = []
    for clade_name, config_path, csv_path, fitting_procedure in results_list:
        res_df = pd.read_csv(csv_path)
        res_dir = os.path.split(config_path)[0]
        res_dir_list.append(
            {
                "fitting procedure": fitting_procedure,
                "clade name": clade_name,
                "res dir": res_dir,
            }
        )
        res_df["fitting procedure"] = fitting_procedure
        res_df["clade name"] = clade_name
        res_df_list.append(res_df)
    big_res_df = pd.concat(res_df_list).reset_index(drop=True)
    res_dir_df = pd.DataFrame(res_dir_list)
    plot_df = big_res_df[
        (big_res_df["species_1"] == which_species)
        | (big_res_df["species_2"] == which_species)
    ]
    nplots = int(plot_df.shape[0] / 2)
    if nplots == 0:
        print(f"Error: not found {which_species}")
        print(
            f"Possible species names = {big_res_df['species_1'].drop_duplicates().tolist()}"
        )
    nrow = -(-nplots // wrap)  # ceil division
    fig, axes = plt.subplots(
        nrow, wrap, layout="constrained", figsize=(wrap * 4, nrow * 4)
    )
    flat_axes = axes.flatten()
    used = 0
    for ax, (pair, pair_df) in zip(
        flat_axes, plot_df.groupby(["species_1", "species_2"])
    ):
        row_global = (
            pair_df[pair_df["fitting procedure"] == "mosaic-tree"].iloc[0].to_dict()
        )
        row_pair = (
            pair_df[pair_df["fitting procedure"] == "mosaic-pair"].iloc[0].to_dict()
        )
        binned_mld_path = res_dir_df[
            (res_dir_df["fitting procedure"] == "mosaic-pair")
            & (res_dir_df["clade name"] == row_global["clade name"])
        ]["res dir"].values[0]

        plot_example_fits(
            row_pair,
            row_global,
            ax,
            binned_mld_path,
            DELTA,
            PAIR_MUC,
            PAIR_MUS,
        )
        used += 1
    for ax in flat_axes[used:]:
        ax.set_visible(False)
    fig.savefig(out_path, dpi=300)
