#!/usr/bin/env python3

"""Main script for running the inference without snakemake."""

import concurrent.futures
from importlib.resources import files
import multiprocessing
import itertools
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"
import shutil
import subprocess as sp
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from mosaic_method.aligning import create_lastz_db
from mosaic_method.fitting import fit_params
from mosaic_method.global_fitting_reparam_kappa_inv import (
    load_mlds_from_config,
    precompute_pairs,
    run_global_fitting,
    scale_mus,
    MUS,
)
from mosaic_method.parsing import get_genome_comp, sum_mlds, bin_mld
from mosaic_method.plot_mld_fit import plot_mld_fit
from mosaic_method.length_analysis import length_analysis
from mosaic_method.get_L0 import L0_calc
from mosaic_method.pretreatment import pretreat
from mosaic_method.overall_inflexion import inflexions
from mosaic_method.compute_ani import compute_ani

PAIRWISE = False


def get_script_path(script_name):
    """Get absolute path to script"""
    scripts_dir = files("mosaic_method").parent / "scripts"
    script_path = scripts_dir / script_name
    return str(script_path)


def run_inference(cfg, genomes_dir=None):
    # pretreatment
    L0_csv = os.path.join(cfg["results_dir"], "L0.csv")
    length_dir = os.path.join(cfg["results_dir"], "length_distributions")
    database_path = os.path.join(cfg["results_dir"], cfg["database_name"])
    if not cfg["alignment"] == "no":
        if genomes_dir:
            cfg["genomes_dir"] = genomes_dir
        if not os.path.exists(cfg["taxon_csv"]):
            cfg["taxon_csv"] = os.path.join(cfg["genomes_dir"], cfg["taxon_csv"])
            if not os.path.exists(cfg["taxon_csv"]):
                print(f"Taxon csv not found at {cfg['taxon_csv']}")
                sys.exit(1)

        if cfg["aligner"] == "lastz":
            print("Pretreating genomes")
            masked_genomes_dir = os.path.join(cfg["genomes_dir"], "masked_genomes")
            os.makedirs(masked_genomes_dir, exist_ok=True)
            # TODO test granularly if already masked
            genomes = [
                fa
                for fa in os.listdir(cfg["genomes_dir"])
                if fa.endswith(("fa", "fasta", "fna"))
            ]
            masked_genomes = [
                fa
                for fa in os.listdir(masked_genomes_dir)
                if fa.endswith(("fa", "fasta", "fna"))
            ]
            if len(masked_genomes) != len(genomes):
                pretreat.pretreat_genomes(
                    cfg["genomes_dir"],
                    cfg["taxon_csv"],
                    masked_genomes_dir,
                    above=2,
                    threads=cfg["max_threads"],
                    transition=True,
                )
        else:
            masked_genomes_dir = cfg["genomes_dir"]

        # length analysis
        print("Analyzing genome lengths")
        os.makedirs(length_dir, exist_ok=True)
        length_analysis(
            cfg["taxon_csv"], cfg["cluster_name"], masked_genomes_dir, length_dir
        )

        # L0s
        print("Calculating L0s")
        L0_calc(length_dir, L0_csv)
        L0_df = pd.read_csv(L0_csv)

        # alignment
        print("Aligning genomes")
        if os.path.exists(database_path):
            update_db = True
        else:
            update_db = False
        con = create_lastz_db(
            cfg["taxon_csv"],
            masked_genomes_dir,
            cfg["cluster_name"],
            database_path,
            cfg["max_threads"],
            update_db,
            cfg["aligner"],
        )
        con.close()
    else:
        L0_df = pd.read_csv(L0_csv)

    # fits
    print("Computing MLDs")
    taxon_df = pd.read_csv(cfg["taxon_csv"], index_col=0)
    level_list = sorted(taxon_df[cfg["cluster_name"]].unique())
    levels = list(itertools.combinations(level_list, 2))
    if PAIRWISE == True:
        # mlds
        binned_mlds = {}
        for level in levels:
            genome_comps = get_genome_comp(level, taxon_df, cfg["cluster_name"])
            summed_mld = sum_mlds(genome_comps, database_path)
            binned_mld = bin_mld(
                summed_mld,
                linear_bin_width=3,
                limit_size=30.5,
                power_increment=0.1,
                ncomp=len(genome_comps),
                censor=float(cfg["censor"]),
            )
            binned_mlds[level] = binned_mld

        print("Fitting MLDs")
        mus = float(cfg["mus"])
        muc = float(cfg["muc"])
        res_opt_full_dict = {}
        res_opt_minus3_dict = {}
        L0s = [
            L0_df.loc[
                (L0_df["bac1"] == level[0]) & (L0_df["bac2"] == level[1]), "L0"
            ].values[0]
            for level in levels
        ]
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=cfg["max_threads"],
            mp_context=multiprocessing.get_context("spawn"),
        ) as executor:
            res_opt_full_list = executor.map(
                fit_params,
                itertools.repeat(cfg["optim"]),
                itertools.repeat(np.array([8, -10])),
                [binned_mlds[level]["freq"] for level in levels],
                itertools.repeat(0.1),
                [binned_mlds[level]["match_length"] for level in levels],
                itertools.repeat(mus),
                itertools.repeat(muc),
                itertools.repeat(cfg["delta"]),
                L0s,
            )
        for level, res_opt_full in zip(levels, res_opt_full_list):
            res_opt_full_dict[level] = res_opt_full[0]
            res_opt_minus3_dict[level] = res_opt_full[1]

        # for level in levels:
        #     current_L0 = L0_df.loc[
        #         (L0_df["bac1"] == level[0]) & (L0_df["bac2"] == level[1])
        #         , "L0"
        #     ].values[0]
        #     res_opt_full, res_opt_minus3 = fit.fit_params(
        #         "dual-annealing",
        #         np.array([8, -8]),
        #         binned_mlds[level]["freq"],
        #         0.1,
        #         binned_mlds[level]["match_length"],
        #         mus,
        #         muc,
        #         cfg["delta"],
        #         current_L0
        #     )
        #     if not res_opt_full.success:
        #         print(f"Optimization failed for {level} full")
        #         sys.exit(1)
        #     res_opt_full_dict[level] = res_opt_full
        #     res_opt_minus3_dict[level] = res_opt_minus3
        # save results -------------------------------------------------------------
        print("Saving results")
        os.makedirs(os.path.join(cfg["results_dir"], "binned_mlds"), exist_ok=True)
        for level in levels:
            binned_mlds[level].to_csv(
                os.path.join(
                    cfg["results_dir"], "binned_mlds", f"{level[0]}_{level[1]}.csv"
                ),
                index=False,
            )

        res_list = []
        for level in levels:
            res_list.append(
                {
                    "species_1": f"{level[0]}",
                    "species_2": f"{level[1]}",
                    "log10tau": res_opt_full_dict[level].x[0],
                    "log10rho": res_opt_full_dict[level].x[1],
                    "minimum": res_opt_full_dict[level].fun,
                    "minimum_minus3": res_opt_minus3_dict[level].fun,
                }
            )
        res_df = pd.DataFrame(res_list)
        res_df["aligner"] = cfg["aligner"]
        res_df["delta"] = cfg["delta"]
    else:
        binned_mlds, L0_df, ncomps = load_mlds_from_config(cfg, 0)
        precomputed = precompute_pairs(binned_mlds, L0_df)

        rows = run_global_fitting(
            precomputed,
            kappa_grid_n=int(cfg["kappa_grid_n"]),
            top_k=int(cfg["top_k"]),
            cond_threshold=float(cfg["cond_threshold"]),
            out_path=os.path.join(cfg["results_dir"], "kappa_tranche.png"),
        )

        out_path = os.path.join(
            cfg["results_dir"], f"censor_{cfg['censor']:.1f}_results_global.csv"
        )
        res_df = pd.DataFrame(rows)
        ncomps_df = pd.DataFrame(ncomps, columns=["species_1", "species_2", "ncomp"])
        res_df = res_df.merge(ncomps_df, "inner", ["species_1", "species_2"])
        scale_mus(res_df, MUS)
        mus = MUS
        muc = res_df["muc"].values[0]
        res_df["censor"] = cfg["censor"]
        res_df.to_csv(out_path)
        print(f"\nResults written to {out_path}")

    # plot fits
    print("Plotting fits")
    os.makedirs(os.path.join(cfg["results_dir"], "fits"), exist_ok=True)
    fig, ax = plt.subplots()
    for level in levels:
        current_L0 = L0_df.loc[
            (L0_df["bac1"] == level[0]) & (L0_df["bac2"] == level[1]), "L0"
        ].values[0]
        plot_mld_fit(
            binned_mlds[level],
            muc,
            mus,
            cfg["delta"],
            res_df.loc[
                (res_df["species_1"] == level[0]) & (res_df["species_2"] == level[1]),
                ["log10tau", "log10rho"],
            ].values[0],
            level,
            current_L0,
            os.path.join(cfg["results_dir"], "fits", f"{level[0]}_{level[1]}_fit.png"),
            ax=ax,
        )
    plt.close(fig)

    # inflexions
    print("Computing inflexions")
    res_df = inflexions(
        res_df, muc, mus, cfg["delta"], L0_df, 0.1, min_r_infl=cfg["min_r_infl"]
    )
    # compute ani values
    res_df = compute_ani(res_df, database_path, taxon_df, cfg["cluster_name"])
    res_df.to_csv(os.path.join(cfg["results_dir"], "results.csv"), index=False)
    # plot tree
    print("Plotting tree")
    if taxon_df[cfg["cluster_name"]].nunique() > 2:
        os.makedirs(os.path.join(cfg["results_dir"], "tree"), exist_ok=True)
        make_tree_args = [
            "Rscript",
            get_script_path("make_trees_cli.R"),
            os.path.join(cfg["results_dir"], "results.csv"),
            cfg["taxon_csv"],
            cfg["cluster_name"],
            cfg["results_dir"],
            length_dir,
            cfg["tree_annotation"],
            str(cfg["min_r_infl"]),
        ]
        try:
            sp.run(make_tree_args, check=True)
        except sp.CalledProcessError as e:
            print(e)
    return res_df
