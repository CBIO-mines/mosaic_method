#!/usr/bin/env python3

"""Main script for running the inference without snakemake."""

import argparse
import concurrent.futures
import itertools
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import subprocess as sp
import sys

import numpy as np
import pandas as pd
import yaml

script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(script_dir)
from parse import fun as parse_fun
from lastz_parallel_db import utils as lastz_utils
from fitting import fun as fit
from plot_mld_fit import plot_mld_fit
from pretreatment import pretreat
from length_analysis import length_analysis
from get_L0 import L0_calc
from overall_inflexion import inflexions


def run_inference(cfg, genomes_dir=None):
    # pretreatment
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
            masked_genomes_dir = cfg["genomes_dir"] + "_masked"
            os.makedirs(masked_genomes_dir, exist_ok=True)
            # TODO test granularly if already masked
            genomes = [fa for fa in os.listdir(cfg["genomes_dir"]) if fa.endswith(("fa", "fasta", "fna"))]
            masked_genomes = [fa for fa in os.listdir(masked_genomes_dir) if fa.endswith(("fa", "fasta", "fna"))]
            if len(masked_genomes) != len(genomes):
                pretreat.pretreat_genomes(
                    cfg["genomes_dir"],
                    cfg["taxon_csv"],
                    masked_genomes_dir,
                    os.path.join(script_dir, "pretreatment/tools"),
                    above=2,
                    threads=cfg["max_threads"],
                    transition=True
                )
        else:
            masked_genomes_dir = cfg["genomes_dir"]


        # length analysis
        print("Analyzing genome lengths")
        length_dir = os.path.join(cfg["results_dir"], "length_distributions")
        os.makedirs(length_dir, exist_ok=True)
        length_analysis(
            cfg["taxon_csv"],
            cfg["cluster_name"],
            masked_genomes_dir,
            length_dir
        )

        # L0s
        print("Calculating L0s")
        L0_csv = os.path.join(cfg["results_dir"], "L0.csv")
        L0_calc(
            length_dir,
            L0_csv
        )
        L0_df = pd.read_csv(L0_csv)


        # alignment
        print("Aligning genomes")
        database_path = os.path.join(cfg["results_dir"], cfg["database_name"])
        if os.path.exists(database_path):
            update_db = True
        else:
            update_db = False
        con = lastz_utils.create_lastz_db(
            cfg["taxon_csv"],
            masked_genomes_dir,
            cfg["cluster_name"],
            database_path,
            cfg["max_threads"],
            update_db,
            cfg["aligner"]
        )
        con.close()
    else:
        database_path = os.path.join(cfg["results_dir"], cfg["database_name"])

    # mlds
    print("Computing MLDs")
    taxon_df = pd.read_csv(cfg["taxon_csv"])
    level_list = sorted(taxon_df[cfg["cluster_name"]].unique())
    levels = list(itertools.combinations(level_list, 2))
    binned_mlds = {}
    for level in levels:
        genome_comps = parse_fun.get_genome_comp(level, cfg["taxon_csv"], "", cfg["cluster_name"], output_csv=False)
        full_mld = parse_fun.get_all_mlds(genome_comps, database_path, threads=cfg["max_threads"])
        summed_mld = parse_fun.sum_mlds(full_mld)
        binned_mld = parse_fun.bin_mld(
            summed_mld,
            linear_bin_width=3,
            limit_size=30.5,
            power_increment=0.1,
            ncomp=full_mld.shape[0]
        )
        binned_mlds[level] = binned_mld

    # fits
    print("Fitting MLDs")
    mus = float(cfg["mus"])
    muc = float(cfg["muc"])
    res_opt_full_dict = {}
    res_opt_minus3_dict = {}
    L0s = [
        L0_df.loc[(L0_df["bac1"] == level[0]) & (L0_df["bac2"] == level[1]), "L0"].values[0]
        for level in levels
    ]
    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg["max_threads"]) as executor:
        res_opt_full_list = executor.map(
            fit.fit_params,
            itertools.repeat("dual-annealing"),
            itertools.repeat(np.array([8, -8])),
            [binned_mlds[level]["freq"] for level in levels],
            itertools.repeat(0.1),
            [binned_mlds[level]["match_length"] for level in levels],
            itertools.repeat(mus),
            itertools.repeat(muc),
            itertools.repeat(cfg["delta"]),
            L0s
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

    # plot fits
    print("Plotting fits")
    os.makedirs(os.path.join(cfg["results_dir"], "fits"), exist_ok=True)
    for level in levels:
        current_L0 = L0_df.loc[
            (L0_df["bac1"] == level[0]) & (L0_df["bac2"] == level[1])
            , "L0"
        ].values[0]
        plot_mld_fit(
            binned_mlds[level],
            muc,
            mus,
            cfg["delta"],
            res_opt_full_dict[level].x,
            level,
            current_L0,
            os.path.join(cfg["results_dir"], "fits", f"{level[0]}_{level[1]}_fit.png")
        )


    # save results -------------------------------------------------------------
    print("Saving results")
    os.makedirs(os.path.join(cfg["results_dir"], "binned_mlds"), exist_ok=True)
    for level in levels:
        binned_mlds[level].to_csv(
            os.path.join(cfg["results_dir"], "binned_mlds", f"{level[0]}_{level[1]}.csv"),
            index=False
        )

    res_list = []
    for level in levels:
        res_list.append({
            "species_1": f"{level[0]}",
            "species_2": f"{level[1]}",
            "log10tau": res_opt_full_dict[level].x[0],
            "log10rho": res_opt_full_dict[level].x[1],
            "minimum": res_opt_full_dict[level].fun,
            "minimum_minus3": res_opt_minus3_dict[level].fun
        })
    res_df = pd.DataFrame(res_list)
    res_df["aligner"] = cfg["aligner"]
    res_df["delta"] = cfg["delta"]


    # inflexions
    print("Computing inflexions")
    res_df = inflexions(
        res_df,
        muc,
        mus,
        cfg["delta"],
        L0_df,
        0.1,
        min_r_infl=cfg["min_r_infl"]
    )
    res_df.to_csv(
        os.path.join(cfg["results_dir"], "results.csv"),
        index=False
    )
    # plot tree
    print("Plotting tree")
    if taxon_df[cfg["cluster_name"]].nunique() > 2:
        os.makedirs(os.path.join(cfg["results_dir"], "tree"), exist_ok=True)
        make_tree_args = [
            "Rscript",
            os.path.join(script_dir, "make_trees_cli.R"),
            os.path.join(cfg["results_dir"], "results.csv"),
            cfg["taxon_csv"],
            cfg["cluster_name"],
            cfg["results_dir"],
            length_dir,
            cfg["tree_annotation"],
            str(cfg["min_r_infl"])
        ]
        try:
            sp.run(make_tree_args, check=True)
        except sp.CalledProcessError as e:
            print(e)
    return res_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        Run the simulation and inference pipeline without snakemake.
        """
    )
    parser.add_argument(
        "config",
        type=str,
        help="The path to the configuration file."
    )
    parser.add_argument(
        "--genomes_dir",
        type=str,
        help="The directory containing the genomes (overwrites the one in the config file)",
        default=None
    )

    args = parser.parse_args()
    with open(args.config, "r") as config_file:
        cfg = yaml.safe_load(config_file)
    run_inference(cfg, args.genomes_dir)
