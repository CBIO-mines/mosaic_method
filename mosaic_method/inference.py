#!/usr/bin/env python3

"""Main script for running the inference without snakemake."""

import concurrent.futures
from importlib.resources import files
import itertools
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import shutil
import subprocess as sp
import sys

import numpy as np
import pandas as pd
import yaml

from mosaic_method.aligning import create_lastz_db
from mosaic_method.fitting import fit_params
from mosaic_method.parsing import get_genome_comp, get_all_mlds, sum_mlds, bin_mld
from mosaic_method.plot_mld_fit import plot_mld_fit
from mosaic_method.length_analysis import length_analysis
from mosaic_method.get_L0 import L0_calc
from mosaic_method.pretreatment import pretreat
from mosaic_method.overall_inflexion import inflexions
from mosaic_method.compute_ani import compute_ani

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
            genomes = [fa for fa in os.listdir(cfg["genomes_dir"]) if fa.endswith(("fa", "fasta", "fna"))]
            masked_genomes = [fa for fa in os.listdir(masked_genomes_dir) if fa.endswith(("fa", "fasta", "fna"))]
            if len(masked_genomes) != len(genomes):
                pretreat.pretreat_genomes(
                    cfg["genomes_dir"],
                    cfg["taxon_csv"],
                    masked_genomes_dir,
                    above=2,
                    threads=cfg["max_threads"],
                    transition=True
                )
        else:
            masked_genomes_dir = cfg["genomes_dir"]


        # length analysis
        print("Analyzing genome lengths")
        os.makedirs(length_dir, exist_ok=True)
        length_analysis(
            cfg["taxon_csv"],
            cfg["cluster_name"],
            masked_genomes_dir,
            length_dir
        )

        # L0s
        print("Calculating L0s")
        L0_calc(
            length_dir,
            L0_csv
        )
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
            cfg["aligner"]
        )
        con.close()
    else:
        L0_df = pd.read_csv(L0_csv)

    # mlds
    print("Computing MLDs")
    taxon_df = pd.read_csv(cfg["taxon_csv"], index_col=0)
    level_list = sorted(taxon_df[cfg["cluster_name"]].unique())
    levels = list(itertools.combinations(level_list, 2))
    binned_mlds = {}
    for level in levels:
        genome_comps = get_genome_comp(level, taxon_df, cfg["cluster_name"])
        full_mld = get_all_mlds(genome_comps, database_path)
        summed_mld = sum_mlds(full_mld)
        binned_mld = bin_mld(
            summed_mld,
            linear_bin_width=3,
            limit_size=30.5,
            power_increment=0.1,
            ncomp=full_mld.shape[0],
            censor=float(cfg["censor"])
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
            fit_params,
            itertools.repeat(cfg["optim"]),
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
    # compute ani values
    res_df = compute_ani(
        res_df,
        database_path
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
            get_script_path("make_trees_cli.R"),
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

