#!/usr/bin/env python3

"""Subsamples the genome pairs MLDs to normalize kappa inference across clade comparisons"""

from sqlalchemy.dialects.oracle.dictionary import DB_LINK_PLACEHOLDER

import argparse
import os
import yaml
import sys
from mosaic_method.parsing import bin_mld, sum_mlds, get_genome_comp
from mosaic_method.inference import run_inference
import tempfile
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr, wilcoxon

CLADE_ANNOTATION_MAPPING = {
    "Enterobacteriales": "misha_annotation",
    "Bacillales": "species.gtdb",
    "Methanobacteriota": "species.gtdb",
}


def genome_per_taxon(taxon_df: pd.DataFrame, clade_annot: str):
    grouped_df = taxon_df.groupby(clade_annot).count()
    avg_genome = grouped_df["genome"].mean()
    return avg_genome


def subsample_target(taxon_df: pd.DataFrame, target: float, clade_annot: str):
    """
    subsamples a taxon_df to a target level.

    Parameters
    ----------

    taxon_df: the taxon_df

    target: the average number of genomes per taxon of the target clade

    clade_annot: the level at which to sample
    """
    nsamp = int(min(taxon_df.groupby(clade_annot).count()["genome"].min(), target))
    subsampled_df = taxon_df.groupby(clade_annot).sample(nsamp)
    return subsampled_df


def merge_subs_ori(
    res_subs: pd.DataFrame,
    res_original: pd.DataFrame,
):
    res_df_subs = res_subs[["species_1", "species_2", "tau", "kappa", "ncomp"]]
    res_df_ori = res_original[["species_1", "species_2", "tau", "kappa", "ncomp"]]
    res_comp = pd.merge(
        res_df_subs,
        res_df_ori,
        "inner",
        on=["species_1", "species_2"],
        suffixes=("_subs", "_ori"),
    )
    return res_comp


def repeat_subsampling(
    ori_dict: dict[str, list[pd.DataFrame | dict]],
    target_taxon: pd.DataFrame,
    repeat: int,
):
    res_subsampled_df_dic = {}
    for ori_name, (ori_taxon, ori_res, model_cfg) in ori_dict.items():
        ori_res_dir = model_cfg["results_dir"]
        res_rep_df_dic = {}
        for i in range(repeat):
            target_gtax = genome_per_taxon(
                target_taxon, CLADE_ANNOTATION_MAPPING["Methanobacteriota"]
            )
            subsampled_df = subsample_target(
                ori_taxon, target_gtax, CLADE_ANNOTATION_MAPPING[ori_name]
            )
            with tempfile.TemporaryDirectory() as temp_dir:
                # symlink needed files
                os.symlink(
                    os.path.join(ori_res_dir, model_cfg["database_name"]),
                    os.path.join(temp_dir, model_cfg["database_name"]),
                )
                os.symlink(
                    os.path.join(ori_res_dir, "L0.csv"),
                    os.path.join(temp_dir, "L0.csv"),
                )
                taxon_csv = tempfile.NamedTemporaryFile(
                    mode="w", suffix=".csv", delete=False
                )
                subsampled_df.to_csv(taxon_csv)
                model_cfg["taxon_csv"] = taxon_csv.name
                model_cfg["results_dir"] = temp_dir
                res_df = run_inference(model_cfg)
                taxon_csv.close()
                os.remove(taxon_csv.name)
            res_rep_df_dic[i] = merge_subs_ori(res_df, ori_res)
        res_subsampled_df_dic[ori_name] = pd.concat(res_rep_df_dic)
    res_subsampled_df = pd.concat(res_subsampled_df_dic)
    return res_subsampled_df


def only_table(res_subsampled_df: pd.DataFrame, out_path: str):
    res_subsampled_df.index.set_names(["Clade", "repeat", "comp number"], inplace=True)
    res_subsampled_df.reset_index(inplace=True)
    means = {
        "kappa_ori": res_subsampled_df.groupby(["Clade"])["kappa_ori"].mean(),
        "kappa_mean": res_subsampled_df.groupby(["Clade"])["kappa_subs"].mean(),
        "kappa_std": res_subsampled_df.groupby(["Clade"])["kappa_subs"].std(),
        "ncomp_subs": res_subsampled_df[res_subsampled_df["repeat"] == 0]
        .groupby(["Clade"])["ncomp_subs"]
        .sum(),
        "ncomp_ori": res_subsampled_df[res_subsampled_df["repeat"] == 0]
        .groupby(["Clade"])["ncomp_ori"]
        .sum(),
    }
    summary_df = pd.concat(means, axis=1)
    summary_df.to_csv(out_path)


def only_plot(res_subsampled_df: pd.DataFrame, out_path: list[str]):
    res_subsampled_df.index.set_names(["Clade", "repeat", "comp number"], inplace=True)
    res_subsampled_df.reset_index(level=["repeat"], inplace=True)
    res_subsampled_df.reset_index(level=["comp number"], inplace=True, drop=True)
    ncols = res_subsampled_df["repeat"].max() + 1
    for out, (k, gdf) in zip(out_path, res_subsampled_df.groupby("Clade")):
        fig, axes1d = plt.subplots(
            1, ncols, figsize=(ncols * 4, 4), layout="constrained"
        )
        for ax, (rep, repdf) in zip(axes1d, gdf.groupby("repeat")):
            pear = pearsonr(repdf["tau_ori"], repdf["tau_subs"]).statistic
            spear = spearmanr(repdf["tau_ori"], repdf["tau_subs"]).statistic

            ax.scatter(repdf["tau_ori"], repdf["tau_subs"])
            lims = [
                min(repdf["tau_ori"].min(), repdf["tau_subs"].min()),
                max(repdf["tau_ori"].max(), repdf["tau_subs"].max()),
            ]
            ax.plot(lims, lims, "k--", linewidth=0.8)
            ax.text(0.2, 0.8, f"r={pear:.2f}\nρ={spear:.2f}", transform=ax.transAxes)
            per_subs = repdf["ncomp_subs"].sum() / repdf["ncomp_ori"].sum() * 100
            ax.set_title(
                f"κ {repdf['kappa_ori'].values[0]:.2f} -> {repdf['kappa_subs'].values[0]:.2f}"
            )
            ax.set_xlabel("Original distance")
            ax.set_ylabel("Subsampled inference distance")
        fig.savefig(out, dpi=300)


def main(
    ori_dict: dict[str, list[str]],
    target_taxon: str,
    rep: int,
    out_figure: list[str],
    out_csv: str,
):
    """
    Computes and plots the subsampling.

    Parameters
    ----------
    ori_dir: dict
        contains all the info for reinferring. Example:
        ori_dict = {"Bacillales": ["/home/paulimer/Data/bacillales_species/species_taxon_w_staph_w_outg.csv", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_1_kappa/results.csv", "/home/paulimer/Documents/results_bacteria_mlds/results_bacillales_species_1_kappa/main_config_bacillales_species_w_staph_global.yaml"],
                    "Enterobacteriales": ["/home/paulimer/Data/latest_entero_genomes_no_plasmid/misha_taxon_w_annotation_v3.csv", "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz_1_kappa/results.csv", "/home/paulimer/Documents/results_bacteria_mlds/results_entero_v3_lastz_1_kappa/main_config_entero.yaml"]}
    target_taxon: str
        path to the target taxon
    rep: int
       number of repetitions
    out: list[str]
       where to plot each clade subsampling
    """
    ori_df_dict = {}
    for ori_name, (ori_taxon, ori_res, model_cfg_path) in ori_dict.items():
        ori_df_dict[ori_name] = [pd.read_csv(ori_taxon)]
        ori_df_dict[ori_name].append(pd.read_csv(ori_res))
        with open(model_cfg_path, "r") as f:
            model_cfg = yaml.safe_load(f)
        ori_df_dict[ori_name].append(model_cfg)
    target_taxon_df = pd.read_csv(target_taxon)
    res_subsampled_df = repeat_subsampling(ori_df_dict, target_taxon_df, rep)
    res_subsampled_df.to_csv(out_csv)
    only_plot(res_subsampled_df, out_figure)
    only_table(
        res_subsampled_df,
        os.path.join(os.path.split(out_csv)[0], "subsampling_summary.csv"),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
        Run the inference subsampling the number of genomes to the level of genomes/species of another example
        """)
    parser.add_argument("conf", type=str, help="The path to the config file.")
    args = parser.parse_args()
    with open(args.conf, "r") as f:
        cfg = yaml.safe_load(f)
    main(
        cfg["subs_tax"],
        cfg["target_taxon"],
        cfg["repeats"],
        cfg["out_figure"],
        cfg["out_csv"],
    )
