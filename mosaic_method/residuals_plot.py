#!/usr/bin/env python3
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import yaml

from mosaic_method.fitting import theoretical_mld
from mosaic_method.parsing import get_all_mlds, get_genome_comp, sum_mlds

"""
Tommaso idea:
this is what I do:

- say you have a vector w of observed counts in each bin, and one w_th of predicted counts
- I assume the observed is Poisson distributed, so for each bin the std dev is sigma[i] = sqrt(w[i])
- the normalized residual I compute is (w[i] - w_th[i]) / sigma[i]
- then, usually my plot are log binned, so when plotting the MLD the values of y shown in the plot are actually w[i] / width[i] where width[i] is the bin width.
  So also the plotted error needs to be divided by width (of course for the residuals it does not matter)

"""

def get_theoretical_and_observed(summed_mld, logtau, r_infl, muc, mus, delta, L0, ncomp):
    # DONE: plot this again to be sure
    summed_mld = summed_mld[summed_mld["match_length"] < r_infl].copy()
    summed_mld["freq"] = summed_mld["freq"] / ncomp
    _, summed_mld["th_freq"] = theoretical_mld(
        [logtau, -20],
        0.1,
        summed_mld["match_length"].values,
        mus,
        muc,
        delta,
        L0
    )
    return summed_mld


def bin_3_mld(summed_mld):
    summed_mld_max = summed_mld["match_length"].max() + 3
    interval_ind = pd.IntervalIndex.from_arrays(
        np.arange(1, (summed_mld_max // 3) * 3, 3),
        np.arange(3, ((summed_mld_max // 3) + 1) * 3, 3),
        closed="both"
    )
    cuts = pd.cut(summed_mld["match_length"], bins=interval_ind)
    res = summed_mld.groupby(cuts, observed=False)[["freq", "th_freq"]].sum().reset_index()
    res["freq"] = res.apply(lambda x: x["freq"]/((x["match_length"].right - x["match_length"].left)), axis=1)
    res["th_freq"] = res.apply(lambda x: x["th_freq"]/((x["match_length"].right - x["match_length"].left)), axis=1)

    res["match_length"] = res["match_length"].apply(lambda x: (x.left+x.right)/2)
    return res



def calc_residuals(mld_db_path, taxon_df, results_df, main_cfg, L0_df=None, chi2=True):
    if L0_df is not None:
        results_df = pd.merge(results_df, L0_df, "inner", left_on=["species_1", "species_2"], right_on=["bac1", "bac2"])

    with open(main_cfg, "r") as f:
        cfg = yaml.safe_load(f)
    delta = float(cfg["delta"])
        
    summed_mld_w_th_list = []
    simulation = False
    # first summed mlds
    for _, row in results_df.iterrows():
        genomes = get_genome_comp(
            (row["species_1"], row["species_2"]),
            taxon_df,
            cfg["cluster_name"],
        )
        ncomp = len(genomes)
        mld_df = get_all_mlds(genomes, mld_db_path)
        summed_mld = sum_mlds(mld_df)
        muc = row["empirical_muc"]
        mus = row["empirical_mus"]
        try:
            summed_mld = get_theoretical_and_observed(summed_mld, row["log10tau"], row["r_infl"], muc, mus, delta, row["L0"], ncomp)
        except KeyError:
            # case where log10tau is called sim_tau in simulations (and is really the theoretical)
            # also r_infl = infinity, let's say 1e4
            # TODO actually the fit is also interesting
            summed_mld = get_theoretical_and_observed(summed_mld, row["sim_tau"], 1e4, muc, mus, delta, row["L0"], ncomp)
            simulation = True
            # summed_mld = get_theoretical_and_observed(summed_mld, row["fit_tau"], 1e4, muc, mus, delta, row["L0"], ncomp)
        bin3_comp_df = bin_3_mld(summed_mld)
        both_mld = pd.concat({"summed": summed_mld, "bin3": bin3_comp_df}).reset_index(names=["type", "drop"]).drop(["drop"], axis=1)
        both_mld = pd.concat([both_mld, pd.concat([pd.DataFrame(row).T] * both_mld.shape[0], axis=0).reset_index(drop=True)], axis=1)
        summed_mld_w_th_list.append(both_mld)

    bdl_mld = pd.concat(summed_mld_w_th_list)
    bdl_mld = bdl_mld[bdl_mld["freq"] != 0]
    bdl_mld["std"] = np.sqrt(bdl_mld["freq"])
    bdl_mld["residuals"] = (bdl_mld["freq"] - bdl_mld["th_freq"]) / bdl_mld["std"]
    if not simulation:
        for k, g in bdl_mld.groupby("type"):
            g.to_csv(os.path.join(os.path.split(mld_db_path)[0], f"{k}_mlds.csv"))
    else:
        for k, g in bdl_mld.groupby("type"):
            mld_dir, mld_name = os.path.split(mld_db_path)
            g.to_csv(os.path.join(mld_dir, f"{mld_name}_{k}_mlds_residuals.csv"))

    return bdl_mld
       

def chi_square(all_summed_mlds, min_r=0, group_by_keys=["exp", "species_1", "species_2"]):
    """
    Computes the 'chi-square' between the fit and the observed data.
    Normalizes by the number of comparisons and each comparison's number of bins
    
    Parameters

    ----------
    all_summed_mlds: pd.DataFrame
    df of summed/binned mlds and experiments

    min_r: float
    the minimum length from which to consider the residuals

    group_by_keys: list
    The keys to group the mld by. Can be more numerous in case of simulations.

    Returns
    -------
    res: pd.DataFrame
    the dataframe with the sums of residuals for each group.
    """
    residual_list = []
    for exp, all_summed_mld in all_summed_mlds.groupby(group_by_keys):
        all_summed_mld_min_r = all_summed_mld[all_summed_mld["match_length"] > min_r].copy()
        res_exp = [*exp, all_summed_mld_min_r["residuals"].abs().sum(), all_summed_mld_min_r.shape[0]]
        residual_list.append(res_exp)

    res = pd.DataFrame(residual_list, columns=[*group_by_keys, "chi2", "nbins"])
    return res
        

def plot_resid(all_summed_mld, outfile, min_r=0, max_r=200, quantile_divergence=3):
    plt.clf()
    all_summed_mld_min_r = all_summed_mld[(all_summed_mld["match_length"] > min_r) & (all_summed_mld["match_length"] < max_r)].copy()
    if quantile_divergence > 1:
        aver_div_df = all_summed_mld_min_r[["species_1", "species_2", "average_divergence"]].drop_duplicates().reset_index(drop=True)
        aver_div_df["average_divergence"] = aver_div_df["average_divergence"].astype(float)
        aver_div_df["divergence_quantile"] = pd.qcut(aver_div_df["average_divergence"], quantile_divergence)

        all_summed_mld_min_r = pd.merge(all_summed_mld_min_r, aver_div_df, how="left", on=["species_1", "species_2"])
        ax = sns.boxplot(all_summed_mld_min_r, x="match_length", y="residuals", hue="divergence_quantile", native_scale=True, showfliers=False)
        ax.set_xscale("log")
        quantile_path = os.path.join(os.path.split(outfile)[0], f"{os.path.splitext(outfile)[0]}_quantile_{quantile_divergence}{os.path.splitext(outfile)[1]}")
        plt.savefig(quantile_path, dpi=300)
    plt.clf()
    ax = sns.boxplot(all_summed_mld_min_r, x="match_length", y="residuals", native_scale=True, showfliers=False)
    ax.set_xscale("log")
    plt.savefig(outfile, dpi=300)



if __name__ == "__main__":
    logb_residuals = pd.read_csv("/home/paulimer/Documents/results_bacteria_mlds/simulated_results/debug_runs_dual_annealing/log_brownian__unbounded___1e-11_1e-09__entero/dbs/2.28e+08_lastz_nu_1.e-5.db_bin3_mlds_residuals.csv", index_col=0)
    linearrw_residuals = pd.read_csv("/home/paulimer/Documents/results_bacteria_mlds/simulated_results/debug_runs_dual_annealing/random_walk__linear___1e-11_1e-09__entero/dbs/2.28e+08_lastz_rw_step_fraction_1e4.db_bin3_mlds_residuals.csv", index_col=0)
    both_residuals = logb_residuals.merge(linearrw_residuals, "outer", ["species_1", "species_2", "match_length", "type", "sim_tau", "tree_height", "L0", "th_freq"], suffixes=("_logb", "_linearrw")).fillna(0)
    test_species = ["Klebsiellapneumoniae", "Salmonellaenterica"]
    test_residuals = both_residuals.query(
        "species_1 in @test_species and species_2 in @test_species"
    )
    fig, ax = plt.subplots()
    ax.plot(test_residuals["match_length"], test_residuals["freq_logb"], label="logbrownian")
    ax.plot(test_residuals["match_length"], test_residuals["freq_linearrw"], label="random walk")
    ax.plot(test_residuals["match_length"], test_residuals["th_freq"], label="theoretical")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(1, None)
    ax.legend()
    plt.show()
