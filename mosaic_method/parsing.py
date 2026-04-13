#!/usr/bin/env python3

import concurrent.futures
import itertools
import os


import pandas as pd
import numpy as np
import sqlite3


def get_genome_comp(species, taxon_df, level):
    """
    Gets a list of genome comparisons from a taxon dataframe and two taxa, and "level".
    """
    genomes_dic = {}
    for sp in sorted(species):
        genomes_dic[sp] = list(taxon_df[taxon_df[level] == sp]["genome"])
    genomes_1 = sorted(genomes_dic[species[0]])
    genomes_2 = sorted(genomes_dic[species[1]])
    res = []
    for g_1, g_2 in itertools.product(genomes_1, genomes_2):
        res += [(g_1, g_2)]
    return res


def get_all_mlds(genome_comps, lastz_db_path):
    """
    Parses a database and returns the resulting mlds concatenated.
    A version without parallel reads of the db but parallel decoding
    """
    mld_comps = {}
    sqlite_con = sqlite3.connect(lastz_db_path)
    cur = sqlite_con.cursor()
    genome_comps = [sorted(pair) for pair in genome_comps]

    cur.execute("CREATE TEMP TABLE temp_pairs (genome1 TEXT, genome2 TEXT);")
    cur.executemany("INSERT INTO temp_pairs(genome1, genome2) VALUES (?, ?)", genome_comps)
    cur.execute("SELECT l.genome1, l.genome2, l.count_array FROM lastz l JOIN temp_pairs t ON l.genome1 = t.genome1 AND l.genome2 = t.genome2;")
    rows = cur.fetchall()
    cur.execute("DROP TABLE temp_pairs;")

    for genome_1, genome_2, byte_array in rows:
       comp_vs = f"{genome_1}_vs_{genome_2}"
       mld_array = np.frombuffer(byte_array, dtype=np.int64)
       mld_comps[comp_vs] = mld_array


    max_len = max([len(l) for l in mld_comps.values()])
    df_mlds = pd.DataFrame.from_dict(mld_comps, orient="index", columns=range(1, max_len + 1))
    df_mlds = df_mlds.fillna(int(0))
    df_mlds = df_mlds.reset_index(names=["comp"])
    return df_mlds



# def parse_florian_mld(path):
#     """
#     Parses Florian's assembly-wise MLD files and merges them into a taxon-wide full_mld.

#     Parameters
#     ----------
#     path: str
#     the path to the directory (named bac1_bac2) where the MLDs are stored

#     Returns
#     -------
#     A pandas.DataFrame containing a taxa comparison's MLD
#     """
#     comp_files = [mld_file for mld_file in os.listdir(path) if mld_file.endswith(".MLD")]
#     res = {}
#     for c_f in comp_files:
#         tmp_dic = {}
#         c_f_full = os.path.join(path, c_f)
#         with open(c_f_full, "r") as filein:
#             for line in filein:
#                 try:
#                     tmp_dic[line.split()[0]] = int(line.split()[1])
#                 except IndexError:
#                     print(line)
#                 except:
#                     print("other error")
#                 finally:
#                     continue
#         res[c_f.split(".")[0].replace("-", "_")] = tmp_dic
#     res_df = pd.DataFrame.from_dict(res, orient="index")
#     res_df = res_df.rename_axis("comp").reset_index()
#     return res_df


def sum_mlds(mld_comp_df):
    """
    Takes a df of mlds by comparison and sums it.
    """
    summed_colmuns = mld_comp_df.drop(labels="comp", axis=1).sum(axis=0, skipna=True)
    summed_df = summed_colmuns.reset_index()
    summed_df.columns = ["match_length", "freq"]
    summed_df = summed_df.astype({"match_length" : "int64"})
    return summed_df.sort_values(by=["match_length"]).reset_index(drop=True)


# def bin_mld(summed_df, linear_bin_width, limit_size, power_increment, ncomp):
#     """
#     Bins and normalizes a summed mld according to a specific pattern.

#     Parameters
#     ----------
#     summed_df: pd.DataFrame
#     an unbinned dataframe with a column "match_length" and a column "freq"
#     linear_bin_width: float
#     the width of the bin in the linear part of bin vector
#     limit_size: float
#     the limit at which the bin vector switches from linear to log
#     power_increment: float
#     the "bin width" of the log part
#     ncomp: float
#     the number of comparison summed here

#     Returns
#     -------
#     a binned pd.Dataframe with columns :
#     - "match_length" containing the geometric mean of the bin
#     - "freq" containing the counts corresponding to the bin
#     """
#     match_bin = list(np.arange(0.5, limit_size, linear_bin_width))
#     initial_len = len(match_bin)
#     cur_power = 0.1
#     while match_bin[-1] < max(summed_df["match_length"]):
#         match_bin += [match_bin[initial_len - 1]*10**(cur_power)]
#         cur_power += power_increment

#     res = pd.DataFrame.from_dict(
#         {"match_length" : match_bin,
#          "freq" : [0.0] * len(match_bin)}
#     )
#     summed_row = 0
#     binned_row = 0
#     while(summed_row < summed_df.shape[0]):
#         if(summed_df.loc[summed_row, "match_length"] <= res.loc[binned_row, "match_length"]):
#             res.loc[binned_row, "freq"] += summed_df.loc[summed_row, "freq"]
#             summed_row += 1
#         else:
#             binned_row += 1
#             if binned_row > res.shape[0] - 1:
#                 break

#     for binned_row in range(res.shape[0] - 1):
#         len_bin = res.loc[binned_row + 1, "match_length"] - res.loc[binned_row, "match_length"]
#         res.loc[binned_row + 1, "freq"] /= len_bin * ncomp

#     gmean_match_length = np.sqrt(np.array(match_bin)[1:]*np.array(match_bin)[:-1])
#     res.drop([0], inplace=True)
#     res["match_length"] = gmean_match_length
#     res.reset_index(drop=True)

#     return res

def bin_mld(summed_df, linear_bin_width, limit_size, power_increment, ncomp, censor=0):
   "I swear this is the last one - it wasn't"
   linear_bin_width = int(linear_bin_width)
   limit_size = np.round(limit_size).astype(int)
   linear_ind = pd.IntervalIndex.from_arrays(
        np.arange(1, (limit_size // linear_bin_width) * linear_bin_width, linear_bin_width),
        np.arange(linear_bin_width, ((limit_size // linear_bin_width) + 1) * linear_bin_width, linear_bin_width),
        closed="both"
   )
   log_breaks = np.power(10, np.arange(np.log10(limit_size), np.log10(max(summed_df["match_length"])) + power_increment, power_increment))
   # it needs to be exact
   log_breaks = np.round(log_breaks).astype(int)
   ind = pd.IntervalIndex.from_arrays(
       log_breaks[0:-1] + 1,
       log_breaks[1:],
       closed="both"
   )
   overall_ind = linear_ind.union(ind)
   cuts = pd.cut(summed_df["match_length"], bins=overall_ind)
   binned_df = summed_df.groupby(cuts, observed=False)["freq"].sum().reset_index()
   binned_df["freq"] = binned_df.apply(lambda x: x["freq"]/((x["match_length"].right - x["match_length"].left)*ncomp), axis=1)
   binned_df["match_length"] = binned_df["match_length"].apply(lambda x: np.sqrt((x.left*x.right))).astype(float)
   censor_index = np.searchsorted(binned_df["match_length"].values, censor, "right")
   res = binned_df.iloc[censor_index:, ].reset_index(drop=True)
   return res
