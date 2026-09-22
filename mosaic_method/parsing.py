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


def sum_mlds(genome_comps, lastz_db_path):
    """
    Parses a database and returns the resulting mlds concatenated.
    A version without parallel reads of the db but parallel decoding
    """
    sqlite_con = sqlite3.connect(lastz_db_path)
    cur = sqlite_con.cursor()
    genome_comps = [sorted(pair) for pair in genome_comps]

    cur.execute("CREATE TEMP TABLE temp_pairs (genome1 TEXT, genome2 TEXT);")
    cur.executemany("INSERT INTO temp_pairs(genome1, genome2) VALUES (?, ?)", genome_comps)
    cur.execute("SELECT l.genome1, l.genome2, l.count_array FROM lastz l JOIN temp_pairs t ON l.genome1 = t.genome1 AND l.genome2 = t.genome2;")
    rows = cur.fetchall()
    cur.execute("DROP TABLE temp_pairs;")

    max_len = max(len(np.frombuffer(l, dtype=np.int64)) for _, _, l in rows)
    summed_mld = np.zeros(max_len)
    for _, _, byte_array in rows:
       mld_array = np.frombuffer(byte_array, dtype=np.int64)
       padded_array = np.pad(mld_array, (0, max_len - len(mld_array)), constant_values=0)
       summed_mld += padded_array

    summed_df = pd.DataFrame(summed_mld,  index=range(1, max_len + 1), columns=["freq"], dtype=np.int64).reset_index(names=["match_length"])
    summed_df["match_length"] = summed_df["match_length"].astype(np.int64)
    return summed_df


def sum_mlds_old(mld_comp_df):
    """
    Takes a df of mlds by comparison and sums it.
    """
    summed_colmuns = mld_comp_df.drop(labels="comp", axis=1).sum(axis=0, skipna=True)
    summed_df = summed_colmuns.reset_index()
    summed_df.columns = ["match_length", "freq"]
    summed_df = summed_df.astype({"match_length" : "int64"})
    return summed_df.sort_values(by=["match_length"]).reset_index(drop=True)

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
