#!/usr/bin/env python3

import collections
import os

import pandas as pd
import numpy as np

def read_lastz_file(lastz_file):
    """
    Reads a line from a lastz file alignement file output.
    """
    with open(lastz_file, "r") as filein:
        for line in filein:
            yield line.strip("\n")


def parse_cigarx_line(line):
    """
    Parses a cigarx line and counts the length of matches.
    """
    li = len(line) - 1
    matches = []
    while li >= 0:
        if line[li] == "=":
            len_match = ""
            li -= 1
            while line[li].isdigit():
                len_match = line[li] + len_match
                li -= 1
                if li == -1:
                    break
            if len_match:
                matches.append(int(len_match))
            else:
                matches.append(1)
        else:
            li -= 1
    res = collections.Counter(matches)
    return res


def read_lastz_output(lastz_file):
    """
    Adds all the cigarx matches length counts together.
    """
    exact_matches = collections.Counter()
    for line in read_lastz_file(lastz_file):
        exact_matches += parse_cigarx_line(line)
    return exact_matches


def parse_cigars(path):
    """
    Parses a directory of lastz output files and returns the resulting length counts.
    """
    lastz_files = [lz_f for lz_f in os.listdir(path) if ".txt" in lz_f]
    matches_dic = {}
    for lz_f in lastz_files:
        matches_dic[lz_f[:-4]] = read_lastz_output(os.path.join(path, lz_f))

    df_mlds = pd.DataFrame.from_dict(matches_dic, orient="index").rename_axis("comp").reset_index()
    return df_mlds


def parse_florian_mld(path):
    """
    Parses Florian's assembly-wise MLD files and merges them into a taxon-wide full_mld.

    Parameters
    ----------
    path: str
    the path to the directory (named bac1_bac2) where the MLDs are stored

    Returns
    -------
    A pandas.DataFrame containing a taxa comparison's MLD
    """
    comp_files = [mld_file for mld_file in os.listdir(path) if mld_file.endswith(".MLD")]
    res = {}
    for c_f in comp_files:
        tmp_dic = {}
        c_f_full = os.path.join(path, c_f)
        with open(c_f_full, "r") as filein:
            for line in filein:
                try:
                    tmp_dic[line.split()[0]] = int(line.split()[1])
                except IndexError:
                    print(line)
                except:
                    print("other error")
                finally:
                    continue
        res[c_f.split(".")[0].replace("-", "_")] = tmp_dic
    res_df = pd.DataFrame.from_dict(res, orient="index")
    res_df = res_df.rename_axis("comp").reset_index()
    return res_df


def sum_mlds(mld_comp_df):
    """
    Takes a df of mlds by comparison and sums it.
    """
    summed_colmuns = mld_comp_df.drop(labels="comp", axis=1).sum(axis=0, skipna=True)
    summed_df = summed_colmuns.reset_index()
    summed_df.columns = ["match_length", "freq"]
    summed_df = summed_df.astype({"match_length" : "int64"})
    return summed_df.sort_values(by=["match_length"]).reset_index(drop=True)


def bin_mld(summed_df, linear_bin_width, limit_size, power_increment, ncomp):
    """
    Bins and normalizes a summed mld according to a specific pattern.

    Parameters
    ----------
    summed_df: pd.DataFrame
    an unbinned dataframe with a column "match_length" and a column "freq"
    linear_bin_width: float
    the width of the bin in the linear part of bin vector
    limit_size: float
    the limit at which the bin vector switches from linear to log
    power_increment: float
    the "bin width" of the log part
    ncomp: float
    the number of comparison summed here

    Returns
    -------
    a binned pd.Dataframe with columns :
    - "match_length" containing the geometric mean of the bin
    - "freq" containing the counts corresponding to the bin
    """
    match_bin = list(np.arange(0.5, limit_size, linear_bin_width))
    initial_len = len(match_bin)
    cur_power = 0.1
    while match_bin[-1] < max(summed_df["match_length"]):
        match_bin += [match_bin[initial_len - 1]*10**(cur_power)]
        cur_power += power_increment

    res = pd.DataFrame.from_dict(
        {"match_length" : match_bin,
         "freq" : [0.0] * len(match_bin)}
    )
    summed_row = 0
    binned_row = 0
    while(summed_row < summed_df.shape[0]):
        if(summed_df.loc[summed_row, "match_length"] <= res.loc[binned_row, "match_length"]):
            res.loc[binned_row, "freq"] += summed_df.loc[summed_row, "freq"]
            summed_row += 1
        else:
            binned_row += 1
            if binned_row > res.shape[0] - 1:
                break

    for binned_row in range(res.shape[0] - 1):
        len_bin = res.loc[binned_row + 1, "match_length"] - res.loc[binned_row, "match_length"]
        res.loc[binned_row + 1, "freq"] /= len_bin * ncomp

    gmean_match_length = np.sqrt(np.array(match_bin)[1:]*np.array(match_bin)[:-1])
    res.drop([0], inplace=True)
    res["match_length"] = gmean_match_length
    res.reset_index(drop=True)

    return res
