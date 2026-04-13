"""Small module that calculates the L0 from two distributions of lengths of sums of contigs"""

import itertools
import os

import numpy as np
import pandas as pd


def L0_calculation(bac1_csv, bac2_csv, full_min=True):
    bac1_df = pd.read_csv(bac1_csv, index_col=False)
    bac2_df = pd.read_csv(bac2_csv, index_col=False)
    res = 0
    if full_min:
        res = np.minimum.outer(bac1_df["length"].values, bac2_df["length"].values).mean()
    else :
        res = (bac1_df["length"].min() + bac2_df["length"].min()) / 2
    return res


def L0_concatenation(distr_dir):
    distr_files = [dfile for dfile in os.listdir(distr_dir) if dfile.endswith("_distribution.csv")]
    # sorting not taking _distribution.csv into account
    # thus Bacillus < Bacillus_A (whereas Bacillus_A_distribution.csv < Bacillus_distribution.csv)
    distr_files.sort(key=lambda x : x.rpartition("_")[0])
    bac1_2 = [(L0_calculation(os.path.join(distr_dir, bac_1), os.path.join(distr_dir, bac_2)), bac_1.rpartition("_")[0], bac_2.rpartition("_")[0]) for (bac_1, bac_2) in itertools.combinations(distr_files, 2)]
    res = pd.DataFrame(bac1_2, columns = ["L0", "bac1", "bac2"])
    return res


def write_L0_df(L0_df, outfile):
    L0_df.to_csv(outfile, index=False)


def L0_calc(distr_dir, L0_csv):
    L0_df = L0_concatenation(distr_dir)
    write_L0_df(L0_df, L0_csv)
