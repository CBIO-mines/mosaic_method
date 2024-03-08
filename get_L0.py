"""Small module that calculates the L0 from two distributions of lengths of sums of contigs"""

import argparse
import itertools
import os

import pandas as pd


def L0_calculation(bac1_csv, bac2_csv, full_min=True):
    bac1_df = pd.read_csv(bac1_csv, index_col=False)
    bac2_df = pd.read_csv(bac2_csv, index_col=False)
    res = 0
    if full_min:
        for bac1_l, bac2_l in itertools.product(bac1_df["length"], bac2_df["length"]):
            res += min(bac1_l, bac2_l)
        res /= bac1_df.shape[0] * bac2_df.shape[0]
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


def main():
    parser = argparse.ArgumentParser(
        description="""
        From the combination of all possible pairs of species, from the sum of contigs lengths distributions (from each assembly for a species)
        stored in csv files (one for each species), get all L0s."""
    )
    parser.add_argument(
        "distr_dir",
        type=str,
        help="The directory containing the sums of contigs distributions"
    )
    parser.add_argument(
        "L0_csv",
        type=str,
        help="The path to where the L0 csv must be stored"
    )
    args = parser.parse_args()
    L0_df = L0_concatenation(args.distr_dir)
    write_L0_df(L0_df, args.L0_csv)


if __name__ == "__main__":
    main()
