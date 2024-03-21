#!/usr/bin/env python3

import argparse
import sys

import numpy as np
import pandas as pd
from fit import *

def main():
    parser = argparse.ArgumentParser(
        description="""
        Fits the parameters to the match length distribution according to Sheinman et al. 2023""")
    parser.add_argument(
        "binned_mld",
        type=str,
        help="The binned mld from the parsing script"
    )
    parser.add_argument(
        "out_params",
        type=str,
        help="fitted params output file path"
    )
    parser.add_argument(
        "--bacs",
        type=str,
        help="Names of the two taxon/bacteria of the comparison, separated by a comma"
    )
    parser.add_argument(
        "--L0",
        type=str,
        help="If given, the L0 used in calculating the theoretical_mld. Otherwise, it is fitted"
    )
    parser.add_argument(
        "--save_surface_plot",
        type=str,
        help="Whether to draw and where to save a surface plot of the fitting of parameters"
    )
    parser.add_argument(
        "--mus",
        type=float,
        help="The least conserved segments' mutation rate",
        default=5.e-9
    )
    parser.add_argument(
        "--muc",
        type=float,
        help="The most conserved segments' mutation rate",
        default=6.e-11
    )
    parser.add_argument(
        "--delta",
        type=float,
        help="The aligner's sensitivity",
        default=0.55
    )

    args = parser.parse_args()
    binned_mld = pd.read_csv(args.binned_mld)
    mus = args.mus
    muc = args.muc
    delta = args.delta
    if args.L0:
        L0_df = pd.read_csv(args.L0, index_col=False)
        species = args.bacs.split(",")
        L0 = L0_df[ (L0_df["bac1"] == species[0]) & (L0_df["bac2"] == species[1]) ]["L0"].squeeze()
        init_params = np.array([5, -5])
    else:
        L0 = None
        init_params = np.array([5, -5, 10e6])


    res_opt_full, res_opt_minus3 = fit_params(
        opt_method="dual-annealing",
        init_pars=init_params,
        empirical_mld=np.array(binned_mld["freq"], dtype=np.float64),
        smal_dif=np.double(0.1),
        match_lengths=np.array(binned_mld["match_length"], dtype=np.float64),
        mus=mus,
        muc=muc,
        delta=delta,
        L0=np.double(L0)
    )
    if not res_opt_full.success:
        sys.exit("Fitting failed")
    else:
        write_results(res_opt_full, args.out_params, L0, res_opt_minus3)

    if args.save_surface_plot:
        plot_surface(
            res_opt_full.x[0] - 1,
            res_opt_full.x[0] + 1,
            res_opt_full.x[1] - 1,
            res_opt_full.x[1] + 1,
            100,
            args.save_surface_plot,
            np.array(binned_mld["freq"]),
            0.1,
            np.array(binned_mld["match_length"]),
            mus,
            muc,
            delta,
            L0,
            res_opt_full.x
        )


if __name__ == "__main__":
    main()
