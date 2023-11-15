#!/usr/bin/env python3

import argparse
import sys

import numpy as np
import pandas as pd
from parsefit import parse, fit


def main():
    parser = argparse.ArgumentParser(
        description="""
        Concatenates cigarx stats from alignements,
        then fits the parameters to the match length distribution according to Sheinman et al. 2023""")
    parser.add_argument(
        "--from_cigarx",
        type=str,
        help="Start from cigarx files. Incompatible with --from_full_mld"
    )
    parser.add_argument(
        "--from_full_mld",
        help="Start from a binned mld csv file. Incompatible with --from_cigarx",
        type=str
    )
    parser.add_argument(
        "out_params",
        type=str,
        help="fitted params output file path"
    )
    parser.add_argument(
        "bin_out_file",
        type=str,
        help="binned mld output file path"
    )
    parser.add_argument(
        "--L0",
        type=str,
        help="If given, the L0 used in calculating the theoretical_mld. Otherwise, it is fitted"
    )
    parser.add_argument(
        "--save_full_mld",
        type=str,
        help="if specified, full mld by comparison output file path. Incompatible with --from_full_mld"
    )

    args = parser.parse_args()
    if args.from_cigarx and args.from_full_mld:
        sys.exit("--from_cigarx and --from_full_mld are incompatible options")
    if args.from_full_mld and args.save_full_mld:
        sys.exit("--from_full_mld and --save_full_mld are incompatible options")


    mus = 5e-9
    muc = 6e-11
    delta = 0.55
    if args.L0:
        L0_df = pd.read_csv(args.L0, index_col=False)
        # Search relevant L0 based on output file species (a bit hacky)
        species = args.bin_out_file.split("/")[1].split("_")[:2]
        L0 = L0_df[ (L0_df["bac1"] == species[0]) & (L0_df["bac2"] == species[1]) ]["L0"].squeeze()
        init_params = np.array([5, -5])
    else:
        L0 = None
        init_params = np.array([5, -5, 10e6])

    if not args.from_full_mld:
        full_mld = parse.parse_cigars(args.from_cigarx, args.save_full_mld)
        summed_mld = fit.sum_mlds(full_mld)
        binned_mld = fit.bin_mld(
            summed_df=summed_mld,
            linear_bin_width=3,
            limit_size=35.5,
            power_increment=0.1,
            ncomp=full_mld.shape[0]
        )
    else:
        full_mld = pd.read_csv(args.from_full_mld, index_col=0)
        summed_mld = fit.sum_mlds(full_mld)
        binned_mld = fit.bin_mld(
            summed_df=summed_mld,
            linear_bin_width=3,
            limit_size=35.5,
            power_increment=0.1,
            ncomp=full_mld.shape[0]
        )

    res_opt = fit.fit_params(
        opt_method="nelder-mead",
        init_pars=init_params,
        empirical_mld=np.array(binned_mld["freq"]),
        smal_dif=0.1,
        match_lengths=np.array(binned_mld["match_length"]),
        mus=mus,
        muc=muc,
        delta=delta,
        L0=L0
    )
    if not res_opt.success:
        sys.exit("Fitting failed")
    else:
        fit.write_results(binned_mld, res_opt.x, args.bin_out_file, args.out_params, L0)


if __name__ == "__main__":
    main()
