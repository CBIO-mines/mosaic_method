#!/usr/bin/env python3
import argparse
import sys

import numpy as np
import pandas as pd
from fun import *

def main():
    parser = argparse.ArgumentParser(
        description="""
        Concatenates single mlds from alignements, saves them. Optionnaly.
        """)
    parser.add_argument(
        "bin_out_file",
        type=str,
        help="binned mld output file path"
    )
    parser.add_argument(
        "--from_csv",
        type=str,
        help="Start from csv files. Incompatible with from_florian_mld"
    )
    parser.add_argument(
        "--from_florian_mld",
        help="Start from assembly-wise MLD files. Incompatible with --from_csv or --from_full_mld",
        type=str
    )
    parser.add_argument(
        "--species",
        help="Comma-separated species whose mlds must be summed and merged",
        type=str
    )
    parser.add_argument(
        "--species_csv",
        help="File defining genome-species relationships",
        type=str
    )
    parser.add_argument(
        "--save_full_mld",
        type=str,
        help="full mld by comparison output file path"
    )
    args = parser.parse_args()

    if args.from_csv and args.from_florian_mld:
        sys.exit("--from_csv and --from_florian_mld are incompatible options")

    if args.from_csv:
        species = args.species.split(",")
        genome_comps = get_genome_comp(species, args.species_csv, args.from_csv)
        full_mld = parse_csv(genome_comps)
        summed_mld = sum_mlds(full_mld)
        binned_mld = bin_mld(
            summed_df=summed_mld,
            linear_bin_width=3,
            limit_size=35.5,
            power_increment=0.1,
            ncomp=full_mld.shape[0]
        )
    elif args.from_florian_mld:
        full_mld = parse_florian_mld(args.from_florian_mld)
        summed_mld = sum_mlds(full_mld)
        binned_mld = bin_mld(
            summed_df=summed_mld,
            linear_bin_width=3,
            limit_size=35.5,
            power_increment=0.1,
            ncomp=full_mld.shape[0]
        )

    full_mld.to_csv(args.save_full_mld)
    binned_mld.to_csv(args.bin_out_file, index=False)


if __name__ == "__main__":
    main()
