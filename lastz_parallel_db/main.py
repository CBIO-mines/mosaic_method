#!/usr/bin/env python3

import argparse
import sys

import numpy as np
import pandas as pd
from utils import *

def main():
    parser = argparse.ArgumentParser(
        description="""
        Runs lastz, saves its output as smaller csv files.
        """)
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=10,
        help="Number of threads to use (default: 10)"
    )
    parser.add_argument(
        "-u",
        "--update",
        action="store_true",
        help="Update the existing database instead of creating a new one"
    )
    parser.add_argument(
        "taxon_csv",
        type=str,
        help="Path to the taxon csv file"
    )
    parser.add_argument(
        "genomes_path",
        type=str,
        help="Path to the genomes directory"
    )
    parser.add_argument(
        "cluster_name",
        type=str,
        help="The name of the cluster level(column in taxon_csv) to use to make comparisons"
    )
    parser.add_argument(
        "output_db",
        type=str,
        help="The output sqlite file"
    )
    args = parser.parse_args()
    if args.update:
        update_lastz_db(args.taxon_csv, args.cluster_name, args.output_db, args.threads)
    con = create_lastz_db(args.taxon_csv, args.genomes_path, args.cluster_name, args.output_db, args.threads)
    con.close()

if __name__ == "__main__":
    main()
