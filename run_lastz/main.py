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
        "query",
        type=str,
        help="The query full path"
    )
    parser.add_argument(
        "target",
        type=str,
        help="The target full path"
    )
    parser.add_argument(
        "output_csv",
        type=str,
        help="The output csv file"
    )
    args = parser.parse_args()
    save_lastz(args.query, args.target, args.output_csv)

if __name__ == "__main__":
    main()
