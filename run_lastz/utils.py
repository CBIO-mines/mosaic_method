#!/usr/bin/env python3


import collections
import os
import subprocess as sp

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


def save_lastz(query, target, output_csv):
    """
    Takes the exact matches dictionnary, saves it as a csv file for further processing.
    """
    dict_lastz = run_lastz(query, target)
    res = pd.DataFrame.from_dict(dict_lastz, orient="index").rename_axis("match_length").reset_index()
    res.rename(columns={0: "freq"}, inplace=True)
    res = res.sort_values(by="match_length")
    res.to_csv(output_csv, index=False)


def run_lastz(query, target):
    """
    Runs lastz, adds cigarx matches length counts together and returns them as a counter.
    """
    query = query+"[multiple]"
    target = target+"[multiple]"

    res_lastz = sp.run(
        ['lastz',query,target,"--format=general:cigarx", "--ambiguous=iupac", "--chain"],
        capture_output=True,
        check=True,
        encoding="utf-8"
    )
    exact_matches = collections.Counter()
    for row in res_lastz.stdout.splitlines()[1:]:
        exact_matches += parse_cigarx_line(row)
    return exact_matches
