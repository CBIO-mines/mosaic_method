#!/usr/bin/env python3

import argparse
import collections
import csv
import os

import pandas as pd

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


def parse_cigars(path, outfile=None):
    """
    Parses a directory of lastz output files and writes the resulting length counts.
    """
    lastz_files = [lz_f for lz_f in os.listdir(path) if ".txt" in lz_f]
    matches_dic = {}
    for lz_f in lastz_files:
        matches_dic[lz_f[:-4]] = read_lastz_output(os.path.join(path, lz_f))

    df_mlds = pd.DataFrame.from_dict(matches_dic, orient="index").rename_axis("comp").reset_index()
    if outfile:
        df_mlds.to_csv(outfile)
    return df_mlds
