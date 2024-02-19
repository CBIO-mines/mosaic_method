#!/usr/bin/env python3

import os

import numpy as np
import pandas as pd


def count_aligns(spec_df):
    spec_list = np.array(spec_df["genome"])
    res = 0
    for i in range(len(spec_list)):
        for j in range(i + 1, len(spec_list)):
            res += spec_list[i] * spec_list[j]
    return res



species_file = "species_list_bacillaceae_florian_genomes.txt"
species_csv = "species_florian.csv"

with open(species_file, "r") as filein:
    species_list = filein.read().splitlines()

species_df = pd.read_csv(species_csv)
species_counts = species_df.groupby("species").count()
print(species_counts)
print(f"The number of alignments necessary for this list is : {count_aligns(species_counts)}")
