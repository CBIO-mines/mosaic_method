#!/usr/bin/env python3

import os

import numpy as np
import pandas as pd


def count_aligns(spec_df):
    spec_list = np.array(spec_df["Genome_count"])
    res = 0
    for i in range(len(spec_list)):
        for j in range(i + 1, len(spec_list)):
            res += spec_list[i] * spec_list[j]
    return res




species_file = "/cluster/CBIO/home/petheimer/test_florian/species_list_bacillaceae_all.txt"
species_dir = "/cluster/CBIO/data1/petheimer/refseq_bacillaceae/"

with open(species_file, "r") as filein:
    species_list = filein.read().splitlines()

res = []
for root, _, files in os.walk(species_dir):
    spec = root.split("/")[-1]
    fasta_files = [f for f in files if f.endswith((".fna", ".fasta", ".fa")) and not f.startswith(".")]
    if spec in species_list:
        res.append({"Species": spec, "Genome_count" : len(fasta_files)})

spec_num_df = pd.DataFrame(res)
print(spec_num_df)

print(f"The number of alignments necessary for this list is : {count_aligns(spec_num_df)}")
