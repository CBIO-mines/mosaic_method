#!/usr/bin/env python3

import os

import numpy as np
import pandas as pd


def count_aligns(genome_count):
    """
    Computes the number of alignments necessary : each genus's genome is aligned with each other genus's genome
    """
    res = 0
    for i in range(len(genome_count)):
        for j in range(i + 1, len(genome_count)):
                res += genome_count[i] * genome_count[j]
    return res




cluster_level = "misha_annotation"
taxon_csv = "/project/bacteria_mlds-data/misha_taxon.csv"

taxon_df = pd.read_csv(taxon_csv)
cluster_counts = taxon_df.groupby(cluster_level).count()["genome"]
print(cluster_counts)
print(f"The number of alignments necessary for this list is : {count_aligns(cluster_counts)}")
