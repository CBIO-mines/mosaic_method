#!/usr/bin/env python3

import pandas as pd
import numpy as np


def count_comp(taxon_path, level="species.gtdb"):
    taxon_df = pd.read_csv(taxon_path)
    counts = taxon_df.groupby(level)["genome"].count().values
    n_comps = (counts.sum() ** 2 - (counts**2).sum()) / 2
    average_genome_count = counts.mean()
    return average_genome_count, n_comps
