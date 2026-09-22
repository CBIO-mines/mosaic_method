#!/usr/bin/env python3

import seaborn as sns
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform


def plot_heatmap(res_df, outfile="distance_matrix.png"):
    species = sorted(set(res_df['species_1']) | set(res_df['species_2']))

    pivot = (res_df
        .pivot(index='species_1', columns='species_2', values='log10tau')
        .reindex(index=species, columns=species)
    )
    matrix_tau = pivot.fillna(pivot.T)
    np.fill_diagonal(matrix_tau.values, 0)
    linkage_matrix = linkage(squareform(matrix_tau.values), method='average')
    n = len(matrix_tau)
    g = sns.clustermap(
        10**matrix_tau,
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix,
        cmap='viridis',
        figsize=(n * 0.3, n * 0.3),  # scale with number of species

    )
    g.savefig(outfile, dpi=300)
