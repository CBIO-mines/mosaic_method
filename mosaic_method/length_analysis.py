"""Small module to gather and represent length of fasta files"""

import os

import matplotlib.pyplot as plt
import pandas as pd


def SimpleFastaParser(handle):
    """Iterate over Fasta records as string

    Arguments:
     - handle - input stream opened in text mode

    For each record a string is returned, the sequence (with any
    whitespace removed).
    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in SimpleFastaParser(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('delta', 'CGCGC')

    Adapted from biopython https://github.com/biopython/biopython (no need for the full biopython project)

    """
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if line[0] == ">":
            break
    else:
        # no break encountered - probably an empty file
        return

    # Main logic
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines = []
    for line in handle:
        if line[0] == ">":
            yield "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            continue
        lines.append(line.rstrip())

    yield "".join(lines).replace(" ", "").replace("\r", "")


def get_fasta_len(fasta_file):
    """Get the sum of the lengths of the contig in a fasta file."""
    res = 0
    with open(fasta_file, "r") as fasta_file:
        for seq in SimpleFastaParser(fasta_file):
            res += len(seq)
    return res


def get_len_distribution(fasta_dir, species_df):
    """For a directory containing all genomes and a dataframe of species to genome mapping, gathers all the lengths of the genomes."""
    species_df["length"] = species_df.apply(
        lambda row: get_fasta_len(os.path.join(fasta_dir, row.genome)), axis=1
    )
    return species_df


def plot_histogram(df_len, output_file=None):
    """Represent the lengths distribution as a histogram."""
    fig = plt.Figure()
    ax = fig.subplots()
    ax.hist(df_len["length"], bins=30)
    ax.set_xlabel("Fasta length")
    ax.set_ylabel("Count")
    if output_file:
        fig.savefig(output_file)
    else:
        fig.show()


def length_analysis(taxon_csv, cluster, fasta_dir, save_dir):
    taxon_df = pd.read_csv(taxon_csv, index_col=0)
    len_df = get_len_distribution(fasta_dir, taxon_df)
    grouped_l = len_df.groupby(cluster)
    for level in grouped_l.groups.keys():
        output_csv = os.path.join(save_dir, str(level) + "_distribution.csv")
        output_png = os.path.join(save_dir, str(level) + "_distribution.png")
        grouped_l.get_group(level).drop(cluster, axis=1).to_csv(output_csv, index=False)
        plot_histogram(grouped_l.get_group(level), output_file=output_png)
