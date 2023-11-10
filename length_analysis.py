"""Small module to gather and represent length of fasta files"""

import argparse
import os

import matplotlib . pyplot as plt
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
    with open(fasta_file, 'r') as fasta_file:
        for seq in SimpleFastaParser(fasta_file):
            res += len(seq)
    return res


def get_len_distribution(fasta_dir):
    """For a directory containing fasta, gathers all the lengths of the files."""
    fasta_files = [fafile for fafile in os.listdir(fasta_dir) if fafile.endswith(".fa") or fafile.endswith(".fasta")]
    res_df = pd.DataFrame.from_dict(data={"Genome" : fasta_files})
    res_df["Length"] = res_df.apply(lambda row: get_fasta_len(os.path.join(fasta_dir, row.Genome)), axis=1)
    return res_df


def plot_histogram(df_len, output_file=None):
    """Represent the lengths distribution as a histogram."""
    plt.hist(df_len["Length"], bins = 10)
    plt.xlabel("Fasta length")
    plt.ylabel("Count")
    if output_file:
        plt.savefig(output_file)
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="""
        From a given directory containing fasta files,
        plots a histogram of the lengths of the sum of the contigs,
        optionnaly writes a csv file of those lengths."""
    )
    parser.add_argument(
        "--save_distr",
        type=str,
        help="The file to save a csv of the distribution"
    )
    parser.add_argument(
        "--save_plot",
        type=str,
        help="The file to save the histogram"
    )
    parser.add_argument(
        "fasta_dir",
        type=str,
        help="The directory containing fasta files"
    )
    args = parser.parse_args()

    len_df = get_len_distribution(args.fasta_dir)
    if args.save_distr:
        len_df.to_csv(args.save_distr, index=False)
    if args.save_plot:
        plot_histogram(len_df, output_file=args.save_plot)
    else:
        output_file = args.fasta_dir.strip("/") + "_distribution.png"
        plot_histogram(len_df, output_file)


if __name__ == "__main__":
    main()
