
import collections
import itertools
import os
import re

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn import preprocessing

def SimpleFastaParser(handle):
    """Iterate over Fasta records as string tuples.

    Arguments:
     - handle - input stream opened in text mode

    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.

    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in SimpleFastaParser(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')

    """
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
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
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, "".join(lines).replace(" ", "").replace("\r", "")


def fasta_writer(tiselist, fasta_out):
    os.makedirs(os.path.dirname(fasta_out), exist_ok=True)
    with open(fasta_out, "w") as fileout:
        for title, seq in tiselist:
            lines = [f">{title}\n"]
            for i in range(0, len(seq), 60):
                lines.append(seq[i : i + 60] + "\n")
            fileout.write("".join(lines))


def write_only_chromosomes(refseq_dir_in, spe_file_list, refseq_dir_out):
    """Filters out plasmids from a fasta file, write a new fasta to a new location."""
    with open(spe_file_list, "r") as filein:
        species_selected = filein.read().splitlines()
    spe_list = [spe for spe in os.listdir(refseq_dir_in) if spe in species_selected]
    for spe in spe_list:
        fasta_list = [fa for fa in os.listdir(os.path.join(refseq_dir_in, spe)) if fa.endswith((".fasta", ".fa", ".fna"))]
        for fasta in fasta_list:
            tiselist = []
            with open(os.path.join(refseq_dir_in, spe, fasta), "r") as filein:
                for title, seq in SimpleFastaParser(filein):
                    if not "plasmid" in title:
                        tiselist.append((title, seq))
            fasta_writer(tiselist, os.path.join(refseq_dir_out, spe, fasta))


def get_fasta_entry_number(fasta_file, fasta_dir):
    """Get the number of > (and number of plasmid?) by fasta."""
    res = {
        "fasta": fasta_file,
        "seq_count": 0,
        "plasmid_count": 0,
        "nc": 0,
        "nz": 0,
        "ac": 0,
        "length" : 0,
        "length_plasmid" : 0,
        "g_count": 0,
        "c_count": 0
    }
    with open(os.path.join(fasta_dir, fasta_file), 'r') as filein:
        for title, seq in SimpleFastaParser(filein):
            res["seq_count"] += 1
            res["length"] += len(seq)
            res["g_count"] = seq.count("G")
            res["c_count"] = seq.count("C")
            if "plasmid" in title:
                res["plasmid_count"] += 1
                res["length_plasmid"] += len(seq)
            if title[:2] == "NZ":
                res["nz"] += 1
            if title[:2] == "NC":
                res["nc"] += 1
            if title[:2] == "AC":
                res["ac"] += 1
    res["gc_per"] = (res["g_count"] + res["c_count"]) / res["length"]
    del res["c_count"]
    del res["g_count"]
    return res


def gather_fasta_stats(species_dir):
    """returns a pd.DataFrame containing fasta stats"""
    genomes = [gen for gen in os.listdir(species_dir) if gen.endswith((".fasta", ".fa", ".fna"))]
    stat_list = [get_fasta_entry_number(fasta, species_dir) for fasta in genomes]
    return pd.DataFrame(stat_list)


def gather_species_stat(refseq_dir, file_list):
    """Concatenates_species_info"""
    with open(file_list, "r") as filein:
        species_list = filein.read().splitlines()
    species_dfs = {spe: gather_fasta_stats(os.path.join(refseq_dir, spe)) for spe in os.listdir(refseq_dir) if spe in species_list}
    return pd.concat(objs=species_dfs, axis=0, join="outer").reset_index(level=0, names=['species']).reset_index(drop=True)


def hist_count(df, column):
   plt.hist(df[column])
   plt.show()


def get_reverse_complement(sequence):
    res = list(sequence)
    for i, letter in enumerate(sequence):
        if letter == "A":
            res[i] = "T"
        elif letter == "T":
            res[i] = "A"
        elif letter == "G":
            res[i] = "C"
        elif letter == "C":
            res[i] = "G"
        else:
            res[i] = "N"
    return "".join(res)[::-1]


def get_unique_list(tetra_list):
    """Given a list of kmers, returns a list without reverse complements"""
    res = []
    for mer in tetra_list:
        if not get_reverse_complement(mer) in res:
            res.append(mer)
    return res


def deduplicate_reverse_complements(tetra_dic):
    """Given a dictionary of frequencies of tetramers, sums reverse complements."""
    unique_mers = get_unique_list(tetra_dic.keys())
    res_dic = {}
    for mer in tetra_dic.keys():
        if mer in unique_mers:
            res_dic[mer] = tetra_dic[get_reverse_complement(mer)] + tetra_dic[mer]
    return res_dic


def get_tetramers_count(sequence):
    """Calculate the tetramer frequency, keeping in mind there are only
    128 tetramers when you account for reverse complementation (4**4/2).
    Does not normalize by sequence length."""
    full_tetra_list = ["".join(seq) for seq in itertools.product("ATCG", repeat=4)]
    res = {}
    for tetra in full_tetra_list:
        start = end = None
        overlap_count = 0
        while True:
            start = sequence.find(tetra, start, end) + 1
            if start != 0:
                overlap_count += 1
            else:
                break
        res[tetra] = overlap_count
    res = deduplicate_reverse_complements(res)
    return res


def fasta_tetramers(fasta_file, species, refseq_dir):
    res = {"species" : species, "fasta" : fasta_file, "tetramers": collections.Counter()}
    with open(os.path.join(refseq_dir, species, fasta_file), 'r') as filein:
        for _, seq in SimpleFastaParser(filein):
            res["tetramers"] += collections.Counter(get_tetramers_count(seq))
    return res


def calculate_tetramers_dir(refseq_dir, species_file):
    list_tetra_spe = []
    with open(species_file, "r") as filein:
        species_admitted = filein.read().splitlines()
    species_list = [spe for spe in os.listdir(refseq_dir) if spe in species_admitted]
    for spe in species_list:
        fasta_files = [gen for gen in os.listdir(os.path.join(refseq_dir, spe)) if gen.endswith((".fasta", ".fa", ".fna"))]
        for genome in fasta_files:
            list_tetra_spe.append(fasta_tetramers(genome, spe, refseq_dir))
    res_df = pd.DataFrame(list_tetra_spe)
    # unpack tetramers dic to columns
    return pd.concat([res_df.drop(["tetramers"], axis=1), res_df["tetramers"].apply(pd.Series)], axis=1)


def single_tetramers_norm(tetramer, tetramer_count, gc_per, seq_len):
    # normalizes a tetramer by a gc_per
    gc_count = tetramer.count("G") + tetramer.count("C")
    res = tetramer_count / ( (gc_per ** gc_count) * ((1 - gc_per) ** (4 - gc_count)) * seq_len )
    return res


def normalize_nucleotides_frequencies(tetramers_row):
    # normalizes frequencies, dividing the counts of a kmer by its expected frequency
    # according to a zero-order Markov chain, keeping G and C frequencies undistinct
    length = tetramers_row["length"]
    gc_per = tetramers_row["gc_per"]
    for column, value in tetramers_row.items():
        if re.match(r"[ATGC][ATGC][ATGC][ATGC]", column):
            tetramers_row[column] = single_tetramers_norm(column, value, gc_per, length)
        else:
            continue
    return tetramers_row



if __name__ == "__main__":
    recalculate = True
    # GC correction does not work
    normalize_frequencies_gc = False
    stats_outfile = "tetra_stats_df_all_recode.csv"
    raw_tetra_outfile = "raw_tetramers.csv"
    genome_dir = "../genomes/refseq_bacillaceae_filtered_all"
    species_list_file = "./species_list_bacillaceae_5.txt"
    stats_df = gather_species_stat(genome_dir, species_list_file)

    # gather data
    if recalculate:
        tetra_df = calculate_tetramers_dir(genome_dir, species_list_file)
        tetra_df.to_csv(raw_tetra_outfile, index=False)
    else:
        tetra_df = pd.read_csv(raw_tetra_outfile)

    if normalize_frequencies_gc:
        tetra_df_norm = pd.merge(stats_df[["fasta", "length", "gc_per"]], tetra_df, on=["fasta"], how="inner", validate="one_to_one").apply(normalize_nucleotides_frequencies, axis=1)
        tetra_df_norm = tetra_df_norm.drop(["length", "gc_per"], axis=1)
    else:
        # normalise length
        tetra_df_norm = pd.concat([tetra_df[["species", "fasta"]].reset_index(drop=True), tetra_df.drop(["species", "fasta"], axis=1).div(stats_df["length"], axis=0).reset_index(drop=True)], axis=1)
    # merge data
    total_df = pd.merge(stats_df, tetra_df_norm, on=["fasta", "species"], how="inner", validate="one_to_one")
    species_map = {
        'Lactiplantibacillusall': "Lactiplantibacillus",
        'Staphylococcusall': "Staphylococcus",
        'Cytobacillus': 'Cytobacillus',
        'Metabacillus' : 'Metabacillus',
        'Parageobacillus': 'Parageobacillus',
        'Lysinibacillus': 'Lysinibacillus',
        'BacillusAll': "Bacillus",
        'Peribacillus': 'Peribacillus' ,
        'Halobacillus': 'Halobacillus',
        'Neobacillus': 'Neobacillus',
        'Virgibacillus': 'Virgibacillus',
        'Heyndrickxia': 'Heyndrickxia',
        'EscherichiaColiall': 'EscherichiaColi',
        'Geobacillus': 'Geobacillus',
        'Priestia': 'Priestia',
        'Rossellomorea': 'Rossellomorea',
        'Anoxybacillus' : 'Anoxybacillus',
        'Niallia': 'Niallia'
    }
    total_df = total_df.assign(species = total_df.species.map(species_map))
    total_df.to_csv(stats_outfile, index=False)
        # calculate PCA
    X = total_df.drop(["species", "fasta", "seq_count", "nc", "plasmid_count", "nz", "ac", "length", "gc_per"], axis=1)
    pca = PCA(n_components=2)
    pca.fit(X)
    print(f"Variance explained by the two componenents : {pca.explained_variance_ratio_}")
    pca_res = pca.transform(X)
    le = preprocessing.LabelEncoder()
    le.fit(total_df["species"])
    y = le.transform(total_df["species"])

    # Plot PCA
    fig, ax = plt.subplots(figsize=(10, 10))
    scatter = ax.scatter(pca_res[:, 0], pca_res[:, 1], c=y, cmap="tab10")
    ax.set(xlabel=f"{pca.explained_variance_ratio_[0]*100:.2f} %", ylabel=f"{pca.explained_variance_ratio_[1]*100:.2f} %")
    _ = ax.legend(
        scatter.legend_elements()[0], le.classes_, loc="upper right", title="Bacteria"
    )
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    fig.savefig("PCA_subsamp_wo_plasmids_python_norm.png", dpi=96)
