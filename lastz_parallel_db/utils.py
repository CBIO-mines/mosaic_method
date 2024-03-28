import collections
import concurrent.futures
import io
import itertools
import os
import subprocess as sp

import sqlite3
import pandas as pd
import numpy as np


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
    count_array = np.zeros(max(matches), dtype=int)
    for match_length in matches:
        # first element is matches of length 1
        count_array[match_length - 1] += 1
    return count_array


def run_lastz(query, target):
    """
    Runs lastz, adds cigarx matches length counts together and returns them as a counter.
    """
    query = query+"[multiple]"
    target = target+"[multiple]"

    res_lastz = sp.run(
        ['lastz',query,target,"--format=general:length1,idfrac,cigarx", "--ambiguous=iupac", "--chain"],
        capture_output=True,
        check=True,
        encoding="utf-8"
    )
    res_lastz_df = pd.read_csv(io.StringIO(res_lastz.stdout), sep="\t")
    summed_count_array = None
    summed_matches = 0
    summed_aligned = 0
    for _, row in res_lastz_df.iterrows():
        # cigarx counting
        count_array = parse_cigarx_line(row["cigarx"])
        # Pad the arrays with zeroes if they have different sizes
        if summed_count_array is None:
            summed_count_array = count_array
            continue
        if len(count_array) < len(summed_count_array):
            count_array = np.pad(count_array, (0, len(summed_count_array) - len(count_array)), mode='constant')
        elif len(count_array) > len(summed_count_array):
            summed_count_array = np.pad(summed_count_array, (0, len(count_array) - len(summed_count_array)), mode='constant')
        summed_count_array += count_array

        # idfrac calculations
        summed_matches += int(row["idfrac"].split("/")[0])
        summed_aligned += int(row["idfrac"].split("/")[1])
    average_divergence = 1 - summed_matches / summed_aligned
    return summed_count_array, average_divergence

def lastz_entry(genome_pair):
    """
    computes a lastz entry for two genomes
    """
    genome_1, genome_2 = genome_pair
    count_array, average_divergence = run_lastz(genome_1, genome_2)
    genome_small = min(genome_1, genome_2)
    genome_large = max(genome_1, genome_2)
    return genome_small, genome_large, count_array, average_divergence


def create_lastz_db(taxon_csv, genomes_path, cluster_name, db_name, threads):
    """Creates a sqlite3 database from a taxon csv file, generating all necessary alignments"""
    sqlite3_conn = sqlite3.connect(db_name)
    taxon_df = pd.read_csv(taxon_csv)
    sqlite3_conn.execute("CREATE TABLE lastz (genome1 STRING, genome2 STRING, count_array blob, average_divergence INT);")
    sqlite3_conn.execute("CREATE TABLE taxon (genome STRING, cluster STRING);")
    for _, row in taxon_df.iterrows():
        sqlite3_conn.execute("INSERT INTO taxon VALUES (?, ?)", (row["genome"], row[cluster_name]))
    sorted_clusters = sorted(taxon_df[cluster_name].unique())

    # gather necessary genome comparisons
    genomes_comps = []
    for cluster_1, cluster_2 in itertools.combinations(sorted_clusters, 2):
        genomes_1 = taxon_df[taxon_df[cluster_name] == cluster_1]["genome"]
        genomes_2 = taxon_df[taxon_df[cluster_name] == cluster_2]["genome"]
        genomes_1 = [os.path.join(genomes_path, genome) for genome in sorted(genomes_1)]
        genomes_2 = [os.path.join(genomes_path, genome) for genome in sorted(genomes_2)]
        genomes_comps += list(itertools.product(genomes_1, genomes_2))

    # run lastz in parallel
    res_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for result in executor.map(lastz_entry, genomes_comps):
            res_list.append(result)
    for res in res_list:
        res = list(res)
        res[2] = res[2].tobytes()
        res[0] = os.path.basename(res[0])
        res[1] = os.path.basename(res[1])
        sqlite3_conn.execute("INSERT INTO lastz VALUES (?, ?, ?, ?)", res)
    sqlite3_conn.commit()
    return sqlite3_conn


def update_lastz_db(taxon_csv, genomes_path, cluster_name, db_name, threads):
    """Updates a sqlite3 database from a taxon csv file, generating all necessary alignments"""
    sqlite3_conn = sqlite3.connect(db_name)
    taxon_df = pd.read_csv(taxon_csv)
    sorted_clusters = sorted(taxon_df[cluster_name].unique())
    already_compared = sqlite3_conn.execute("SELECT genome1, genome2 FROM lastz").fetchall()
    # gather necessary genome comparisons
    genomes_comps = []
    for cluster_1, cluster_2 in itertools.combinations(sorted_clusters, 2):
        genomes_1 = sorted(taxon_df[taxon_df[cluster_name] == cluster_1]["genome"])
        genomes_2 = sorted(taxon_df[taxon_df[cluster_name] == cluster_2]["genome"])
        genomes_comps += list(itertools.product(genomes_1, genomes_2))
    sorted_genomes_comps = [(g1, g2) if g1 < g2 else (g2, g1) for g1, g2 in genomes_comps]
    genomes_comps = [genome_pair for genome_pair in sorted_genomes_comps if (genome_pair[0], genome_pair[1]) not in already_compared]
    genomes_comps = [(os.path.join(genomes_path, genome_1), os.path.join(genomes_path, genome_2)) for genome_1, genome_2 in genomes_comps]

    # run lastz in parallel
    res_list = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for result in executor.map(lastz_entry, genomes_comps):
            res_list.append(result)
    for res in res_list:
        res = list(res)
        res[2] = res[2].tobytes()
        res[0] = os.path.basename(res[0])
        res[1] = os.path.basename(res[1])
        sqlite3_conn.execute("INSERT INTO lastz VALUES (?, ?, ?, ?)", res)
    sqlite3_conn.commit()
    return res_list
