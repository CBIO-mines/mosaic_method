import concurrent.futures
import io
import itertools
import multiprocessing
import os
import re
import subprocess as sp
import tempfile
import sqlite3
import pandas as pd
import numpy as np

from mosaic_method.lastz_alignments import run_lastz
from mosaic_method.mummer_alignments import run_mummer
from mosaic_method.no_aligner_alignments import run_no_aligner


def align_exec(genome_pair, align="lastz", prefix=""):
    """
    computes a lastz entry for two genomes
    """
    genome_1, genome_2 = genome_pair
    if align == "lastz":
        align_dic = run_lastz(genome_1, genome_2)
    elif align == "mummer":
        align_dic = run_mummer(genome_1, genome_2, prefix)
    else:
        align_dic = run_no_aligner(genome_1, genome_2)
    genome_small = min(genome_1, genome_2)
    genome_large = max(genome_1, genome_2)
    return genome_small, genome_large, align_dic


def align_entry(align_res, cur):
    """
    inserts a result in database
    """
    res = []
    res.append(os.path.basename(align_res[0]))
    res.append(os.path.basename(align_res[1]))
    res.append(align_res[2]["result"][0].tobytes())
    res.append(align_res[2]["result"][1])
    res.append(align_res[2]["result"][2])
    res.append(align_res[2]["result"][3])
    cur.execute("INSERT INTO lastz VALUES (?, ?, ?, ?, ?, ?)", res)


def create_lastz_db(
    taxon_csv,
    genomes_path,
    cluster_name,
    db_name,
    threads,
    update=False,
    aligner="lastz",
):
    """Creates or updates a sqlite3 database from a taxon csv file, generating all necessary alignments"""
    # connect to database
    if not update:
        if os.path.exists(db_name):
            os.remove(db_name)
    sqlite3_conn = sqlite3.connect(db_name)
    taxon_df = pd.read_csv(taxon_csv)
    cur = sqlite3_conn.cursor()
    sorted_clusters = sorted(taxon_df[cluster_name].unique())

    # Gather genome comps
    if update:
        cur.executemany(
            "INSERT OR IGNORE INTO taxon VALUES (?, ?)",
            taxon_df[["genome", cluster_name]].itertuples(index=False),
        )
        already_compared = cur.execute("SELECT genome1, genome2 FROM lastz").fetchall()
        # gather necessary genome comparisons
        genomes_comps = []
        for cluster_1, cluster_2 in itertools.combinations(sorted_clusters, 2):
            genomes_1 = sorted(taxon_df[taxon_df[cluster_name] == cluster_1]["genome"])
            genomes_2 = sorted(taxon_df[taxon_df[cluster_name] == cluster_2]["genome"])
            genomes_comps += list(itertools.product(genomes_1, genomes_2))

        sorted_genomes_comps = [
            (g1, g2) if g1 < g2 else (g2, g1) for g1, g2 in genomes_comps
        ]
        genomes_comps = [
            genome_pair
            for genome_pair in sorted_genomes_comps
            if (genome_pair[0], genome_pair[1]) not in already_compared
        ]
        genomes_comps = [
            (os.path.join(genomes_path, genome_1), os.path.join(genomes_path, genome_2))
            for genome_1, genome_2 in genomes_comps
        ]
        # print genome comps and already compared
        print(
            f"Already in database : {len(already_compared)} comparisons, {len(genomes_comps)} to align, total : {len(sorted_genomes_comps)} comparisons."
        )
    else:
        cur.execute(
            "CREATE TABLE lastz (genome1 STRING, genome2 STRING, count_array blob, average_divergence REAL, total_matches INT, total_aligned INT);"
        )
        cur.execute("CREATE TABLE taxon (genome STRING, cluster STRING);")
        cur.executemany(
            "INSERT INTO taxon VALUES (?, ?)",
            taxon_df[["genome", cluster_name]].itertuples(index=False),
        )
        # gather necessary genome comparisons
        genomes_comps = []
        for cluster_1, cluster_2 in itertools.combinations(sorted_clusters, 2):
            genomes_1 = taxon_df[taxon_df[cluster_name] == cluster_1]["genome"]
            genomes_2 = taxon_df[taxon_df[cluster_name] == cluster_2]["genome"]
            genomes_1 = [
                os.path.join(genomes_path, genome) for genome in sorted(genomes_1)
            ]
            genomes_2 = [
                os.path.join(genomes_path, genome) for genome in sorted(genomes_2)
            ]
            genomes_comps += list(itertools.product(genomes_1, genomes_2))

    if genomes_comps == []:
        print("No genomes to align")
        return sqlite3_conn

    # run aligner
    batch_size = 200
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=threads, mp_context=multiprocessing.get_context("spawn")
    ) as executor:
        for i in range(0, len(genomes_comps), batch_size):
            num_genomes = len(genomes_comps[i : i + batch_size])
            prefixes = [f"mummer_{j}" for j in range(num_genomes)]
            # this dictionary is supposed to tell me in which comp an error occured
            future_res = {
                executor.submit(align_exec, genome_comp, aligner, prefix): genome_comp
                for genome_comp, aligner, prefix in zip(
                    genomes_comps[i : i + batch_size], [aligner] * num_genomes, prefixes
                )
            }
            for future in concurrent.futures.as_completed(future_res):
                try:
                    res = future.result()
                except Exception as e:
                    print(f"{future_res[future]} raised an error : {e}")
                    raise
                else:
                    align_entry(res, cur)
            sqlite3_conn.commit()

    cur.execute(
        "CREATE INDEX IF NOT EXISTS idx_lastz_genome1_genome2 ON lastz(genome1, genome2);"
    )
    return sqlite3_conn
