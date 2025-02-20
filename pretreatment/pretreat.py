#!/usr/bin/env python3
import argparse
import concurrent.futures
import itertools
import os
import re
import subprocess as sp
import tempfile

import Bio.Seq
import Bio.SeqIO
import pandas as pd


def remove_AmbiguousIUPAC(genome_path, output_dir):
    """Replaces letters RYWSMKHBVD by N"""
    records = []
    with open(genome_path, 'r') as genome_handle:
        for record in Bio.SeqIO.parse(genome_handle, 'fasta'):
            seq = re.sub(r'[RYWSMKHBVD]', 'N', str(record.seq))
            record.seq = Bio.Seq.Seq(seq)
            records.append(record)
    for record in records:
        non_valid_chars = set(record.seq) - set("ATCGN")
        if non_valid_chars:
            raise ValueError(f"Invalid characters in {record.id}: {non_valid_chars}")

    output_path = os.path.join(output_dir, os.path.basename(genome_path))
    with open(output_path, 'w') as output_handle:
        Bio.SeqIO.write(records, output_handle, 'fasta')


def mask_repeats(genome_path, output_dir, lastz_tools_dir, above=2, transition=True):
    """Masks repeats in a genome"""
    if transition:
        transition_param = ""
    else:
        transition_param = "--notransition"

    genome_name = os.path.basename(genome_path)
    output_dir = output_dir.rstrip("/")
    output_path = os.path.join(output_dir, genome_name)
    masked_intervals_path = os.path.join(output_dir, genome_name + ".masked_intervals.dat")
    fasta_fragments_py = os.path.join(lastz_tools_dir, "fasta_fragments.py")
    fasta_softmask_intervals_py = os.path.join(lastz_tools_dir, "fasta_softmask_intervals.py")
    detect_repeats_cmd = f"""cat {genome_path} | python {fasta_fragments_py} --fragment=200 --step=100 | \
    lastz {genome_path}[multiple,unmask,nameparse=darkspace] /dev/stdin --masking={above+1} \
    --progress+masking=10K --format=none --outputmasking+:soft={masked_intervals_path} {transition_param}"""
    mask_repeats_cmd = f"cat  {genome_path} | python {fasta_softmask_intervals_py} --origin=1 {masked_intervals_path} > {output_path}"

    # TODO import lastz tools instead of using subprocess
    sp.run(detect_repeats_cmd, shell=True, check=True)
    sp.run(mask_repeats_cmd, shell=True, check=True)



def pretreat_genomes(directory, taxon_csv, output, lastz_tools_dir, above=2, threads=1, transition=True):
    """Masks repeats and replaces ambiguous IUPAC letters by N"""
    taxon_df = pd.read_csv(taxon_csv)
    genomes = taxon_df["genome"].tolist()
    genomes_path = [os.path.join(directory, g) for g in genomes if os.path.exists(os.path.join(directory, g))]
    n_genomes = len(genomes_path)
    if len(genomes) != n_genomes:
        raise ValueError(f"Genomes not found: {set(genomes) - set([os.path.basename(g) for g in genomes_path])}")
    print(f"Found {n_genomes} genomes to pretreat.")
    os.makedirs(output, exist_ok=True)
    with tempfile.TemporaryDirectory() as temp_dir:
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            try:
                executor.map(remove_AmbiguousIUPAC, genomes_path, itertools.repeat(temp_dir))
            except ValueError as e:
                print(e)
                return
        temp_genomes = [os.path.join(temp_dir, g) for g in genomes]
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            executor.map(mask_repeats, temp_genomes, itertools.repeat(output), itertools.repeat(lastz_tools_dir), itertools.repeat(above), itertools.repeat(transition))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Masks repeats and replaces ambiguous IUPAC letters by N")
    parser.add_argument(
        "--threads",
        "-t",
        type=int,
        default=1,
        help="Number of threads to use"
    )
    parser.add_argument(
        "--above",
        "-a",
        type=int,
        default=2,
        help="Mask repeats present in at least this number of copies"
    )
    parser.add_argument(
        "lastz_tools_dir",
        type=str,
        help="Directory containing lastz tools"
    )
    parser.add_argument(
        "directory",
        type=str,
        help="Directory containing genomes to pretreat"
    )
    parser.add_argument(
        "taxon_csv",
        type=str,
        help="CSV file containing genomes and their taxon"
    )
    parser.add_argument(
        "output",
        type=str,
        help="Output directory"
    )
    args = parser.parse_args()
    pretreat_genomes(args.directory, args.taxon_csv, args.output, args.lastz_tools_dir, args.threads)
