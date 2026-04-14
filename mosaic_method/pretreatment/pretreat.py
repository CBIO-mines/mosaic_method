#!/usr/bin/env python3
import concurrent.futures
from importlib.resources import files
import itertools
import os
import re
import subprocess as sp
import tempfile

import Bio.Seq
import Bio.SeqIO
import pandas as pd

LASTZ_TOOL_DIR = files("mosaic_method")

def get_lastz_tool_dir():
    return str(files("mosaic_method") / "pretreatment" / "tools")

def remove_AmbiguousIUPAC(genome_path, output_dir):
    """Replaces letters RYWSMKHBVD by N"""
    records = []
    with open(genome_path, 'r') as genome_handle:
        for record in Bio.SeqIO.parse(genome_handle, 'fasta'):
            seq = re.sub(r'[RYWSMKHBVDrywsmkhbvd]', 'N', str(record.seq))
            record.seq = Bio.Seq.Seq(seq)
            records.append(record)
    for record in records:
        non_valid_chars = set(record.seq) - set("ATCGN")
        if non_valid_chars:
            raise ValueError(f"Invalid characters in {record.id}: {non_valid_chars}")

    output_path = os.path.join(output_dir, os.path.basename(genome_path))
    with open(output_path, 'w') as output_handle:
        Bio.SeqIO.write(records, output_handle, 'fasta')


def mask_repeats(genome_path, output_dir, above=2, transition=True, softmask=False):
    """Masks repeats in a genome"""
    genome_name = os.path.basename(genome_path)
    output_path = os.path.join(output_dir, genome_name)
    masked_intervals_path = os.path.join(output_dir, genome_name + ".masked_intervals.dat")
    fasta_fragments_py = os.path.join(get_lastz_tool_dir(), "fasta_fragments.py")
    fasta_softmask_intervals_py = os.path.join(get_lastz_tool_dir(), "fasta_softmask_intervals.py")

    # Step 1: generate fragments into a temp file, pass to LASTZ as a regular path
    lastz_cmd = [
        'lastz',
        f'{genome_path}[multiple,unmask,nameparse=darkspace]',
        # placeholder, filled in below once we have frag_path
        f'--masking={above + 1}',
        '--progress+masking=10K',
        '--format=none',
        f'--outputmasking+:soft={masked_intervals_path}',
    ]
    if not transition:
        lastz_cmd.append('--notransition')

    frag_f = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
    frag_path = frag_f.name
    try:
        with open(genome_path, 'r') as genome_f:
            sp.run(
                ['python', fasta_fragments_py, '--fragment=200', '--step=100'],
                stdin=genome_f, stdout=frag_f, check=True
            )
        frag_f.close()
        sp.run(lastz_cmd + [frag_path], check=True)
    finally:
        frag_f.close()
        os.remove(frag_path)

    # Step 2: apply mask intervals, redirect stdout to output file directly
    mask_cmd = [
        'python', fasta_softmask_intervals_py,
        '--origin=1',
        masked_intervals_path,
    ]
    if not softmask:
        mask_cmd += ['--mask=N']

    with open(genome_path, 'r') as genome_f, open(output_path, 'w') as output_f:
        sp.run(mask_cmd, stdin=genome_f, stdout=output_f, check=True)

    os.remove(masked_intervals_path)


def pretreat_genomes(directory, taxon_csv, output, above=2, threads=1, transition=True):
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
                list(executor.map(remove_AmbiguousIUPAC, genomes_path, itertools.repeat(temp_dir)))
            except ValueError as e:
                print(e)
                return

        temp_genomes = [os.path.join(temp_dir, g) for g in genomes]
        with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
            list(executor.map(mask_repeats, temp_genomes, itertools.repeat(output), itertools.repeat(above), itertools.repeat(transition)))
