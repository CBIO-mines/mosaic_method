#!/usr/bin/env python3
import itertools
import os
from pathlib import Path

genomes = [p for p in Path('../genomes/refseq_bacillaceae_filtered_all/BacillusAll').rglob("*.fna")]

query_file = "test_ani_query.txt"
ref_file = "test_ani_ref.txt"
n_genomes = len(genomes)

with open(query_file, "w") as query:
    for path in genomes[:n_genomes//2]:
        query.write(f"{str(path)}\n")

with open(ref_file, "w") as ref:
    for path in genomes[n_genomes//2:]:
        ref.write(f"{str(path)}\n")
