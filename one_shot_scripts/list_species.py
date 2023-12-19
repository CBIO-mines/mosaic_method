#!/usr/bin/env python3

import os

def species_enough_fasta(base_path, x):
    result = []
    for root, dirs, files in os.walk(base_path):
        if len(files) > x:
            result.append(root.split("/")[-1])
    return result

# Example usage:
base_directory = "/cluster/CBIO/data1/petheimer/refseq_bacillaceae/"
x_value = 4  # Change this to your desired threshold
result_directories = species_enough_fasta(base_directory, x_value)

with open("species_list_bacillaceae_5.txt", "w") as fileout:
    for directory in result_directories:
        fileout.write(f"{directory}\n")
