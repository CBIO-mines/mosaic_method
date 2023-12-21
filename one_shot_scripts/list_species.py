#!/usr/bin/env python3

import os

def species_enough_fasta(base_path, x):
    result = []
    for root, dirs, files in os.walk(base_path):
        if root != base_path and len(files) > x and not any(directory.endswith("all") for directory in dirs):#not root.endswith("all"):
            result.append(root.split("/")[-1])
    return result

# Example usage:
base_directory = "/cluster/CBIO/data1/petheimer/refseq_bacillaceae/"
x_value = 0  # Change this to your desired threshold
result_directories = species_enough_fasta(base_directory, x_value)


with open("species_list_bacillaceae_all.txt", "w") as fileout:
    for directory in result_directories:
        fileout.write(f"{directory}\n")
