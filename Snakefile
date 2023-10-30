
import itertools
import os

ESCH_COL = [fa_file[:-6] for fa_file in os.listdir("EscherichiaColi") if "fasta" in fa_file]
SALMO = [fa_file[:-6] for fa_file in os.listdir("Salmonella") if "fasta" in fa_file]

onstart:
    print("##### Creating profile pipeline #####\n")
    print("\t Creating jobs output subfolders...\n")
    shell("mkdir -p jobs/lastz")
    shell("mkdir -p jobs/concat")
    shell("mkdir -p jobs/plot")

rule all:
    input:
        "matches_matrix_2.csv",
        "plot_fig2_2.png"


rule lastz:
    input:
        salmo_fa="Salmonella/{salmo_genome_fa}.fasta",
        esch_fa="EscherichiaColi/{esch_genome_fa}.fasta"
    output:
        "lastz_out_2/{esch_genome_fa}_{salmo_genome_fa}.txt"
    shell:
        "lastz {input.salmo_fa}[multiple] {input.esch_fa} "
        "--format=general:cigarx > {output}"

rule concat:
    input:
        expand("lastz_out_2/{esch}_{salmo}.txt", esch=ESCH_COL, salmo=SALMO)
    output:
        "matches_matrix_2.csv"
    params:
        lastz_dir="lastz_out_2"
    shell:
        "python parse_cigars.py {params.lastz_dir} {output}"

rule plot:
    input:
        "matches_matrix_2.csv"
    output:
        "plot_fig2_2.png"
    script:
        "analyse_matches.R"
