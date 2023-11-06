
import itertools
import os

ESCH_COL = [fa_file[:-6] for fa_file in os.listdir("EscherichiaColi") if "fasta" in fa_file]
SALMO = [fa_file[:-6] for fa_file in os.listdir("Salmonella") if "fasta" in fa_file]

onstart:
    print("##### Creating profile pipeline #####\n")
    print("\t Creating jobs output subfolders...\n")
    shell("mkdir -p jobs/lastz")
    shell("mkdir -p jobs/merge_fit")
    shell("mkdir -p jobs/plot")

rule all:
    input:
        plot="plot_fig2.png",
        full_mld="full_mld_comp.csv",
        fitted_params="fitted_params.csv",
        binned_mld="binned_mld.csv"


rule lastz:
    input:
        salmo_fa="Salmonella/{salmo_genome_fa}.fasta",
        esch_fa="EscherichiaColi/{esch_genome_fa}.fasta"
    output:
        "lastz_out_2/{esch_genome_fa}_{salmo_genome_fa}.txt"
    shell:
        "lastz {input.salmo_fa}[multiple] {input.esch_fa} "
        "--format=general:cigarx > {output}"

rule merge_fit:
    input:
        expand("lastz_out_2/{esch}_{salmo}.txt", esch=ESCH_COL, salmo=SALMO)
    output:
        full_mld="full_mld_comp.csv",
        fitted_params="fitted_params.csv",
        binned_mld="binned_mld.csv"
    params:
        lastz_dir="lastz_out_2",
        L0=4903888.5
    shell:
        "python parsefit.py --from_cigarx {params.lastz_dir} --L0 {params.L0} "
        "--full_mld {output.full_mld} {output.fitted_params} {output.binned_mld}"

rule plot:
    input:
        binned_mld="binned_mld.csv",
        fitted_params="fitted_params.csv"
    output:
        "plot_fig2.png"
    script:
        "plot_mld_fit.R"
