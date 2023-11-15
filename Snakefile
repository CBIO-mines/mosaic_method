
import itertools
import os

ESCH_DIR = "EscherichiaColi"
SALM_DIR = "Salmonella"
ESCH_COL = [fa_file[:-6] for fa_file in os.listdir(ESCH_DIR) if "fasta" in fa_file]
SALMO = [fa_file[:-6] for fa_file in os.listdir(SALM_DIR) if "fasta" in fa_file]

onstart:
    print("##### Creating profile pipeline #####\n")
    print("\t Creating jobs output subfolders...\n")
    shell("mkdir -p jobs/lastz")
    shell("mkdir -p jobs/merge_fit")
    shell("mkdir -p jobs/plot")
    shell("mkdir -p jobs/lengths")
    shell("mkdir -p jobs/L0")

rule all:
    input:
        plot="plot_fig2.png",
        full_mld="full_mld_comp.csv",
        fitted_params="fitted_params.csv",
        binned_mld="binned_mld.csv",
        lengths_a=expand("{fasta_dir}_distribution.{ext}", fasta_dir = [ESCH_DIR, SALM_DIR], ext = ["png", "csv"])
        L0s="all_L0s.csv"


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
        all_L0s="all_L0s.csv"
    output:
        full_mld="full_mld_comp.csv",
        fitted_params="fitted_params.csv",
        binned_mld="binned_mld.csv"
    params:
        lastz_dir="lastz_out_2",
    shell:
        "python parsefit.py --from_cigarx {params.lastz_dir} --L0 {input.all_L0s} "
        "--save_full_mld {output.full_mld} {output.fitted_params} {output.binned_mld}"


rule plot:
    input:
        binned_mld="binned_mld.csv",
        fitted_params="fitted_params.csv"
    output:
        "plot_fig2.png"
    script:
        "plot_mld_fit.R"


rule lengths:
    input:
        config["species_dir"] + "{fasta_directory}"
    output:
        csv_distr="lengths_distributions/{fasta_directory}_distribution.csv",
        histo="lengths_distributions/{fasta_directory}_distribution.png"
    shell:
        "python length_analysis.py --save_distr {output.csv_distr} --save_plot {output.histo} {input}"


rule L0:
    input:
        distribs=expand("lengths_distributions/{fasta_dir}_distribution.csv", fasta_dir = SPECIES_LIST)
    output:
        "all_L0s.csv"
    params:
        distr_dir="lengths_distributions/"
    shell:
        "python get_L0.py {params.distr_dir} {output}"
