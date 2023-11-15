configfile: "./smk_config.yml"
import itertools
import os


with open(config["species_file"], "r") as filein:
    SPECIES_LIST = sorted(filein.read().splitlines())

def all_lastz_align(wildcards):
    res = []
    for fa_file_1, fa_file_2 in itertools.product(sorted(os.listdir(wildcards.species_1)), sorted(os.listdir(wildcards.species_2))):
        if fa_file_1.endswith(".fa") and fa_file_2.endswith(".fa"):
            res += [f'{config["lastz"]}/{wildcards.species_1}_{wildcards.species_2}/{fa_file_1[:-3]}_{fa_file_2[:-3]}.txt']
        if fa_file_1.endswith(".fasta") and fa_file_2.endswith(".fasta"):
            res += [f'{config["lastz"]}/{wildcards.species_1}_{wildcards.species_2}/{fa_file_1[:-6]}_{fa_file_2[:-6]}.txt']
    return res

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
        plot=[f"{bac1}_{bac2}_plot_fig2.png" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        full_mld=[f"full_mlds/{bac1}_{bac2}_full_mld_comp.csv" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        fitted_params=[f"fitted_params/{bac1}_{bac2}_fitted_params.csv" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        binned_mld=[f"binned_mlds/{bac1}_{bac2}_binned_mld.csv" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        lengths_a=expand("lengths_distributions/{fasta_dir}_distribution.{ext}", fasta_dir = SPECIES_LIST, ext = ["png", "csv"]),
        L0s="all_L0s.csv"


rule lastz:
    input:
        spec_1_fa=config["species_dir"] + "{species_1}/{fasta_1}.fasta",
        spec_2_fa=config["species_dir"] + "{species_2}/{fasta_2}.fasta"
    output:
        config["lastz"] + "{species_1}_{species_2}/{fasta_1}_{fasta_2}.txt"
    shell:
        "lastz {input.spec_1_fa}[multiple] {input.spec_2_fa} "
        "--format=general:cigarx > {output}"

rule merge_fit:
    input:
        all_aligns=all_lastz_align,
        all_L0s="all_L0s.csv"
    output:
        full_mld="full_mlds/{species_1}_{species_2}_full_mld_comp.csv",
        fitted_params="fitted_params/{species_1}_{species_2}_fitted_params.csv",
        binned_mld="binned_mlds/{species_1}_{species_2}_binned_mld.csv"
    params:
        lastz_dir=lambda wildcards: f'{config["lastz"]}/{wildcards.species_1}_{wildcards.species_2}',
    shell:
        "python parsefit.py --from_cigarx {params.lastz_dir} --L0 {input.all_L0s} "
        "--save_full_mld {output.full_mld} {output.fitted_params} {output.binned_mld}"


rule plot:
    input:
        binned_mld="binned_mlds/{species_1}_{species_2}_binned_mld.csv",
        fitted_params="fitted_params/{species_1}_{species_2}_fitted_params.csv"
    output:
        "{species_1}_{species_2}_plot_fig2.png"
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
