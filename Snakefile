configfile: "./smk_config.yml"
import itertools
import os


with open(config["species_file"], "r") as filein:
    SPECIES_LIST = sorted(filein.read().splitlines())

if config["from_f_mld"]:
    include: "florian_mld.smk"
else:
    include: "lastz_rule.smk"

onstart:
    print("##### Creating profile pipeline #####\n")
    print("\t Creating jobs output subfolders...\n")
    shell("mkdir -p jobs/merge_fit")
    shell("mkdir -p jobs/plot")
    shell("mkdir -p jobs/lengths")
    shell("mkdir -p jobs/L0")
    if not config["from_f_mld"] == "yes":
        shell("mkdir -p jobs/lastz")

rule all:
    input:
        plot=[f"{bac1}_{bac2}_plot_fig2.png" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        full_mld=[f"full_mlds/{bac1}_{bac2}_full_mld_comp.csv" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        fitted_params=[f"fitted_params/{bac1}_{bac2}_fitted_params.csv" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        binned_mld=[f"binned_mlds/{bac1}_{bac2}_binned_mld.csv" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        lengths_a=expand("lengths_distributions/{fasta_dir}_distribution.{ext}", fasta_dir = SPECIES_LIST, ext = ["png", "csv"]),
        surfaces=[f"surfaces/{bac1}_{bac2}_surface_plot.png" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        L0s="all_L0s.csv"


rule plot:
    input:
        binned_mld="binned_mlds/{species_1}_{species_2}_binned_mld.csv",
        fitted_params="fitted_params/{species_1}_{species_2}_fitted_params.csv"
    output:
        "{species_1}_{species_2}_plot_fig2.png"
    params:
        species=lambda w: f"{w.species_1},{w.species_2}",
        mus=config["mus"],
        muc=config["muc"],
        delta=config["delta"]
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
