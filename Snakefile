configfile: "./smk_config.yml"
import itertools
import os


with open(config["species_file"], "r") as filein:
    SPECIES_LIST = sorted(filein.read().splitlines())

if config["from_f_mld"] == "yes":
    include: "florian_mld.smk"
else:
    include: "lastz_rule.smk"

onstart:
    print("##### Creating profile pipeline #####\n")
    print("\t Creating jobs output subfolders...\n")
    shell("mkdir -p jobs/fit")
    shell("mkdir -p jobs/merge")
    shell("mkdir -p jobs/plot")
    shell("mkdir -p jobs/lengths")
    shell("mkdir -p jobs/L0")
    shell("mkdir -p jobs/trees")
    shell("mkdir -p jobs/analyse_comparisons")
    shell("mkdir -p jobs/gather_comparisons")
    if not config["from_f_mld"] == "yes":
        shell("mkdir -p jobs/lastz")

rule all:
    input:
        plot=[f"{config['results_dir']}fig2_plots/{bac1}_vs_{bac2}_plot_fig2.png" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        full_mld=[f"{config['results_dir']}full_mlds/{bac1}_vs_{bac2}_full_mld_comp.csv" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        fitted_params=[f"{config['results_dir']}fitted_params/{bac1}_vs_{bac2}_fitted_params.csv" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        binned_mld=[f"{config['results_dir']}binned_mlds/{bac1}_vs_{bac2}_binned_mld.csv" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        lengths_a=[config['results_dir'] + len_distr for len_distr in expand("lengths_distributions/{fasta_dir}_distribution.{ext}", fasta_dir = SPECIES_LIST, ext = ["png", "csv"])],
        surfaces=[f"{config['results_dir']}surfaces/{bac1}_vs_{bac2}_surface_plot.png" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        L0s=f"{config['results_dir']}all_L0s.csv",
        tree=config["results_dir"] + "family_tree.svg",
        comparisons=[f"{config['results_dir']}analyse_comparisons/{bac1}_vs_{bac2}_inflexion_res.csv" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)]


rule plot:
    input:
        binned_mld=config["results_dir"] + "binned_mlds/{species_1}_vs_{species_2}_binned_mld.csv",
        fitted_params=config["results_dir"] + "fitted_params/{species_1}_vs_{species_2}_fitted_params.csv"
    output:
        config["results_dir"] + "fig2_plots/{species_1}_vs_{species_2}_plot_fig2.png"
    params:
        species=lambda w: f"{w.species_1},{w.species_2}",
        results_dir=config["results_dir"],
        mus=config["mus"],
        muc=config["muc"],
        delta=config["delta"]
    script:
        "plot_mld_fit.R"


rule lengths:
    input:
        config["genomes_dir"]
    output:
        csv_distr=expand(config["results_dir"] + "lengths_distributions/{species}_distribution.csv", species = SPECIES_LIST),
        histo=expand(config["results_dir"] + "lengths_distributions/{species}_distribution.png", species = SPECIES_LIST)
    params:
        species_csv=config["species_csv"],
        output_dir=config["results_dir"] + "lengths_distributions/"
    shell:
        "python length_analysis.py --species {params.species_csv} --save_dir {params.output_dir} {input}"


rule L0:
    input:
        distribs=[config['results_dir'] + len_distr for len_distr in expand("lengths_distributions/{fasta_dir}_distribution.csv", fasta_dir = SPECIES_LIST)]
    output:
        config["results_dir"] + "all_L0s.csv"
    params:
        distr_dir=config["results_dir"] + "lengths_distributions/"
    shell:
        "python get_L0.py {params.distr_dir} {output}"


rule fit:
    input:
        binned_mld=config["results_dir"] + "binned_mlds/{species_1}_vs_{species_2}_binned_mld.csv",
        all_L0s=config["results_dir"] + "all_L0s.csv"
    output:
        fitted_params=config["results_dir"] + "fitted_params/{species_1}_vs_{species_2}_fitted_params.csv",
        surface_plot=config["results_dir"] + "surfaces/{species_1}_vs_{species_2}_surface_plot.png"
    params:
        species=lambda w: f"{w.species_1},{w.species_2}",
        mus=config["mus"],
        muc=config["muc"],
        delta=config["delta"]
    shell:
        "python fit/main.py --bacs {params.species}  --L0 {input.all_L0s} "
        "--mus {params.mus} --muc {params.muc} --delta {params.delta} "
        "--save_surface_plot {output.surface_plot} {input.binned_mld} {output.fitted_params}"


rule trees:
    input:
        fitted_params=[f"{config['results_dir']}fitted_params/{bac1}_vs_{bac2}_fitted_params.csv" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)],
        lengths=[f"{config['results_dir']}lengths_distributions/{bac}_distribution.csv" for bac in SPECIES_LIST],
        inflexions=config["results_dir"] + "inflexion_exists.csv"
    output:
        config["results_dir"] + "family_tree.svg",
        config["results_dir"] + "fitteddistance_vs_founddistance.png",
        config["results_dir"] + "hist_fitteddistance.png"
    params:
        metadata=config["metadata"],
        species_list=SPECIES_LIST,
        fitted_params_dir=config["results_dir"] + "fitted_params/",
        mock="yes",
        results_dir=config["results_dir"],
        use_inflexion="no",
        genome_lengths=config["results_dir"] + "lengths_distributions/"
    script:
        "make_trees.R"


rule analyse_comparisons:
    input:
        config['results_dir'] + "fitted_params/{bac1}_vs_{bac2}_fitted_params.csv"
    output:
        config['results_dir'] + "analyse_comparisons/{bac1}_vs_{bac2}_fitted_single_mlds_plot.png",
        config['results_dir'] + "analyse_comparisons/{bac1}_vs_{bac2}_inflexion_res.csv",
        config['results_dir'] + "analyse_comparisons/{bac1}_vs_{bac2}_single_comp_r_infl.csv"
    params:
        species=lambda w: f"{w.bac1},{w.bac2}",
        fitted_params_dir=config["results_dir"] + "fitted_params/",
        full_mlds_dir=config["results_dir"] + "full_mlds/",
        binned_mld_dir=config["results_dir"] + "binned_mlds/",
        results_dir=config["results_dir"],
        min_r_infl=50
    script:
        "analyse_comparisons.R"


rule gather_comparisons:
    input:
        res_analyse = [f"{config['results_dir']}analyse_comparisons/{bac1}_vs_{bac2}_fitted_single_mlds_plot.png" for bac1, bac2 in itertools.combinations(SPECIES_LIST, 2)]
    output:
        config['results_dir'] + "inflexion_exists.csv"
    params:
        fitted_params_dir=config["results_dir"] + "fitted_params/",
        analyse_dir=config["results_dir"] + "analyse_comparisons/",
        results_dir=config["results_dir"],
        min_r_infl=30
    script:
        "gather_comparisons.R"
