def all_lastz_align(wildcards):
    res = []
    for fa_file_1, fa_file_2 in itertools.product(sorted(os.listdir(config["species_dir"] + wildcards.species_1)), sorted(os.listdir(config["species_dir"] + wildcards.species_2))):
        if not any(ext in fa_file_1 for ext in ["fna", "fasta", "fa"]) or not any(ext in fa_file_2 for ext in ["fna", "fasta", "fa"]):
            continue
        if fa_file_1.startswith(".") or fa_file_2.startswith("."):
            continue
        res += [f'{config["results_dir"]}{config["lastz"]}{wildcards.species_1}_{wildcards.species_2}/{fa_file_1}_vs_{fa_file_2}.txt']
    return res


rule lastz:
    input:
        spec_1_fa=config["species_dir"] + "{species_1}/{fasta_1}",
        spec_2_fa=config["species_dir"] + "{species_2}/{fasta_2}"
    output:
        temp(config["results_dir"] + config["lastz"] + "{species_1}_{species_2}/{fasta_1}_vs_{fasta_2}.txt")
    shell:
        "lastz {input.spec_1_fa}[multiple] {input.spec_2_fa}[multiple] "
        "--format=general:cigarx --ambiguous=iupac > {output}"


rule merge:
    input:
        all_aligns=all_lastz_align
    output:
        full_mld=config["results_dir"] + "full_mlds/{species_1}_vs_{species_2}_full_mld_comp.csv",
        binned_mld=config["results_dir"] + "binned_mlds/{species_1}_vs_{species_2}_binned_mld.csv"
    params:
        lastz_dir=lambda w: f'{config["results_dir"]}{config["lastz"]}{w.species_1}_{w.species_2}',
    shell:
        "python parse/main.py --from_cigarx {params.lastz_dir}  "
        "--save_full_mld {output.full_mld}  {output.binned_mld}"


rule fit:
    input:
        binned_mld=config["results_dir"] + "binned_mlds/{species_1}_vs_{species_2}_binned_mld.csv"
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
