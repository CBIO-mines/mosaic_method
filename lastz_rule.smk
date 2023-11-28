def all_lastz_align(wildcards):
    res = []
    for fa_file_1, fa_file_2 in itertools.product(sorted(os.listdir(config["species_dir"] + wildcards.species_1)), sorted(os.listdir(config["species_dir"] + wildcards.species_2))):
        if fa_file_1.endswith(".fa") and fa_file_2.endswith(".fa"):
            res += [f'{config["lastz"]}{wildcards.species_1}_{wildcards.species_2}/{fa_file_1[:-3]}_{fa_file_2[:-3]}.txt']
        if fa_file_1.endswith(".fasta") and fa_file_2.endswith(".fasta"):
            res += [f'{config["lastz"]}{wildcards.species_1}_{wildcards.species_2}/{fa_file_1[:-6]}_{fa_file_2[:-6]}.txt']
    return res


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
        binned_mld="binned_mlds/{species_1}_{species_2}_binned_mld.csv",
        surface_plot="surfaces/{species_1}_{species_2}_surface_plot.png"
    params:
        lastz_dir=lambda w: f'{config["lastz"]}{w.species_1}_{w.species_2}',
        species=lambda w: f"{w.species_1},{w.species_2}"?
        mus=config["mus"],
        muc=config["muc"],
        delta=config["delta"]
    shell:
        "python parsefit.py --bacs {params.species} --from_cigarx {params.lastz_dir} --L0 {input.all_L0s} "
        "--mus {params.mus} --muc {params.muc} --delta {params.delta} "
        "--save_full_mld {output.full_mld} --save_surface_plot {output.surface_plot} {output.fitted_params} {output.binned_mld}"
