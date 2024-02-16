from parse.fun import get_genome_comp

def all_lastz_align(wildcards):
    species = [wildcards.species_1, wildcards.species_2]
    lastz_res_path = f'{config["results_dir"]}{config["lastz"]}'
    res = get_genome_comp(species, config["species_csv"], lastz_res_path)
    return res


rule lastz:
    input:
        spec_1_fa=ancient(config['genomes_dir'] + "{fasta_1}"),
        spec_2_fa=ancient(config['genomes_dir'] + "{fasta_2}")
    output:
        config["results_dir"] + config["lastz"] + "{fasta_1}_vs_{fasta_2}.csv"
    shell:
        "python run_lastz/main.py {input.spec_1_fa} {input.spec_2_fa} {output}"


rule merge:
    input:
        all_aligns=all_lastz_align
    output:
        full_mld=config["results_dir"] + "full_mlds/{species_1}_vs_{species_2}_full_mld_comp.csv",
        binned_mld=config["results_dir"] + "binned_mlds/{species_1}_vs_{species_2}_binned_mld.csv"
    params:
        csv_dir=f'{config["results_dir"]}{config["lastz"]}',
        species_csv=config["species_csv"],
        species=lambda w: f'{w.species_1},{w.species_2}'
    shell:
        "python parse/main.py --from_csv {params.csv_dir} "
        "--save_full_mld {output.full_mld} --species {params.species} "
        "--species_csv {params.species_csv} {output.binned_mld}"
