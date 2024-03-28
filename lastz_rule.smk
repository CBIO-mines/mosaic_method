
rule lastz:
    output:
        config["results_dir"] + "lastz_sqlite_database.db"
    params:
        genomes_dir=config["genomes_dir"],
        taxon_csv=config["taxon_csv"],
        cluster_name=config["cluster_name"]
    threads: config["lastz_threads"]
    shell:
        "python lastz_parallel_db/main.py -t {threads} {params.taxon_csv} "
        "{params.genomes_dir} {params.cluster_name} {output}"



rule merge:
    input:
        config["results_dir"] + "lastz_sqlite_database.db"
    output:
        full_mld=[f"{config['results_dir']}full_mlds/{bac1}_vs_{bac2}_full_mld_comp.csv" for bac1, bac2 in itertools.combinations(CLUSTER_LIST, 2)],
        binned_mld=[f"{config['results_dir']}binned_mlds/{bac1}_vs_{bac2}_binned_mld.csv" for bac1, bac2 in itertools.combinations(CLUSTER_LIST, 2)]
    params:
        taxon_csv=config["taxon_csv"],
        cluster_name=config["cluster_name"],
        full_mld_dir=config["results_dir"] + "full_mlds/",
        binned_mld_dir=config["results_dir"] + "binned_mlds/"
    shell:
        "python parse/main.py --from_sqlite_db {input} --full_mld {params.full_mld_dir} "
        "--binned_mld {params.binned_mld_dir} --taxon_csv {params.taxon_csv} --cluster_name {params.cluster_name}"
