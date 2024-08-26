rule database:
    input:
        genomes=[os.path.join(config["pretreatment_dir"], g) for g in GENOME_LIST]
    params:
        up=lambda wc: "--update" if config["database_state"] == "update" else "",
        taxon_csv=config["taxon_csv"],
        genomes_dir=config["pretreatment_dir"],
        cluster_name=config["cluster_name"]
    threads: config["max_threads"]
    output:
        database
    shell:
        "python lastz_parallel_db/main.py --threads {threads} {params.up}  {params.taxon_csv} {params.genomes_dir} "
        "{params.cluster_name}  {output}"
