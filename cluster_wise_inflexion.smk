rule overall_inflexion:
    input:
        fitted_params=[f"{config['results_dir']}fitted_params/{bac1}_vs_{bac2}_fitted_params.csv" for bac1, bac2 in itertools.combinations(CLUSTER_LIST, 2)]
    output:
        config['results_dir'] + "inflexion_exists.csv",
        config['results_dir'] + "inflexion_by_cluster.csv"
    params:
        fitted_params_dir=config["results_dir"] + "fitted_params/",
        results_dir=config["results_dir"],
        min_r_infl=15,
        taxon_csv=config["taxon_csv"],
        cluster_name=config["cluster_name"]
    script:
        "overall_inflexion.R"
