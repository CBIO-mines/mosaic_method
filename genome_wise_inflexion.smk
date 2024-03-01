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
        res_analyse = [f"{config['results_dir']}analyse_comparisons/{bac1}_vs_{bac2}_fitted_single_mlds_plot.png" for bac1, bac2 in itertools.combinations(CLUSTER_LIST, 2)]
    output:
        config['results_dir'] + "inflexion_exists.csv"
    params:
        fitted_params_dir=config["results_dir"] + "fitted_params/",
        analyse_dir=config["results_dir"] + "analyse_comparisons/",
        results_dir=config["results_dir"],
        min_r_infl=30
    script:
        "gather_comparisons.R"
