rule merge_fit:
    input:
        f_mlds=config["f_mlds_dir"] + "{species_1}_{species_2}/",
        all_L0s=config["results_dir"] + "all_L0s.csv"
    output:
        full_mld=config["results_dir"] + "full_mlds/{species_1}_vs_{species_2}_full_mld_comp.csv",
        fitted_params=config["results_dir"] + "fitted_params/{species_1}_vs_{species_2}_fitted_params.csv",
        binned_mld=config["results_dir"] + "binned_mlds/{species_1}_vs_{species_2}_binned_mld.csv",
        surface_plot=config["results_dir"] + "surfaces/{species_1}_vs_{species_2}_surface_plot.png"
    params:
        species=lambda w: f"{w.species_1},{w.species_2}",
        mus=config["mus"],
        muc=config["muc"],
        delta=config["delta"]
    shell:
        "python parsefit.py --bacs {params.species} --from_florian_mld {input.f_mlds} --L0 {input.all_L0s} "
        "--mus {params.mus} --muc {params.muc} --delta {params.delta} "
        "--save_full_mld {output.full_mld} --save_surface_plot {output.surface_plot} {output.fitted_params} {output.binned_mld}"
