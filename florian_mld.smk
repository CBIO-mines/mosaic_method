rule merge:
    input:
        f_mlds=config["f_mlds_dir"] + "{species_1}_{species_2}/"
    output:
        full_mld=config["results_dir"] + "full_mlds/{species_1}_vs_{species_2}_full_mld_comp.csv",
        binned_mld=config["results_dir"] + "binned_mlds/{species_1}_vs_{species_2}_binned_mld.csv"
    shell:
        "python parse/main.py --from_florian_mld {input.f_mlds} "
        "--save_full_mld {output.full_mld} {output.binned_mld}"
