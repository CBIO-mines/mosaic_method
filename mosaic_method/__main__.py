import argparse
import os
import shutil
import yaml

from .inference import run_inference


def main():
    parser = argparse.ArgumentParser(description="""
        Run the inference pipeline without snakemake.
        """)
    parser.add_argument("config", type=str, help="The path to the configuration file.")
    parser.add_argument(
        "--genomes_dir",
        type=str,
        help="The directory containing the genomes (overwrites the one in the config file)",
        default=None,
    )

    args = parser.parse_args()
    with open(args.config, "r") as config_file:
        cfg = yaml.safe_load(config_file)
    os.makedirs(cfg["results_dir"], exist_ok=True)
    if not os.path.exists(
        os.path.join(cfg["results_dir"], os.path.basename(args.config))
    ):
        shutil.copy(args.config, cfg["results_dir"])
    run_inference(cfg, args.genomes_dir)
