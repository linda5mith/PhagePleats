import argparse
import subprocess
import sys
from importlib.resources import files
import warnings
import os

os.environ["PYTHONWARNINGS"] = "ignore"
warnings.simplefilter("ignore")

def main():
    parser = argparse.ArgumentParser(
        prog="phagepleats",
        description="PhagePleats: structural phage taxonomy prediction",
    )

    subparsers = parser.add_subparsers(dest="cmd", required=True)

    # --------------------
    # run subcommand
    # --------------------
    run = subparsers.add_parser(
        "run",
        help="Run the PhagePleats pipeline",
    )

    run.add_argument(
        "--path_to_pdbs",
        required=True,
        help="Directory containing input PDB files",
    )
    run.add_argument(
        "--metadata",
        required=True,
        help="Genome–protein mapping metadata CSV",
    )
    run.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )
    run.add_argument(
        "--cores",
        type=int,
        default=8,
        help="Number of CPU cores (default: 8)",
    )

    args = parser.parse_args()

    if args.cmd == "run":
        snakefile = files("phagepleats") / "resources" / "Snakefile"

        cmd = [
            "snakemake",
            "-s", str(snakefile),
            "--cores", str(args.cores),
            "--config",
            f"PDBs={args.path_to_pdbs}",
            f"metadata={args.metadata}",
            f"outdir={args.outdir}",
        ]

        sys.exit(subprocess.call(cmd))
