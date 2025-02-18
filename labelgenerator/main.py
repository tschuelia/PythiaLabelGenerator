import argparse
import pathlib
import shutil
import sys
import time

from labelgenerator.logger import get_header, logger, log_runtime_information
from labelgenerator.raxmlng import infer_ml_trees, rf_distance
from pypythia.msa import parse


DEFAULT_RAXMLNG_EXE = (
    pathlib.Path(shutil.which("raxml-ng")) if shutil.which("raxml-ng") else None
)

DEFAULT_IQTREE_EXE = (
    pathlib.Path(shutil.which("iqtree2")) if shutil.which("iqtree2") else None
)


def _parse_cli():
    parser = argparse.ArgumentParser(description="Generate the ground truth difficulty for the given MSA.")
    parser.add_argument(
        "-m",
        "--msa",
        type=str,
        required=True,
        help="Multiple Sequence Alignment to compute the ground truth difficulty for. "
             "Must be in either phylip or fasta format.",
    )

    parser.add_argument(
        "-r",
        "--raxmlng",
        type=str,
        default=DEFAULT_RAXMLNG_EXE,
        required=DEFAULT_RAXMLNG_EXE is None,
        help="Path to the binary of RAxML-NG. For install instructions see https://github.com/amkozlov/raxml-ng."
             "(default: 'raxml-ng' if in $PATH, otherwise this option is mandatory).",
    )

    parser.add_argument(
        "-i",
        "--iqtree",
        type=str,
        default=DEFAULT_IQTREE_EXE,
        required=DEFAULT_IQTREE_EXE is None,
        help="Path to the binary of IQ-TREE2. For install instructions see http://www.iqtree.org."
             "(default: 'iqtree2' if in $PATH, otherwise this option is mandatory).",
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        required=False,
        help="Number of threads to use for the RAxML-NG tree inference and IQ-TREE statistical tests (default: all available).",
    )

    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        default=0,
        required=False,
        help="Seed for the RAxML-NG tree inference (default: 0).",
    )

    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=False,
        help="Prefix of the RAxML-NG and IQ-TREE log and result files (default: MSA file name).",
    )

    parser.add_argument(
        "--model",
        type=str,
        required=False,
        help="Model to use for the RAxML-NG tree inference (default: 'GTR+G' for DNA, 'LG+G' for AA, "
             "and 'MULTIx_GTR' for morphological data).",
    )

    parser.add_argument(
        "--ntrees",
        type=int,
        required=False,
        default=100,
        help="Number of ML trees to infer (default: 100)",
    )

    parser.add_argument(
        "--redo",
        action="store_true",
        help="Redo all computations, even if the results already exist.",
    )


    return parser.parse_args()


def main():
    logger.info(get_header())
    args = _parse_cli()

    msa_file = pathlib.Path(args.msa)

    prefix = pathlib.Path(args.prefix) if args.prefix else msa_file

    log_file = pathlib.Path(f"{prefix}.labelGen.log")
    logger.add(log_file, format="{message}")
    log_file.write_text(get_header() + "\n")

    logger.info(
        f"LabelGenerator was called at {time.strftime('%d-%b-%Y %H:%M:%S')} as follows:\n"
    )
    logger.info(" ".join(sys.argv))
    logger.info("")

    msa_file = pathlib.Path(args.msa)
    msa_obj = parse(msa_file)

    if (model := args.model) is None:
        model = msa_obj.get_raxmlng_model()

    raxmlng = pathlib.Path(args.raxmlng)

    log_runtime_information("Starting label computation.")

    # 1. Infer 100 ML trees for the given MSA using RAxML-NG
    log_runtime_information(f"Inferring {args.ntrees} ML trees using RAxML-NG.")
    infer_ml_trees(msa_file, raxmlng, model, prefix, args.ntrees, args.seed, args.threads)

    # 2. RF-Distance
    log_runtime_information("Computing RF-Distance between ML trees.")
    prop_unique_topos, rel_rfdist = rf_distance(prefix, raxmlng)
