import argparse
import pathlib
import shutil
import sys
import time

from labelgenerator import __version__
from labelgenerator.label import compute_label
from labelgenerator.logger import (
    SCRIPT_START,
    get_header,
    log_runtime_information,
    logger,
)

DEFAULT_RAXMLNG_EXE = (
    pathlib.Path(shutil.which("raxml-ng")) if shutil.which("raxml-ng") else None
)

DEFAULT_IQTREE_EXE = (
    pathlib.Path(shutil.which("iqtree2")) if shutil.which("iqtree2") else None
)


def _parse_cli():
    parser = argparse.ArgumentParser(
        description="Generate the ground truth difficulty for the given MSA."
    )
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
        help="Number of threads to use for the RAxML-NG tree inference and IQ-TREE statistical tests (default: autoconfig in RAxML-NG and IQ-TREE).",
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

    parser.add_argument(
        "-V",
        "--version",
        action="version",
        version=__version__,
        help="Print the version number and exit.",
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

    log_runtime_information("Starting label computation.")

    difficulty = compute_label(
        msa_file=msa_file,
        raxmlng=pathlib.Path(args.raxmlng),
        iqtree=pathlib.Path(args.iqtree),
        prefix=prefix,
        model=args.model,
        n_trees=args.ntrees,
        seed=args.seed,
        threads=args.threads,
        redo=args.redo,
        log_info=True,
    )

    script_end = time.perf_counter()

    logger.info("")
    logger.info(f"Ground Truth Difficulty for {msa_file}: {difficulty:.3f}")

    if args.ntrees < 100:
        logger.info(
            "WARNING: The number of inferred ML trees is less than 100. The computed label may be less reliable."
        )

    logger.info("")
    total_runtime = script_end - SCRIPT_START
    hours, remainder = divmod(total_runtime, 3600)
    minutes, seconds = divmod(remainder, 60)

    if hours > 0:
        logger.info(
            f"Total runtime: {int(hours):02d}:{int(minutes):02d}:{seconds:02d} hours ({round(total_runtime)} seconds)."
        )
    elif minutes > 0:
        logger.info(
            f"Total runtime: {int(minutes):02d}:{int(seconds):02d} minutes ({round(total_runtime)} seconds)."
        )
    else:
        logger.info(f"Total runtime: {seconds:.2f} seconds.")


if __name__ == "__main__":
    main()
