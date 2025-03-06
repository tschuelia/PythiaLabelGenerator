import math
import pathlib
from typing import Optional

from pypythia.raxmlng import RAxMLNG, get_raxmlng_rfdist_results, run_raxmlng_command


def _inference_results_exist_and_correct(prefix: pathlib.Path, n_trees: int) -> bool:
    ml_trees = pathlib.Path(f"{prefix}.raxml.mlTrees")
    best_tree = pathlib.Path(f"{prefix}.raxml.bestTree")

    if n_trees == 1:
        ml_trees = best_tree

    logfile = pathlib.Path(f"{prefix}.raxml.log")

    files_exist = ml_trees.exists() and best_tree.exists() and logfile.exists()
    if not files_exist:
        # Files don't exist yet, nothing to check.
        return False

    # Check if the number of ML trees is correct, if there is a mismatch, raise an error
    n_trees_in_file = sum(1 for _ in ml_trees.open())
    if n_trees_in_file != n_trees:
        raise ValueError(
            f"Number of trees in {ml_trees} ({n_trees_in_file}) does not match the expected number of trees ({n_trees})."
            f"Please set the `redo` flag to recompute the trees."
        )

    # Finally, check if the previous RAxML-NG run completed successfully
    return "Elapsed time:" in logfile.read_text()


def infer_ml_trees(
    msa: pathlib.Path,
    raxmlng: pathlib.Path,
    model: str,
    prefix: pathlib.Path,
    n_trees: int = 100,
    seed: int = 0,
    threads: Optional[int] = None,
    redo: bool = False,
) -> None:
    """
    Infers ML trees using RAxML-NG.

    If the results for the given prefix and number of trees already exist, the function will return without doing anything.
    If the number of trees in the existing file does not match the expected number of trees, a ValueError is raised.
    If you want to redo the computation, set redo=True.

    Args:
        msa (pathlib.Path): Path to the MSA file.
        raxmlng (pathlib.Path): Path to the RAxML-NG executable.
        model (str): Model to use for the RAxML-NG tree inference.
        prefix (pathlib.Path): Prefix to use for the RAxML-NG output files.
        n_trees (int): Number of ML trees to infer.
        seed (int): Seed to use for the RAxML-NG inference.
        threads (Optional[int]): Number of threads to use for the RAxML-NG inference.
            Per default, uses the automatic setting of RAxML-NG which is likely to use all cores of your machine.
        redo (bool): Flag to redo the computation even if the results already exist.

    Returns:
        None

    Raises:
        ValueError:
        - If the number of trees is less than 1.
        - If the number of trees in a previously computed file does not match the expected number of trees.

    """
    if n_trees < 1:
        raise ValueError("Number of trees needs to be at least 1.")

    if not redo and _inference_results_exist_and_correct(prefix, n_trees):
        return

    n_pars_trees = math.ceil(n_trees / 2)
    n_rand_trees = n_trees - n_pars_trees
    rand_string = f",rand{{{n_rand_trees}}}" if n_rand_trees > 0 else ""

    cmd = [
        raxmlng,
        "--msa",
        msa,
        "--model",
        model,
        "--seed",
        seed,
        "--prefix",
        prefix,
        "--tree",
        f"pars{{{n_pars_trees}}}{rand_string}",
    ]

    if threads is not None:
        cmd.extend(["--threads", threads])

    if redo:
        cmd.append("--redo")

    run_raxmlng_command(list(map(str, cmd)))


def _rfdist_results_exists_and_correct(prefix: pathlib.Path, n_trees: int) -> bool:
    rfdist = pathlib.Path(f"{prefix}.raxml.rfDistances")
    logfile = pathlib.Path(f"{prefix}.raxml.log")

    # 1. Check if all RAxML-NG files already exist
    files_exist = rfdist.exists() and logfile.exists()
    if not files_exist:
        return False

    # 2. Run is complete
    run_complete = "Elapsed time:" in logfile.read_text()

    # 3. Check if the number of pairwise RF-Distance results is correct
    expected_number_of_pairs = n_trees * (n_trees - 1) // 2
    n_pairs_in_file = sum(1 for _ in rfdist.open())
    if n_pairs_in_file != expected_number_of_pairs:
        raise ValueError(
            f"Number of pairwise RF-Distances in {rfdist} ({n_pairs_in_file}) does not match "
            f"the expected number of pairs ({expected_number_of_pairs})."
            f"Please set the `redo` flag to recompute the distances."
        )

    return files_exist and run_complete


def rf_distance(
    ml_trees: pathlib.Path,
    prefix: pathlib.Path,
    raxmlng: pathlib.Path,
    n_trees: Optional[int] = None,
    redo: bool = False,
) -> tuple[int, float]:
    """
    Compute the number of unique topologies and the average relative RF distance for a set of ML trees.
    If the results already exist, the function will return the results without recomputing them.

    Args:
        ml_trees (pathlib.Path): Path to the file containing the ML trees.
        prefix (pathlib.Path): Prefix to use for the RAxML-NG output files.
        raxmlng (pathlib.Path): Path to the RAxML-NG executable.
        n_trees (Optional[int]): Number of trees that were inferred.
            If not provided, the number of trees will be inferred from the file.
            Explicitly provide the number of trees if you want to check if existing results for the given
            prefix contain the results for the correct number of trees.
        redo (bool): Flag to redo the computation if the results already exist.

    Returns:
        tuple[int, float]: The number of unique topologies and the average relative RF distance.

    """
    if n_trees is None:
        n_trees = sum(1 for _ in ml_trees.open())

    if n_trees == 0:
        raise ValueError("At least 1 tree is required.")

    if n_trees == 1:
        return 1, 0.0

    if not redo and _rfdist_results_exists_and_correct(prefix, n_trees):
        num_topos, rel_rfdist = get_raxmlng_rfdist_results(
            pathlib.Path(f"{prefix}.raxml.log")
        )
    else:
        raxmlng = RAxMLNG(raxmlng)
        num_topos, rel_rfdist = raxmlng.get_rfdistance_results(
            ml_trees, prefix, **{"redo": None} if redo else {}
        )

    return num_topos, rel_rfdist
