import pathlib
from typing import Optional

from pypythia.raxmlng import RAxMLNG, get_raxmlng_rfdist_results, run_raxmlng_command


def _raxmlng_ml_inference_done(prefix: pathlib.Path, n_trees: int) -> bool:
    ml_trees = pathlib.Path(f"{prefix}.raxml.mlTrees")
    best_tree = pathlib.Path(f"{prefix}.raxml.bestTree")
    logfile = pathlib.Path(f"{prefix}.raxml.log")

    # 1. Check if all RAxML-NG files already exist
    files_exist = ml_trees.exists() and best_tree.exists() and logfile.exists()
    if not files_exist:
        return False

    # 2. Number of ML trees is correct
    n_trees_in_file = sum(1 for _ in ml_trees.open())
    n_trees_correct = n_trees_in_file == n_trees

    # 3. Run is complete
    run_complete = "Elapsed time:" in logfile.read_text()

    return files_exist and n_trees_correct and run_complete


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
    if _raxmlng_ml_inference_done(prefix, n_trees) and not redo:
        return

    n_pars_trees = n_trees // 2
    n_rand_trees = n_trees - n_pars_trees

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
        f"pars{{{n_pars_trees}}},rand{{{n_rand_trees}}}",
    ]

    if threads is not None:
        cmd.extend(["--threads", threads])

    if redo:
        cmd.append("--redo")

    run_raxmlng_command(list(map(str, cmd)))


def _raxmlng_rfdist_done(prefix: pathlib.Path) -> bool:
    rfdist = pathlib.Path(f"{prefix}.rfdist")
    logfile = pathlib.Path(f"{prefix}.rfdist.raxml.log")

    # 1. Check if all RAxML-NG files already exist
    files_exist = rfdist.exists() and logfile.exists()
    if not files_exist:
        return False

    # 2. Run is complete
    run_complete = "Elapsed time:" in logfile.read_text()

    return files_exist and run_complete


def rf_distance(
    prefix: pathlib.Path, raxmlng: pathlib.Path, redo: bool = False
) -> tuple[float, float]:
    ml_trees = pathlib.Path(f"{prefix}.raxml.mlTrees")

    if _raxmlng_rfdist_done(prefix) and not redo:
        num_topos, rel_rfdist, _ = get_raxmlng_rfdist_results(
            pathlib.Path(f"{prefix}.rfdist.raxml.log")
        )
    else:
        raxmlng = RAxMLNG(raxmlng)
        num_topos, rel_rfdist, _ = raxmlng.get_rfdistance_results(
            ml_trees, pathlib.Path(f"{prefix}.rfdist"), **{redo: None} if redo else {}
        )

    n_trees = sum(1 for _ in ml_trees.open())
    return num_topos / n_trees, rel_rfdist
