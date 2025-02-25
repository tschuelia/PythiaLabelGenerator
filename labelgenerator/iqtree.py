import pathlib
import subprocess
from typing import Optional

from pypythia.custom_types import DataType

from labelgenerator.iqtree_parser import get_iqtree_results


def _iqtree_results_exist_and_done(prefix: pathlib.Path) -> bool:
    iqtree_file = pathlib.Path(f"{prefix}.iqtree")
    logfile = pathlib.Path(f"{prefix}.log")

    if not iqtree_file.exists() or not logfile.exists():
        return False

    return (
        "Total wall-clock time used" in iqtree_file.read_text()
        and "Date and Time:" in logfile.read_text()
    )


def get_iqtree_model(data_type: DataType) -> str:
    return {
        DataType.DNA: "GTR+G4+FO",
        DataType.AA: "LG+G4+FO",
        DataType.MORPH: "MK",
    }[data_type]


def run_statstests(
    msa: pathlib.Path,
    ml_trees: pathlib.Path,
    best_tree: pathlib.Path,
    iqtree: pathlib.Path,
    model: str,
    prefix: pathlib.Path,
    seed: int = 0,
    threads: Optional[int] = None,
    is_morph: bool = False,
    redo: bool = False,
) -> None:
    if not redo and _iqtree_results_exist_and_done(prefix):
        return

    cmd = [
        iqtree,
        "-s",
        msa,
        "-m",
        model,
        "-pre",
        prefix,
        "-z",
        ml_trees,
        "-treediff",
        "-te",
        best_tree,
        "-n",
        0,
        "-zb",
        10000,
        "-zw",
        "-au",
        "-seed",
        seed,
    ]

    if is_morph:
        cmd.extend(["-st", "MORPH"])

    if threads is not None:
        cmd.extend(["-nt", threads])

    if redo:
        cmd.append("-redo")

    try:
        subprocess.check_output(list(map(str, cmd)), encoding="utf-8")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Running IQ-TREE command failed: {e.stdout}")
    except Exception as e:
        raise RuntimeError("Running IQ-TREE command failed.") from e


def filter_plausible_trees(
    ml_trees: pathlib.Path,
    iqtree_results: pathlib.Path,
    plausible_ml_trees: pathlib.Path,
):
    iqtree_results = get_iqtree_results(iqtree_results)
    newick_trees = [t.strip() for t in ml_trees.read_text().splitlines()]

    if not len(iqtree_results) == len(newick_trees):
        raise ValueError(
            "Number of IQ-TREE results does not match the number of ML trees."
        )

    plausible_indices = [
        i
        for i, x in enumerate(iqtree_results)
        if x["plausible"] and not x["duplicate_topology"]
    ]
    plausible_trees = [newick_trees[i] for i in plausible_indices]

    plausible_ml_trees.write_text("\n".join(plausible_trees))
