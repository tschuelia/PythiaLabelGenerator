import pathlib
from typing import Optional

from pypythia.msa import parse_msa

from labelgenerator.iqtree import (
    filter_plausible_trees,
    get_iqtree_model,
    run_statstests,
)
from labelgenerator.logger import log_runtime_information
from labelgenerator.raxmlng import infer_ml_trees, rf_distance


def get_label(
    n_all: int,
    rf_all: float,
    n_unique_all: int,
    n_plausible: int,
    rf_plausible: float,
    n_unique_plausible: int,
) -> float:
    if n_unique_all > n_all:
        raise ValueError(
            "Number of unique trees cannot be higher than the total number of trees."
        )
    if n_unique_plausible > n_plausible:
        raise ValueError(
            "Number of unique plausible trees cannot be higher than the number of plausible trees."
        )
    if n_plausible > n_all:
        raise ValueError(
            "Number of plausible trees cannot be higher than the total number of trees."
        )
    if not 0 <= rf_all <= 1:
        raise ValueError("RF distance for all trees must be between 0 and 1.")
    if not 0 <= rf_plausible <= 1:
        raise ValueError("RF distance for plausible trees must be between 0 and 1.")

    proportion_unique_all = n_unique_all / n_all
    proportion_unique_plausible = n_unique_plausible / n_plausible
    proportion_plausible = n_plausible / n_all

    total = (
        rf_all
        + proportion_unique_all
        + rf_plausible
        + proportion_unique_plausible
        + (1 - proportion_plausible)
    )
    label = total / 5

    eps = 1e-9
    assert -eps <= label <= 1 + eps, (
        f"Label {label} is not between 0 and 1. Check the input values."
    )

    return label


def compute_label(
    msa_file: pathlib.Path,
    raxmlng: pathlib.Path,
    iqtree: pathlib.Path,
    prefix: pathlib.Path,
    model: Optional[str] = None,
    n_trees: int = 100,
    seed: int = 0,
    threads: Optional[int] = None,
    redo: bool = False,
    log_info: bool = True,
):
    msa_obj = parse_msa(msa_file)
    model = model or msa_obj.get_raxmlng_model()

    # 1. Infer 100 ML trees for the given MSA using RAxML-NG
    if log_info:
        log_runtime_information(f"Inferring {n_trees} ML trees using RAxML-NG.")

    infer_ml_trees(
        msa=msa_file,
        raxmlng=raxmlng,
        model=model,
        prefix=prefix,
        n_trees=n_trees,
        seed=seed,
        threads=threads,
        redo=redo,
    )

    # 2. RF-Distance ML trees
    ml_trees = pathlib.Path(f"{prefix}.raxml.mlTrees")
    rfdistance_prefix = pathlib.Path(f"{prefix}.rfdist")

    if log_info:
        log_runtime_information("Computing RF-Distance between ML trees.")

    n_unique_all, rf_all = rf_distance(
        ml_trees=ml_trees, prefix=rfdistance_prefix, raxmlng=raxmlng, redo=redo
    )

    if log_info:
        log_runtime_information(f"> RF-Distance ML trees: {round(rf_all, 2)}")
        log_runtime_information(f"> Unique topologies ML trees: {n_unique_all}")

    # 3. IQ-TREE statistical tests
    best_tree = pathlib.Path(f"{prefix}.raxml.bestTree")
    iqtree_prefix = pathlib.Path(f"{prefix}.iqtree")
    if log_info:
        log_runtime_information("Running IQ-TREE statistical tests.")
    run_statstests(
        msa=msa_file,
        ml_trees=ml_trees,
        best_tree=best_tree,
        iqtree=iqtree,
        model=get_iqtree_model(msa_obj.data_type),
        prefix=iqtree_prefix,
        seed=seed,
        threads=threads,
        redo=redo,
    )

    # 4. Filter the plausible trees
    iqtree_results = pathlib.Path(f"{iqtree_prefix}.iqtree")
    plausible_ml_trees = pathlib.Path(f"{prefix}.raxml.plausibleTrees")
    if log_info:
        log_runtime_information("Filtering plausible ML trees.")
    filter_plausible_trees(
        ml_trees=ml_trees,
        iqtree_results=iqtree_results,
        plausible_ml_trees=plausible_ml_trees,
    )

    n_plausible_trees = sum(1 for _ in plausible_ml_trees.open())

    if log_info:
        log_runtime_information(f"> Found {n_plausible_trees} plausible trees.")

    # 5. RF-Distance plausible trees
    rfdistance_plausible_prefix = pathlib.Path(f"{prefix}.rfdist.plausible")
    if log_info:
        log_runtime_information("Computing RF-Distance between plausible ML trees.")
    n_unique_plausible, rf_plausible = rf_distance(
        ml_trees=plausible_ml_trees,
        prefix=rfdistance_plausible_prefix,
        raxmlng=raxmlng,
        redo=redo,
    )

    if log_info:
        log_runtime_information(
            f"> RF-Distance plausible trees: {round(rf_plausible, 2)}"
        )
        log_runtime_information(
            f"> Unique topologies plausible trees: {n_unique_plausible}"
        )

    # 6. Compute the ground truth difficulty
    difficulty = get_label(
        n_all=n_trees,
        rf_all=rf_all,
        n_unique_all=n_unique_all,
        n_plausible=n_plausible_trees,
        rf_plausible=rf_plausible,
        n_unique_plausible=n_unique_plausible,
    )
    return difficulty
