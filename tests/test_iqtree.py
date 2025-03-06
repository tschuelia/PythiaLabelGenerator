import datetime
import pathlib
import re
import tempfile

import pytest
from pypythia.custom_types import DataType

from labelgenerator.iqtree import (
    _iqtree_results_exist_and_done,
    filter_plausible_trees,
    get_iqtree_model,
    run_statstests,
)


@pytest.mark.parametrize("data_type", [DataType.DNA, DataType.AA, DataType.MORPH])
def test_iqtree_results_exist_and_done(log_dir, data_type):
    prefix = log_dir / f"{data_type.name}.iqtree"
    assert _iqtree_results_exist_and_done(prefix)


def test_iqtree_results_exist_and_done_non_existing_results(log_dir):
    prefix = log_dir / "nonexisting.iqtree"
    assert not _iqtree_results_exist_and_done(prefix)


@pytest.mark.parametrize("data_type", [DataType.DNA, DataType.AA, DataType.MORPH])
def test_get_iqtree_model(data_type):
    model = get_iqtree_model(data_type)

    if data_type == DataType.DNA:
        assert model == "GTR+G4+FO"
    elif data_type == DataType.AA:
        assert model == "LG+G4+FO"
    elif data_type == DataType.MORPH:
        assert model == "MK"
    else:
        raise ValueError(f"Unknown data type: {data_type}")


@pytest.mark.parametrize("data_type", [DataType.DNA, DataType.AA, DataType.MORPH])
def test_run_statstests(iqtree_command, data_dir, ml_tree_dir, data_type):
    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = pathlib.Path(tmpdir) / "test"
        run_statstests(
            msa=data_dir / f"{data_type.name}.phy",
            ml_trees=ml_tree_dir / f"{data_type.name}.raxml.mlTrees",
            best_tree=ml_tree_dir / f"{data_type.name}.raxml.bestTree",
            iqtree=iqtree_command,
            model=get_iqtree_model(data_type),
            prefix=prefix,
            seed=42,
            is_morph=data_type == DataType.MORPH,
            redo=True,
        )
        assert _iqtree_results_exist_and_done(prefix)


@pytest.mark.parametrize("data_type", [DataType.DNA, DataType.AA, DataType.MORPH])
def test_run_statstests_existing_files(
    iqtree_command, log_dir, data_dir, ml_tree_dir, data_type
):
    test_start = datetime.datetime.now()
    prefix = log_dir / f"{data_type.name}.iqtree"

    # sanity check
    assert _iqtree_results_exist_and_done(prefix)

    run_statstests(
        msa=data_dir / f"{data_type.name}.phy",
        ml_trees=ml_tree_dir / f"{data_type.name}.raxml.mlTrees",
        best_tree=ml_tree_dir / f"{data_type.name}.raxml.bestTree",
        iqtree=iqtree_command,
        model=get_iqtree_model(data_type),
        prefix=prefix,
        seed=42,
        is_morph=data_type == DataType.MORPH,
        redo=False,
    )
    assert _iqtree_results_exist_and_done(prefix)
    # The results should have existed already, so the logged date in the IQ-Tree log file should be before test_start
    # Logline looks like this: Date and Time: Tue Feb 25 07:56:49 2025
    logfile = pathlib.Path(f"{prefix}.log")
    logged_iqtree_end = datetime.datetime.strptime(
        re.search(r"Date and Time: (.+)", (logfile.read_text())).group(1),
        "%a %b %d %H:%M:%S %Y",
    )
    assert logged_iqtree_end < test_start


def test_run_statstests_fails():
    with pytest.raises(RuntimeError, match="Running IQ-TREE command failed"):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = pathlib.Path(tmpdir) / "test"
            run_statstests(
                msa=pathlib.Path("nonexistent_file.phy"),
                ml_trees=pathlib.Path("nonexistent_file.mlTrees"),
                best_tree=pathlib.Path("nonexistent_file.bestTree"),
                iqtree=pathlib.Path("nonexistent_command"),
                model="GTR+G4+FO",
                prefix=prefix,
                seed=42,
                redo=True,
            )


def test_run_statstests_fails_with_CalledProcessError(ml_tree_dir, iqtree_command):
    with pytest.raises(
        RuntimeError,
        match=re.compile(
            r"Running IQ-TREE command failed:.*?Reading alignment file", re.DOTALL
        ),
    ):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = pathlib.Path(tmpdir) / "test"
            run_statstests(
                msa=pathlib.Path("nonexistent_file.phy"),
                ml_trees=ml_tree_dir / "DNA.raxml.mlTrees",
                best_tree=ml_tree_dir / "DNA.raxml.bestTree",
                iqtree=iqtree_command,
                model="GTR+G4+FO",
                prefix=prefix,
                seed=42,
                redo=True,
            )


@pytest.mark.parametrize(
    "data_type, plausible_tree_indices",
    [
        (DataType.DNA, [2, 6, 8, 9, 10, 11, 12, 13, 14, 18, 19]),
        (DataType.AA, list(range(20))),
        (DataType.MORPH, list(range(20))),
    ],
)
def test_filter_plausible_trees_filters_correct_trees(
    ml_tree_dir, log_dir, data_type, plausible_tree_indices
):
    ml_trees = ml_tree_dir / f"{data_type.name}.raxml.mlTrees"
    iqtree_results = log_dir / f"{data_type.name}.iqtree.iqtree"

    with tempfile.NamedTemporaryFile(mode="w") as tmpfile:
        filter_plausible_trees(
            ml_trees=ml_trees,
            iqtree_results=iqtree_results,
            plausible_ml_trees=pathlib.Path(tmpfile.name),
        )

        ml_trees = ml_trees.read_text().strip().splitlines()
        expected_plausible_ml_trees = [ml_trees[i] for i in plausible_tree_indices]
        actual_plausible_ml_trees = (
            pathlib.Path(tmpfile.name).read_text().strip().splitlines()
        )

        assert len(expected_plausible_ml_trees) == len(actual_plausible_ml_trees)
        assert expected_plausible_ml_trees == actual_plausible_ml_trees


def test_filter_plausible_trees_fails_for_incorrect_n_trees(log_dir):
    iqtree_results = log_dir / "DNA.iqtree.iqtree"

    with pytest.raises(
        ValueError,
        match="Number of IQ-TREE results does not match the number of ML trees.",
    ):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = pathlib.Path(tmpdir)
            ml_trees = tmpdir / "incorrect.mlTrees"
            ml_trees.write_text("A,B,(C,D));")
            filter_plausible_trees(
                ml_trees=ml_trees,
                iqtree_results=iqtree_results,
                plausible_ml_trees=tmpdir / "plausible.mlTrees",
            )
