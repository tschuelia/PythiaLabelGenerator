import pathlib
import re
import subprocess
import tempfile

import pytest
from pypythia.custom_types import DataType

from labelgenerator.iqtree import (
    _iqtree_results_exist_and_done,
    get_iqtree_model,
    run_statstests,
)


@pytest.mark.parametrize("data_type", [DataType.DNA, DataType.AA, DataType.MORPH])
def test_iqtree_results_exist_and_done(log_dir, data_type):
    prefix = log_dir / f"{data_type.name}.iqtree"
    assert _iqtree_results_exist_and_done(prefix)


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
            threads=4,
            is_morph=data_type == DataType.MORPH,
            redo=True,
        )
        assert _iqtree_results_exist_and_done(prefix)


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
                threads=4,
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
                threads=4,
                redo=True,
            )
