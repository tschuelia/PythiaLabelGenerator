import pytest
import tempfile
import pathlib

from labelgenerator.iqtree import get_iqtree_model, run_statstests, _iqtree_results_exist_and_done
from pypythia.custom_types import DataType


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
    model = get_iqtree_model(data_type)
    msa = data_dir / f"{data_type}.phy"
    ml_trees = ml_tree_dir / f"{data_type}.raxml.mlTrees"
    best_tree = data_dir / f"{data_type}.raxml.bestTree"

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        prefix = tmpdir / "test"
        run_statstests(
            msa,
            ml_trees,
            best_tree,
            iqtree_command,
            model,
            prefix,
            42,
            threads=4,
            redo=True,
        )
        assert _iqtree_results_exist_and_done(prefix)