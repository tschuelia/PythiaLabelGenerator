import pathlib
import tempfile

import pytest

from labelgenerator.label import get_label, compute_label
from pypythia.custom_types import DataType


def test_get_label_fails_for_invalid_input():
    with pytest.raises(ValueError, match="Number of unique trees cannot be higher than the total number of trees."):
        get_label(n_all=10, rf_all=0.5, n_unique_all=20, n_plausible=5, rf_plausible=0.5, n_unique_plausible=5)

    with pytest.raises(ValueError, match="Number of unique plausible trees cannot be higher than the number of plausible trees."):
        get_label(n_all=10, rf_all=0.5, n_unique_all=5, n_plausible=5, rf_plausible=0.5, n_unique_plausible=10)

    with pytest.raises(ValueError, match="Number of plausible trees cannot be higher than the total number of trees."):
        get_label(n_all=5, rf_all=0.5, n_unique_all=5, n_plausible=10, rf_plausible=0.5, n_unique_plausible=5)

    with pytest.raises(ValueError, match="RF distance for all trees must be between 0 and 1."):
        get_label(n_all=5, rf_all=1.5, n_unique_all=5, n_plausible=5, rf_plausible=0.5, n_unique_plausible=5)

    with pytest.raises(ValueError, match="RF distance for plausible trees must be between 0 and 1."):
        get_label(n_all=5, rf_all=0.5, n_unique_all=5, n_plausible=5, rf_plausible=1.5, n_unique_plausible=5)


@pytest.mark.parametrize("values, expected", [
    ([100, 0, 1, 100, 0, 1], 0.004),
    ([100, 1, 100, 100, 1, 100], 0.8),
    ([100, 1, 100, 50, 1, 50], 0.9)
])
def test_get_label(values, expected):
    computed = get_label(*values)
    assert computed == pytest.approx(expected)


@pytest.mark.parametrize("data_type, expected_label", [(DataType.DNA, 0.729), (DataType.AA, 0.04), (DataType.MORPH, 0.214)])
def test_compute_label(raxmlng_command, iqtree_command, data_dir, data_type, expected_label):
    msa_file = data_dir / f"{data_type.name}.phy"

    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = pathlib.Path(tmpdir) / "test"

        label = compute_label(
            msa_file=msa_file,
            raxmlng=raxmlng_command,
            iqtree=iqtree_command,
            prefix=prefix,
            n_trees=10,
            threads=4,
            log_info=False
        )
        assert label == pytest.approx(expected_label, abs=0.01)