import datetime
import math
import pathlib
import re
import tempfile

import pytest
from pypythia.custom_types import DataType
from pypythia.msa import parse_msa

from labelgenerator.raxmlng import (
    _inference_results_exist_and_correct,
    _rfdist_results_exists_and_correct,
    infer_ml_trees,
    rf_distance,
)


def test_inference_results_exist_and_correct(done_raxml_inference_prefix):
    assert _inference_results_exist_and_correct(done_raxml_inference_prefix, 6)


def test_inference_results_exist_and_correct_nonexisting_files():
    prefix = pathlib.Path("nonexisting")
    assert not _inference_results_exist_and_correct(prefix, 6)


def test_inference_results_exist_and_correct_incorrect_ntrees(
    done_raxml_inference_prefix,
):
    with pytest.raises(ValueError, match="does not match the expected number of trees"):
        _inference_results_exist_and_correct(done_raxml_inference_prefix, 4)


def test_rfdist_results_exists_and_correct(done_raxml_rfdist_prefix):
    assert _rfdist_results_exists_and_correct(done_raxml_rfdist_prefix, 6)


def test_rfdist_results_exists_and_correct_incorrect_ntrees(done_raxml_rfdist_prefix):
    with pytest.raises(ValueError, match="does not match the expected number of pairs"):
        _rfdist_results_exists_and_correct(done_raxml_rfdist_prefix, 2)


@pytest.mark.parametrize("n_trees", [1, 2, 5, 10])
def test_infer_ml_trees(raxmlng_command, dna_msa, n_trees):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        prefix = tmpdir / "test"
        model = "GTR+G"

        infer_ml_trees(
            dna_msa,
            raxmlng_command,
            model,
            prefix,
            n_trees,
            seed=42,
            redo=True,
        )
        assert _inference_results_exist_and_correct(prefix, n_trees)

        # Check the RAxML-NG log file:
        # for even number of trees, the number of parsimony and random trees is equal
        # for odd number of trees, the number of parsimony trees is one more than the random trees
        n_pars_expected = math.ceil(n_trees / 2)
        n_rand_expected = n_trees - n_pars_expected
        expected_random = (
            f"random ({n_rand_expected}) + " if n_rand_expected > 0 else ""
        )
        expected_parsimony = (
            f"parsimony ({n_pars_expected})" if n_pars_expected > 0 else ""
        )
        expected_log = f"start tree(s): {expected_random}{expected_parsimony}"

        log_file = prefix.with_suffix(".raxml.log")
        assert expected_log in log_file.read_text()


@pytest.mark.parametrize("data_type", [DataType.DNA, DataType.AA, DataType.MORPH])
def test_infer_ml_trees_for_dtypes(raxmlng_command, data_dir, data_type):
    msa = data_dir / f"{data_type.name}.phy"
    model = parse_msa(msa).get_raxmlng_model()

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        prefix = tmpdir / "test"
        n_trees = 4

        infer_ml_trees(
            msa,
            raxmlng_command,
            model,
            prefix,
            n_trees,
            42,
            threads=4,
            redo=True,
        )
        assert _inference_results_exist_and_correct(prefix, n_trees)


@pytest.mark.parametrize("data_type", [DataType.DNA, DataType.AA, DataType.MORPH])
def test_infer_ml_trees_existing_files(
    raxmlng_command, ml_tree_dir, data_dir, data_type
):
    test_start = datetime.datetime.now()

    msa = data_dir / f"{data_type.name}.phy"
    model = parse_msa(msa).get_raxmlng_model()
    prefix = ml_tree_dir / f"{data_type.name}"
    n_trees = 20

    # sanity check
    assert _inference_results_exist_and_correct(prefix, n_trees)

    infer_ml_trees(
        msa,
        raxmlng_command,
        model,
        prefix,
        n_trees,
        42,
        threads=4,
        redo=False,
    )
    assert _inference_results_exist_and_correct(prefix, n_trees)

    # The results should have existed already, so the end date in the RAxML-NG log file should be before test_start
    # Logline looks like this: Analysis started: 24-Feb-2025 12:46:35 / finished: 24-Feb-2025 12:46:40
    logfile = pathlib.Path(f"{prefix}.raxml.log")
    logged_raxmlng_end = datetime.datetime.strptime(
        re.search(r"finished: (.+)$", logfile.read_text(), re.MULTILINE).group(1),
        "%d-%b-%Y %H:%M:%S",
    )
    assert logged_raxmlng_end < test_start


def test_infer_ml_trees_fails_for_n_trees_less_than_1(raxmlng_command, dna_msa):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        prefix = tmpdir / "test"
        model = "GTR+G"

        with pytest.raises(ValueError, match="Number of trees needs to be at least 1."):
            infer_ml_trees(
                dna_msa,
                raxmlng_command,
                model,
                prefix,
                0,
                42,
                threads=4,
                redo=True,
            )


def test_rf_distance(raxmlng_command, done_raxml_inference_prefix):
    ml_trees = pathlib.Path(f"{done_raxml_inference_prefix}.raxml.mlTrees")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        prefix = tmpdir / "test"
        num_topos, rfdist = rf_distance(ml_trees, prefix, raxmlng_command)

        assert _rfdist_results_exists_and_correct(
            pathlib.Path(f"{done_raxml_inference_prefix}.rfdist"), 6
        )
        assert num_topos == 2
        assert 0.15 == pytest.approx(rfdist, 0.0)


def test_rf_distance_single_tree(raxmlng_command):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        prefix = tmpdir / "test"

        ml_trees = tmpdir / "single_tree.newick"
        ml_trees.write_text("((A,B),C);")

        num_topos, rfdist = rf_distance(ml_trees, prefix, raxmlng_command)

        assert not _rfdist_results_exists_and_correct(prefix, 1)

        assert num_topos == 1
        assert rfdist == 0.0


def test_rf_distance_existing_results(
    raxmlng_command, done_raxml_inference_prefix, done_raxml_rfdist_prefix
):
    ml_trees = pathlib.Path(f"{done_raxml_inference_prefix}.raxml.mlTrees")

    # make sure the files actually exist
    assert _rfdist_results_exists_and_correct(done_raxml_rfdist_prefix, 6)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        num_topos, rfdist = rf_distance(
            ml_trees, done_raxml_rfdist_prefix, raxmlng_command, n_trees=6, redo=False
        )
        assert num_topos == 2
        assert 0.15 == pytest.approx(rfdist, 0.0)


def test_rf_distance_existing_fails_for_zero_trees(raxmlng_command):
    with pytest.raises(ValueError, match="At least 1 tree is required."):
        # n_trees passed explicitly
        rf_distance(
            pathlib.Path("nonexisting"),
            pathlib.Path("nonexisting"),
            raxmlng_command,
            n_trees=0,
            redo=False,
        )

    with pytest.raises(ValueError, match="At least 1 tree is required."):
        # n_trees not set, but no trees in file
        with tempfile.NamedTemporaryFile(mode="w") as tmpfile:
            rf_distance(
                pathlib.Path(tmpfile.name),
                pathlib.Path("nonexisting"),
                raxmlng_command,
                redo=False,
            )
