import math
import pathlib
import tempfile

import pytest
from pypythia.custom_types import DataType
from pypythia.msa import parse

from labelgenerator.raxmlng import (
    _inference_results_exist_and_correct,
    _rfdist_results_exists_and_correct,
    infer_ml_trees,
    rf_distance,
)


def test_inference_results_exist_and_correct(done_raxml_inference_prefix):
    assert _inference_results_exist_and_correct(done_raxml_inference_prefix, 6)


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
            42,
            threads=4,
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
    model = parse(msa).get_raxmlng_model()

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
