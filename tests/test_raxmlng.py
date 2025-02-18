import math
import pathlib
import tempfile

import pytest

from labelgenerator.raxmlng import infer_ml_trees, rf_distance, _check_existing_inference_results, _raxmlng_rfdist_done


def test_check_existing_inference_results(done_raxml_inference_prefix):
    assert _check_existing_inference_results(done_raxml_inference_prefix, 6)


def test_check_existing_inference_results_incorrect_ntrees(done_raxml_inference_prefix):
    with pytest.raises(ValueError, match="does not match the expected number of trees"):
        _check_existing_inference_results(done_raxml_inference_prefix, 4)


def test_raxmlng_rfdist_done(done_raxml_rfdist_prefix):
    assert _raxmlng_rfdist_done(done_raxml_rfdist_prefix)


@pytest.mark.parametrize("n_trees", [1, 2, 5, 10])
def test_infer_ml_trees(raxmlng_command, phylip_msa_file, n_trees):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        prefix = tmpdir / "test"
        model = "GTR+G"

        infer_ml_trees(phylip_msa_file, raxmlng_command, model, prefix, n_trees, 42, threads=4, redo=True)
        assert _check_existing_inference_results(prefix, n_trees)

        # Check the RAxML-NG log file:
        # for even number of trees, the number of parsimony and random trees is equal
        # for odd number of trees, the number of parsimony trees is one more than the random trees
        n_pars_expected = math.ceil(n_trees / 2)
        n_rand_expected = n_trees - n_pars_expected
        expected_random = f"random ({n_rand_expected}) + " if n_rand_expected > 0 else ""
        expected_parsimony = f"parsimony ({n_pars_expected})" if n_pars_expected > 0 else ""
        expected_log = f"start tree(s): {expected_random}{expected_parsimony}"

        log_file = prefix.with_suffix(".raxml.log")
        assert expected_log in log_file.read_text()


def test_infer_ml_trees_fails_for_n_trees_less_than_1(raxmlng_command, phylip_msa_file):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        prefix = tmpdir / "test"
        model = "GTR+G"

        with pytest.raises(ValueError, match="Number of trees needs to be at least 1."):
            infer_ml_trees(phylip_msa_file, raxmlng_command, model, prefix, 0, 42, threads=4, redo=True)


def test_rf_distance(raxmlng_command, done_raxml_inference_prefix):
    ml_trees = pathlib.Path(f"{done_raxml_inference_prefix}.raxml.mlTrees")
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = pathlib.Path(tmpdir)
        prefix = tmpdir / "test"
        proportion_unique, rfdist = rf_distance(ml_trees, prefix, raxmlng_command)

        assert _raxmlng_rfdist_done(pathlib.Path(f"{done_raxml_inference_prefix}.rfdist"))
        assert 1/3 == pytest.approx(proportion_unique, 0.0)
        assert 0.15 == pytest.approx(rfdist, 0.0)


