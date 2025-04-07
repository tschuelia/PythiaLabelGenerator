import tempfile
import pytest

import pandas as pd

from labelgenerator.main import main
from labelgenerator import __version__
import pathlib

from pypythia.custom_types import DataType
from pypythia.msa import parse_msa


@pytest.mark.parametrize(
    "data_type, expected_label",
    [(DataType.DNA, 0.772), (DataType.AA, 0.04), (DataType.MORPH, 0.173)],
)
def test_main(data_type, expected_label, data_dir, raxmlng_command, iqtree_command):
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a temporary directory for the test
        prefix = pathlib.Path(tmpdir) / "test"
        seed = 42
        n_trees = 10

        logfile = prefix.with_suffix(".labelGen.log")
        features_file = prefix.with_suffix(".csv")

        msa_file = data_dir / f"{data_type.name}.phy"
        msa_obj = parse_msa(msa_file)

        args = [
            "--msa",
            str(msa_file),
            "--raxmlng",
            str(raxmlng_command),
            "--iqtree",
            str(iqtree_command),
            "--seed",
            str(seed),
            "--prefix",
            str(prefix),
            "--ntrees",
            str(n_trees),
            "--redo",
        ]

        main(args)

        # Check if all output files exist
        expected_files = [
            logfile,
            features_file,
            # RAxML-NG output files
            prefix.with_suffix(".raxml.bestTree"),
            prefix.with_suffix(".raxml.log"),
            prefix.with_suffix(".raxml.mlTrees"),
            prefix.with_suffix(".raxml.plausibleTrees"),
            # IQ-TREE output files
            prefix.with_suffix(".iqtree.iqtree"),
            prefix.with_suffix(".iqtree.log"),
        ]

        for file in expected_files:
            assert file.exists(), f"Expected file {file} does not exist."

        # Check if the features csv file is correct
        features_content = pd.read_csv(features_file)
        assert features_content.shape[0] == 1, "Expected one row in the features file."
        pd.testing.assert_index_equal(
            features_content.columns,
            pd.Index(
                [
                    "num_taxa",
                    "num_sites",
                    "num_patterns",
                    "num_patterns/num_taxa",
                    "num_sites/num_taxa",
                    "num_patterns/num_sites",
                    "proportion_gaps",
                    "proportion_invariant",
                    "entropy",
                    "bollback",
                    "pattern_entropy",
                    "avg_rfdist_parsimony",
                    "proportion_unique_topos_parsimony",
                    "difficulty",
                ]
            ),
        )

        # Check if the label is correct
        label = features_content["difficulty"].values[0]
        assert label == pytest.approx(expected_label, abs=0.01)

        # Check if the log file is correct and contains the expected output
        expected_lines = [
            f"PyDLG version {__version__} released by The Exelixis Lab",
            "PyDLG was called at",
            "Starting label computation.",
            "Computing RF-Distance between ML trees.",
            "Inferring 10 ML trees using RAxML-NG.",
            "RF-Distance ML trees:",
            "Unique topologies ML trees:",
            "Running IQ-TREE statistical tests.",
            "Computing RF-Distance between plausible ML trees.",
            "RF-Distance plausible trees:",
            "Unique topologies plausible trees:",
            "Computing corresponding Pythia features.",
            f"Ground Truth Difficulty for {msa_file}:",
            "WARNING: The number of inferred ML trees is less than 100. The computed label may be less reliable.",
            "Label computation runtime",
            "Total runtime",
        ]

        if msa_obj.contains_duplicate_sequences():
            expected_lines.append(
                "WARNING: The input MSA contains duplicate sequences. Consider deduplicating the MSA before running the label generator."
            )
        if msa_obj.contains_full_gap_sequences():
            expected_lines.append(
                "WARNING: The input MSA contains sequences that are only gaps. Consider removing these sequences before running the label generator."
            )

        log_content = logfile.read_text()

        for line in expected_lines:
            assert line in log_content
