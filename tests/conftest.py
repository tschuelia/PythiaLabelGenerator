import pathlib

import pytest

from .test_config import RAXMLNG_COMMAND


@pytest.fixture
def raxmlng_command():
    return pathlib.Path(RAXMLNG_COMMAND)


@pytest.fixture
def phylip_msa_file():
    return pathlib.Path.cwd() / "tests" / "data" / "DNA" / "1.phy"


@pytest.fixture
def data_dir():
    return pathlib.Path.cwd() / "tests" / "data"


@pytest.fixture
def ml_tree_dir(data_dir):
    return data_dir / "mltrees"


@pytest.fixture
def dna_msa(data_dir):
    return data_dir / "DNA.phy"


@pytest.fixture
def aa_msa(data_dir):
    return data_dir / "AA.phy"


@pytest.fixture
def morph_msa(data_dir):
    return data_dir / "MORPH.phy"


@pytest.fixture
def done_raxml_inference_prefix(data_dir):
    return data_dir / "logs" / "test"


@pytest.fixture
def done_raxml_rfdist_prefix(data_dir):
    return data_dir / "logs" / "test.rfdist"
