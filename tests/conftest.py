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
def done_raxml_inference_prefix():
    return pathlib.Path.cwd() / "tests" / "data" / "logs" / "test"


@pytest.fixture
def done_raxml_rfdist_prefix():
    return pathlib.Path.cwd() / "tests" / "data" / "logs" / "test.rfdist"
