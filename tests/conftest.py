import os
import pathlib
import shutil

import pytest


def _command_path(environment_variable, executable):
    configured_command = os.environ.get(environment_variable)
    command = configured_command or shutil.which(executable)
    if command is None:
        pytest.fail(
            f"Could not find {executable!r}. Add it to PATH or set "
            f"{environment_variable}."
        )
    return pathlib.Path(command)


@pytest.fixture
def raxmlng_command():
    return _command_path("RAXMLNG_COMMAND", "raxml-ng")


@pytest.fixture
def iqtree_command():
    return _command_path("IQTREE_COMMAND", "iqtree2")


@pytest.fixture
def data_dir():
    return pathlib.Path.cwd() / "tests" / "data"


@pytest.fixture
def ml_tree_dir(data_dir):
    return data_dir / "mltrees"


@pytest.fixture
def log_dir(data_dir):
    return data_dir / "logs"


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
def done_raxml_inference_prefix(log_dir):
    return log_dir / "test"


@pytest.fixture
def done_raxml_rfdist_prefix(log_dir):
    return log_dir / "test.rfdist"
