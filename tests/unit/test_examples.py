"""Every example script in ``examples/`` runs without error (SPEC §8 code-only parity, §11)."""

import runpy
from pathlib import Path

import pytest

EXAMPLES_DIR = Path(__file__).parents[2] / "examples"
EXAMPLE_SCRIPTS = sorted(EXAMPLES_DIR.glob("*.py"))


def test_examples_exist():
    assert EXAMPLE_SCRIPTS, f"no example scripts found in {EXAMPLES_DIR}"


@pytest.mark.slow
@pytest.mark.parametrize("script", EXAMPLE_SCRIPTS, ids=lambda path: path.name)
def test_example_runs(script):
    runpy.run_path(str(script), run_name="__main__")
