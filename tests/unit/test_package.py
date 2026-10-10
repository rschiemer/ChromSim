"""Basic package health: imports, version, and the optional-GUI rule (SPEC §1, §9)."""

import importlib
import importlib.util
import os
import subprocess
import sys

import pytest

import chromsim

SUBPACKAGES = [
    "core",
    "isotherms",
    "backends",
    "data",
    "estimation",
    "optimizers",
    "analysis",
    "io",
    "jobs",
    "storage",
]


def test_version_is_set():
    assert isinstance(chromsim.__version__, str)
    assert chromsim.__version__.count(".") == 2  # SemVer MAJOR.MINOR.PATCH


@pytest.mark.parametrize("name", SUBPACKAGES)
def test_subpackages_import(name):
    importlib.import_module(f"chromsim.{name}")


def test_importing_chromsim_does_not_import_nicegui():
    # Run in a fresh interpreter so modules loaded by other tests or pytest plugins don't count.
    code = (
        "import sys, chromsim, chromsim.api, chromsim.cli; "
        "sys.exit(1 if any(m.split('.')[0] == 'nicegui' for m in sys.modules) else 0)"
    )
    completed = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert completed.returncode == 0, "importing chromsim pulled in nicegui: " + completed.stderr


@pytest.mark.skipif(
    os.environ.get("PIXI_ENVIRONMENT_NAME") != "noui",
    reason="only meaningful in the pixi `noui` environment",
)
def test_noui_environment_really_has_no_nicegui():
    # Guards the no-gui test tier: if NiceGUI leaks into `noui`, that tier proves nothing.
    assert importlib.util.find_spec("nicegui") is None
