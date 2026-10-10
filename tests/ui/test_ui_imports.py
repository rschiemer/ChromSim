"""The UI may import nothing from chromsim except ``chromsim.api`` (SPEC §8, CLAUDE.md rule 1).

This test reads the source with ``ast`` (it never imports the UI), so it also runs where NiceGUI
is not installed.
"""

import ast
from pathlib import Path

import pytest

UI_DIR = Path(__file__).parents[2] / "src" / "chromsim" / "ui"


def forbidden_chromsim_imports(source: str, module_name: str) -> list[str]:
    """Return the chromsim modules imported by ``source`` that a UI module may not use.

    Parameters
    ----------
    source : str
        Python source code of one module.
    module_name : str
        Dotted name of that module, e.g. ``"chromsim.ui.pages.simulate"``; used to resolve
        relative imports.

    Returns
    -------
    list of str
        Offending fully qualified module names (empty if the module follows the rule).
    """
    package_parts = module_name.split(".")[:-1]
    imported = []
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            imported += [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            if node.level == 0:
                base = node.module or ""
            else:
                parent = package_parts[: len(package_parts) - (node.level - 1)]
                base = ".".join(parent + ([node.module] if node.module else []))
            # `from chromsim import api` imports chromsim.api, so check each imported name.
            imported += [f"{base}.{alias.name}" for alias in node.names] + [base]

    def allowed(name: str) -> bool:
        if name.split(".")[0] != "chromsim":
            return True  # third-party and standard library
        return (
            name == "chromsim"  # bare package, the imported names are checked separately
            or name == "chromsim.api"
            or name.startswith("chromsim.api.")
            or name == "chromsim.ui"
            or name.startswith("chromsim.ui.")
        )

    return sorted({name for name in imported if not allowed(name)})


def test_checker_accepts_api_and_ui_internal_imports():
    source = (
        "import nicegui\n"
        "from chromsim import api\n"
        "import chromsim.api\n"
        "from chromsim.api import doctor\n"
        "from . import widgets\n"
        "from .widgets import plot\n"
    )
    assert forbidden_chromsim_imports(source, "chromsim.ui.pages.simulate") == []


@pytest.mark.parametrize(
    "source",
    [
        "from chromsim.core import column\n",
        "import chromsim.backends.cadet_backend\n",
        "from chromsim import units\n",
        "from ... import diagnostics\n",  # relative import escaping the ui package
    ],
)
def test_checker_rejects_science_imports(source):
    assert forbidden_chromsim_imports(source, "chromsim.ui.pages.simulate") != []


def test_ui_package_only_imports_chromsim_api():
    violations = {}
    for path in sorted(UI_DIR.rglob("*.py")):
        relative = path.relative_to(UI_DIR.parents[1]).with_suffix("")
        module_name = ".".join(relative.parts)
        bad = forbidden_chromsim_imports(path.read_text(encoding="utf-8"), module_name)
        if bad:
            violations[str(path)] = bad
    assert not violations, f"UI modules import chromsim internals: {violations}"
