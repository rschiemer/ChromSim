"""No module-level mutable state anywhere in the package (SPEC §8.1, CLAUDE.md rule 1).

Read with ``ast`` (nothing is imported). A module-level name bound to a mutable container is
state shared by every user of a server process; only the plugin registries and a few constant
lookup tables may do that, and they are listed (with the reason) in ``ALLOWED`` below.
"""

import ast
from pathlib import Path

import pytest

PACKAGE = Path(__file__).parents[2] / "src" / "chromsim"

MUTABLE_CALLS = {"list", "dict", "set", "defaultdict", "OrderedDict", "Counter", "deque"}

# (module, name): why it is acceptable.
# Constant tables are made read-only instead (tuple, frozenset, types.MappingProxyType).
ALLOWED = {
    ("isotherms/registry.py", "_ISOTHERMS"): "plugin registry, filled at import (SPEC §5)",
    ("data/importers/registry.py", "_IMPORTERS"): "plugin registry, filled at import",
    ("optimizers/registry.py", "_OPTIMIZERS"): "plugin registry, filled at import (SPEC §6)",
    ("estimation/objectives.py", "_RESIDUALS"): "residual registry, filled at import",
    ("io/project_file.py", "MIGRATIONS"): "schema-migration registry, filled at import",
}
# ``__all__`` is the conventional list of public names, read by tools and never changed.
IGNORED_NAMES = {"__all__"}


def mutable_module_names(source: str) -> list[str]:
    """Names bound at module level to a mutable container (list/dict/set literal or call)."""
    names = []
    for node in ast.parse(source).body:
        if isinstance(node, ast.Assign):
            targets, value = node.targets, node.value
        elif isinstance(node, ast.AnnAssign) and node.value is not None:
            targets, value = [node.target], node.value
        else:
            continue
        mutable = isinstance(
            value, ast.List | ast.Dict | ast.Set | ast.ListComp | ast.DictComp | ast.SetComp
        ) or (
            isinstance(value, ast.Call)
            and getattr(value.func, "id", getattr(value.func, "attr", None)) in MUTABLE_CALLS
        )
        if mutable:
            names += [t.id for t in targets if isinstance(t, ast.Name)]
    return names


def test_checker_finds_mutable_bindings():
    source = "A = []\nB: dict = {}\nC = set()\nD = (1, 2)\nE = frozenset()\n"
    assert mutable_module_names(source) == ["A", "B", "C"]


@pytest.mark.parametrize(
    "path", sorted(PACKAGE.rglob("*.py")), ids=lambda p: p.relative_to(PACKAGE).as_posix()
)
def test_no_module_level_mutable_state(path):
    module = path.relative_to(PACKAGE).as_posix()  # ALLOWED uses "/" on every OS
    offending = [
        name
        for name in mutable_module_names(path.read_text(encoding="utf-8"))
        if (module, name) not in ALLOWED and name not in IGNORED_NAMES
    ]
    assert not offending, (
        f"{module} binds {offending} to a mutable container at module level. Make it a tuple/"
        "frozenset/MappingProxyType, move it into a function or class instance, or (for a plugin "
        "registry) add it to ALLOWED with the reason."
    )
