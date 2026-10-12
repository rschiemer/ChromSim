"""Code-only parity (SPEC §8): every UI page has an example script doing the same with the API.

Read with ``ast`` (NiceGUI is not needed). For each page, every ``api.<name>`` the page uses must
appear in its example scripts, unless it is UI plumbing (``UI_ONLY``) or a convenience wrapper
whose equivalence to the tutorial code is tested elsewhere (``EQUIVALENT``). The examples
themselves run in CI (tests/unit/test_examples.py).
"""

import ast
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[2]
PAGES_DIR = ROOT / "src" / "chromsim" / "ui" / "pages"
EXAMPLES_DIR = ROOT / "examples"

PAGE_EXAMPLES = {
    "project": ["00_project_files_and_units", "01_simulate_lge"],
    "data": ["02_compare_with_data"],
    "simulate": ["01_simulate_lge", "05_results_and_workspace"],
    "calibrate": ["03_calibrate"],
    "optimize": ["04_optimize_process"],
    "results": ["05_results_and_workspace"],
}

UI_ONLY = {
    # forms, errors and session plumbing
    "ChromSimError": "error display",
    "Workspace": "the session's workspace (Tutorial 05 shows it)",
    "codegen": "'Show Python code' itself",
    "Column": "form generated from the model",
    "Component": "form generated from the model",
    "Resin": "form generated from the model",
    "SolverSettings": "form generated from the model",
    "Step": "form generated from the model",
    "System": "form generated from the model",
    "isotherm_parameter_table": "isotherm parameter form",
    "project_from_dict": "validates the form (load_project in the tutorials)",
    "project_from_yaml": "YAML upload (load_project in the tutorials)",
    "fit_parameter_candidates": "check boxes offered in the fit form",
    "process_variable_candidates": "check boxes offered in the optimization form",
    "tracer_summary": "tracer result table (the tutorial prints the same values)",
    "plot_correlation_matrix": "correlation plot from a stored record (plot_correlation)",
    "Pool": "rebuilds the stored best pool for plotting",
}

EQUIVALENT = {
    # wrapper: tutorial calls; equivalence tested in tests/unit/test_ui_helpers.py and
    # tests/workflows/test_job_setups.py
    "tracer_run_from_file": {"read_chromatogram", "subtract_baseline", "TracerRun"},
    "lge_run_from_file": {"read_chromatogram", "conductivity_to_salt", "elution_salt_at_peak"},
    "FitSetup": {"inverse_fit", "parameter_uncertainty"},
    "AssayInput": {"read_assay"},
    "ProcessSetup": {"optimize_process"},
}


def api_names(path: Path) -> set[str]:
    """Names used as ``api.<name>`` in a Python file."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return {
        node.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute)
        and isinstance(node.value, ast.Name)
        and node.value.id == "api"
    }


def test_every_page_has_examples():
    pages = {path.stem for path in PAGES_DIR.glob("*.py") if path.stem != "__init__"}
    assert pages == set(PAGE_EXAMPLES)
    for examples in PAGE_EXAMPLES.values():
        for name in examples:
            assert (EXAMPLES_DIR / f"{name}.py").is_file(), name


@pytest.mark.parametrize("page", PAGE_EXAMPLES)
def test_page_uses_only_what_its_examples_show(page):
    shown = set().union(*(api_names(EXAMPLES_DIR / f"{e}.py") for e in PAGE_EXAMPLES[page]))
    missing = []
    for name in sorted(api_names(PAGES_DIR / f"{page}.py") - shown - set(UI_ONLY)):
        if name in EQUIVALENT and EQUIVALENT[name] <= shown:
            continue
        missing.append(name)
    assert not missing, (
        f"The {page} page uses api.{', api.'.join(missing)}, which its example scripts "
        f"{PAGE_EXAMPLES[page]} do not show. Add it to an example (SPEC §8 code-only parity)."
    )
