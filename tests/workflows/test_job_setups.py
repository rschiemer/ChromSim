"""Fit and process-optimization setups run from stored workspace data (what the UI submits)."""

from pathlib import Path

import pytest
from project_examples import make_binding_project

from chromsim import api
from chromsim.estimation.parameters import FitParameter
from chromsim.jobs import AssayInput, FitSetup, ProcessSetup, run_fit, run_process_optimization
from chromsim.storage import FileStorage

pytestmark = pytest.mark.slow
ROOT = Path(__file__).parents[2]


def stored_workspace(tmp_path) -> FileStorage:
    storage = FileStorage(tmp_path / "ws")
    project = api.load_project(ROOT / "examples/projects/cex_mab_hmw_lge.yaml")
    project.solver.n_col, project.solver.n_par = 20, 4
    storage.save_project("first guess", project)
    for path in (
        ROOT / "tests/data/akta/cex_mab_hmw_run01.csv",
        ROOT / "tests/data/assays/sec_hplc.xlsx",
        ROOT / "examples/assay_profiles/sec_hplc_long.yaml",
    ):
        storage.save_raw_file(path.name, path.read_bytes())
    return storage


def test_fit_setup_recovers_k_eq_from_stored_data(tmp_path):
    storage = stored_workspace(tmp_path)
    setup = FitSetup(
        project="first guess",
        result_name="fit1",
        chromatogram="cex_mab_hmw_run01.csv",
        fit_uv=False,
        assays=[AssayInput(file="sec_hplc.xlsx", profile="sec_hplc_long.yaml")],
        parameters=[
            FitParameter(
                path="isotherm.parameters.mAb.k_eq",
                lower=1e-3,
                upper=1e-1,
                scale="log",
                name="k_eq",
            )
        ],
        max_evaluations=30,
    )
    calls = []
    assert run_fit(storage, setup, callback=lambda n, cost: calls.append(n)) == "fit1"
    record = storage.load_record("fits", "fit1")
    # The data were made with k_eq = 1.5e-2 (other parameters as in the project); the coarse grid
    # moves the estimate slightly.
    assert record["values"]["k_eq"] == pytest.approx(1.5e-2, rel=0.1)
    assert record["uncertainty"][0]["parameter"] == "k_eq"
    assert storage.load_project("fit1-fitted").isotherm.parameters["mAb"]["k_eq"] == pytest.approx(
        record["values"]["k_eq"]
    )
    assert "fit1" in storage.list_simulation_results()
    assert calls


def test_process_setup_stores_best_method(tmp_path):
    storage = FileStorage(tmp_path / "ws")
    components = [{"name": "salt"}, {"name": "A", "mw_kDa": 150.0}, {"name": "B", "mw_kDa": 150.0}]
    isotherm = {
        "model": "SMA",
        "parameters": {
            "A": {"k_eq": 1e-2, "nu": 5.0, "sigma": 20.0, "k_kin": 10.0},
            "B": {"k_eq": 3e-2, "nu": 5.0, "sigma": 20.0, "k_kin": 10.0},
        },
    }
    project = make_binding_project(
        components,
        isotherm,
        {"salt": "50 mM", "A": "0.1 mM", "B": "0.02 mM"},
        elute={"salt": "1000 mM"},
        n_col=40,
    )
    project.method[2].duration = "10 min"
    storage.save_project("separation", project)
    setup = ProcessSetup(
        project="separation",
        result_name="opt1",
        target="A",
        purity=0.99,
        variables=[
            FitParameter(
                path="method.2.duration", lower=5, upper=120, unit="min", name="gradient_min"
            )
        ],
        pool_steps=["elute"],
        options={"max_evaluations": 12},
    )
    assert run_process_optimization(storage, setup) == "opt1"
    record = storage.load_record("processes", "opt1")
    assert record["best_kpis"]["purity"] >= 0.99 - 1e-9
    assert record["best_kpis"]["recovery"] > 0.41  # better than the 10 min start
    assert "opt1-best" in storage.list_projects()
    assert storage.load_optimization_result("opt1").n_evaluations == len(record["table"])
