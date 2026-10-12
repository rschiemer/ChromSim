"""Generated "equivalent Python code" (SPEC §8) is valid and uses only ``chromsim.api``."""

import ast

import pytest
from project_examples import make_physics_project

from chromsim import api
from chromsim.estimation.parameters import FitParameter

ROOT = "my-workspace"


def fit_setup():
    return api.FitSetup(
        project="lge",
        result_name="fit",
        chromatogram="run.csv",
        assays=[api.AssayInput(file="assay.xlsx", profile="profile.yaml")],
        parameters=[FitParameter(path="isotherm.parameters.mAb.k_eq", lower=1e-3, upper=1.0)],
    )


def process_setup():
    return api.ProcessSetup(
        project="lge",
        result_name="opt",
        target="mAb",
        purity=0.99,
        variables=[FitParameter(path="method.3.duration", lower=5, upper=40, unit="CV")],
        options={"max_evaluations": 10},
    )


SCRIPTS = {
    "project": lambda: api.codegen.code_for_project(ROOT, "lge"),
    "simulation": lambda: api.codegen.code_for_simulation(ROOT, "lge", "result"),
    "data": lambda: api.codegen.code_for_data(ROOT, "run.csv", [("assay.xlsx", "profile.yaml")]),
    "fit": lambda: api.codegen.code_for_fit(ROOT, fit_setup()),
    "process": lambda: api.codegen.code_for_process(ROOT, process_setup()),
    "results": lambda: api.codegen.code_for_results(ROOT, "result"),
}


@pytest.mark.parametrize("name", SCRIPTS)
def test_generated_code_compiles_and_uses_existing_api_names(name):
    tree = ast.parse(SCRIPTS[name]())
    imports = [node for node in ast.walk(tree) if isinstance(node, ast.Import | ast.ImportFrom)]
    assert [ast.unparse(node) for node in imports] == ["from chromsim import api"]
    used = {
        node.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute) and getattr(node.value, "id", None) == "api"
    }
    missing = sorted(name for name in used if not hasattr(api, name))
    assert not missing, f"generated code uses names chromsim.api does not export: {missing}"


@pytest.mark.parametrize(
    ("write_code", "build", "model"),
    [
        (api.codegen.code_for_fit, fit_setup, api.FitSetup),
        (api.codegen.code_for_process, process_setup, api.ProcessSetup),
    ],
)
def test_fit_and_process_scripts_rebuild_the_same_setup(write_code, build, model):
    script = write_code(ROOT, build())
    call = next(
        node
        for node in ast.walk(ast.parse(script))
        if isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "model_validate"
    )
    assert model.model_validate(ast.literal_eval(call.args[0])) == build()


def test_simulation_script_runs_and_stores_the_result(tmp_path, monkeypatch):
    storage = api.FileStorage(tmp_path / "ws")
    storage.save_project("tracer", make_physics_project("LRM", n_col=20))
    script = api.codegen.code_for_simulation(str(storage.root), "tracer", "from-script")
    monkeypatch.chdir(tmp_path)
    exec(compile(script, "generated", "exec"), {"__name__": "__main__"})  # our own template
    assert storage.list_simulation_results() == ["from-script"]
