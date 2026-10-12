"""Helpers added for the web UI (SPEC §8) give the same answers as the tutorial code."""

from pathlib import Path

import pytest
from project_examples import make_project

from chromsim import api
from chromsim.core.column import Column

ROOT = Path(__file__).parents[2]
AKTA = ROOT / "tests" / "data" / "akta"
EXAMPLES = ROOT / "examples"
ML = 1e-6


def test_unit_choices_follow_field_dimensions():
    length = api.field_dimensions(Column.model_fields["length"])
    assert length == ("length",)
    assert api.unit_choices(length)[0] == "cm"
    assert api.field_dimensions(Column.model_fields["model"]) is None


def test_isotherm_parameter_table_lists_sma_parameters():
    rows = {row["name"]: row for row in api.isotherm_parameter_table("SMA")}
    assert {"k_eq", "nu", "sigma"} <= set(rows)
    assert rows["k_eq"]["log_scale"]
    assert rows["nu"]["per_component"]


def test_project_from_dict_reports_the_location_of_an_error():
    data = make_project().model_dump(mode="json")
    data["column"]["length"] = "20 kg"
    with pytest.raises(api.ChromSimError, match="column.length"):
        api.project_from_dict(data)
    data["column"]["length"] = "25 cm"
    assert api.project_from_dict(data).column.length == make_project().column.length.model_copy(
        update={"magnitude": 25.0}
    )


def test_fit_candidates_cover_every_isotherm_parameter_and_contain_the_current_value():
    project = api.load_project(EXAMPLES / "projects" / "cex_mab_hmw_lge.yaml")
    candidates = api.fit_parameter_candidates(project)
    paths = {candidate.path for candidate in candidates}
    for component, values in project.isotherm.parameters.items():
        for name in values:
            assert f"isotherm.parameters.{component}.{name}" in paths
    for candidate in candidates:
        value = api.get_value(project, candidate.path, candidate.unit)
        assert candidate.lower <= value <= candidate.upper, candidate.path


def test_process_candidates_are_step_durations_and_gradient_ends():
    project = make_project()
    candidates = api.process_variable_candidates(project)
    assert sum(c.path.endswith(".duration") for c in candidates) == len(project.method)
    for candidate in candidates:
        value = api.get_value(project, candidate.path, candidate.unit)
        assert candidate.lower == pytest.approx(value / 2)
        assert candidate.upper == pytest.approx(value * 2)


def test_tracer_helpers_reproduce_tutorial_03():
    project = api.load_project(EXAMPLES / "projects" / "cex_mab_hmw_lge.yaml")
    flow = project.resolved_method()[0].flow_rate_m3_per_s
    runs = [
        api.tracer_run_from_file(AKTA / f"{file}.csv", kind)
        for file, kind in (
            ("tracer_system", "system"),
            ("tracer_dextran", "non_penetrating"),
            ("tracer_acetone", "pore_penetrating"),
        )
    ]
    analysis = api.analyse_tracers(runs, project.column, flow)
    assert analysis.col_porosity == pytest.approx(0.356, abs=0.005)
    assert analysis.total_porosity == pytest.approx(0.901, abs=0.005)
    assert analysis.par_porosity == pytest.approx(0.846, abs=0.005)
    rows = {row["quantity"]: row for row in api.tracer_summary(analysis)}
    assert any(row["unit"] == "mL" for row in rows.values())


def test_lge_helper_reproduces_tutorial_03_yamamoto():
    project = api.load_project(EXAMPLES / "projects" / "cex_mab_hmw_lge.yaml")
    calibration = api.load_conductivity_calibration(
        EXAMPLES / "calibrations" / "nacl_conductivity.yaml"
    )
    runs = [
        api.lge_run_from_file(AKTA / f"lge_{cv}cv.csv", project, cv, calibration)
        for cv in (10, 20, 30, 40)
    ]
    result = api.yamamoto_fit(
        runs,
        project.column.volume_m3,
        project.resin.total_porosity,
        project.resin.ionic_capacity_mol_per_m3,
    )
    assert result.nu == pytest.approx(6.20, abs=0.05)
