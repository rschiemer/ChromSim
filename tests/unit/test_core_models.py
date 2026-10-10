"""Domain model: components, column, resin, system, method steps and project (SPEC §4.2)."""

import math

import pytest
from project_examples import make_project
from pydantic import ValidationError

from chromsim.core.column import Column, ColumnModel
from chromsim.core.components import Component
from chromsim.core.method import Step, resolve_method
from chromsim.core.project import IsothermConfig, SolverSettings
from chromsim.core.resin import Resin
from chromsim.core.system import System

# --- Component ------------------------------------------------------------------------------------


def test_component_si_properties():
    mab = Component(name="mAb", mw_kDa=148, extinction_coeff_L_per_g_cm=1.4)
    assert mab.molar_mass_kg_per_mol == pytest.approx(148.0)
    # 1.4 L/(g·cm) = 1.4 × 1e-3 m³ / (1e-3 kg × 1e-2 m) = 140 m²/kg
    assert mab.extinction_coeff_m2_per_kg == pytest.approx(140.0)
    assert Component(name="salt").molar_mass_kg_per_mol is None


@pytest.mark.parametrize("field", ["mw_kDa", "extinction_coeff_L_per_g_cm"])
def test_component_rejects_non_positive_values(field):
    with pytest.raises(ValidationError, match=field):
        Component(name="x", **{field: 0})


# --- Column ---------------------------------------------------------------------------------------


def test_column_geometry_in_si():
    column = Column(model="GRM", length="20 cm", diameter="0.5 cm")
    assert column.model is ColumnModel.GRM
    assert column.length_m == pytest.approx(0.20)
    assert column.diameter_m == pytest.approx(0.005)
    assert column.cross_section_m2 == pytest.approx(math.pi * 0.0025**2)
    assert column.volume_m3 == pytest.approx(math.pi * 0.0025**2 * 0.20)


def test_column_defaults_to_grm():
    assert Column(length="10 cm", diameter="1 cm").model is ColumnModel.GRM


@pytest.mark.parametrize(
    ("field", "value", "hint"),
    [
        ("length", "20 mL", "not a length"),
        ("length", "0 cm", "zero"),
        ("diameter", "-1 cm", "negative"),
        ("model", "XYZ", "LRM"),
    ],
)
def test_column_rejects_bad_values(field, value, hint):
    data = {"length": "20 cm", "diameter": "0.5 cm", field: value}
    with pytest.raises(ValidationError, match=hint):
        Column(**data)


def test_column_transport_parameters_are_per_component():
    column = Column(
        length="20 cm",
        diameter="0.5 cm",
        axial_dispersion="1e-7 m^2/s",
        film_diffusion={"mAb": "1e-5 m/s"},
        pore_diffusion={"mAb": "1e-11 m^2/s"},
    )
    assert column.film_diffusion["mAb"].unit == "m/s"


# --- Resin ----------------------------------------------------------------------------------------


def test_resin_total_porosity():
    resin = Resin(name="r", particle_radius="45 um", col_porosity=0.36, par_porosity=0.85)
    # ε_t = ε_c + (1 − ε_c) ε_p = 0.36 + 0.64 × 0.85 = 0.904
    assert resin.total_porosity == pytest.approx(0.904)
    assert resin.particle_radius_m == pytest.approx(45e-6)
    assert resin.ionic_capacity_mol_per_m3 is None


@pytest.mark.parametrize("porosity", [0.0, 1.0, -0.1, 1.2])
def test_resin_rejects_porosity_outside_unit_interval(porosity):
    with pytest.raises(ValidationError, match="col_porosity"):
        Resin(name="r", particle_radius="45 um", col_porosity=porosity, par_porosity=0.85)


def test_resin_ionic_capacity_must_be_molar():
    with pytest.raises(ValidationError, match="not a molar concentration"):
        Resin(
            name="r",
            particle_radius="45 um",
            ionic_capacity="50 g/L",
            col_porosity=0.36,
            par_porosity=0.85,
        )


# --- System ---------------------------------------------------------------------------------------


def test_system_defaults_to_zero_dead_volume():
    system = System()
    assert system.tubing_volume_m3 == 0.0
    assert system.mixer_volume_m3 == 0.0


def test_system_volumes_in_si():
    system = System(tubing_volume="0.25 mL", mixer_volume="0.2 mL")
    assert system.tubing_volume_m3 == pytest.approx(0.25e-6)
    assert system.mixer_volume_m3 == pytest.approx(0.2e-6)


# --- Step -----------------------------------------------------------------------------------------


def test_step_rejects_gradient_end_without_start():
    with pytest.raises(ValidationError, match="inlet_end.*'mAb'.*inlet"):
        Step(
            name="elute",
            duration="20 CV",
            flow="150 cm/h",
            inlet={"salt": "50 mM"},
            inlet_end={"mAb": "1 g/L"},
        )


def test_step_rejects_wrong_duration_unit():
    with pytest.raises(ValidationError, match="not a duration"):
        Step(name="wash", duration="5 mM", flow="150 cm/h", inlet={})


def test_step_is_gradient():
    step = Step(
        name="e",
        duration="1 CV",
        flow="1 mL/min",
        inlet={"salt": "1 mM"},
        inlet_end={"salt": "2 mM"},
    )
    assert step.is_gradient
    assert not Step(name="w", duration="1 CV", flow="1 mL/min", inlet={}).is_gradient


# --- Project validation ---------------------------------------------------------------------------


def test_example_project_builds():
    project = make_project()
    assert project.component_names == ["salt", "mAb"]
    assert project.isotherm == IsothermConfig(
        model="SMA", parameters={"mAb": {"k_eq": 0.01, "nu": 6.0, "sigma": 30.0, "k_kin": 100.0}}
    )
    assert project.solver == SolverSettings(backend="cadet", n_col=100, n_par=25, abstol=1e-8)


def test_project_rejects_duplicate_component_names():
    with pytest.raises(ValidationError, match="Component names must be unique.*'salt'"):
        make_project(components=[{"name": "salt"}, {"name": "salt"}])


def test_project_rejects_unknown_component_in_inlet():
    method = [{"name": "load", "duration": "1 CV", "flow": "1 mL/min", "inlet": {"hcp": "1 mM"}}]
    with pytest.raises(ValidationError, match="Step 'load'.*unknown component 'hcp'.*salt, mAb"):
        make_project(method=method)


def test_project_rejects_unknown_component_in_isotherm():
    isotherm = {"model": "SMA", "parameters": {"IgG": {"k_eq": 1.0}}}
    with pytest.raises(ValidationError, match="isotherm parameters.*unknown component 'IgG'"):
        make_project(isotherm=isotherm)


def test_project_rejects_unknown_component_in_transport():
    column = {"length": "20 cm", "diameter": "0.5 cm", "film_diffusion": {"IgG": "1e-5 m/s"}}
    with pytest.raises(ValidationError, match="film_diffusion.*unknown component 'IgG'"):
        make_project(column=column)


def test_project_requires_molecular_weight_for_mass_concentration():
    components = [{"name": "salt"}, {"name": "mAb"}]  # no mw_kDa
    with pytest.raises(ValidationError, match="'mAb'.*1 g/L.*mw_kDa"):
        make_project(components=components)


def test_project_rejects_empty_method():
    with pytest.raises(ValidationError, match="method"):
        make_project(method=[])


def test_project_rejects_misspelled_keys():
    # extra="forbid": a typo must not be silently ignored
    with pytest.raises(ValidationError, match="lenght"):
        make_project(column={"lenght": "20 cm", "diameter": "0.5 cm"})


def test_project_validates_assignment():
    project = make_project()
    with pytest.raises(ValidationError):
        project.column.length = "20 mL"


# --- method timeline in SI ------------------------------------------------------------------------


def test_resolve_method_for_example_project():
    project = make_project()
    steps = resolve_method(project.method, project.components, project.column)

    # Q = 150 cm/h × π (0.25 cm)² = 8.18123e-9 m³/s; 1 CV = L/u = 480 s.
    flow = 1.5 / 3600 * math.pi * 0.0025**2
    expected_durations = [2400.0, 2e-6 / flow, 2400.0, 9600.0]  # 5 CV, 2 mL, 5 CV, 20 CV
    assert [step.name for step in steps] == ["equilibrate", "load", "wash", "elute"]

    start = 0.0
    for step, duration in zip(steps, expected_durations, strict=True):
        assert step.start_s == pytest.approx(start)
        assert step.end_s == pytest.approx(start + duration)
        assert step.flow_rate_m3_per_s == pytest.approx(flow)
        start += duration

    # Concentrations in mol/m³, ordered like project.components = [salt, mAb].
    # 1 g/L mAb = 1 kg/m³ ÷ 148 kg/mol = 6.75676e-3 mol/m³; unspecified components are 0.
    assert steps[0].inlet_start_mol_per_m3 == pytest.approx((50.0, 0.0))
    assert steps[1].inlet_start_mol_per_m3 == pytest.approx((50.0, 1.0 / 148.0))
    assert steps[1].inlet_end_mol_per_m3 == steps[1].inlet_start_mol_per_m3
    assert steps[3].inlet_start_mol_per_m3 == pytest.approx((50.0, 0.0))
    assert steps[3].inlet_end_mol_per_m3 == pytest.approx((500.0, 0.0))
