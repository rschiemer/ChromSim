"""Contract every registered isotherm must meet (SPEC §5, §9 plugins tier; CLAUDE.md rule 6).

The tests are generated for ``list_isotherms()``, so a new isotherm is checked automatically:
documentation (equation, reference, units, docs page), a valid example configuration, the
mapping onto CADET-Process names, and a load–elute simulation that conserves mass.
"""

import re
from pathlib import Path

import CADETProcess.processModel as process_model
import numpy as np
import pytest

from chromsim import api
from chromsim.isotherms.cadet_native import CadetNativeIsotherm
from chromsim.isotherms.python_isotherm import PythonIsotherm, binding_parameters
from chromsim.isotherms.registry import get_isotherm, list_isotherms

DOCS = Path(__file__).parents[2] / "docs" / "theory" / "isotherms"
ISOTHERMS = list_isotherms()


def docs_slug(name: str) -> str:
    """Docs page name: CamelCase → snake_case ("BiLangmuir" → "bi_langmuir", "SMA" → "sma")."""
    return re.sub(r"(?<=[a-z])(?=[A-Z])|(?<=[A-Z])(?=[A-Z][a-z])", "_", name).lower()


def example_project(isotherm, n_col: int = 30) -> api.Project:
    """Small LRM column: equilibrate, load 2 min, then 60 min elution (salt gradient if any)."""
    example = isotherm.example
    names = [c["name"] for c in example["components"]]
    modifiers = {k: v for k, v in example["load"].items() if k not in example["parameters"]}
    elute_start = dict(modifiers)
    elute_step = {"name": "elute", "duration": "60 min", "flow": "1 mL/min", "inlet": elute_start}
    if example["elute"]:
        elute_step["inlet_end"] = {k: v for k, v in example["elute"].items() if k in elute_start}
    data = {
        "name": f"{isotherm.name} contract",
        "components": example["components"],
        "column": {
            "model": "LRM",
            "length": "5 cm",
            "diameter": "1 cm",
            "axial_dispersion": "1e-7 m^2/s",
        },
        "resin": {
            "name": "r",
            "particle_radius": "30 um",
            "ionic_capacity": "500 mM",
            "col_porosity": 0.4,
            "par_porosity": 0.5,
        },
        "isotherm": {"model": isotherm.name, "parameters": example["parameters"]}
        | (
            {"global_parameters": example["global_parameters"]}
            if example.get("global_parameters")
            else {}
        ),
        "method": [
            {"name": "equilibrate", "duration": "1 min", "flow": "1 mL/min", "inlet": modifiers},
            {"name": "load", "duration": "2 min", "flow": "1 mL/min", "inlet": example["load"]},
            elute_step,
        ],
        "solver": {"n_col": n_col},
    }
    assert set(names) >= set(example["parameters"])
    return api.Project.model_validate(data)


@pytest.mark.parametrize("name", ISOTHERMS)
def test_documentation(name):
    isotherm = get_isotherm(name)
    doc = type(isotherm).__doc__ or ""
    assert ".. math::" in doc, f"{name}: docstring needs the model equation (.. math::)"
    assert "Reference" in doc, f"{name}: docstring needs a literature reference"
    for parameter_name, parameter in isotherm.parameters.items():
        assert parameter.unit and parameter.description, (
            f"{name}.{parameter_name}: unit/description"
        )
    page = DOCS / f"{docs_slug(name)}.md"
    assert page.is_file(), f"{name}: add the theory page {page.relative_to(DOCS.parents[2])}"


@pytest.mark.parametrize("name", ISOTHERMS)
def test_example_is_a_valid_project(name):
    project = example_project(get_isotherm(name))
    assert project.isotherm.model == name


@pytest.mark.parametrize(
    "name", [n for n in ISOTHERMS if isinstance(get_isotherm(n), CadetNativeIsotherm)]
)
def test_cadet_mapping_names_exist(name):
    isotherm = get_isotherm(name)
    project = example_project(isotherm)
    component_system = process_model.ComponentSystem(project.component_names)
    binding = getattr(process_model, isotherm.cadet_class)(component_system, name="binding")
    values = isotherm.cadet_parameters(
        project.component_names,
        project.isotherm.parameters,
        500.0,
        project.isotherm.global_parameters,
    )
    for attribute in values:
        assert attribute in binding.parameters, (
            f"{isotherm.cadet_class} has no parameter '{attribute}'"
        )


@pytest.mark.parametrize(
    "name", [n for n in ISOTHERMS if isinstance(get_isotherm(n), CadetNativeIsotherm)]
)
def test_cadet_simulation_conserves_mass(name):
    isotherm = get_isotherm(name)
    project = example_project(isotherm)
    result = api.simulate(project)
    load = project.resolved_method()[1]
    for component in isotherm.example["parameters"]:
        index = project.component_names.index(component)
        injected = load.inlet_start_mol_per_m3[index] * load.duration_s  # constant flow
        eluted = np.trapezoid(result.outlet(component), result.time_s)
        # Everything elutes within the 60 min elution of the example; CADET conserves mass to
        # ~1e-6, so 1 % only allows for material still leaving at the end.
        assert eluted == pytest.approx(injected, rel=0.01), f"{name}: {component} mass balance"


PYTHON_ISOTHERMS = [n for n in ISOTHERMS if isinstance(get_isotherm(n), PythonIsotherm)]
FLUX_ISOTHERMS = [n for n in ISOTHERMS if hasattr(get_isotherm(n), "flux")]


@pytest.mark.parametrize("name", FLUX_ISOTHERMS)
def test_flux_shape_and_non_binding_rows(name):
    isotherm = get_isotherm(name)
    project = example_project(isotherm)
    names = project.component_names + ["inert"]
    parameters = binding_parameters(
        names, project.isotherm.parameters, project.isotherm.global_parameters, 500.0
    )
    rng = np.random.default_rng(0)
    c = rng.uniform(0.1, 1.0, (len(names), 7))
    q = np.zeros_like(c)
    q[0] = 500.0  # bound salt where an isotherm uses it; ignored otherwise
    rate = isotherm.flux(c, q, parameters)
    assert rate.shape == c.shape
    assert np.all(rate[-1] == 0.0), "non-binding components must have zero flux"
    assert np.all(np.isfinite(rate))


@pytest.mark.parametrize("name", PYTHON_ISOTHERMS)
def test_python_isotherm_runs_on_reference_backend_and_conserves_mass(name):
    isotherm = get_isotherm(name)
    project = example_project(isotherm)
    project.solver.backend = "reference"
    result = api.simulate(project)
    load = project.resolved_method()[1]
    for component in isotherm.example["parameters"]:
        index = project.component_names.index(component)
        injected = load.inlet_start_mol_per_m3[index] * load.duration_s
        eluted = np.trapezoid(result.outlet(component), result.time_s)
        assert eluted == pytest.approx(injected, rel=0.01), f"{name}: {component} mass balance"
