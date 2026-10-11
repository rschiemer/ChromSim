"""Shared example projects for tests (importable because pytest adds tests/ to sys.path)."""

from chromsim.core.project import Project


def make_project(**overrides) -> Project:
    """The SPEC §4.3 example project, built in code."""
    data = {
        "schema_version": 1,
        "name": "CEX_mAb_LGE",
        "components": [
            {"name": "salt"},
            {"name": "mAb", "mw_kDa": 148, "extinction_coeff_L_per_g_cm": 1.4},
        ],
        "column": {"model": "GRM", "length": "20 cm", "diameter": "0.5 cm"},
        "resin": {
            "name": "Generic CEX",
            "particle_radius": "45 um",
            "ionic_capacity": "1200 mM",
            "col_porosity": 0.36,
            "par_porosity": 0.85,
        },
        "system": {"tubing_volume": "0.25 mL", "mixer_volume": "0.2 mL"},
        "isotherm": {
            "model": "SMA",
            "parameters": {"mAb": {"k_eq": 1.0e-2, "nu": 6.0, "sigma": 30.0, "k_kin": 1.0e2}},
        },
        "method": [
            {
                "name": "equilibrate",
                "duration": "5 CV",
                "flow": "150 cm/h",
                "inlet": {"salt": "50 mM"},
            },
            {
                "name": "load",
                "duration": "2 mL",
                "flow": "150 cm/h",
                "inlet": {"salt": "50 mM", "mAb": "1 g/L"},
            },
            {"name": "wash", "duration": "5 CV", "flow": "150 cm/h", "inlet": {"salt": "50 mM"}},
            {
                "name": "elute",
                "duration": "20 CV",
                "flow": "150 cm/h",
                "inlet": {"salt": "50 mM"},
                "inlet_end": {"salt": "500 mM"},
            },
        ],
        "solver": {"backend": "cadet", "n_col": 100, "n_par": 25, "abstol": 1.0e-8},
    }
    data.update(overrides)
    return Project.model_validate(data)


# Physics-test column: 10 cm × 1 cm, ε_c = 0.4, ε_p = 0.5 → ε_t = 0.4 + 0.6 × 0.5 = 0.7.
# V_c = π (0.5 cm)² × 10 cm = 7.85398 mL; at 1 mL/min the unretained time is
# t_0 = ε_t V_c / Q = 5.49779 min = 329.867 s.
PHYSICS_COLUMN_VOLUME_ML = 7.853981633974483
PHYSICS_FLOW_ML_PER_MIN = 1.0
PHYSICS_COL_POROSITY = 0.4
PHYSICS_PAR_POROSITY = 0.5
PHYSICS_TOTAL_POROSITY = 0.7


def make_physics_project(
    column_model: str = "LRM",
    components: list[str] | None = None,
    isotherm: dict | None = None,
    method: list[dict] | None = None,
    system: dict | None = None,
    n_col: int = 100,
    n_par: int = 10,
) -> Project:
    """A small, fast column for physics tests. Defaults: non-binding tracer, LRM, no system."""
    components = components or ["tracer"]
    data = {
        "name": "physics",
        "components": [{"name": name} for name in components],
        "column": {
            "model": column_model,
            "length": "10 cm",
            "diameter": "1 cm",
            "axial_dispersion": "1e-7 m^2/s",
            "film_diffusion": {name: "1e-5 m/s" for name in components},
            "pore_diffusion": {name: "1e-10 m^2/s" for name in components},
        },
        "resin": {
            "name": "test resin",
            "particle_radius": "30 um",
            "ionic_capacity": "500 mM",
            "col_porosity": PHYSICS_COL_POROSITY,
            "par_porosity": PHYSICS_PAR_POROSITY,
        },
        "system": system or {},
        "isotherm": isotherm or {"model": "Linear", "parameters": {}},
        "method": method
        or [
            {
                "name": "pulse",
                "duration": "6 s",
                "flow": "1 mL/min",
                "inlet": {components[0]: "1 mM"},
            },
            {"name": "elute", "duration": "30 min", "flow": "1 mL/min", "inlet": {}},
        ],
        "solver": {"n_col": n_col, "n_par": n_par},
    }
    return Project.model_validate(data)


def make_binding_project(
    components: list[dict],
    isotherm: dict,
    load: dict,
    elute: dict | None = None,
    column_model: str = "LRM",
    n_col: int = 30,
) -> Project:
    """Small column (5 cm × 1 cm, ε_c 0.4, ε_p 0.5, Λ 500 mM): equilibrate, load, elute 60 min.

    Modifiers (components in ``load`` without isotherm parameters) are present from the start;
    ``elute`` gives their concentrations at the end of a linear elution gradient.
    """
    modifiers = {k: v for k, v in load.items() if k not in isotherm.get("parameters", {})}
    elute_step = {
        "name": "elute",
        "duration": "60 min",
        "flow": "1 mL/min",
        "inlet": dict(modifiers),
    }
    if elute:
        elute_step["inlet_end"] = elute
    data = {
        "name": "binding",
        "components": components,
        "column": {
            "model": column_model,
            "length": "5 cm",
            "diameter": "1 cm",
            "axial_dispersion": "1e-7 m^2/s",
            "film_diffusion": {c["name"]: "1e-5 m/s" for c in components},
            "pore_diffusion": {c["name"]: "1e-10 m^2/s" for c in components},
        },
        "resin": {
            "name": "r",
            "particle_radius": "30 um",
            "ionic_capacity": "500 mM",
            "col_porosity": 0.4,
            "par_porosity": 0.5,
        },
        "isotherm": isotherm,
        "method": [
            {"name": "equilibrate", "duration": "1 min", "flow": "1 mL/min", "inlet": modifiers},
            {"name": "load", "duration": "2 min", "flow": "1 mL/min", "inlet": load},
            elute_step,
        ],
        "solver": {"n_col": n_col},
    }
    return Project.model_validate(data)
