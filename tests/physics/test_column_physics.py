"""Physics of simulated columns against analytical results (SPEC §9, physics tier).

Every test uses the small column from ``project_examples.make_physics_project``:
10 cm × 1 cm, ε_c = 0.4, ε_p = 0.5 (ε_t = 0.7), V_c = 7.854 mL, Q = 1 mL/min.

Moment results hold exactly for the model equations (closed Danckwerts boundaries), so the
remaining differences are numerical: the 1 s output grid, trapezoidal integration and the finite
cut-off of the simulated time. Tolerances are set for those.
"""

import numpy as np
import pytest
from project_examples import (
    PHYSICS_COL_POROSITY,
    PHYSICS_COLUMN_VOLUME_ML,
    PHYSICS_FLOW_ML_PER_MIN,
    PHYSICS_PAR_POROSITY,
    PHYSICS_TOTAL_POROSITY,
    make_physics_project,
)

from chromsim.backends import get_backend

FLOW_ML_PER_S = PHYSICS_FLOW_ML_PER_MIN / 60.0
PULSE_S = 6.0
# Unretained time t_0 = ε_t V_c / Q = 329.867 s.
T0_S = PHYSICS_TOTAL_POROSITY * PHYSICS_COLUMN_VOLUME_ML / FLOW_ML_PER_S

# Relative tolerance on first moments and breakthrough times: the 1 s output grid and trapezoidal
# integration give errors of order (1 s)² / σ² relative to peak widths σ ≥ 10 s, i.e. < 1e-3.
MOMENT_RTOL = 1e-3
# Mass balance: CADET integrates to abstol 1e-8 mol/m³; trapezoidal integration of smooth peaks
# on a 1 s grid adds < 1e-4 relative error.
MASS_RTOL = 1e-3


def simulate(project):
    return get_backend("cadet").simulate(project)


def first_moment(time_s, concentration):
    """μ_1 = ∫ t c dt / ∫ c dt, s."""
    return np.trapezoid(time_s * concentration, time_s) / np.trapezoid(concentration, time_s)


def recovered_fraction(result, component, injected_mol_per_m3_times_s):
    """Eluted amount ∫ c dt divided by the injected amount c_feed · t_pulse (constant flow)."""
    return np.trapezoid(result.outlet(component), result.time_s) / injected_mol_per_m3_times_s


# --- non-binding tracer -------------------------------------------------------------------------


@pytest.mark.parametrize("column_model", ["LRM", "LRMP", "GRM"])
def test_tracer_first_moment_is_residence_time(column_model):
    # Pore-penetrating, non-binding tracer: μ_1 = t_pulse/2 + ε_t V_c / Q
    # (Kucera, J. Chromatogr. 19 (1965) 237; also Guiochon et al. 2006, ch. 6).
    result = simulate(make_physics_project(column_model))
    expected = PULSE_S / 2 + T0_S
    assert first_moment(result.time_s, result.outlet("tracer")) == pytest.approx(
        expected, rel=MOMENT_RTOL
    )


@pytest.mark.parametrize("column_model", ["LRM", "GRM"])
def test_tracer_mass_balance(column_model):
    result = simulate(make_physics_project(column_model))
    assert recovered_fraction(result, "tracer", 1.0 * PULSE_S) == pytest.approx(1.0, rel=MASS_RTOL)


def test_system_dead_volume_delays_tracer_by_v_over_q():
    # Mixer (CSTR) and tubing (DPFR) each add their mean residence time V/Q to the first moment:
    # (0.5 mL + 0.5 mL) / (1 mL/min) = 60 s.
    without = simulate(make_physics_project("LRM"))
    with_system = simulate(
        make_physics_project("LRM", system={"mixer_volume": "0.5 mL", "tubing_volume": "0.5 mL"})
    )
    shift_s = first_moment(with_system.time_s, with_system.outlet("tracer")) - first_moment(
        without.time_s, without.outlet("tracer")
    )
    assert shift_s == pytest.approx(60.0, rel=5 * MOMENT_RTOL)  # difference of two moments


def test_linear_salt_gradient_is_delayed_by_residence_time():
    # A linear ramp through a linear, non-binding system leaves it shifted by the mean residence
    # time: c_out(t) = c_in(t − t_0) once the transient from the ramp start has passed.
    method = [
        {
            "name": "equilibrate",
            "duration": "5 min",
            "flow": "1 mL/min",
            "inlet": {"salt": "50 mM"},
        },
        {
            "name": "gradient",
            "duration": "40 min",
            "flow": "1 mL/min",
            "inlet": {"salt": "50 mM"},
            "inlet_end": {"salt": "450 mM"},
        },
    ]
    result = simulate(make_physics_project("LRM", components=["salt"], method=method))
    gradient_start_s, slope_mM_per_s = 300.0, 400.0 / 2400.0
    t_check = gradient_start_s + 1200.0 + T0_S  # middle of the ramp, after the delay
    salt_out = np.interp(t_check, result.time_s, result.outlet("salt"))
    assert salt_out == pytest.approx(50.0 + slope_mM_per_s * 1200.0, rel=1e-3)


# --- linear isotherm ----------------------------------------------------------------------------


def test_linear_isotherm_first_moment_grm():
    # GRM, rapid-equilibrium linear isotherm q = K c (q per solid volume):
    #   μ_1 = t_pulse/2 + (ε_c V_c / Q) [1 + F (ε_p + (1 − ε_p) K)],   F = (1 − ε_c)/ε_c
    # (Kucera 1965; Schmidt-Traub et al., Preparative Chromatography, 2nd ed. 2012, ch. 6).
    henry = 2.0
    project = make_physics_project(
        "GRM", isotherm={"model": "Linear", "parameters": {"tracer": {"k_eq": henry}}}
    )
    result = simulate(project)
    eps_c, eps_p = PHYSICS_COL_POROSITY, PHYSICS_PAR_POROSITY
    phase_ratio = (1 - eps_c) / eps_c
    interstitial_time_s = eps_c * PHYSICS_COLUMN_VOLUME_ML / FLOW_ML_PER_S
    expected = PULSE_S / 2 + interstitial_time_s * (1 + phase_ratio * (eps_p + (1 - eps_p) * henry))
    tracer = result.outlet("tracer")
    assert first_moment(result.time_s, tracer) == pytest.approx(expected, rel=MOMENT_RTOL)
    assert recovered_fraction(result, "tracer", 1.0 * PULSE_S) == pytest.approx(1.0, rel=MASS_RTOL)


def test_kinetic_linear_isotherm_has_same_first_moment_as_equilibrium():
    # Finite binding kinetics broaden the peak (second moment) but do not move its first moment.
    equilibrium = simulate(
        make_physics_project(
            "LRM", isotherm={"model": "Linear", "parameters": {"tracer": {"k_eq": 2.0}}}
        )
    )
    kinetic = simulate(
        make_physics_project(
            "LRM",
            isotherm={"model": "Linear", "parameters": {"tracer": {"k_eq": 2.0, "k_kin": 0.5}}},
        )
    )
    assert first_moment(kinetic.time_s, kinetic.outlet("tracer")) == pytest.approx(
        first_moment(equilibrium.time_s, equilibrium.outlet("tracer")), rel=MOMENT_RTOL
    )


# --- Langmuir breakthrough ----------------------------------------------------------------------

LANGMUIR_K_EQ = 1.0  # m³/mol (= 1/mM)
LANGMUIR_Q_MAX = 10.0  # mol/m³ solid
FEED_MM = 1.0
# q* = q_max K c / (1 + K c) = 10 × 1 × 1 / 2 = 5 mol/m³ solid
Q_STAR = LANGMUIR_Q_MAX * LANGMUIR_K_EQ * FEED_MM / (1 + LANGMUIR_K_EQ * FEED_MM)
LOAD_METHOD = [
    {"name": "load", "duration": "60 min", "flow": "1 mL/min", "inlet": {"protein": "1 mM"}}
]


def langmuir_project(column_model):
    return make_physics_project(
        column_model,
        components=["protein"],
        isotherm={
            "model": "Langmuir",
            "parameters": {"protein": {"k_eq": LANGMUIR_K_EQ, "q_max": LANGMUIR_Q_MAX}},
        },
        method=LOAD_METHOD,
    )


@pytest.mark.parametrize(
    ("column_model", "solid_fraction"),
    [
        ("LRM", 1 - PHYSICS_TOTAL_POROSITY),  # LRM: q per (1 − ε_t) of the column
        ("GRM", (1 - PHYSICS_COL_POROSITY) * (1 - PHYSICS_PAR_POROSITY)),  # particle skeleton
    ],
)
def test_langmuir_breakthrough(column_model, solid_fraction):
    result = simulate(langmuir_project(column_model))
    c_out = result.outlet("protein")

    # Plateau: once saturated, the outlet equals the feed.
    assert c_out[-1] == pytest.approx(FEED_MM, rel=1e-6)

    # Stoichiometric breakthrough time from a mass balance on the saturated column:
    #   t_s = ∫ (1 − c_out/c_feed) dt = (V_c / Q) [ε_t + solid_fraction · q*/c_feed]
    # (Guiochon et al. 2006, ch. 4).
    t_s = np.trapezoid(1 - c_out / FEED_MM, result.time_s)
    expected = (
        PHYSICS_COLUMN_VOLUME_ML
        / FLOW_ML_PER_S
        * (PHYSICS_TOTAL_POROSITY + solid_fraction * Q_STAR / FEED_MM)
    )
    assert t_s == pytest.approx(expected, rel=MOMENT_RTOL)


# --- SMA load–wash–elute ------------------------------------------------------------------------


def test_sma_gradient_elution_mass_balance():
    # SPEC §4.3 project with illustrative transport parameters on a coarse grid: all loaded
    # protein elutes in the 50 → 500 mM gradient.
    from project_examples import make_project

    project = make_project()
    project.column.axial_dispersion = "2e-7 m^2/s"
    project.column.film_diffusion = {"salt": "1e-4 m/s", "mAb": "2e-6 m/s"}
    project.column.pore_diffusion = {"salt": "1e-10 m^2/s", "mAb": "5e-12 m^2/s"}
    project.solver.n_col, project.solver.n_par = 40, 10
    result = simulate(project)

    load = project.resolved_method()[1]
    injected = load.inlet_start_mol_per_m3[1] * load.duration_s  # mol/m³ · s (constant flow)
    assert recovered_fraction(result, "mAb", injected) == pytest.approx(1.0, rel=MASS_RTOL)
    # No significant breakthrough during load and wash: at 50 mM salt the protein binds strongly.
    # A small numerical/kinetic tail (~5e-5 of the feed on this coarse grid) is expected; leakage
    # above 0.1 % of the feed would mean the salt dependence or the capacity is mapped wrongly.
    wash_end_s = project.resolved_method()[2].end_s
    leakage = result.outlet("mAb")[result.time_s < wash_end_s].max()
    assert leakage < 1e-3 * load.inlet_start_mol_per_m3[1]
