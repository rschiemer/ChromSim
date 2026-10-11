"""Tracer analysis on exact synthetic peaks (no CADET): every parameter is recovered exactly."""

import numpy as np
import pytest
from project_examples import make_project
from synthetic import rectangular_pulse, through_system

from chromsim.data.experiment import Signal
from chromsim.errors import ChromSimError
from chromsim.estimation.workflows.tracer import (
    TracerRun,
    analyse_tracers,
    peclet_from_variance_ratio,
    variance_ratio,
)

ML = 1e-6
COLUMN = make_project().column  # 20 cm × 0.5 cm, V_c = 3.92699 mL
FLOW = 8.181231e-9  # m³/s (150 cm/h)
VOLUME = np.linspace(0.0, 15.0, 30001) * ML
INJECTION = 0.05 * ML
# The synthetic system pulse resolves its edges to half a grid step (0.00025 mL), which shifts the
# system volume and hence the porosities by ~2e-4 relative; 1e-3 covers that with margin.
POROSITY_RTOL = 1e-3
TRUE = {"mixer": 0.2 * ML, "tubing": 0.25 * ML, "eps_c": 0.36, "eps_t": 0.904, "peclet": 800.0}


def gaussian_signal(name, mean, variance):
    values = np.exp(-0.5 * (VOLUME - mean) ** 2 / variance)
    return Signal(name, "uv", "volume", VOLUME, values, "AU", "mAU")


def runs():
    """System pulse (exact) and two column peaks with exactly known moments."""
    system_signal = through_system(
        VOLUME, rectangular_pulse(VOLUME, 1 * ML, INJECTION), TRUE["mixer"], TRUE["tubing"]
    )
    system_variance = TRUE["mixer"] ** 2
    v_sys = TRUE["mixer"] + TRUE["tubing"]
    retention_np = TRUE["eps_c"] * COLUMN.volume_m3
    variance_np = variance_ratio(TRUE["peclet"]) * retention_np**2
    retention_p = TRUE["eps_t"] * COLUMN.volume_m3
    injection_terms = INJECTION / 2, INJECTION**2 / 12
    return [
        TracerRun(
            "system",
            "system",
            Signal("sys", "uv", "volume", VOLUME, system_signal, "AU", "mAU"),
            1 * ML,
            INJECTION,
        ),
        TracerRun(
            "dextran",
            "non_penetrating",
            gaussian_signal(
                "dex",
                1 * ML + v_sys + retention_np + injection_terms[0],
                system_variance + variance_np + injection_terms[1],
            ),
            1 * ML,
            INJECTION,
        ),
        TracerRun(
            "acetone",
            "pore_penetrating",
            gaussian_signal(
                "ace",
                1 * ML + v_sys + retention_p + injection_terms[0],
                system_variance + 4e-15 + injection_terms[1],
            ),
            1 * ML,
            INJECTION,
        ),
    ]


def test_peclet_inversion():
    for peclet in (0.5, 10.0, 800.0, 1e5):
        assert peclet_from_variance_ratio(variance_ratio(peclet)) == pytest.approx(peclet, rel=1e-8)
    with pytest.raises(ChromSimError, match="outside"):
        peclet_from_variance_ratio(1.5)


def test_all_parameters_recovered():
    analysis = analyse_tracers(runs(), COLUMN, FLOW)
    assert analysis.system_volume_m3 == pytest.approx(0.45 * ML, rel=1e-3)
    assert analysis.mixer_volume_m3 == pytest.approx(0.2 * ML, rel=5e-3)
    assert analysis.tubing_volume_m3 == pytest.approx(0.25 * ML, rel=5e-3)
    assert analysis.col_porosity == pytest.approx(0.36, rel=POROSITY_RTOL)
    assert analysis.total_porosity == pytest.approx(0.904, rel=POROSITY_RTOL)
    # ε_p = (0.904 − 0.36)/(1 − 0.36) = 0.85
    assert analysis.par_porosity == pytest.approx(0.85, rel=1e-3)
    assert analysis.peclet == pytest.approx(800.0, rel=1e-2)
    # D_ax = u_int L / Pe = Q/(A ε_c) × 0.2 m / 800
    u_int = FLOW / (COLUMN.cross_section_m2 * 0.36)
    assert analysis.axial_dispersion_m2_per_s == pytest.approx(u_int * 0.2 / 800.0, rel=1e-2)
    dextran = next(peak for peak in analysis.peaks if peak.name == "dextran")
    assert dextran.asymmetry == pytest.approx(1.0, rel=1e-2)
    assert dextran.hetp_m == pytest.approx(0.2 / (1 / variance_ratio(800.0)), rel=1e-2)


def test_apply_to_project():
    project = analyse_tracers(runs(), COLUMN, FLOW).apply_to(make_project())
    assert project.resin.col_porosity == pytest.approx(0.36, rel=POROSITY_RTOL)
    assert project.resin.par_porosity == pytest.approx(0.85, rel=1e-3)
    assert project.system.mixer_volume.unit == "mL"
    assert project.system.mixer_volume.magnitude == pytest.approx(0.2, rel=5e-3)
    assert project.column.axial_dispersion.unit == "m^2/s"


def test_column_runs_need_system_volume():
    with pytest.raises(ChromSimError, match="need the system volume"):
        analyse_tracers(runs()[1:], COLUMN, FLOW)
    analysis = analyse_tracers(runs()[1:], COLUMN, FLOW, system_volume_m3=0.45 * ML)
    assert analysis.col_porosity == pytest.approx(0.36, rel=1e-3)
    assert analysis.mixer_volume_m3 is None


def test_system_dominated_variance_is_flagged():
    # TRUE peclet 800 on this small column: the column variance is ~6 % of the system variance.
    analysis = analyse_tracers(runs(), COLUMN, FLOW)
    assert any("D_ax from moments is therefore unreliable" in note for note in analysis.notes)


def test_negative_column_variance_gives_note_not_crash():
    broken = runs()
    dextran = broken[1]
    narrow = gaussian_signal("dex", dextran.signal.x[np.argmax(dextran.signal.values)], 1e-18)
    broken[1] = TracerRun(
        "dextran", "non_penetrating", narrow, dextran.injection_start_m3, INJECTION
    )
    analysis = analyse_tracers(broken, COLUMN, FLOW)
    assert analysis.axial_dispersion_m2_per_s is None
    assert any("No axial dispersion from moments" in note for note in analysis.notes)


def test_apply_to_subtracts_uv_delay_from_tubing():
    project = make_project()
    project.system.uv_delay_volume = "0.05 mL"
    updated = analyse_tracers(runs(), COLUMN, FLOW).apply_to(project)
    # measured: mixer 0.2 + "tubing" 0.25 (incl. the 0.05 mL UV delay) → tubing 0.20 mL
    assert updated.system.tubing_volume.magnitude == pytest.approx(0.20, rel=1e-2)
    project.system.uv_delay_volume = "0.5 mL"
    with pytest.raises(ChromSimError, match="check system.uv_delay_volume"):
        analyse_tracers(runs(), COLUMN, FLOW).apply_to(project)
