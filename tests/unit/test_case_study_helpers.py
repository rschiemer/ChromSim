"""Helpers added for published data sets: time → volume axis, method inlet profiles and
concentration targets for inverse fits."""

from pathlib import Path

import numpy as np
import pytest

from chromsim import api
from chromsim.errors import ChromSimError

EXAMPLES = Path(__file__).parents[2] / "examples"
ML = 1e-6


@pytest.fixture(scope="module")
def project():
    return api.load_project(EXAMPLES / "projects" / "cex_mab_lge.yaml")


def test_time_to_volume():
    signal = api.Signal(
        "c", "concentration", "time", np.array([0.0, 60.0]), np.ones(2), "mol/m**3", "mM"
    )
    flow = 0.5 * ML / 60.0  # 0.5 mL/min
    on_volume = api.time_to_volume(signal, flow, start_volume_m3=1 * ML)
    assert on_volume.x_axis == "volume"
    assert np.allclose(on_volume.x, [1 * ML, 1.5 * ML])
    with pytest.raises(ChromSimError, match="already on the volume axis"):
        api.time_to_volume(on_volume, flow)
    with pytest.raises(ChromSimError, match="must be positive"):
        api.time_to_volume(signal, 0.0)


def test_inlet_profile_follows_the_method(project):
    profile = api.inlet_profile(project, "salt", delay_volume_m3=2 * ML)
    steps = project.resolved_method()
    gradient = steps[3]
    start = sum(s.flow_rate_m3_per_s * s.duration_s for s in steps[:3]) + 2 * ML
    end = start + gradient.flow_rate_m3_per_s * gradient.duration_s
    assert profile.kind == "salt" and profile.unit == "mol/m**3"
    assert np.interp(start - 0.1 * ML, profile.x, profile.values) == pytest.approx(50.0)
    assert np.interp(0.5 * (start + end), profile.x, profile.values) == pytest.approx(275.0)
    assert np.interp(end, profile.x, profile.values) == pytest.approx(500.0)
    with pytest.raises(ChromSimError, match="Unknown component"):
        api.inlet_profile(project, "HCP")


def test_concentration_target_compares_component_outlet():
    coarse = api.load_project(EXAMPLES / "projects" / "cex_mab_hmw_lge.yaml")
    coarse.solver.n_col, coarse.solver.n_par = 20, 4
    result = api.simulate(coarse)
    measured = api.Signal(
        "mAb",
        "concentration",
        "volume",
        result.volume_m3 + coarse.system.uv_delay_volume_m3,
        result.outlet("mAb", "mM"),
        "mol/m**3",
        "mM",
    )
    target = api.SignalTarget(measured, component="mAb")
    assert np.allclose(target.residuals(result, coarse), 0.0)
    in_g_per_l = api.Signal(
        "mAb",
        "concentration",
        "volume",
        measured.x,
        result.outlet("mAb", "g/L"),
        "kg/m**3",
        "g/L",
    )
    assert np.allclose(api.SignalTarget(in_g_per_l, component="mAb").residuals(result, coarse), 0)
    with pytest.raises(ChromSimError, match="need component"):
        api.SignalTarget(measured).residuals(result, coarse)
