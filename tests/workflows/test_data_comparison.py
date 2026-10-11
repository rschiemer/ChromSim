"""End-to-end: imported fixtures agree with a simulation using the parameters that generated them.

The fixtures in tests/data/ were generated (tests/data/generate_fixtures.py) from the tutorial
project with mAb k_eq = 1.5e-2, plus 3–4 % assay noise and UV/conductivity noise. Simulating with
those parameters must therefore reproduce every readout to within the noise, and the first-guess
project (k_eq = 1.0e-2) must not. This checks importers, unit handling, fraction matching, fraction
averaging and detector predictions together.
"""

from pathlib import Path

import numpy as np
import pytest

from chromsim import api

ROOT = Path(__file__).parents[2]
DATA = ROOT / "tests" / "data"
TRUE_MAB_K_EQ = 1.5e-2


def rms_relative_deviation(comparison, analyte):
    table = comparison.table[
        (comparison.table["analyte"] == analyte) & ~comparison.table["below_loq"]
    ]
    return float(
        np.sqrt(np.mean(((table["simulated"] - table["measured"]) / table["measured"].max()) ** 2))
    )


@pytest.fixture(scope="module")
def data():
    experiment = api.read_chromatogram(DATA / "akta" / "cex_mab_hmw_run01.csv")
    total = api.read_assay(
        DATA / "assays" / "total_protein_a280.csv",
        ROOT / "examples/assay_profiles/total_protein_a280.yaml",
    )
    sec = api.read_assay(
        DATA / "assays" / "sec_hplc.xlsx", ROOT / "examples/assay_profiles/sec_hplc_long.yaml"
    )
    return experiment, total, sec


def simulate(mab_k_eq):
    project = api.load_project(ROOT / "examples/projects/cex_mab_hmw_lge.yaml")
    project.isotherm.parameters["mAb"]["k_eq"] = mab_k_eq
    return project, api.simulate(project)


@pytest.mark.slow
def test_true_parameters_reproduce_all_data_within_noise(data):
    experiment, total, sec = data
    project, result = simulate(TRUE_MAB_K_EQ)

    # Assays: generated with 3 % (A280) and 4 % (SEC) relative noise; the coarser project grid
    # (50×12 vs 100×25 for the fixtures) adds well under 1 %. 5 % of the peak value is a safe bound.
    total_cmp = api.compare_with_assay(result, total)
    sec_cmp = api.compare_with_assay(result, sec, fractions=experiment.fractions)
    assert rms_relative_deviation(total_cmp, "total_protein") < 0.05
    assert rms_relative_deviation(sec_cmp, "monomer") < 0.05
    assert rms_relative_deviation(sec_cmp, "hmw") < 0.05
    assert sec_cmp.assay_fractions_not_collected == ["Load", "A99"]

    # UV: noise 0.2 mAU on a ~390 mAU peak, baseline removed; peak positions agree within 0.3 mL.
    uv = api.subtract_baseline(experiment.signal("uv"), (2e-6, 15e-6), (140e-6, 147e-6))
    uv_sim = api.predict_uv(result, project)
    assert uv.x[np.argmax(uv.values)] == pytest.approx(
        uv_sim.x[np.argmax(uv_sim.values)], abs=0.3e-6
    )
    assert uv.values.max() == pytest.approx(uv_sim.values.max(), rel=0.02)

    # Conductivity → salt through the calibration reproduces the simulated salt within 2 mM
    # (noise 0.03 mS/cm ≈ 0.4 mM, plus interpolation of the calibration).
    calibration = api.load_conductivity_calibration(
        ROOT / "examples/calibrations/nacl_conductivity.yaml"
    )
    salt = api.conductivity_to_salt(experiment.signal("conductivity"), calibration)
    simulated_salt = np.interp(
        salt.x,
        result.volume_m3 + project.system.conductivity_delay_volume_m3,
        result.outlet("salt"),
    )
    gradient = (salt.x > 60e-6) & (salt.x < 125e-6)
    assert np.max(np.abs(salt.values[gradient] - simulated_salt[gradient])) < 2.0


@pytest.mark.slow
def test_first_guess_parameters_are_detectably_wrong(data):
    _, total, _ = data
    _, result = simulate(1.0e-2)
    assert rms_relative_deviation(api.compare_with_assay(result, total), "total_protein") > 0.15
