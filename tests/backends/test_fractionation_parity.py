"""ChromSim's pooling KPIs equal CADET-Process's Fractionator on the same CADET run (ADR 0010).

ChromSim fractionates its backend-independent SimulationResult itself; this test proves the
numbers agree with CADET-Process for identical cut points. Conventions differ only by fixed
factors, which the test converts:

- productivity: ChromSim per packed-bed volume V_c, CADET-Process per solid volume
  V_solid = (1 − ε_t) V_c;
- buffer: ChromSim V_pumped / n_pool, CADET-Process n_pool / V_eluent (the inverse);
- purity: CADET-Process by amount of substance; ChromSim by mass when molecular weights are
  known, so purity is compared from the pooled amounts.
"""

import tempfile
from pathlib import Path

import numpy as np
import pytest
from CADETProcess.fractionation import Fractionator
from CADETProcess.simulator import Cadet

from chromsim import api
from chromsim.analysis.fractionation import Pool, pool_performance
from chromsim.backends.cadet_backend import build_process

# Both sides integrate the same 1 s outlet samples with the trapezoidal rule; CADET-Process may
# interpolate cut points slightly differently, so 0.1 % covers the remaining difference.
RTOL = 1e-3


@pytest.fixture(scope="module")
def runs():
    project = api.load_project(Path(__file__).parents[2] / "examples/projects/cex_mab_hmw_lge.yaml")
    project.solver.n_col, project.solver.n_par = 20, 4
    result = api.simulate(project)

    process = build_process(project)
    flow_sheet = process.flow_sheet
    flow_sheet.add_feed_inlet("inlet")
    flow_sheet.add_eluent_inlet("inlet")
    flow_sheet.add_product_outlet("outlet")
    with tempfile.TemporaryDirectory() as folder:
        cadet_results = Cadet(temp_dir=Path(folder)).simulate(process)
    return project, result, cadet_results


def test_pool_kpis_equal_cadet_process(runs):
    project, result, cadet_results = runs
    start_s, end_s = 11500.0, 13300.0  # pool across the mAb elution peak (max at ~12400 s)
    ours = pool_performance(result, Pool("mAb", start_s, end_s), impurities=["HMW"])

    fractionator = Fractionator(cadet_results, components=["mAb", "HMW"])
    fractionator.add_fractionation_event("start", "mAb", start_s)
    fractionator.add_fractionation_event("end", -1, end_s)
    mab = 0  # index in the fractionator's component list

    assert ours.amounts_mol["mAb"] == pytest.approx(fractionator.mass[mab], rel=RTOL)
    assert ours.recovery == pytest.approx(fractionator.recovery[mab], rel=RTOL)
    v_c, eps_t = project.column.volume_m3, project.resin.total_porosity
    assert ours.productivity_mol_per_m3_s * v_c == pytest.approx(
        fractionator.productivity[mab] * (1 - eps_t) * v_c, rel=RTOL
    )
    assert 1.0 / ours.buffer_consumption_m3_per_mol == pytest.approx(
        fractionator.eluent_consumption[mab], rel=RTOL
    )
    molar_purity = ours.amounts_mol["mAb"] / (ours.amounts_mol["mAb"] + ours.amounts_mol["HMW"])
    assert molar_purity == pytest.approx(fractionator.purity[mab], rel=RTOL)
    assert np.isfinite(ours.purity) and 0 < ours.purity <= 1
