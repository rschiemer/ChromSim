"""The U1–U6 workflow (SPEC §2) can be done entirely in the web UI (Phase 7 acceptance).

Drives the pages with NiceGUI's simulated user on the Tutorial 01–04 data: set up and simulate a
project (U1), import a chromatogram (U2), tracer analysis (U3), Yamamoto and an inverse fit (U4,
U5), and a process optimization (U6). Fit and optimization use very few evaluations; their
numerical quality is tested in tests/workflows.
"""

from pathlib import Path

import pytest

pytest.importorskip("nicegui")

from nicegui.testing import user_simulation  # noqa: E402

from chromsim import api  # noqa: E402

pytestmark = [pytest.mark.ui, pytest.mark.anyio, pytest.mark.slow]

ROOT = Path(__file__).parents[2]
AKTA = ROOT / "tests" / "data" / "akta"
EXAMPLES = ROOT / "examples"
JOB_WAIT = 3000  # × 0.1 s


@pytest.fixture
def anyio_backend():
    return "asyncio"


@pytest.fixture
def workspace_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("NICEGUI_STORAGE_PATH", str(tmp_path / "nicegui"))
    root = tmp_path / "workspace"
    storage = api.FileStorage(root)
    project = api.load_project(EXAMPLES / "projects" / "cex_mab_hmw_lge.yaml")
    project.solver.n_col, project.solver.n_par = 20, 4  # coarse: this tests the UI, not physics
    storage.save_project("hmw", project)
    files = [AKTA / f"{name}.csv" for name in ("tracer_system", "tracer_dextran", "tracer_acetone")]
    files += [AKTA / f"lge_{cv}cv.csv" for cv in (10, 20, 30, 40)]
    files += [AKTA / "cex_mab_hmw_run01.csv", EXAMPLES / "calibrations" / "nacl_conductivity.yaml"]
    for file in files:
        storage.save_raw_file(file.name, file.read_bytes())
    return root


def set_value(user, marker, value):
    element = min(user.find(marker=marker).elements, key=lambda e: e.id)
    with user.client:
        element.value = value


async def test_u1_to_u6_in_the_ui(workspace_dir):
    from chromsim.ui.app import create_root

    storage = api.FileStorage(workspace_dir)
    runner = api.LocalProcessRunner(storage, max_workers=1)
    try:
        async with user_simulation(root=create_root(workspace_dir, runner=runner)) as user:
            # U1: check and save the project, then simulate it
            await user.open("/")
            user.find(marker="save-project").click()
            await user.should_see("Saved project")
            await user.open("/simulate")
            set_value(user, "result-name", "hmw-sim")
            user.find(marker="run-simulation").click()
            await user.should_see(marker="chromatogram", retries=JOB_WAIT)

            # U2: import an ÄKTA chromatogram
            await user.open("/data")
            set_value(user, "chromatogram-select", "cex_mab_hmw_run01.csv")
            await user.should_see("fractions")

            # U3: tracer analysis, applied to a new project
            await user.open("/calibrate")
            set_value(user, "tracer-system", "tracer_system.csv")
            set_value(user, "tracer-excluded", "tracer_dextran.csv")
            set_value(user, "tracer-penetrating", "tracer_acetone.csv")
            user.find(marker="analyse-tracers").click()
            await user.should_see(marker="tracer-summary")
            user.find("Apply and save").click()
            await user.should_see("Saved project 'hmw-tracer'")

            # U4: Yamamoto from four LGE runs
            set_value(user, "lge-calibration", "nacl_conductivity.yaml")
            for index, cv in enumerate((10, 20, 30, 40)):
                set_value(user, f"lge-file-{index}", f"lge_{cv}cv.csv")
                set_value(user, f"lge-cv-{index}", cv)
            user.find(marker="yamamoto-fit").click()
            await user.should_see(marker="yamamoto-result")

            # U4/U5: inverse fit of k_eq to the UV trace, with uncertainty
            set_value(user, "fit-chromatogram", "cex_mab_hmw_run01.csv")
            set_value(user, "fit-uv-from", 60.0)
            set_value(user, "fit-uv-to", 140.0)
            set_value(user, "fit-isotherm.parameters.mAb.k_eq", True)
            set_value(user, "fit-max-evaluations", 4)
            user.find(marker="start-fit").click()
            await user.should_see(marker="fit-done", retries=JOB_WAIT)

            # U6: optimise the gradient length for recovery at 99 % purity
            await user.open("/optimize")
            set_value(user, "vary-method.3.duration", True)
            set_value(user, "process-max-evaluations", 3)
            user.find(marker="start-optimization").click()
            await user.should_see(marker="process-done", retries=JOB_WAIT)
    finally:
        runner.shutdown()

    assert {"hmw-tracer", "hmw-fit-fitted", "hmw-opt-best"} <= set(storage.list_projects())
    assert storage.list_records("fits") == ["hmw-fit"]
    assert storage.list_records("processes") == ["hmw-opt"]
