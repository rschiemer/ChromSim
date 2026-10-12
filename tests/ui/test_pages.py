"""Smoke tests for the web UI (SPEC §8): every page opens, and a simulation runs as a job.

Uses NiceGUI's simulated user (no browser). Skipped where NiceGUI is not installed.
"""

import pytest

pytest.importorskip("nicegui")

from nicegui.testing import User, user_simulation  # noqa: E402
from project_examples import make_physics_project, make_project  # noqa: E402

from chromsim import api  # noqa: E402

pytestmark = [pytest.mark.ui, pytest.mark.anyio]

PAGES = [
    ("/", "Check and save"),
    ("/data", "Upload data files"),
    ("/simulate", "Simulate"),
    ("/calibrate", "Start fit"),
    ("/optimize", "Start optimization"),
    ("/results", "Jobs"),
]


@pytest.fixture
def anyio_backend():
    return "asyncio"


@pytest.fixture
def workspace_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("NICEGUI_STORAGE_PATH", str(tmp_path / "nicegui"))
    root = tmp_path / "workspace"
    storage = api.FileStorage(root)
    storage.save_project("lge", make_project())
    storage.save_project("tracer", make_physics_project("LRM", n_col=20))
    return root


@pytest.fixture
def runner(workspace_dir):
    runner = api.LocalProcessRunner(api.FileStorage(workspace_dir), max_workers=1)
    yield runner
    runner.shutdown()


def set_value(user, marker, value):
    """Set the value of the (first) element with ``marker``, as a user selecting or typing it."""
    element = min(user.find(marker=marker).elements, key=lambda e: e.id)
    with user.client:
        element.value = value


async def simulated_user(workspace_dir, runner):
    from chromsim.ui.app import create_root

    return user_simulation(root=create_root(workspace_dir, runner=runner))


@pytest.mark.parametrize(("path", "expected"), PAGES)
async def test_page_opens(workspace_dir, runner, path, expected):
    async with await simulated_user(workspace_dir, runner) as user:
        user: User
        await user.open(path)
        await user.should_see(expected)
        await user.should_see(marker="show-code")
        await user.should_see(marker="download-yaml")


async def test_show_code_dialog_on_every_page(workspace_dir, runner):
    async with await simulated_user(workspace_dir, runner) as user:
        for path, _ in PAGES:
            await user.open(path)
            user.find(marker="show-code").click()
            await user.should_see("Equivalent Python code")


async def test_project_page_checks_and_saves(workspace_dir, runner):
    async with await simulated_user(workspace_dir, runner) as user:
        await user.open("/")
        user.find(marker="save-project").click()
        await user.should_see("Saved project")


@pytest.mark.slow
async def test_simulation_job_from_the_simulate_page(workspace_dir, runner):
    async with await simulated_user(workspace_dir, runner) as user:
        await user.open("/simulate")
        set_value(user, "project-select", "tracer")
        set_value(user, "result-name", "ui-result")
        user.find(marker="run-simulation").click()
        await user.should_see(marker="chromatogram", retries=600)
    assert "ui-result" in api.FileStorage(workspace_dir).list_simulation_results()
