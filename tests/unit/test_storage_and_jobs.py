"""Workspace storage, the local job runner and the Workspace (SPEC §8.1)."""

import time

import numpy as np
import pytest
from job_tasks import cancellable_task, failing_task, quick_task
from project_examples import make_physics_project

from chromsim.errors import ChromSimError
from chromsim.jobs import LocalProcessRunner
from chromsim.storage import FileStorage
from chromsim.workspace import Workspace


def wait_until_finished(runner, job_id, timeout_s=120):
    deadline = time.monotonic() + timeout_s
    while time.monotonic() < deadline:
        status = runner.status(job_id)
        if status.finished:
            return status
        time.sleep(0.2)
    raise TimeoutError(job_id)


# --- storage --------------------------------------------------------------------------------------


def test_projects_and_records_round_trip(tmp_path):
    storage = FileStorage(tmp_path / "ws")
    project = make_physics_project("LRM", n_col=20)
    storage.save_project("tracer run", project)
    assert storage.list_projects() == ["tracer run"]
    assert storage.load_project("tracer run") == project
    storage.save_record("notes", "a", {"x": 1})
    assert storage.load_record("notes", "a") == {"x": 1}
    assert storage.list_records("notes") == ["a"]
    assert storage.load_record("provenance", "project tracer run")["created_by"]


def test_raw_files_are_never_overwritten(tmp_path):
    storage = FileStorage(tmp_path)
    path = storage.save_raw_file("run.csv", b"original")
    assert storage.save_raw_file("run.csv", b"original") == path  # same content: fine
    with pytest.raises(ChromSimError, match="never overwritten"):
        storage.save_raw_file("run.csv", b"changed")
    assert path.read_bytes() == b"original"
    assert storage.list_raw_files() == ["run.csv"]


@pytest.mark.parametrize("name", ["../escape", "", "a/b", ".hidden"])
def test_unsafe_names_are_rejected(tmp_path, name):
    with pytest.raises(ChromSimError, match="Invalid name"):
        FileStorage(tmp_path).save_record("notes", name, {})


def test_missing_items_have_clear_errors(tmp_path):
    storage = FileStorage(tmp_path)
    with pytest.raises(ChromSimError, match="No project 'x'"):
        storage.load_project("x")
    with pytest.raises(ChromSimError, match="No simulation result 'x'"):
        storage.load_simulation_result("x")


# --- job runner (real worker processes) -----------------------------------------------------------


@pytest.fixture
def runner(tmp_path):
    runner = LocalProcessRunner(FileStorage(tmp_path / "ws"), max_workers=2)
    yield runner
    runner.shutdown()


@pytest.mark.slow  # spawns worker processes (~5 s start-up each)
def test_job_runs_and_reports(runner):
    job_id = runner.submit("demo", quick_task, "out-1")
    status = wait_until_finished(runner, job_id)
    assert status.state == "done" and status.result == "out-1" and status.progress == 1.0
    assert runner.storage.load_record("outputs", "out-1") == {"ok": True}
    assert status.created_by and job_id.startswith("demo-")
    assert [s.job_id for s in runner.list_jobs()] == [job_id]
    assert runner.result(job_id) == "out-1"


@pytest.mark.slow  # spawns worker processes (~5 s start-up each)
def test_failed_job_records_the_error(runner):
    status = wait_until_finished(runner, runner.submit("demo", failing_task))
    assert status.state == "failed"
    assert status.message == "ValueError: deliberate failure"
    with pytest.raises(ChromSimError, match="failed"):
        runner.result(status.job_id)


@pytest.mark.slow  # spawns worker processes (~5 s start-up each)
def test_running_job_can_be_cancelled(runner):
    job_id = runner.submit("demo", cancellable_task, 60.0)
    deadline = time.monotonic() + 60
    while runner.status(job_id).state != "running" and time.monotonic() < deadline:
        time.sleep(0.1)
    runner.cancel(job_id)
    status = wait_until_finished(runner, job_id, timeout_s=30)
    assert status.state == "cancelled"


@pytest.mark.slow  # spawns worker processes (~5 s start-up each)
def test_job_status_survives_a_new_runner(runner, tmp_path):
    # State lives in storage: a new runner (e.g. after a browser refresh / app restart) sees it.
    job_id = runner.submit("demo", quick_task, "out-2")
    wait_until_finished(runner, job_id)
    other = LocalProcessRunner(FileStorage(tmp_path / "ws"))
    assert other.status(job_id).state == "done"


# --- workspace ------------------------------------------------------------------------------------


@pytest.mark.slow
def test_workspace_simulation_job(tmp_path):
    workspace = Workspace(tmp_path / "ws")
    try:
        workspace.storage.save_project("tracer", make_physics_project("LRM", n_col=20))
        job_id = workspace.submit_simulation("tracer", "tracer-result")
        status = wait_until_finished(workspace.runner, job_id)
        assert status.state == "done", status.message
        result = workspace.storage.load_simulation_result("tracer-result")
        assert np.trapezoid(result.outlet("tracer"), result.time_s) == pytest.approx(6.0, rel=1e-3)
    finally:
        workspace.close()


def test_fit_progress_counts_every_start():
    from chromsim.estimation.parameters import FitParameter
    from chromsim.jobs import FitSetup

    parameter = FitParameter(path="resin.col_porosity", lower=0.2, upper=0.6)
    setup = FitSetup(project="p", result_name="r", parameters=[parameter], max_evaluations=50)
    assert setup.expected_evaluations() == 50
    assert setup.model_copy(update={"multi_start": 3}).expected_evaluations() == 200
