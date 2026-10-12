"""A ``Workspace``: everything one user works on, passed explicitly (SPEC §8.1).

There is no global state: the web UI creates one Workspace per browser session (they may share a
folder and a job runner), notebooks create their own. All data live in the workspace's storage,
so nothing is lost when a browser is refreshed or a session ends.
"""

from pathlib import Path

from chromsim.jobs.base import JobRunner, JobStatus
from chromsim.jobs.local import LocalProcessRunner
from chromsim.jobs.tasks import FitSetup, ProcessSetup, fit_task, process_task, simulation_task
from chromsim.storage.file_storage import FileStorage


class Workspace:
    """A workspace folder with its storage and a job runner.

    Parameters
    ----------
    root : str or Path
        Workspace folder (created if missing).
    runner : JobRunner, optional
        Shared runner; by default a :class:`LocalProcessRunner` for this workspace.
    max_workers : int
        Parallel jobs of the default runner.

    Attributes
    ----------
    storage : FileStorage
    runner : JobRunner
    current_project : str or None
        The project selected in this session.
    """

    def __init__(self, root: str | Path, runner: JobRunner | None = None, max_workers: int = 2):
        self.storage = FileStorage(root)
        self.runner = (
            runner if runner is not None else LocalProcessRunner(self.storage, max_workers)
        )
        self.current_project: str | None = None

    def submit_simulation(self, project_name: str, result_name: str | None = None) -> str:
        """Simulate a stored project in the background; returns the job id."""
        return self.runner.submit(
            "simulation", simulation_task, project_name, result_name or project_name
        )

    def submit_fit(self, setup: FitSetup) -> str:
        """Run an inverse fit in the background; returns the job id."""
        return self.runner.submit("fit", fit_task, setup.model_dump(mode="json"))

    def submit_process_optimization(self, setup: ProcessSetup) -> str:
        """Run a process optimization in the background; returns the job id."""
        return self.runner.submit("process", process_task, setup.model_dump(mode="json"))

    def job(self, job_id: str) -> JobStatus:
        """Current status of a job."""
        return self.runner.status(job_id)

    def result(self, job_id: str) -> str:
        """Name of the stored output of a finished job (raises JobError otherwise)."""
        return self.runner.result(job_id)

    def cancel(self, job_id: str) -> None:
        """Ask a job to stop."""
        self.runner.cancel(job_id)

    def jobs(self) -> list[JobStatus]:
        """All jobs of this workspace, newest first."""
        return self.runner.list_jobs()

    def close(self) -> None:
        """Stop the default runner's worker processes."""
        self.runner.shutdown()
