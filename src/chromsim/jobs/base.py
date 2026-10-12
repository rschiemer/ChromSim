"""The ``JobRunner`` protocol for long jobs (fits, optimizations, simulations) (SPEC §8.1).

A job is a module-level function ``task(context, *args) -> str`` that stores its outputs in the
workspace :class:`~chromsim.storage.FileStorage` and returns the name under which it stored them.
It reports progress with ``context.report(fraction, message)``; that call raises
:class:`JobCancelled` once a cancel was requested, so long jobs stop at their next checkpoint.

Job state lives in the storage (``records/jobs/<id>.json``), not in memory: any process can read
it, a browser refresh loses nothing, and a queue-backed runner (Celery/RQ/arq) can later replace
:class:`~chromsim.jobs.local.LocalProcessRunner` without changing callers.
"""

from collections.abc import Callable
from dataclasses import asdict, dataclass
from typing import Literal, Protocol

from chromsim.errors import ChromSimError
from chromsim.storage.file_storage import FileStorage, stamp

JobState = Literal["queued", "running", "done", "failed", "cancelled"]
FINISHED_STATES = ("done", "failed", "cancelled")
JOB_RECORDS = "jobs"
CANCEL_RECORDS = "job-cancel"


class JobCancelled(ChromSimError):
    """Raised inside a job when a cancel was requested."""


@dataclass
class JobStatus:
    """State of one job.

    Attributes
    ----------
    job_id, kind : str
    state : {"queued", "running", "done", "failed", "cancelled"}
    progress : float
        Fraction done, 0–1 (best estimate).
    message : str
        Latest progress or error message.
    result : str or None
        Name of the stored output (when done).
    created_at, created_by : str
    """

    job_id: str
    kind: str
    state: JobState
    progress: float = 0.0
    message: str = ""
    result: str | None = None
    created_at: str = ""
    created_by: str = ""

    @property
    def finished(self) -> bool:
        """``True`` when done, failed or cancelled."""
        return self.state in FINISHED_STATES


def save_status(storage: FileStorage, status: JobStatus) -> None:
    storage.save_record(JOB_RECORDS, status.job_id, asdict(status))


def load_status(storage: FileStorage, job_id: str) -> JobStatus:
    return JobStatus(**storage.load_record(JOB_RECORDS, job_id))


@dataclass(frozen=True)
class JobContext:
    """Handed to every task: its job id and workspace, and the way to report progress."""

    job_id: str
    storage_root: str

    @property
    def storage(self) -> FileStorage:
        return FileStorage(self.storage_root)

    def cancel_requested(self) -> bool:
        return self.storage.has_record(CANCEL_RECORDS, self.job_id)

    def report(self, progress: float, message: str = "") -> None:
        """Record progress (0–1) and a message; raises :class:`JobCancelled` if cancelled."""
        if self.cancel_requested():
            raise JobCancelled(f"Job {self.job_id} was cancelled.")
        storage = self.storage
        status = load_status(storage, self.job_id)
        status.progress = float(min(max(progress, 0.0), 1.0))
        status.message = message
        save_status(storage, status)


def new_status(job_id: str, kind: str) -> JobStatus:
    provenance = stamp()
    return JobStatus(
        job_id=job_id,
        kind=kind,
        state="queued",
        message="waiting",
        created_at=provenance["created_at"],
        created_by=provenance["created_by"],
    )


class JobRunner(Protocol):
    """Runs jobs in the background."""

    def submit(self, kind: str, task: Callable[..., str], *args) -> str:
        """Start ``task(context, *args)``; returns the job id immediately."""
        ...

    def status(self, job_id: str) -> JobStatus: ...

    def cancel(self, job_id: str) -> None: ...

    def result(self, job_id: str) -> str:
        """Name of the stored output of a finished job (raises JobError otherwise)."""
        ...

    def list_jobs(self) -> list[JobStatus]: ...

    def shutdown(self) -> None: ...
