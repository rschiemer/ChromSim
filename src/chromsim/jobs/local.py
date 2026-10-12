"""``LocalProcessRunner``: jobs in a pool of local worker processes (SPEC §8.1, v1)."""

import logging
import multiprocessing
import traceback
import uuid
from collections.abc import Callable
from concurrent.futures import Future, ProcessPoolExecutor

from chromsim.errors import JobError
from chromsim.jobs.base import (
    CANCEL_RECORDS,
    JOB_RECORDS,
    JobCancelled,
    JobContext,
    JobStatus,
    load_status,
    new_status,
    save_status,
)
from chromsim.storage.file_storage import FileStorage

logger = logging.getLogger(__name__)


def execute(task: Callable[..., str], context: JobContext, args: tuple) -> str | None:
    """Run one job in a worker process and record how it ended."""
    storage = context.storage
    status = load_status(storage, context.job_id)
    status.state, status.message = "running", "started"
    save_status(storage, status)
    try:
        result = task(context, *args)
    except JobCancelled:
        status = load_status(storage, context.job_id)
        status.state, status.message = "cancelled", "cancelled"
        save_status(storage, status)
        return None
    except Exception as error:  # report every failure in the job record instead of losing it
        logger.error("Job %s failed:\n%s", context.job_id, traceback.format_exc())
        status = load_status(storage, context.job_id)
        status.state, status.message = "failed", f"{type(error).__name__}: {error}"
        save_status(storage, status)
        return None
    status = load_status(storage, context.job_id)
    status.state, status.progress, status.result = "done", 1.0, result
    status.message = "finished"
    save_status(storage, status)
    return result


class LocalProcessRunner:
    """Runs jobs in ``max_workers`` local processes (spawn context: Windows/macOS safe).

    Parameters
    ----------
    storage : FileStorage
        Workspace where job records and outputs are kept.
    max_workers : int
        Jobs running at the same time.
    """

    def __init__(self, storage: FileStorage, max_workers: int = 2):
        self.storage = storage
        self.max_workers = max_workers
        self._executor: ProcessPoolExecutor | None = None
        self._futures: dict[str, Future] = {}

    def submit(self, kind: str, task: Callable[..., str], *args) -> str:
        job_id = f"{kind}-{uuid.uuid4().hex[:8]}"
        save_status(self.storage, new_status(job_id, kind))
        if self._executor is None:
            context = multiprocessing.get_context("spawn")
            self._executor = ProcessPoolExecutor(self.max_workers, mp_context=context)
        job_context = JobContext(job_id=job_id, storage_root=str(self.storage.root))
        self._futures[job_id] = self._executor.submit(execute, task, job_context, args)
        logger.info("Submitted job %s", job_id)
        return job_id

    def status(self, job_id: str) -> JobStatus:
        status = load_status(self.storage, job_id)
        future = self._futures.get(job_id)
        if future is not None and future.done() and not status.finished:
            # The worker died before it could record the outcome (e.g. killed process).
            error = future.exception() if not future.cancelled() else None
            status.state = "failed" if error else "cancelled"
            status.message = f"{type(error).__name__}: {error}" if error else "cancelled"
            save_status(self.storage, status)
        return status

    def cancel(self, job_id: str) -> None:
        """Request cancellation: queued jobs stop at once, running ones at their next report."""
        self.storage.save_record(CANCEL_RECORDS, job_id, {"requested": True})
        future = self._futures.get(job_id)
        if future is not None and future.cancel():
            status = load_status(self.storage, job_id)
            status.state, status.message = "cancelled", "cancelled before start"
            save_status(self.storage, status)

    def result(self, job_id: str) -> str:
        """Name of the stored output of a job that is done (e.g. a simulation result name).

        Raises
        ------
        JobError
            If the job is still queued or running, failed, or was cancelled.
        """
        status = self.status(job_id)
        if status.state != "done" or status.result is None:
            raise JobError(
                f"Job {job_id} has no result: it is {status.state} ({status.message}). "
                "Check `status(job_id)` until it is finished."
            )
        return status.result

    def list_jobs(self) -> list[JobStatus]:
        statuses = [self.status(job_id) for job_id in self.storage.list_records(JOB_RECORDS)]
        return sorted(statuses, key=lambda s: s.created_at, reverse=True)

    def shutdown(self) -> None:
        if self._executor is not None:
            self._executor.shutdown(wait=False, cancel_futures=True)
            self._executor = None
