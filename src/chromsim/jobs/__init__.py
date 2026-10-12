"""``JobRunner`` protocol and the local process runner for long jobs (SPEC §8.1)."""

from chromsim.jobs.base import JobCancelled, JobContext, JobRunner, JobStatus
from chromsim.jobs.local import LocalProcessRunner
from chromsim.jobs.tasks import (
    AssayInput,
    FitSetup,
    ProcessSetup,
    run_fit,
    run_process_optimization,
    run_simulation,
)

__all__ = [
    "AssayInput",
    "FitSetup",
    "JobCancelled",
    "JobContext",
    "JobRunner",
    "JobStatus",
    "LocalProcessRunner",
    "ProcessSetup",
    "run_fit",
    "run_process_optimization",
    "run_simulation",
]
