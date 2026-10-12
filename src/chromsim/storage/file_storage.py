"""File-based storage of a workspace folder (SPEC §8.1, v1).

Layout::

    <root>/projects/<name>.yaml          project files
    <root>/results/<name>.h5             simulation results (HDF5 with provenance)
    <root>/optimizations/<name>.h5       optimization results with full history
    <root>/data/<name>                   uploaded raw files (written once, never changed)
    <root>/records/<kind>/<name>.json    small JSON records (jobs, fit summaries, settings)

Every project, result and record stores ``created_by`` (the OS user) and a timestamp for future
multi-user provenance. JSON records are written atomically (temporary file + rename), so a
reader in another process never sees half a file.
"""

import datetime
import getpass
import json
import os
import re
import tempfile
from pathlib import Path

from chromsim.backends.base import SimulationResult
from chromsim.core.project import Project
from chromsim.errors import ChromSimError
from chromsim.io.project_file import load_project, save_project
from chromsim.io.results import (
    read_optimization_result,
    read_simulation_result,
    write_optimization_result,
    write_simulation_result,
)
from chromsim.optimizers.base import OptimizationResult

_VALID_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\- ]{0,100}$")


def _check_name(name: str) -> str:
    if not _VALID_NAME.match(name) or ".." in name:
        raise ChromSimError(
            f"Invalid name {name!r}: use letters, digits, spaces, '.', '_' or '-' (max. 100), "
            "starting with a letter or digit."
        )
    return name


def stamp() -> dict:
    """Provenance stamp: creation time (UTC) and OS user."""
    return {
        "created_at": datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds"),
        "created_by": getpass.getuser(),
    }


class FileStorage:
    """A workspace folder.

    Parameters
    ----------
    root : str or Path
        Folder; created if missing.
    """

    def __init__(self, root: str | Path):
        self.root = Path(root).expanduser().resolve()
        for folder in ("projects", "results", "optimizations", "data", "records"):
            (self.root / folder).mkdir(parents=True, exist_ok=True)

    # --- projects -----------------------------------------------------------------------------

    def save_project(self, name: str, project: Project) -> None:
        save_project(project, self.root / "projects" / f"{_check_name(name)}.yaml")
        self.save_record("provenance", f"project {name}", stamp())

    def load_project(self, name: str) -> Project:
        path = self.root / "projects" / f"{_check_name(name)}.yaml"
        if not path.is_file():
            raise ChromSimError(f"No project '{name}' in workspace {self.root}.")
        return load_project(path)

    def list_projects(self) -> list[str]:
        return sorted(path.stem for path in (self.root / "projects").glob("*.yaml"))

    # --- results ------------------------------------------------------------------------------

    def save_simulation_result(self, name: str, result: SimulationResult) -> None:
        write_simulation_result(result, self.root / "results" / f"{_check_name(name)}.h5")

    def load_simulation_result(self, name: str) -> SimulationResult:
        path = self.root / "results" / f"{_check_name(name)}.h5"
        if not path.is_file():
            raise ChromSimError(f"No simulation result '{name}' in workspace {self.root}.")
        return read_simulation_result(path)

    def result_file(self, name: str) -> Path:
        """Path of a stored simulation result's HDF5 file (e.g. for downloads)."""
        path = self.root / "results" / f"{_check_name(name)}.h5"
        if not path.is_file():
            raise ChromSimError(f"No simulation result '{name}' in workspace {self.root}.")
        return path

    def list_simulation_results(self) -> list[str]:
        return sorted(path.stem for path in (self.root / "results").glob("*.h5"))

    def save_optimization_result(self, name: str, result: OptimizationResult) -> None:
        write_optimization_result(result, self.root / "optimizations" / f"{_check_name(name)}.h5")

    def load_optimization_result(self, name: str) -> OptimizationResult:
        path = self.root / "optimizations" / f"{_check_name(name)}.h5"
        if not path.is_file():
            raise ChromSimError(f"No optimization result '{name}' in workspace {self.root}.")
        return read_optimization_result(path)

    # --- raw data (read-only once stored) -----------------------------------------------------

    def save_raw_file(self, name: str, content: bytes) -> Path:
        """Store an uploaded file. Raw data are never overwritten (CLAUDE.md)."""
        path = self.root / "data" / _check_name(name)
        if path.exists():
            if path.read_bytes() == content:
                return path
            raise ChromSimError(
                f"A different file named '{name}' is already stored; raw data are never "
                "overwritten. Upload it under another name."
            )
        path.write_bytes(content)
        path.chmod(0o444)  # read-only, as a reminder
        self.save_record("provenance", f"data {name}", stamp())
        return path

    def raw_file_path(self, name: str) -> Path:
        path = self.root / "data" / _check_name(name)
        if not path.is_file():
            raise ChromSimError(f"No data file '{name}' in workspace {self.root}.")
        return path

    def list_raw_files(self) -> list[str]:
        return sorted(path.name for path in (self.root / "data").iterdir() if path.is_file())

    # --- JSON records -------------------------------------------------------------------------

    def _record_path(self, kind: str, name: str) -> Path:
        folder = self.root / "records" / _check_name(kind)
        folder.mkdir(exist_ok=True)
        return folder / f"{_check_name(name)}.json"

    def save_record(self, kind: str, name: str, record: dict) -> None:
        path = self._record_path(kind, name)
        handle, temporary = tempfile.mkstemp(dir=path.parent, suffix=".tmp")
        with os.fdopen(handle, "w", encoding="utf-8") as file:
            json.dump(record, file, indent=1, default=str)
        os.replace(temporary, path)  # atomic on all platforms

    def load_record(self, kind: str, name: str) -> dict:
        path = self._record_path(kind, name)
        if not path.is_file():
            raise ChromSimError(f"No {kind} record '{name}' in workspace {self.root}.")
        return json.loads(path.read_text(encoding="utf-8"))

    def list_records(self, kind: str) -> list[str]:
        folder = self.root / "records" / _check_name(kind)
        return sorted(path.stem for path in folder.glob("*.json")) if folder.is_dir() else []

    def has_record(self, kind: str, name: str) -> bool:
        return self._record_path(kind, name).is_file()
