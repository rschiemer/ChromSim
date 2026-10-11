"""What every chromatogram importer provides, and helpers they share."""

import datetime
import getpass
import hashlib
from pathlib import Path
from typing import Any, Protocol

from chromsim.data.experiment import Experiment
from chromsim.errors import ChromSimError


class ChromatogramImporter(Protocol):
    """A chromatogram importer.

    Importers only ever *read* files: raw data is read-only (CLAUDE.md).
    """

    name: str
    description: str
    file_extensions: tuple[str, ...]

    def can_read(self, path: Path) -> bool:
        """``True`` if this importer recognises the file (used for automatic selection)."""
        ...

    def read(self, path: Path, **options: Any) -> Experiment:
        """Read the file into an :class:`Experiment`; options are importer-specific."""
        ...


def file_metadata(path: Path, importer_name: str) -> dict:
    """Provenance of an imported file: path, importer, SHA-256, time and user."""
    return {
        "source_file": str(path.resolve()),
        "importer": importer_name,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "imported_at": datetime.datetime.now(datetime.UTC).isoformat(timespec="seconds"),
        "created_by": getpass.getuser(),
    }


def guess_kind(curve_name: str, unit: str) -> str:
    """Classify a curve by name and unit: uv, conductivity, ph, concentration_b, pressure, other."""
    name = curve_name.lower().replace(" ", "")
    if name.startswith("uv") or unit.strip().endswith("AU"):
        return "uv"
    if "cond" in name and "%" not in unit:
        return "conductivity"
    if name.startswith("ph"):
        return "ph"
    if name.startswith("conc") or name in ("%b", "concb"):
        return "concentration_b"
    if "press" in name:
        return "pressure"
    return "other"


def parse_number(text) -> float | None:
    """A float from a cell, accepting decimal commas; ``None`` for empty or non-numeric cells."""
    if text is None:
        return None
    if isinstance(text, int | float):
        return float(text)
    cleaned = str(text).strip().strip('"')
    if not cleaned:
        return None
    try:
        return float(cleaned)
    except ValueError:
        pass
    try:
        return float(cleaned.replace(",", "."))
    except ValueError:
        return None


def check_options(importer_name: str, options: dict, allowed: tuple[str, ...]) -> None:
    """Reject options an importer does not understand (so typos are not silently ignored)."""
    unknown = sorted(set(options) - set(allowed))
    if unknown:
        accepted = ", ".join(allowed) if allowed else "none"
        raise ChromSimError(
            f"{importer_name} does not accept option(s) {', '.join(unknown)}; accepted: {accepted}."
        )
