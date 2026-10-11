"""Chromatogram importer registry.

A new importer is a new file with a class decorated with ``@register_importer``, plus a fixture
file, contract tests and a docs section (CLAUDE.md rule 6; see docs/how_to/add_importer.md).
The registry is filled at import time only.
"""

from pathlib import Path

from chromsim.data.experiment import Experiment
from chromsim.data.importers.base import ChromatogramImporter
from chromsim.errors import ChromSimError, PluginRegistrationError

_IMPORTERS: dict[str, ChromatogramImporter] = {}


def register_importer(importer_class: type) -> type:
    """Class decorator that adds a chromatogram importer under its ``name``."""
    importer: ChromatogramImporter = importer_class()
    if importer.name in _IMPORTERS:
        raise PluginRegistrationError(
            f"An importer named '{importer.name}' is already registered; choose a different name."
        )
    _IMPORTERS[importer.name] = importer
    return importer_class


def get_importer(name: str) -> ChromatogramImporter:
    """The importer registered under ``name``."""
    if name not in _IMPORTERS:
        raise ChromSimError(f"Unknown importer '{name}'. Available: {', '.join(list_importers())}.")
    return _IMPORTERS[name]


def list_importers() -> list[str]:
    """Names of all registered importers, sorted."""
    return sorted(_IMPORTERS)


def read_chromatogram(path: str | Path, importer: str | None = None, **options) -> Experiment:
    """Import a chromatogram file.

    Parameters
    ----------
    path : str or Path
        File to read. It is never modified.
    importer : str, optional
        Importer name (see :func:`list_importers`). By default the first importer that
        recognises the file is used.
    **options
        Passed to the importer, e.g. ``mapping=`` for ``generic_csv``.

    Returns
    -------
    Experiment
        Signals in SI, fraction marks, injections, run log and provenance.

    Raises
    ------
    ChromSimError
        If the file does not exist or no importer recognises it.
    """
    path = Path(path)
    if not path.is_file():
        raise ChromSimError(f"Chromatogram file not found: {path}")
    if importer is not None:
        return get_importer(importer).read(path, **options)
    for candidate in _IMPORTERS.values():
        if candidate.can_read(path):
            return candidate.read(path, **options)
    raise ChromSimError(
        f"No importer recognises {path.name}. Choose one with importer=..., e.g. 'generic_csv' "
        f"with a column mapping. Available importers: {', '.join(list_importers())}."
    )
