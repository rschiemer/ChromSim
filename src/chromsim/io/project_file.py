"""Read and write project files (YAML ⇄ :class:`~chromsim.core.project.Project`).

Project files are plain YAML in the layout of SPEC §4.3. Quantities are written exactly as the
user gave them (``"20 cm"`` stays ``"20 cm"``), so a load–save cycle loses nothing.

Every file carries a ``schema_version``. When the Project structure changes incompatibly, bump
``CURRENT_SCHEMA_VERSION`` in ``chromsim/core/project.py`` and add a function to ``MIGRATIONS``
that turns version *n* data into version *n + 1* data. Old files are then upgraded on load.
"""

import copy
import logging
from collections.abc import Callable
from pathlib import Path

import yaml
from pydantic import ValidationError

from chromsim.core.project import CURRENT_SCHEMA_VERSION, Project
from chromsim.errors import ProjectFileError

logger = logging.getLogger(__name__)

# from-version -> function converting the data dict to from-version + 1. Constant: never modified
# at runtime. Empty while only schema version 1 exists.
MIGRATIONS: dict[int, Callable[[dict], dict]] = {}


def migrate(data: dict, source: str) -> dict:
    """Upgrade project data to the current schema version.

    Parameters
    ----------
    data : dict
        Project data as read from YAML. Not modified.
    source : str
        Where the data came from (file name), for error messages.

    Returns
    -------
    dict
        Data at ``CURRENT_SCHEMA_VERSION``.

    Raises
    ------
    ProjectFileError
        If the version is missing, not an integer, newer than this ChromSim, or has no migration.
    """
    version = data.get("schema_version")
    if version is None:
        raise ProjectFileError(
            f"{source} has no schema_version; add `schema_version: {CURRENT_SCHEMA_VERSION}` "
            "as the first line."
        )
    if not isinstance(version, int) or isinstance(version, bool):
        raise ProjectFileError(f"{source}: schema_version must be a whole number, got {version!r}.")
    if version > CURRENT_SCHEMA_VERSION:
        raise ProjectFileError(
            f"{source} uses schema_version {version}, written by a newer ChromSim (this one reads "
            f"up to {CURRENT_SCHEMA_VERSION}). Update ChromSim."
        )

    migrated = copy.deepcopy(data)
    while version < CURRENT_SCHEMA_VERSION:
        if version not in MIGRATIONS:
            raise ProjectFileError(
                f"{source}: no migration from schema_version {version} to {version + 1}."
            )
        logger.info("Migrating %s from schema_version %d to %d", source, version, version + 1)
        migrated = MIGRATIONS[version](migrated)
        version += 1
        migrated["schema_version"] = version
    return migrated


def _format_validation_error(error: ValidationError, source: str) -> str:
    """Turn a pydantic error into one readable line per problem: ``location: message``."""
    lines = [f"{source} has {error.error_count()} problem(s):"]
    for problem in error.errors():
        location = ".".join(str(part) for part in problem["loc"]) or "project"
        message = problem["msg"].removeprefix("Value error, ")
        lines.append(f"  - {location}: {message}")
    return "\n".join(lines)


def project_from_yaml(text: str, source: str = "project text") -> Project:
    """Create a Project from YAML text.

    Parameters
    ----------
    text : str
        YAML in the project-file layout (SPEC §4.3).
    source : str, optional
        Name used in error messages, e.g. the file name.

    Returns
    -------
    Project
        The validated project.

    Raises
    ------
    ProjectFileError
        If the text is not valid YAML, has an unsupported schema version, or fails validation.
        The message lists every problem with its location, e.g. ``column.length``.
    """
    try:
        data = yaml.safe_load(text)
    except yaml.YAMLError as error:
        raise ProjectFileError(f"{source} is not valid YAML: {error}") from error
    if not isinstance(data, dict):
        raise ProjectFileError(
            f"{source} must be a YAML mapping (`key: value` lines) at the top level."
        )

    data = migrate(data, source)
    try:
        return Project.model_validate(data)
    except ValidationError as error:
        raise ProjectFileError(_format_validation_error(error, source)) from None


def project_to_yaml(project: Project) -> str:
    """Write a Project as YAML text.

    Optional fields that are not set are left out; every other field is written explicitly,
    including defaults, so the file fully documents the run.

    Parameters
    ----------
    project : Project
        The project to write.

    Returns
    -------
    str
        YAML text that :func:`project_from_yaml` reads back into an equal Project.
    """
    data = project.model_dump(mode="json", exclude_none=True)
    return yaml.safe_dump(data, sort_keys=False, allow_unicode=True, width=100)


def load_project(path: str | Path) -> Project:
    """Read a project file.

    Parameters
    ----------
    path : str or Path
        Path to a ``.yaml`` project file.

    Returns
    -------
    Project
        The validated project.

    Raises
    ------
    ProjectFileError
        If the file is missing, unreadable or invalid (see :func:`project_from_yaml`).
    """
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        raise ProjectFileError(f"Project file not found: {path}") from None
    except OSError as error:
        raise ProjectFileError(f"Cannot read project file {path}: {error}") from error
    return project_from_yaml(text, source=str(path))


def save_project(project: Project, path: str | Path) -> Path:
    """Write a project file, replacing an existing file at ``path``.

    Parameters
    ----------
    project : Project
        The project to write.
    path : str or Path
        Destination ``.yaml`` file. Parent folders must exist.

    Returns
    -------
    Path
        The path written.
    """
    path = Path(path)
    path.write_text(project_to_yaml(project), encoding="utf-8")
    logger.info("Saved project '%s' to %s", project.name, path)
    return path
