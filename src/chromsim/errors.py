"""Exception classes raised by ChromSim.

Every error raised on purpose by the package is a subclass of :class:`ChromSimError`, so callers
(notebooks, CLI, UI) can catch ChromSim problems separately from programming errors. Messages must
say what went wrong *and* what the user can do about it (CLAUDE.md rule 7).
"""


class ChromSimError(Exception):
    """Base class for all errors raised deliberately by ChromSim."""


class EnvironmentCheckError(ChromSimError):
    """The installation is broken, e.g. CADET-Core cannot be found or a test simulation fails."""


class UnitError(ChromSimError, ValueError):
    """A quantity has a missing, unknown or wrong unit, or cannot be converted.

    Also a ``ValueError``, so pydantic reports it as a normal validation error when it is raised
    while checking a model field.
    """


class ProjectFileError(ChromSimError):
    """A project file cannot be read, has an unsupported schema version, or fails validation."""
