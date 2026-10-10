"""Domain model: components, column, resin, system, method and project (SPEC §4.2).

All models share ``MODEL_CONFIG``:

- ``extra="forbid"``: a misspelled key in a project file is an error, never silently ignored.
- ``validate_assignment=True``: changing a field later (e.g. from a notebook) is checked too.
"""

from pydantic import ConfigDict

MODEL_CONFIG = ConfigDict(extra="forbid", validate_assignment=True)
