"""Chromatogram importers and their registry. Importing this package registers the built-ins."""

from chromsim.data.importers import (  # noqa: F401  (registration)
    akta_csv,
    akta_xls,
    akta_xlsx,
    generic_csv,
    generic_xlsx,
    unicorn_res,
)
from chromsim.data.importers.registry import get_importer, list_importers, read_chromatogram

__all__ = ["get_importer", "list_importers", "read_chromatogram"]
