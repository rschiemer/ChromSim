"""ÄKTA UNICORN chromatogram export, Excel format (``.xlsx``).

Same layout as the text export (see :mod:`chromsim.data.importers.akta_csv`): two columns per
curve, a row of curve names and a row of units. The first worksheet is read unless ``sheet=`` is
given. Files are opened read-only.
"""

from pathlib import Path

import openpyxl

from chromsim.data.experiment import Experiment
from chromsim.data.importers.akta_csv import looks_like_unicorn, parse_unicorn_table
from chromsim.data.importers.base import check_options, file_metadata
from chromsim.data.importers.registry import register_importer
from chromsim.errors import ChromSimError


def read_sheet_rows(path: Path, sheet: str | None = None) -> list[list[str]]:
    """All cells of a worksheet as text (numbers via ``repr`` so no precision is lost)."""
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        worksheet = workbook[sheet] if sheet else workbook.worksheets[0]
        rows = []
        for row in worksheet.iter_rows(values_only=True):
            rows.append(
                [
                    "" if cell is None else (repr(cell) if isinstance(cell, float) else str(cell))
                    for cell in row
                ]
            )
        return rows
    finally:
        workbook.close()


@register_importer
class AktaXlsxImporter:
    """ÄKTA UNICORN Excel export (.xlsx)."""

    name = "akta_xlsx"
    description = "ÄKTA UNICORN chromatogram export as Excel (.xlsx)"
    file_extensions = (".xlsx",)

    def can_read(self, path: Path) -> bool:
        if path.suffix.lower() not in self.file_extensions:
            return False
        try:
            return looks_like_unicorn(read_sheet_rows(path)[:6])
        except Exception:  # unreadable or non-Excel content: just "not mine"
            return False

    def read(self, path: Path, **options) -> Experiment:
        """Options: ``sheet`` (worksheet name; default the first sheet)."""
        check_options(self.name, options, allowed=("sheet",))
        sheet = options.get("sheet")
        try:
            rows = read_sheet_rows(path, sheet)
        except KeyError as error:
            raise ChromSimError(f"{path.name} has no worksheet {sheet!r}.") from error
        return parse_unicorn_table(rows, name=path.stem, metadata=file_metadata(path, self.name))
