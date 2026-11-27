"""ÄKTA UNICORN chromatogram export, legacy Excel 97–2003 format (``.xls``).

UNICORN 5 (and the Excel export of older UNICORN versions) writes ``.xls`` workbooks with the
same layout as the text export (see :mod:`chromsim.data.importers.akta_csv`): a row of curve
names such as ``Protein A:10_UV``, a row of units (``ml``, `` mAU``, `` %B`` …), then two columns
per curve. Above them there is often a row with the run number (``10``). The first worksheet is
read unless ``sheet=`` is given.

The file is read with ``xlrd`` (ADR 0012), which supports only this legacy format; ``.xlsx``
files go to ``akta_xlsx``.
"""

from pathlib import Path

import xlrd

from chromsim.data.experiment import Experiment
from chromsim.data.importers.akta_csv import looks_like_unicorn, parse_unicorn_table
from chromsim.data.importers.base import check_options, file_metadata
from chromsim.data.importers.registry import register_importer
from chromsim.errors import ChromSimError


def read_xls_rows(
    path: Path, sheet: str | None = None, max_rows: int | None = None
) -> list[list[str]]:
    """Cells of a legacy Excel worksheet as text (numbers via ``repr``, so no precision is lost).

    ``max_rows`` reads only the first rows (enough to recognise the header).
    """
    workbook = xlrd.open_workbook(str(path), on_demand=True)
    try:
        try:
            worksheet = workbook.sheet_by_name(sheet) if sheet else workbook.sheet_by_index(0)
        except xlrd.biffh.XLRDError as error:
            sheets = ", ".join(workbook.sheet_names())
            raise ChromSimError(
                f"{path.name}: no worksheet {sheet!r}; worksheets: {sheets}."
            ) from error
        n_rows = worksheet.nrows if max_rows is None else min(max_rows, worksheet.nrows)
        rows = []
        for row_index in range(n_rows):
            row = []
            for cell in worksheet.row(row_index):
                if cell.ctype == xlrd.XL_CELL_NUMBER:
                    row.append(repr(float(cell.value)))
                elif cell.ctype in (xlrd.XL_CELL_EMPTY, xlrd.XL_CELL_BLANK):
                    row.append("")
                else:
                    row.append(str(cell.value))
            rows.append(row)
        return rows
    finally:
        workbook.release_resources()


@register_importer
class AktaXlsImporter:
    """ÄKTA UNICORN legacy Excel export (.xls)."""

    name = "akta_xls"
    description = "ÄKTA UNICORN chromatogram export as legacy Excel 97–2003 (.xls)"
    file_extensions = (".xls",)

    def can_read(self, path: Path) -> bool:
        if path.suffix.lower() not in self.file_extensions:
            return False
        try:
            return looks_like_unicorn(read_xls_rows(path, max_rows=6))
        except Exception:  # unreadable or not a legacy Excel file: just "not mine"
            return False

    def read(self, path: Path, **options) -> Experiment:
        """Options: ``sheet`` (worksheet name; default the first sheet)."""
        check_options(self.name, options, allowed=("sheet",))
        rows = read_xls_rows(path, options.get("sheet"))
        return parse_unicorn_table(rows, name=path.stem, metadata=file_metadata(path, self.name))
