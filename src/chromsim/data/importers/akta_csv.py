"""ÄKTA UNICORN chromatogram export, text format (``.csv``, ``.asc``, ``.txt``).

UNICORN exports curves side by side, two columns per curve (x and value)::

    [optional title row, e.g. the result name]
    Chrom.1_UV1_280nm |        | Chrom.1_Cond |       | Chrom.1_Fractions |             | ...
    ml                | mAU    | ml           | mS/cm | ml                | (Fractions) | ...
    0.000             | 0.120  | 0.000        | 5.02  | 0.000             | Waste       | ...

Every curve has its own x values and length. ``Fractions``, ``Injection`` and ``Run Log`` curves
contain text marks instead of numbers. UNICORN 6/7 often writes these files as UTF-16 with tabs
even when the extension is ``.csv``. Encoding (UTF-16, UTF-8, Windows-1252), delimiter (tab,
comma, semicolon) and decimal commas are detected.

The parser was written against the layout above and synthetic files in ``tests/data/akta/``;
validate it with an anonymised export from your own system before relying on it.
"""

import csv
import io
import logging
import re
from pathlib import Path

import numpy as np

from chromsim import units
from chromsim.data.experiment import Experiment, Signal, fractions_from_marks
from chromsim.data.importers.base import check_options, file_metadata, guess_kind, parse_number
from chromsim.data.importers.registry import register_importer
from chromsim.errors import ChromSimError

logger = logging.getLogger(__name__)

AXIS_UNITS = {"ml": "volume", "min": "time"}
_CHROMATOGRAM_PREFIX = re.compile(r"^Chrom\.\d+_")


def _find_unit_row(rows: list[list[str]], max_header_rows: int = 6) -> int | None:
    """Index of the row with axis units (``ml``/``min``) in every other column, or None."""
    for index, row in enumerate(rows[:max_header_rows]):
        x_cells = [cell.strip().lower() for cell in row[0::2] if cell.strip()]
        if x_cells and all(cell in AXIS_UNITS for cell in x_cells):
            return index
    return None


def looks_like_unicorn(rows: list[list[str]]) -> bool:
    """``True`` if the first rows have the UNICORN curve-name / unit header."""
    unit_row = _find_unit_row(rows)
    return unit_row is not None and unit_row >= 1


def _text_marks(cells, x_factor: float, curve: str, name: str) -> list[tuple[float, str]]:
    """(x in SI, text) for a mark curve (fractions, injections, run log); skips non-numeric x."""
    marks = []
    for x_text, text in cells:
        x = parse_number(x_text)
        if x is None:
            logger.warning(
                "%s: curve '%s' has a mark at non-numeric x %r; skipped", name, curve, x_text
            )
            continue
        marks.append((x * x_factor, text.strip()))
    return marks


def parse_unicorn_table(rows: list[list[str]], name: str, metadata: dict) -> Experiment:
    """Build an Experiment from the cells of a UNICORN export (text or spreadsheet).

    Parameters
    ----------
    rows : list of list of str
        All cells as text, row by row.
    name : str
        Experiment name.
    metadata : dict
        Provenance to attach (see :func:`file_metadata`).

    Raises
    ------
    ChromSimError
        If the header cannot be found or a curve uses an axis unit other than ml/min.
    """
    unit_row_index = _find_unit_row(rows)
    if unit_row_index is None or unit_row_index == 0:
        raise ChromSimError(
            f"{name}: no UNICORN header found (a row of curve names followed by a row of units "
            "such as 'ml | mAU | ml | mS/cm'). Use the generic_csv importer with a column mapping."
        )
    name_row = rows[unit_row_index - 1]
    unit_row = rows[unit_row_index]
    data_rows = rows[unit_row_index + 1 :]

    signals: dict[str, Signal] = {}
    fraction_marks: list[tuple[float, str]] = []
    injections: list[float] = []
    log: list[tuple[float, str]] = []

    for column in range(0, len(unit_row) - 1, 2):
        curve = (
            _CHROMATOGRAM_PREFIX.sub("", name_row[column].strip().strip('"'))
            if column < len(name_row)
            else ""
        )
        x_unit = unit_row[column].strip()
        value_unit = unit_row[column + 1].strip()
        if not curve or not x_unit:
            continue
        axis = AXIS_UNITS[x_unit.lower()]
        x_si_unit = {"volume": "m**3", "time": "s"}[axis]
        x_factor = (
            units.ureg.Quantity(1.0, x_unit.lower().replace("ml", "mL")).to(x_si_unit).magnitude
        )

        cells = [
            (row[column], row[column + 1] if column + 1 < len(row) else "")
            for row in data_rows
            if column < len(row) and row[column].strip()
        ]
        lowered = curve.lower().replace(" ", "")
        if "fraction" in lowered or "inject" in lowered or "log" in lowered:
            marks = _text_marks(cells, x_factor, curve, name)
            if "fraction" in lowered:
                fraction_marks += marks
            elif "inject" in lowered:
                injections += [x for x, _ in marks]
            else:
                log += marks
        else:
            pairs = [(parse_number(x), parse_number(y)) for x, y in cells]
            pairs = [(x, y) for x, y in pairs if x is not None and y is not None]
            if not pairs:
                logger.warning("%s: curve '%s' has no numeric data; skipped", name, curve)
                continue
            x_values = np.array([x for x, _ in pairs]) * x_factor
            values, si_unit = units.signal_to_si(np.array([y for _, y in pairs]), value_unit)
            signals[curve] = Signal(
                name=curve,
                kind=guess_kind(curve, value_unit),
                x_axis=axis,
                x=x_values,
                values=values,
                unit=si_unit,
                original_unit=value_unit,
            )

    if not signals:
        raise ChromSimError(f"{name}: the UNICORN export contains no numeric curves.")
    experiment = Experiment(
        name=name, signals=signals, injections_m3=injections, log=log, metadata=metadata
    )
    experiment.fractions = fractions_from_marks(fraction_marks, experiment.end_m3)
    return experiment


def read_text_rows(path: Path) -> list[list[str]]:
    """Read a delimited text file with unknown encoding and delimiter into rows of cells."""
    raw = path.read_bytes()
    if raw.startswith((b"\xff\xfe", b"\xfe\xff")):
        text = raw.decode("utf-16")
    elif raw.startswith(b"\xef\xbb\xbf"):
        text = raw.decode("utf-8-sig")
    else:
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError:
            text = raw.decode("cp1252")
    sample = "\n".join(text.splitlines()[:10])
    try:
        delimiter = csv.Sniffer().sniff(sample, delimiters="\t;,").delimiter
    except csv.Error:
        delimiter = "\t"
    return [row for row in csv.reader(io.StringIO(text), delimiter=delimiter)]


@register_importer
class AktaCsvImporter:
    """ÄKTA UNICORN text export (CSV/ASC/TXT; UTF-16 or UTF-8; tab, comma or semicolon)."""

    name = "akta_csv"
    description = "ÄKTA UNICORN chromatogram export as text (.csv, .asc, .txt)"
    file_extensions = (".csv", ".asc", ".txt")

    def can_read(self, path: Path) -> bool:
        if path.suffix.lower() not in self.file_extensions:
            return False
        try:
            return looks_like_unicorn(read_text_rows(path)[:6])
        except (OSError, UnicodeError, csv.Error):
            return False

    def read(self, path: Path, **options) -> Experiment:
        check_options(self.name, options, allowed=())
        rows = read_text_rows(path)
        return parse_unicorn_table(rows, name=path.stem, metadata=file_metadata(path, self.name))
