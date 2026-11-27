"""Generic delimited-text chromatogram importer driven by a column mapping (SPEC §7.2).

Use it for exports from any instrument. The mapping says which column is the x axis, in which
unit, and which columns are signals::

    x_column: Volume (mL)
    x_unit: mL
    signals:
      - {column: "UV 280 (mAU)", name: UV280, unit: mAU}
      - {column: "Conductivity (mS/cm)", name: Cond, unit: mS/cm, kind: conductivity}

Files without a header row (e.g. two numeric columns from MATLAB or another simulator) are read
with ``header: false``; columns are then given by position, counting from 0::

    header: false
    delimiter: whitespace
    x_column: 0
    x_unit: s
    signals:
      - {column: 1, name: lysozyme, unit: mM}

Concentrations (``mM``, ``g/L`` …) become signals of kind ``concentration`` (or ``salt`` if the
name starts with "salt"). The mapping can be given as a dict, a :class:`GenericCsvMapping` or a
path to a YAML file. :mod:`chromsim.data.importers.generic_xlsx` uses the same mapping for Excel
tables.
"""

from pathlib import Path

import pandas as pd
import yaml
from pydantic import BaseModel, ConfigDict, ValidationError

from chromsim import units
from chromsim.data.experiment import Experiment, Signal
from chromsim.data.importers.base import check_options, file_metadata, guess_kind
from chromsim.data.importers.registry import register_importer
from chromsim.errors import ChromSimError


class SignalColumn(BaseModel):
    """One signal column: ``column`` header (or 0-based position), display ``name``, ``unit``
    and optional ``kind``."""

    model_config = ConfigDict(extra="forbid")

    column: str | int
    name: str | None = None
    unit: str
    kind: str | None = None


class GenericCsvMapping(BaseModel):
    """Column mapping for the generic importers.

    Attributes
    ----------
    x_column : str or int
        Header of the x-axis column, or its position (0 = first column).
    x_unit : str
        Its unit: a volume (``mL``) or a time (``min``, ``s``).
    signals : list of SignalColumn
        Signal columns to import.
    header : bool
        ``True`` (default) if the first row after ``skip_rows`` holds column names; ``False``
        for purely numeric files, whose columns are then addressed by position.
    delimiter : str, optional
        Column separator; ``"whitespace"`` for space/tab-aligned columns; detected if omitted.
        Text files only.
    decimal : str
        Decimal mark, ``"."`` (default) or ``","``. Text files only.
    skip_rows : int
        Lines to skip before the header row (or before the data, without a header).
    """

    model_config = ConfigDict(extra="forbid")

    x_column: str | int
    x_unit: str
    signals: list[SignalColumn]
    header: bool = True
    delimiter: str | None = None
    decimal: str = "."
    skip_rows: int = 0


def load_mapping(mapping) -> GenericCsvMapping:
    """Accept a GenericCsvMapping, a dict or a YAML file path."""
    if isinstance(mapping, GenericCsvMapping):
        return mapping
    try:
        if isinstance(mapping, str | Path):
            mapping = yaml.safe_load(Path(mapping).read_text(encoding="utf-8"))
        return GenericCsvMapping.model_validate(mapping)
    except (OSError, yaml.YAMLError, ValidationError) as error:
        raise ChromSimError(f"Invalid generic column mapping: {error}") from error


def _column_label(table: pd.DataFrame, column: str | int, file_name: str):
    """The table's label for a mapped column (a header text, or a 0-based position)."""
    if isinstance(column, int):
        if not 0 <= column < table.shape[1]:
            raise ChromSimError(
                f"{file_name}: column position {column} does not exist; the file has "
                f"{table.shape[1]} columns (positions 0–{table.shape[1] - 1})."
            )
        return table.columns[column]
    if column not in table.columns:
        return None
    return column


def experiment_from_table(
    table: pd.DataFrame, mapping: GenericCsvMapping, path: Path, importer_name: str
) -> Experiment:
    """Build an Experiment from a table read by a generic importer.

    Parameters
    ----------
    table : pandas.DataFrame
        All cells of the file; with ``header: false`` its column labels are positions.
    mapping : GenericCsvMapping
        Which columns to import, in which units.
    path : Path
        The source file (for the name and provenance).
    importer_name : str
        Recorded in the provenance.

    Raises
    ------
    ChromSimError
        If a mapped column is missing or the x unit is not a volume or time.
    """
    wanted = [mapping.x_column, *(s.column for s in mapping.signals)]
    labels = [_column_label(table, column, path.name) for column in wanted]
    missing = [column for column, label in zip(wanted, labels, strict=True) if label is None]
    if missing:
        raise ChromSimError(
            f"{path.name}: column(s) {', '.join(repr(c) for c in missing)} not found. "
            f"Columns in the file: {', '.join(repr(c) for c in table.columns)}."
        )
    x_label, signal_labels = labels[0], labels[1:]
    x_quantity = units.Quantity(magnitude=1.0, unit=mapping.x_unit)
    axis = {"volume": "volume", "time": "time"}[units.check_dimension(x_quantity, "volume", "time")]
    x_factor = units.to_si(x_quantity, axis)
    signals = {}
    for column, label in zip(mapping.signals, signal_labels, strict=True):
        frame = table[[x_label, label]].apply(pd.to_numeric, errors="coerce").dropna()
        values, si_unit = units.signal_to_si(frame[label].to_numpy(dtype=float), column.unit)
        name = column.name or str(column.column)
        signals[name] = Signal(
            name=name,
            kind=column.kind or guess_kind(name, column.unit),
            x_axis=axis,
            x=frame[x_label].to_numpy(dtype=float) * x_factor,
            values=values,
            unit=si_unit,
            original_unit=column.unit,
        )
    return Experiment(name=path.stem, signals=signals, metadata=file_metadata(path, importer_name))


@register_importer
class GenericCsvImporter:
    """Any delimited text file, with a column mapping (``mapping=``)."""

    name = "generic_csv"
    description = "Any delimited text chromatogram, with a column mapping"
    file_extensions = (".csv", ".txt", ".tsv", ".dat")

    def can_read(self, path: Path) -> bool:
        return False  # needs a mapping, so never chosen automatically

    def read(self, path: Path, **options) -> Experiment:
        """Options: ``mapping`` (required): dict, GenericCsvMapping or path to a YAML file."""
        check_options(self.name, options, allowed=("mapping",))
        mapping = options.get("mapping")
        if mapping is None:
            raise ChromSimError(
                "generic_csv needs mapping=... (dict, GenericCsvMapping or YAML path)."
            )
        mapping = load_mapping(mapping)
        separator = r"\s+" if mapping.delimiter == "whitespace" else mapping.delimiter
        table = pd.read_csv(
            path,
            sep=separator,
            engine="python",
            decimal=mapping.decimal,
            skiprows=mapping.skip_rows,
            header=0 if mapping.header else None,
            encoding_errors="replace",
        )
        return experiment_from_table(table, mapping, path, self.name)
