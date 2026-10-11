"""Generic delimited-text chromatogram importer driven by a column mapping (SPEC §7.2).

Use it for exports from any instrument. The mapping says which column is the x axis, in which
unit, and which columns are signals::

    x_column: Volume (mL)
    x_unit: mL
    signals:
      - {column: "UV 280 (mAU)", name: UV280, unit: mAU}
      - {column: "Conductivity (mS/cm)", name: Cond, unit: mS/cm, kind: conductivity}

The mapping can be given as a dict, a :class:`GenericCsvMapping` or a path to a YAML file.
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
    """One signal column: ``column`` header, display ``name``, ``unit`` and optional ``kind``."""

    model_config = ConfigDict(extra="forbid")

    column: str
    name: str | None = None
    unit: str
    kind: str | None = None


class GenericCsvMapping(BaseModel):
    """Column mapping for :class:`GenericCsvImporter`.

    Attributes
    ----------
    x_column : str
        Header of the x-axis column.
    x_unit : str
        Its unit: a volume (``mL``) or a time (``min``, ``s``).
    signals : list of SignalColumn
        Signal columns to import.
    delimiter : str, optional
        Column separator; detected if omitted.
    decimal : str
        Decimal mark, ``"."`` (default) or ``","``.
    skip_rows : int
        Lines to skip before the header row.
    """

    model_config = ConfigDict(extra="forbid")

    x_column: str
    x_unit: str
    signals: list[SignalColumn]
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
        raise ChromSimError(f"Invalid generic_csv column mapping: {error}") from error


@register_importer
class GenericCsvImporter:
    """Any delimited text file, with a column mapping (``mapping=``)."""

    name = "generic_csv"
    description = "Any delimited text chromatogram, with a column mapping"
    file_extensions = (".csv", ".txt", ".tsv")

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
        table = pd.read_csv(
            path,
            sep=mapping.delimiter,
            engine="python",
            decimal=mapping.decimal,
            skiprows=mapping.skip_rows,
            encoding_errors="replace",
        )
        missing = [
            c
            for c in [mapping.x_column, *(s.column for s in mapping.signals)]
            if c not in table.columns
        ]
        if missing:
            raise ChromSimError(
                f"{path.name}: column(s) {', '.join(repr(c) for c in missing)} not found. "
                f"Columns in the file: {', '.join(repr(c) for c in table.columns)}."
            )
        x_quantity = units.Quantity(magnitude=1.0, unit=mapping.x_unit)
        axis = {"volume": "volume", "time": "time"}[
            units.check_dimension(x_quantity, "volume", "time")
        ]
        x_factor = units.to_si(x_quantity, axis)
        signals = {}
        for column in mapping.signals:
            frame = (
                table[[mapping.x_column, column.column]]
                .apply(pd.to_numeric, errors="coerce")
                .dropna()
            )
            values, si_unit = units.signal_to_si(
                frame[column.column].to_numpy(dtype=float), column.unit
            )
            name = column.name or column.column
            signals[name] = Signal(
                name=name,
                kind=column.kind or guess_kind(name, column.unit),
                x_axis=axis,
                x=frame[mapping.x_column].to_numpy(dtype=float) * x_factor,
                values=values,
                unit=si_unit,
                original_unit=column.unit,
            )
        return Experiment(name=path.stem, signals=signals, metadata=file_metadata(path, self.name))
