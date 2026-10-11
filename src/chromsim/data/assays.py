"""Offline fraction assays (SEC/CEX/HIC-HPLC, ELISA …) read through mapping profiles (SPEC §7.2b).

Labs lay out assay results in many ways. A *mapping profile* (YAML) describes one layout once,
and :func:`read_assay` then reads any file in that layout.

**Wide** layout: one row per fraction, one column per analyte::

    Fraction | Start [ml] | End [ml] | Total protein [mg/ml] | SD [mg/ml]

**Long** layout: one row per fraction and analyte::

    Sample ID | Peak    | Concentration | Unit | SD
    A3        | Monomer | 0.82          | g/L  | 0.02

Each profile also maps analytes (readouts) to model components. A readout is normally one
component, but it can be a sum, e.g. ``total_protein: [mAb, HMW]``.

Values below the limit of quantification may be written as ``<LOQ``, ``< 0.05``, ``BLQ`` or
``n.d.``; they are kept as missing values flagged ``below_loq``.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import yaml
from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator

from chromsim import units
from chromsim.data.importers.base import file_metadata, parse_number
from chromsim.errors import ChromSimError, UnitError
from chromsim.units import Quantity

logger = logging.getLogger(__name__)

BELOW_LOQ_TEXTS = {"<loq", "blq", "n.d.", "nd", "<lod", "below loq"}


class AnalyteColumn(BaseModel):
    """Wide layout: the column holding one analyte, its unit and an optional SD column."""

    model_config = ConfigDict(extra="forbid")

    column: str
    analyte: str
    unit: str
    sd_column: str | None = None


class AssayMappingProfile(BaseModel):
    """How to read one lab's assay result files.

    Attributes
    ----------
    name : str
        Profile name, e.g. ``"SEC-HPLC lab 42"``.
    layout : {"wide", "long"}
        One row per fraction (wide) or per fraction and analyte (long).
    sheet : str or int, optional
        Worksheet for Excel files (name or 0-based index); default the first sheet.
    skip_rows : int
        Lines before the header row (titles, metadata).
    delimiter : str, optional
        Column separator for text files; detected if omitted.
    decimal : str
        Decimal mark in text files: ``"."`` or ``","``.
    fraction_column : str
        Column with fraction ids (matched to the ÄKTA fraction labels).
    start_volume_column, end_volume_column : str, optional
        Columns with fraction start/end volumes; otherwise volumes come from the chromatogram.
    volume_unit : str
        Unit of those volume columns, default ``mL``.
    analytes : list of AnalyteColumn
        Wide layout: one entry per analyte column.
    analyte_column, value_column : str, optional
        Long layout: columns with the analyte name and the value.
    unit_column : str, optional
        Long layout: column with the unit of each value; or give a fixed ``unit``.
    unit : str, optional
        Long layout: unit of all values when there is no unit column.
    sd_column : str, optional
        Long layout: column with the replicate standard deviation.
    analyte_names : dict[str, str]
        Long layout: rename analyte labels in the file to analyte ids, e.g. ``Monomer: monomer``.
    loq : dict[str, str]
        Limit of quantification per analyte, e.g. ``{monomer: "0.02 g/L"}``.
    readouts : dict[str, list[str]]
        Analyte → model components whose concentrations add up to it.
    """

    model_config = ConfigDict(extra="forbid")

    name: str
    layout: Literal["wide", "long"]
    sheet: str | int | None = None
    skip_rows: int = 0
    delimiter: str | None = None
    decimal: str = "."
    fraction_column: str
    start_volume_column: str | None = None
    end_volume_column: str | None = None
    volume_unit: str = "mL"
    analytes: list[AnalyteColumn] = Field(default_factory=list)
    analyte_column: str | None = None
    value_column: str | None = None
    unit_column: str | None = None
    unit: str | None = None
    sd_column: str | None = None
    analyte_names: dict[str, str] = Field(default_factory=dict)
    loq: dict[str, str] = Field(default_factory=dict)
    readouts: dict[str, list[str]]

    @model_validator(mode="after")
    def _layout_fields(self) -> "AssayMappingProfile":
        if self.layout == "wide" and not self.analytes:
            raise ValueError("layout 'wide' needs `analytes` (one entry per analyte column).")
        if self.layout == "long":
            missing = [f for f in ("analyte_column", "value_column") if getattr(self, f) is None]
            if missing:
                raise ValueError(f"layout 'long' needs {', '.join(missing)}.")
            if self.unit_column is None and self.unit is None:
                raise ValueError("layout 'long' needs `unit_column` or a fixed `unit`.")
        if (self.start_volume_column is None) != (self.end_volume_column is None):
            raise ValueError("give both start_volume_column and end_volume_column, or neither.")
        return self


def load_mapping_profile(path: str | Path) -> AssayMappingProfile:
    """Read an assay mapping profile from YAML.

    Raises
    ------
    ChromSimError
        If the file is missing or the profile is invalid (message lists the problems).
    """
    path = Path(path)
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
        return AssayMappingProfile.model_validate(data)
    except FileNotFoundError:
        raise ChromSimError(f"Mapping profile not found: {path}") from None
    except (yaml.YAMLError, ValidationError) as error:
        raise ChromSimError(f"Invalid mapping profile {path}: {error}") from error


@dataclass
class FractionAssay:
    """Assay results per fraction.

    Attributes
    ----------
    name : str
        Assay name (default: file name).
    results : pandas.DataFrame
        One row per fraction and analyte, columns ``fraction``, ``analyte``, ``value``, ``unit``,
        ``sd``, ``below_loq``, ``loq``. Values are in the unit given in ``unit`` (as measured).
    fraction_volumes : dict[str, tuple[float, float]]
        Fraction id → (start, end) volume in m³, if the file contains volumes.
    readouts : dict[str, list[str]]
        Analyte → model components (from the profile).
    metadata : dict
        Source file, profile name, SHA-256, import time and user.
    """

    name: str
    results: pd.DataFrame
    fraction_volumes: dict[str, tuple[float, float]] = field(default_factory=dict)
    readouts: dict[str, list[str]] = field(default_factory=dict)
    metadata: dict = field(default_factory=dict)

    @property
    def fractions(self) -> list[str]:
        """Fraction ids in file order."""
        return list(dict.fromkeys(self.results["fraction"]))

    @property
    def analytes(self) -> list[str]:
        """Analyte ids in file order."""
        return list(dict.fromkeys(self.results["analyte"]))


def _read_table(path: Path, profile: AssayMappingProfile) -> pd.DataFrame:
    """Read the raw table as text cells, so values like '<LOQ' survive."""
    if path.suffix.lower() in (".xlsx", ".xlsm"):
        return pd.read_excel(
            path,
            sheet_name=profile.sheet or 0,
            skiprows=profile.skip_rows,
            dtype=object,
            engine="openpyxl",
        )
    return pd.read_csv(
        path,
        sep=profile.delimiter,
        engine="python",
        skiprows=profile.skip_rows,
        dtype=str,
        keep_default_na=False,
        encoding="utf-8-sig",
    )


def _parse_value(cell, decimal: str) -> tuple[float, bool]:
    """(value, below_loq) from one cell; raises ValueError for text that is neither."""
    if cell is None or (isinstance(cell, float) and np.isnan(cell)):
        return np.nan, False
    text = str(cell).strip()
    if not text:
        return np.nan, False
    if text.lower() in BELOW_LOQ_TEXTS or text.startswith("<"):
        return np.nan, True
    if decimal == ",":
        text = text.replace(",", ".")
    value = parse_number(text)
    if value is None:
        raise ValueError(text)
    return value, False


def _check_concentration_unit(unit: str, analyte: str) -> None:
    try:
        units.check_dimension(
            Quantity(magnitude=1.0, unit=unit),
            "molar_concentration",
            "mass_concentration",
            what="concentration",
        )
    except UnitError as error:
        raise ChromSimError(
            f"Analyte '{analyte}' is reported in '{unit}', which is not a concentration. ChromSim "
            "compares concentrations (e.g. g/L, mg/mL, mM); convert percentages using a total "
            "concentration first."
        ) from error


def read_assay(path: str | Path, profile: AssayMappingProfile | str | Path) -> FractionAssay:
    """Read an assay result file using a mapping profile.

    Parameters
    ----------
    path : str or Path
        CSV/TXT or XLSX file. It is never modified.
    profile : AssayMappingProfile, or path to its YAML file
        Describes the layout.

    Returns
    -------
    FractionAssay
        Results in long form with units, SDs and LOQ flags.

    Raises
    ------
    ChromSimError
        If columns are missing, a value cannot be read, or a unit is not a concentration.
    """
    path = Path(path)
    if not isinstance(profile, AssayMappingProfile):
        profile = load_mapping_profile(profile)
    if not path.is_file():
        raise ChromSimError(f"Assay file not found: {path}")
    try:
        table = _read_table(path, profile)
    except (pd.errors.ParserError, ValueError) as error:
        raise ChromSimError(
            f"{path.name} cannot be read with profile '{profile.name}' ({error}). "
            "Check skip_rows, delimiter and sheet in the profile."
        ) from error
    table.columns = [str(column).strip() for column in table.columns]

    needed = [profile.fraction_column]
    if profile.start_volume_column and profile.end_volume_column:
        needed += [profile.start_volume_column, profile.end_volume_column]
    if profile.layout == "wide":
        needed += [a.column for a in profile.analytes] + [
            a.sd_column for a in profile.analytes if a.sd_column
        ]
    else:
        needed += [
            c
            for c in (
                profile.analyte_column,
                profile.value_column,
                profile.unit_column,
                profile.sd_column,
            )
            if c
        ]
    missing = [column for column in needed if column not in table.columns]
    if missing:
        raise ChromSimError(
            f"{path.name}: column(s) {', '.join(repr(c) for c in missing)} not found with profile "
            f"'{profile.name}'. Columns in the file: {', '.join(repr(c) for c in table.columns)}. "
            "Check skip_rows and the column names in the profile."
        )

    records = []
    fraction_volumes = {}
    volume_factor = units.to_si(Quantity(magnitude=1.0, unit=profile.volume_unit), "volume")
    for row_number, row in table.iterrows():
        fraction = str(row[profile.fraction_column]).strip()
        if not fraction or fraction.lower() == "nan":
            continue
        if profile.start_volume_column:
            start, _ = _parse_value(row[profile.start_volume_column], profile.decimal)
            end, _ = _parse_value(row[profile.end_volume_column], profile.decimal)
            fraction_volumes[fraction] = (start * volume_factor, end * volume_factor)

        if profile.layout == "wide":
            entries = [
                (a.analyte, row[a.column], a.unit, row[a.sd_column] if a.sd_column else None)
                for a in profile.analytes
            ]
        else:
            label = str(row[profile.analyte_column]).strip()
            analyte = profile.analyte_names.get(label, label)
            unit = (
                str(row[profile.unit_column]).strip() if profile.unit_column else str(profile.unit)
            )
            entries = [
                (
                    analyte,
                    row[profile.value_column],
                    unit,
                    row[profile.sd_column] if profile.sd_column else None,
                )
            ]

        for analyte, cell, unit, sd_cell in entries:
            try:
                value, below_loq = _parse_value(cell, profile.decimal)
                sd, _ = _parse_value(sd_cell, profile.decimal)
            except ValueError as error:
                raise ChromSimError(
                    f"{path.name}, row {row_number + profile.skip_rows + 2}: cannot read value "
                    f"'{error}' for analyte '{analyte}' in fraction '{fraction}'."
                ) from None
            records.append(
                {
                    "fraction": fraction,
                    "analyte": analyte,
                    "value": value,
                    "unit": unit,
                    "sd": sd,
                    "below_loq": below_loq,
                    "loq": profile.loq.get(analyte),
                }
            )

    results = pd.DataFrame.from_records(
        records, columns=["fraction", "analyte", "value", "unit", "sd", "below_loq", "loq"]
    )
    for analyte, unit in results[["analyte", "unit"]].drop_duplicates().itertuples(index=False):
        _check_concentration_unit(unit, analyte)
    unmapped = sorted(set(results["analyte"]) - set(profile.readouts))
    if unmapped:
        logger.warning(
            "%s: analyte(s) %s have no `readouts` entry in profile '%s' and cannot be compared "
            "with simulations.",
            path.name,
            ", ".join(unmapped),
            profile.name,
        )
    metadata = file_metadata(path, importer_name="assay") | {"profile": profile.name}
    return FractionAssay(
        name=path.stem,
        results=results,
        fraction_volumes=fraction_volumes,
        readouts=dict(profile.readouts),
        metadata=metadata,
    )
