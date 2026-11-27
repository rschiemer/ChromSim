"""Generic Excel chromatogram importer (``.xlsx``) driven by a column mapping (SPEC §7.2).

For spreadsheets that are not UNICORN exports, e.g. a table with ``time [s] | Salt [mM] |
A [mM]`` prepared by a colleague or exported from another program. It uses the same mapping as
:mod:`chromsim.data.importers.generic_csv` (``delimiter`` and ``decimal`` do not apply)::

    x_column: time [s]
    x_unit: s
    signals:
      - {column: "Salt [mM]", name: salt, unit: mM}
      - {column: "A [mM]", name: A, unit: mM}

The first worksheet is read unless ``sheet=`` is given. The file is opened read-only.
"""

from pathlib import Path

import pandas as pd

from chromsim.data.experiment import Experiment
from chromsim.data.importers.base import check_options
from chromsim.data.importers.generic_csv import experiment_from_table, load_mapping
from chromsim.data.importers.registry import register_importer
from chromsim.errors import ChromSimError


@register_importer
class GenericXlsxImporter:
    """Any Excel table, with a column mapping (``mapping=``, optional ``sheet=``)."""

    name = "generic_xlsx"
    description = "Any Excel (.xlsx) chromatogram table, with a column mapping"
    file_extensions = (".xlsx",)

    def can_read(self, path: Path) -> bool:
        return False  # needs a mapping, so never chosen automatically

    def read(self, path: Path, **options) -> Experiment:
        """Options: ``mapping`` (required) and ``sheet`` (worksheet name; default the first)."""
        check_options(self.name, options, allowed=("mapping", "sheet"))
        if options.get("mapping") is None:
            raise ChromSimError(
                "generic_xlsx needs mapping=... (dict, GenericCsvMapping or YAML path)."
            )
        mapping = load_mapping(options["mapping"])
        sheet = options.get("sheet") or 0
        try:
            table = pd.read_excel(
                path,
                sheet_name=sheet,
                header=0 if mapping.header else None,
                skiprows=mapping.skip_rows,
                engine="openpyxl",
            )
        except (ValueError, KeyError) as error:  # e.g. unknown sheet name
            raise ChromSimError(f"{path.name}: cannot read sheet {sheet!r}: {error}") from error
        return experiment_from_table(table, mapping, path, self.name)
