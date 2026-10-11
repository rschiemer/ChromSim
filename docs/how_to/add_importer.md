# How to add a chromatogram importer

A new importer is **one new file + a fixture + a test entry + a docs section**. No core module
changes (CLAUDE.md rule 6).

## 1. Write the importer

Create `src/chromsim/data/importers/my_instrument.py`:

```python
"""MyInstrument export (.dat): one header line, then volume and UV columns."""

from pathlib import Path

import numpy as np

from chromsim import units
from chromsim.data.experiment import Experiment, Signal
from chromsim.data.importers.base import file_metadata
from chromsim.data.importers.registry import register_importer


@register_importer
class MyInstrumentImporter:
    """MyInstrument text export."""

    name = "my_instrument"
    description = "MyInstrument export (.dat)"
    file_extensions = (".dat",)

    def can_read(self, path: Path) -> bool:
        return path.suffix.lower() == ".dat" and path.read_text().startswith("MyInstrument")

    def read(self, path: Path) -> Experiment:
        table = np.loadtxt(path, skiprows=1)
        values, si_unit = units.signal_to_si(table[:, 1], "mAU")   # always convert via units.py
        uv = Signal("UV", "uv", "volume", table[:, 0] * 1e-6, values, si_unit, "mAU")
        return Experiment(path.stem, {"UV": uv}, metadata=file_metadata(path, self.name))
```

Rules:

- **Only read** the file: never write, rename or touch it.
- Convert all units through `chromsim.units` (`signal_to_si`, `to_si`). x must be in m³ (volume)
  or s (time), and increasing.
- Raise `ChromSimError` with a message that says how to fix the problem.

## 2. Register it

Add the module to the import list in `src/chromsim/data/importers/__init__.py`.

## 3. Add a fixture and the contract test

Put a small, anonymised example file in `tests/data/my_instrument/`. Then add it to
`FIXTURES` in `tests/unit/test_importers.py`:

```python
"my_instrument": (DATA / "my_instrument" / "example.dat", {}),
```

The shared contract test then checks that the file is not modified, that provenance is recorded,
and that every signal has an increasing axis and SI units. Add format-specific tests next to it.

## 4. Document it

Add a row to the importer table, and a short section, in `docs/user_guide/experimental_data.md`.
