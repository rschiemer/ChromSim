"""Native UNICORN 3–5 result files (``.res``) from ÄKTA systems (ÄKTA prime, purifier, explorer).

A ``.res`` file is binary. Its layout (little-endian) is:

- **Start**: 16 identification bytes (``11 47 11 47 18 00 00 00 B0 02 00 00 20 6C 03 00``),
  the file size (int32) at byte 16 and the format version (e.g. ``UNICORN 3.10``) at byte 24.
- **Directory** from byte 686: entries of 344 bytes, the last one named ``LogBook``. Each entry
  has an 8-byte type id, the entry name (cp1252, zero-terminated, 296 bytes) and four int32:
  block size, offset to the next block, block address and offset from the block start to the data.
  Run curves are named ``<result name>:<n>_<curve>``, e.g. ``Manual run 3:1_UV``.
- **Sensor blocks** (marker ``06 00 02 00 01 00 4E 00``): a 240-byte header followed by int32
  pairs (volume, value). The header describes itself: x unit (``ml``) at byte 128 and its scale
  factor (float64) at byte 146; value unit (e.g. ``mAU``) at byte 206 and its scale factor at
  byte 224. Volume in ml = int × x scale, value = int × value scale.
- **Mark blocks** (marker ``06 00 06 00 01 00 4E 00``: logbook, fractions, injection marks):
  records of 180 bytes: time (float64, min), volume (float64, ml), text (164 bytes, cp1252).

The layout was reverse-engineered by the PyCORN project (``docs/RES_files_layout.txt`` in
github.com/pyahmed/PyCORN, which contains no code we reuse); scale factors and units here are read
from each block instead of being assumed. UNICORN 6/7 do not write ``.res`` files: export those
results as CSV/XLSX (``akta_csv``/``akta_xlsx``). Tested with files from ÄKTA prime and
ÄKTA purifier (UNICORN 3.10 format, written by UNICORN 5).

The volume axis is the one recorded by the system (not shifted to the injection); injection
marks are returned in ``Experiment.injections_m3``.
"""

import logging
import struct
from pathlib import Path
from types import MappingProxyType

import numpy as np

from chromsim import units
from chromsim.data.experiment import Experiment, Signal, fractions_from_marks
from chromsim.data.importers.base import check_options, file_metadata, guess_kind
from chromsim.data.importers.registry import register_importer
from chromsim.errors import ChromSimError

logger = logging.getLogger(__name__)

RES_SIGNATURE = bytes.fromhex("1147114718000000b0020000206c0300")
SENSOR_MARKER = bytes.fromhex("0600020001004e00")
MARK_MARKER = bytes.fromhex("0600060001004e00")

DIRECTORY_START = 686
ENTRY_SIZE = 344
ENTRY_NAME_BYTES = slice(8, 304)
ENTRY_NUMBERS_BYTES = slice(304, 320)
LAST_ENTRY_NAME = "LogBook"

SENSOR_X_UNIT = slice(128, 144)
SENSOR_X_SCALE_AT = 146
SENSOR_VALUE_UNIT = slice(206, 222)
SENSOR_VALUE_SCALE_AT = 224

MARK_RECORD_SIZE = 180
MARK_TEXT = slice(16, 180)

# How UNICORN 5 writes some units → how ChromSim's unit handling expects them.
UNIT_SPELLINGS = MappingProxyType({"mAu": "mAU", "%B": "%", "ml/min": "mL/min", "": "pH"})


def _text(raw: bytes) -> str:
    return raw.split(b"\x00", 1)[0].decode("cp1252", errors="replace").strip()


def curve_name(entry_name: str) -> str:
    """Curve name without the ``<result name>:<n>_`` prefix, e.g. ``"UV1_280nm"``."""
    if ":" in entry_name:
        after_colon = entry_name.rsplit(":", 1)[1]
        return after_colon.split("_", 1)[1] if "_" in after_colon else after_colon
    return entry_name


def read_directory(data: bytes, name: str) -> list[tuple[str, int, int, int]]:
    """Directory entries as (entry name, block size, block address, data offset)."""
    if not data.startswith(RES_SIGNATURE):
        raise ChromSimError(f"{name}: not a UNICORN .res file (unknown first bytes).")
    entries = []
    position = DIRECTORY_START
    while position + ENTRY_SIZE <= len(data):
        entry = data[position : position + ENTRY_SIZE]
        entry_name = _text(entry[ENTRY_NAME_BYTES])
        size, _next_offset, address, data_offset = struct.unpack("<4i", entry[ENTRY_NUMBERS_BYTES])
        if not entry_name:
            break
        entries.append((entry_name, size, address, data_offset))
        if entry_name == LAST_ENTRY_NAME:
            break
        position += ENTRY_SIZE
    if not entries:
        raise ChromSimError(f"{name}: the .res file has no directory entries.")
    return entries


def read_sensor(data: bytes, size: int, address: int, data_offset: int):
    """(volume in mL, values, value unit) of a sensor block."""
    header = data[address : address + data_offset]
    x_unit = _text(header[SENSOR_X_UNIT])
    if x_unit.lower() != "ml":
        raise ChromSimError(f"Unexpected x unit '{x_unit}' in a .res sensor block (expected ml).")
    x_scale = struct.unpack_from("<d", header, SENSOR_X_SCALE_AT)[0]
    value_scale = struct.unpack_from("<d", header, SENSOR_VALUE_SCALE_AT)[0]
    value_unit = _text(header[SENSOR_VALUE_UNIT])
    value_unit = UNIT_SPELLINGS.get(value_unit, value_unit)
    raw = data[address + data_offset : address + size]
    pairs = np.frombuffer(raw[: len(raw) // 8 * 8], dtype="<i4").reshape(-1, 2)
    return pairs[:, 0] * x_scale, pairs[:, 1] * value_scale, value_unit


def read_marks(data: bytes, size: int, address: int, data_offset: int) -> list[tuple[float, str]]:
    """(volume in mL, text) of a logbook, fraction or injection block."""
    raw = data[address + data_offset : address + size]
    marks = []
    for start in range(0, len(raw) - MARK_RECORD_SIZE + 1, MARK_RECORD_SIZE):
        record = raw[start : start + MARK_RECORD_SIZE]
        _time_min, volume_ml = struct.unpack_from("<2d", record, 0)
        marks.append((volume_ml, _text(record[MARK_TEXT])))
    return marks


@register_importer
class UnicornResImporter:
    """Native UNICORN 3–5 result file (.res) from ÄKTA prime/purifier/explorer systems."""

    name = "unicorn_res"
    description = "ÄKTA UNICORN 3–5 native result file (.res)"
    file_extensions = (".res",)

    def can_read(self, path: Path) -> bool:
        if path.suffix.lower() not in self.file_extensions:
            return False
        try:
            with path.open("rb") as handle:
                return handle.read(len(RES_SIGNATURE)) == RES_SIGNATURE
        except OSError:
            return False

    def read(self, path: Path, **options) -> Experiment:
        check_options(self.name, options, allowed=())
        data = path.read_bytes()
        ml = units.to_si(units.Quantity(magnitude=1.0, unit="mL"), "volume")
        signals: dict[str, Signal] = {}
        fraction_marks, injections, log = [], [], []
        for entry_name, size, address, data_offset in read_directory(data, path.name):
            if size <= data_offset or address + size > len(data):
                continue  # empty or truncated block (e.g. METHODINFO with size 0)
            marker = data[address : address + len(SENSOR_MARKER)]
            curve = curve_name(entry_name)
            if marker == SENSOR_MARKER:
                volume_ml, values, value_unit = read_sensor(data, size, address, data_offset)
                if volume_ml.size == 0:
                    logger.warning("%s: curve '%s' has no data points; skipped", path.name, curve)
                    continue
                values_si, si_unit = units.signal_to_si(values, value_unit)
                signals[curve] = Signal(
                    name=curve,
                    kind=guess_kind(curve, value_unit),
                    x_axis="volume",
                    x=volume_ml * ml,
                    values=values_si,
                    unit=si_unit,
                    original_unit=value_unit,
                )
            elif marker == MARK_MARKER:
                marks = [(v * ml, text) for v, text in read_marks(data, size, address, data_offset)]
                lowered = curve.lower()
                if "fraction" in lowered:
                    fraction_marks += marks
                elif "inject" in lowered:
                    injections += [volume for volume, _ in marks]
                elif "logbook" in lowered:
                    log += marks
        if not signals:
            raise ChromSimError(f"{path.name}: the .res file contains no sensor curves.")
        experiment = Experiment(
            name=path.stem,
            signals=signals,
            injections_m3=injections,
            log=log,
            metadata=file_metadata(path, self.name)
            | {"format_version": _text(data[24:40]) or "unknown"},
        )
        experiment.fractions = fractions_from_marks(fraction_marks, experiment.end_m3)
        return experiment
