"""Generate the synthetic raw-data fixtures in tests/data/ (run manually, not by pytest).

    pixi run python tests/data/generate_fixtures.py

The files imitate real exports of a CEX run of examples/projects/cex_mab_hmw_lge.yaml, simulated
with "true" parameters (mAb k_eq = 1.5e-2 instead of the project's first guess 1.0e-2) plus
baseline drift and noise:

- akta/cex_mab_hmw_run01.csv  ÄKTA UNICORN 7 text export: UTF-16, tab-separated, title row
- akta/cex_mab_hmw_run01.xlsx the same run as UNICORN Excel export (no title row)
- akta/cex_mab_hmw_run01.xls  the same workbook in legacy Excel 97–2003 format (via LibreOffice)
- assays/total_protein_a280.csv  wide, semicolon + decimal comma, with fraction volumes, <LOQ
- assays/sec_hplc.xlsx           long format, sheet "SEC results", no volumes, includes a feed
                                 sample "Load" and a sample "A99" that was never collected

Calibration data for Tutorial 03 (same column and system; separate random generator, so the files
above are unchanged by adding these):

- akta/tracer_system.csv     acetone pulse with the column bypassed (mixer, tubing, UV delay)
- akta/tracer_dextran.csv    pore-excluded dextran pulse (ε_c and axial dispersion)
- akta/tracer_acetone.csv    pore-penetrating acetone pulse (ε_t)
- akta/lge_{10,20,30,40}cv.csv  low-load mAb runs with 10–40 CV gradients (Yamamoto method)

Native UNICORN 5 format (re-encodes run 01, no new random numbers):

- unicorn_res/cex_mab_hmw_run01.res  UV, conductivity, fractions, injection and logbook as .res

Generic tables (copies of generic/plain_export.csv):

- generic/plain_export.xlsx  the same table as an Excel sheet "Run"
- generic/headerless.dat     volume and UV only, no header, whitespace-aligned

Raw data are read-only for ChromSim (CLAUDE.md); this script is the only thing that writes them.
Re-running it with the same seed reproduces the files.
"""

import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import openpyxl
import yaml

from chromsim import api
from chromsim.analysis.detectors import predict_conductivity, predict_uv
from chromsim.analysis.fractionation import fraction_averages
from chromsim.data.experiment import Fraction
from chromsim.data.preprocessing import ConductivityCalibration

ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT / "tests"))
from synthetic import rectangular_pulse, through_system  # noqa: E402  (needs the path above)

OUT = Path(__file__).parent
TRUE_MAB_K_EQ = 1.5e-2
RNG = np.random.default_rng(20261011)
ML = 1e-6  # m³ per mL

FRACTION_START_ML, FRACTION_SIZE_ML, N_FRACTIONS = 85.0, 2.0, 25
FRACTION_LABELS = [f"{row}{col}" for row in "ABC" for col in range(1, 13)][:N_FRACTIONS]


def simulate_truth():
    project = api.load_project(ROOT / "examples/projects/cex_mab_hmw_lge.yaml")
    project.isotherm.parameters["mAb"]["k_eq"] = TRUE_MAB_K_EQ
    project.solver.n_col, project.solver.n_par = 100, 25
    return project, api.simulate(project)


def curves(project, result):
    """(name, x unit, value unit, x in mL, values) for every exported curve."""
    calibration_file = ROOT / "examples/calibrations/nacl_conductivity.yaml"
    calibration = ConductivityCalibration.model_validate(
        yaml.safe_load(calibration_file.read_text())
    )
    uv = predict_uv(result, project)
    cond = predict_conductivity(result, project, calibration)
    end_ml = result.volume_m3[-1] / ML

    x_uv = np.arange(0.0, end_ml, 0.025)
    uv_mau = np.interp(x_uv, uv.x / ML, uv.values * 1000.0, left=0.0)
    uv_mau += 1.5 + 0.01 * x_uv + RNG.normal(0.0, 0.2, x_uv.size)  # offset, drift, noise

    x_cond = np.arange(0.0, end_ml, 0.05)
    cond_ms = np.interp(x_cond, cond.x / ML, cond.values * 10.0, left=cond.values[0] * 10.0)
    cond_ms += RNG.normal(0.0, 0.03, x_cond.size)

    # Conc B is the pump set point (no delay): 50 mM = 0 %B, 1000 mM = 100 %B.
    x_b = np.arange(0.0, end_ml, 0.1)
    salt_in = np.zeros_like(x_b)
    volume = 0.0
    for step in result.steps:
        step_ml = step.flow_rate_m3_per_s * step.duration_s / ML
        inside = (x_b >= volume) & (x_b <= volume + step_ml)
        fraction = (x_b[inside] - volume) / step_ml
        salt_in[inside] = step.inlet_start_mol_per_m3[0] + fraction * (
            step.inlet_end_mol_per_m3[0] - step.inlet_start_mol_per_m3[0]
        )
        volume += step_ml
    conc_b = (salt_in - 50.0) / 950.0 * 100.0

    x_ph = np.arange(0.0, end_ml, 0.1)
    ph = 5.50 + RNG.normal(0.0, 0.01, x_ph.size)
    return [
        ("UV1_280nm", "ml", "mAU", x_uv, uv_mau),
        ("Cond", "ml", "mS/cm", x_cond, cond_ms),
        ("Conc B", "ml", "%", x_b, conc_b),
        ("pH", "ml", "pH", x_ph, ph),
    ]


def marks(result):
    fraction_starts = FRACTION_START_ML + FRACTION_SIZE_ML * np.arange(N_FRACTIONS)
    fractions_marks = [(0.0, "Waste")] + list(zip(fraction_starts, FRACTION_LABELS, strict=True))
    fractions_marks.append((FRACTION_START_ML + FRACTION_SIZE_ML * N_FRACTIONS, "Waste"))
    load = result.steps[1]
    injection = [(load.start_s * load.flow_rate_m3_per_s / ML, "")]
    log, volume = [(0.0, "Method Run 10/10/2026")], 0.0
    for step in result.steps:
        log.append((volume, f"Block {step.name.capitalize()}"))
        volume += step.flow_rate_m3_per_s * step.duration_s / ML
    log.append((volume, "End"))
    return fractions_marks, injection, log


def unicorn_columns(curve_list, fraction_marks, injection, log):
    """Columns of the UNICORN layout: header cells + data cells, two columns per curve."""
    columns = []
    for name, x_unit, unit, x, values in curve_list:
        columns.append((f"Chrom.1_{name}", x_unit, [f"{v:.3f}" for v in x]))
        columns.append(("", unit, [f"{v:.4f}" for v in values]))
    for name, label, entries in (
        ("Fractions", "(Fractions)", fraction_marks),
        ("Injection", "(Injections)", injection),
        ("Run Log", "(Logbook)", log),
    ):
        columns.append((f"Chrom.1_{name}", "ml", [f"{v:.3f}" for v, _ in entries]))
        columns.append(("", label, [text for _, text in entries]))
    return columns


def _table(columns):
    n_rows = max(len(c[2]) for c in columns)
    header = [c[0] for c in columns]
    units_row = [c[1] for c in columns]
    data = [[c[2][i] if i < len(c[2]) else "" for c in columns] for i in range(n_rows)]
    return header, units_row, data


def write_unicorn_csv(columns, csv_path, title):
    header, units_row, data = _table(columns)
    lines = [title, "\t".join(header), "\t".join(units_row)]
    lines += ["\t".join(row) for row in data]
    csv_path.write_text("\r\n".join(lines) + "\r\n", encoding="utf-16")


def write_unicorn(columns, csv_path, xlsx_path):
    header, units_row, data = _table(columns)
    write_unicorn_csv(columns, csv_path, "cex_mab_hmw_run01 001")

    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = "cex_mab_hmw_run01"
    sheet.append(header)
    sheet.append(units_row)
    for row in data:
        sheet.append([_number_or_text(cell) for cell in row])
    workbook.save(xlsx_path)


def _number_or_text(cell):
    try:
        return float(cell)
    except ValueError:
        return cell or None


def _decimal_comma(value: float, digits: int) -> str:
    return f"{value:.{digits}f}".replace(".", ",")


def write_assays(result, fraction_marks):
    collected = [
        Fraction(label, start * ML, end * ML)
        for (start, label), (end, _) in zip(fraction_marks[1:-1], fraction_marks[2:], strict=True)
    ]
    averages = fraction_averages(result, collected)
    g_per_l = {}
    for row in averages.itertuples():
        mass = (
            row.concentration_mol_per_m3 * {"salt": 0.0, "mAb": 148.0, "HMW": 296.0}[row.component]
        )
        g_per_l[(row.fraction, row.component)] = mass

    # Total protein, wide, European CSV.
    lines = [
        "Total protein (A280) - synthetic data generated by tests/data/generate_fixtures.py",
        "Analyst: anonymised;Instrument: NanoDrop;Date: 2026-10-10",
        "Fraction;Start volume [ml];End volume [ml];Total protein [mg/ml];SD [mg/ml]",
    ]
    for fraction in collected:
        total = g_per_l[(fraction.label, "mAb")] + g_per_l[(fraction.label, "HMW")]
        measured = total * (1 + RNG.normal(0.0, 0.03))
        sd = 0.02 * measured + 0.005
        value = "<LOQ" if measured < 0.05 else _decimal_comma(measured, 3)
        fields = [
            fraction.label,
            _decimal_comma(fraction.start_m3 / ML, 2),
            _decimal_comma(fraction.end_m3 / ML, 2),
            value,
            _decimal_comma(sd, 3),
        ]
        lines.append(";".join(fields))
    (OUT / "assays/total_protein_a280.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # SEC-HPLC, long format Excel, every second fraction with protein, plus feed and a stray id.
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = "SEC results"
    sheet.append(["SEC-HPLC report (synthetic)"])
    sheet.append(["Column: generic SEC, 0.5 mL/min"])
    sheet.append([])
    sheet.append(["Sample ID", "Peak", "Concentration", "Unit", "SD (n=3)"])
    sheet.append(["Load", "Monomer", 1.9, "g/L", 0.02])
    sheet.append(["Load", "HMW", 0.1, "g/L", 0.004])
    for fraction in collected[::2]:
        for peak, component in (("Monomer", "mAb"), ("HMW", "HMW")):
            measured = g_per_l[(fraction.label, component)] * (1 + RNG.normal(0.0, 0.04))
            value = round(measured, 4) if measured >= 0.005 else "< 0.005"
            sheet.append(
                [fraction.label, peak, value, "g/L", round(0.03 * max(measured, 0.0) + 0.001, 4)]
            )
    sheet.append(["A99", "Monomer", 0.5, "g/L", 0.01])
    workbook.save(OUT / "assays/sec_hplc.xlsx")


# --- calibration runs for Tutorial 03 ------------------------------------------------------

CALIBRATION_RNG = np.random.default_rng(20261012)
TRACER_INJECTION_START_ML, TRACER_INJECTION_ML = 1.0, 0.1
TRACER_RUN_ML = 9.0
UV_PER_MM = 400.0  # detector response of the tracers, mAU per mM


def _tutorial_project():
    return api.load_project(ROOT / "examples/projects/cex_mab_hmw_lge.yaml")


def _tracer_project(column_model, par_porosity):
    """Non-binding tracer pulse on the tutorial column, without system volumes."""
    base = _tutorial_project()
    data = {
        "name": "tracer",
        "components": [{"name": "tracer"}],
        "column": {
            "model": column_model,
            "length": str(base.column.length),
            "diameter": str(base.column.diameter),
            "axial_dispersion": str(base.column.axial_dispersion),
            "film_diffusion": {"tracer": "1e-4 m/s"},
            "pore_diffusion": {"tracer": "1e-10 m^2/s"},
        },
        "resin": base.resin.model_dump(mode="json") | {"par_porosity": par_porosity},
        "isotherm": {"model": "Linear", "parameters": {}},
        "method": [
            {
                "name": "equilibrate",
                "duration": f"{TRACER_INJECTION_START_ML} mL",
                "flow": "150 cm/h",
                "inlet": {},
            },
            {
                "name": "inject",
                "duration": f"{TRACER_INJECTION_ML} mL",
                "flow": "150 cm/h",
                "inlet": {"tracer": "1 mM"},
            },
            {
                "name": "elute",
                "duration": f"{TRACER_RUN_ML - 1.1} mL",
                "flow": "150 cm/h",
                "inlet": {},
            },
        ],
        "solver": {"n_col": 400 if column_model == "LRM" else 100, "n_par": 20},
    }
    return api.Project.model_validate(data)


def _uv_curve(volume_m3, concentration_mM, uv_delay_m3, step_ml=0.005):
    """UV trace at the detector, which sees the outlet ``uv_delay_m3`` later."""
    x = np.arange(0.0, volume_m3[-1] / ML, step_ml)
    values = np.interp(x * ML - uv_delay_m3, volume_m3, concentration_mM, left=0.0) * UV_PER_MM
    values += 0.8 + CALIBRATION_RNG.normal(0.0, 0.3, x.size)  # offset and noise, mAU
    return ("UV1_280nm", "ml", "mAU", x, values)


def _injection_columns():
    return [
        ("Chrom.1_Injection", "ml", [f"{TRACER_INJECTION_START_ML:.3f}"]),
        ("", "(Injections)", [""]),
    ]


def _numeric_columns(curve_list):
    columns = []
    for name, x_unit, unit, x, values in curve_list:
        columns.append((f"Chrom.1_{name}", x_unit, [f"{v:.3f}" for v in x]))
        columns.append(("", unit, [f"{v:.4f}" for v in values]))
    return columns


def write_tracer_runs():
    system = _tutorial_project().system
    mixer, tubing = system.mixer_volume_m3, system.tubing_volume_m3
    volume = np.linspace(0.0, TRACER_RUN_ML * ML, 18001)
    pulse = rectangular_pulse(volume, TRACER_INJECTION_START_ML * ML, TRACER_INJECTION_ML * ML)
    signals = {"system": (volume, through_system(volume, pulse, mixer, tubing))}
    for name, model, par_porosity in (("dextran", "LRM", 1e-9), ("acetone", "GRM", None)):
        project = _tracer_project(model, par_porosity or _tutorial_project().resin.par_porosity)
        result = api.simulate(project)
        signals[name] = (
            result.volume_m3,
            through_system(result.volume_m3, result.outlet("tracer"), mixer, tubing),
        )
    for name, (x, concentration) in signals.items():
        uv = _uv_curve(x, concentration, system.uv_delay_volume_m3)
        columns = _numeric_columns([uv]) + _injection_columns()
        write_unicorn_csv(columns, OUT / f"akta/tracer_{name}.csv", f"tracer_{name} 001")


def write_lge_runs():
    calibration_file = ROOT / "examples/calibrations/nacl_conductivity.yaml"
    calibration = ConductivityCalibration.model_validate(
        yaml.safe_load(calibration_file.read_text())
    )
    for gradient_cv in (10, 20, 30, 40):
        project = _tutorial_project()
        data = project.model_dump(mode="json")
        data["components"] = data["components"][:2]  # salt and mAb only
        data["column"]["film_diffusion"].pop("HMW")
        data["column"]["pore_diffusion"].pop("HMW")
        data["isotherm"]["parameters"] = {
            "mAb": data["isotherm"]["parameters"]["mAb"] | {"k_eq": TRUE_MAB_K_EQ}
        }
        data["method"][1] = {
            "name": "load",
            "duration": "0.5 mL",
            "flow": "150 cm/h",
            "inlet": {"salt": "50 mM", "mAb": "1 g/L"},
        }
        data["method"][3]["duration"] = f"{gradient_cv} CV"
        data["solver"] = {"n_col": 100, "n_par": 25}
        project = api.Project.model_validate(data)
        result = api.simulate(project)
        uv = predict_uv(result, project)
        cond = predict_conductivity(result, project, calibration)
        x_uv = np.arange(0.0, result.volume_m3[-1] / ML, 0.05)
        uv_mau = (
            np.interp(x_uv * ML, uv.x, uv.values * 1000.0, left=0.0)
            + 0.5
            + CALIBRATION_RNG.normal(0.0, 0.1, x_uv.size)
        )
        x_cond = np.arange(0.0, result.volume_m3[-1] / ML, 0.1)
        cond_ms = np.interp(x_cond * ML, cond.x, cond.values * 10.0, left=cond.values[0] * 10.0)
        cond_ms += CALIBRATION_RNG.normal(0.0, 0.03, x_cond.size)
        load = result.steps[1]
        injection = load.start_s * load.flow_rate_m3_per_s / ML
        columns = _numeric_columns(
            [("UV1_280nm", "ml", "mAU", x_uv, uv_mau), ("Cond", "ml", "mS/cm", x_cond, cond_ms)]
        )
        columns += [("Chrom.1_Injection", "ml", [f"{injection:.3f}"]), ("", "(Injections)", [""])]
        write_unicorn_csv(
            columns, OUT / f"akta/lge_{gradient_cv}cv.csv", f"lge_{gradient_cv}cv 001"
        )


def _res_entry(name: str, type_id: str, size: int, address: int, data_offset: int) -> bytes:
    entry = bytearray(344)
    entry[0:8] = bytes.fromhex(type_id)
    encoded = name.encode("cp1252")
    entry[8 : 8 + len(encoded)] = encoded
    block_span = -(-size // 512) * 512  # blocks are padded to 512 bytes
    entry[304:320] = struct.pack("<4i", size, block_span, address, data_offset)
    return bytes(entry)


def _res_sensor_block(volume_ml, values, unit: str, value_scale: float) -> bytes:
    header = bytearray(240)
    header[0:8] = bytes.fromhex("0600020001004e00")
    header[128:130] = b"ml"
    header[144:146] = b"\x04\x01"
    struct.pack_into("<d", header, 146, 0.01)
    header[206 : 206 + len(unit)] = unit.encode("cp1252")
    header[222:224] = b"\x04\x01"
    struct.pack_into("<d", header, 224, value_scale)
    pairs = np.column_stack(
        [np.round(np.asarray(volume_ml) / 0.01), np.round(np.asarray(values) / value_scale)]
    ).astype("<i4")
    return bytes(header) + pairs.tobytes()


def _res_mark_block(marks_ml) -> bytes:
    header = bytearray(552)
    header[0:8] = bytes.fromhex("0600060001004e00")
    records = bytearray()
    for volume, text in marks_ml:
        record = bytearray(180)
        struct.pack_into("<2d", record, 0, volume / 0.5, volume)  # time at 0.5 mL/min, unused
        encoded = text.encode("cp1252")[:163]
        record[16 : 16 + len(encoded)] = encoded
        records += record
    return bytes(header) + bytes(records)


def write_unicorn_res(source_csv: Path, res_path: Path) -> None:
    """Re-encode a UNICORN text export as a native UNICORN 5 result file (.res, format 3.10).

    Follows the layout documented in chromsim.data.importers.unicorn_res, so the importer can be
    tested against a file we are free to publish. Volumes are stored in 0.01 mL steps, UV and
    conductivity in 0.001 steps, as UNICORN 5 does.
    """
    experiment = api.read_chromatogram(source_csv)
    run = "cex_mab_hmw_run01"
    sensors = [
        ("UV1_280nm", experiment.signal("uv"), "mAU", 1000.0, 0.001),
        ("Cond", experiment.signal("conductivity"), "mS/cm", 10.0, 0.001),
    ]
    blocks = [
        (
            f"{run}:1_{name}",
            "0000010004000114",
            _res_sensor_block(signal.x / ML, signal.values * to_unit, unit, scale),
            240,
        )
        for name, signal, unit, to_unit, scale in sensors
    ]
    fraction_marks = [(f.start_m3 / ML, f.label) for f in experiment.fractions]
    fraction_marks.append((experiment.fractions[-1].end_m3 / ML, "Waste"))
    blocks += [
        (
            f"{run}:1_Logbook",
            "0000010004004804",
            _res_mark_block([(v / ML, t) for v, t in experiment.log]),
            552,
        ),
        (f"{run}:1_Fractions", "0000010004004404", _res_mark_block(fraction_marks), 552),
        (
            f"{run}:1_Inject",
            "0000010004004604",
            _res_mark_block([(v / ML, "") for v in experiment.injections_m3]),
            552,
        ),
        ("LogBook", "0000010002000113", b"Method Run 10/10/2026 (system logbook)\x00", 0),
    ]
    first_block = 686 + 344 * len(blocks)
    first_block = -(-first_block // 512) * 512
    directory, body, address = b"", b"", first_block
    for name, type_id, block, data_offset in blocks:
        directory += _res_entry(name, type_id, len(block), address, data_offset)
        padded = block + b"\x00" * (-len(block) % 512)
        body += padded
        address += len(padded)
    start = bytearray(686)
    start[0:16] = bytes.fromhex("1147114718000000b0020000206c0300")
    start[24:36] = b"UNICORN 3.10"
    head = bytes(start) + directory
    head += b"\x00" * (first_block - len(head))
    content = bytearray(head + body)
    struct.pack_into("<i", content, 16, len(content))
    res_path.parent.mkdir(exist_ok=True)
    res_path.write_bytes(bytes(content))


def write_unicorn_xls(xlsx_path: Path, xls_path: Path) -> None:
    """Convert the UNICORN Excel fixture to the legacy Excel 97–2003 format (.xls).

    Uses LibreOffice (``soffice``) in headless mode, so no Python package for writing .xls is
    needed; ChromSim itself only reads .xls (with xlrd).
    """
    soffice = shutil.which("soffice") or shutil.which("libreoffice")
    if soffice is None:
        raise SystemExit("LibreOffice (soffice) is needed to write the .xls fixture.")
    with tempfile.TemporaryDirectory() as folder:
        subprocess.run(
            [soffice, "--headless", "--convert-to", "xls", "--outdir", folder, str(xlsx_path)],
            check=True,
            capture_output=True,
        )
        converted = Path(folder) / f"{xlsx_path.stem}.xls"
        xls_path.parent.mkdir(exist_ok=True)
        shutil.copyfile(converted, xls_path)


def write_generic_variants():
    """The generic plain export as an Excel table and as a headerless whitespace file."""
    rows = [line.split(",") for line in (OUT / "generic/plain_export.csv").read_text().splitlines()]
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.title = "Run"
    sheet.append(rows[0])
    for row in rows[1:]:
        sheet.append([_number_or_text(cell) for cell in row])
    workbook.save(OUT / "generic/plain_export.xlsx")
    # Volume in mL, UV in mAU as two right-aligned columns (MATLAB `save -ascii` style).
    lines = [f"{float(row[0]):15.7e}{float(row[1]):15.7e}" for row in rows[1:]]
    (OUT / "generic/headerless.dat").write_text("\n".join(lines) + "\n")


def main():
    (OUT / "akta").mkdir(exist_ok=True)
    (OUT / "assays").mkdir(exist_ok=True)
    project, result = simulate_truth()
    fraction_marks, injection, log = marks(result)
    columns = unicorn_columns(curves(project, result), fraction_marks, injection, log)
    write_unicorn(columns, OUT / "akta/cex_mab_hmw_run01.csv", OUT / "akta/cex_mab_hmw_run01.xlsx")
    write_assays(result, fraction_marks)
    write_tracer_runs()
    write_lge_runs()
    write_unicorn_res(OUT / "akta/cex_mab_hmw_run01.csv", OUT / "unicorn_res/cex_mab_hmw_run01.res")
    write_generic_variants()
    write_unicorn_xls(OUT / "akta/cex_mab_hmw_run01.xlsx", OUT / "akta/cex_mab_hmw_run01.xls")
    print("fixtures written to", OUT)


if __name__ == "__main__":
    main()
