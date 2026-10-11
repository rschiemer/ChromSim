"""Generate the synthetic raw-data fixtures in tests/data/ (run manually, not by pytest).

    pixi run python tests/data/generate_fixtures.py

The files imitate real exports of a CEX run of examples/projects/cex_mab_hmw_lge.yaml, simulated
with "true" parameters (mAb k_eq = 1.5e-2 instead of the project's first guess 1.0e-2) plus
baseline drift and noise:

- akta/cex_mab_hmw_run01.csv  ÄKTA UNICORN 7 text export: UTF-16, tab-separated, title row
- akta/cex_mab_hmw_run01.xlsx the same run as UNICORN Excel export (no title row)
- assays/total_protein_a280.csv  wide, semicolon + decimal comma, with fraction volumes, <LOQ
- assays/sec_hplc.xlsx           long format, sheet "SEC results", no volumes, includes a feed
                                 sample "Load" and a sample "A99" that was never collected

Raw data are read-only for ChromSim (CLAUDE.md); this script is the only thing that writes them.
Re-running it with the same seed reproduces the files.
"""

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


def write_unicorn(columns, csv_path, xlsx_path):
    n_rows = max(len(c[2]) for c in columns)
    header = [c[0] for c in columns]
    units_row = [c[1] for c in columns]
    data = [[c[2][i] if i < len(c[2]) else "" for c in columns] for i in range(n_rows)]

    lines = ["cex_mab_hmw_run01 001", "\t".join(header), "\t".join(units_row)]
    lines += ["\t".join(row) for row in data]
    csv_path.write_text("\r\n".join(lines) + "\r\n", encoding="utf-16")

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


def main():
    (OUT / "akta").mkdir(exist_ok=True)
    (OUT / "assays").mkdir(exist_ok=True)
    project, result = simulate_truth()
    fraction_marks, injection, log = marks(result)
    columns = unicorn_columns(curves(project, result), fraction_marks, injection, log)
    write_unicorn(columns, OUT / "akta/cex_mab_hmw_run01.csv", OUT / "akta/cex_mab_hmw_run01.xlsx")
    write_assays(result, fraction_marks)
    print("fixtures written to", OUT)


if __name__ == "__main__":
    main()
