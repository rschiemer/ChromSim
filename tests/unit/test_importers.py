"""Chromatogram importers: the shared contract every importer must meet, plus format details.

Contract (CLAUDE.md rule 6, SPEC §7.2): every registered importer has a fixture here, reads it
without modifying the file, returns signals with increasing axes in SI units, and records
provenance.
"""

import hashlib
from pathlib import Path

import numpy as np
import pytest

from chromsim.data.importers import get_importer, list_importers, read_chromatogram
from chromsim.data.importers.akta_csv import parse_unicorn_table, read_text_rows
from chromsim.errors import ChromSimError

DATA = Path(__file__).parents[1] / "data"
AKTA_CSV = DATA / "akta" / "cex_mab_hmw_run01.csv"
AKTA_XLSX = DATA / "akta" / "cex_mab_hmw_run01.xlsx"
GENERIC_CSV = DATA / "generic" / "plain_export.csv"
GENERIC_MAPPING = {
    "x_column": "Volume (mL)",
    "x_unit": "mL",
    "signals": [
        {"column": "UV 280 (mAU)", "name": "UV280", "unit": "mAU"},
        {"column": "Conductivity (mS/cm)", "name": "Cond", "unit": "mS/cm"},
    ],
}

# Every registered importer needs an entry here: (fixture file, read options).
FIXTURES = {
    "akta_csv": (AKTA_CSV, {}),
    "akta_xlsx": (AKTA_XLSX, {}),
    "generic_csv": (GENERIC_CSV, {"mapping": GENERIC_MAPPING}),
    "generic_xlsx": (DATA / "generic" / "plain_export.xlsx", {"mapping": GENERIC_MAPPING}),
    "unicorn_res": (DATA / "unicorn_res" / "cex_mab_hmw_run01.res", {}),
}
INTERNAL_UNITS = {"AU", "S/m", "pH", "1", "Pa", "mol/m**3", "kg/m**3"}


def test_every_importer_has_a_fixture():
    assert set(list_importers()) == set(FIXTURES), "add a fixture for each new importer"


@pytest.mark.parametrize("importer_name", sorted(FIXTURES))
def test_importer_contract(importer_name):
    path, options = FIXTURES[importer_name]
    digest_before, mtime_before = (
        hashlib.sha256(path.read_bytes()).hexdigest(),
        path.stat().st_mtime,
    )

    experiment = get_importer(importer_name).read(path, **options)

    # raw data are read-only
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest_before
    assert path.stat().st_mtime == mtime_before
    # provenance
    assert experiment.metadata["importer"] == importer_name
    assert experiment.metadata["sha256"] == digest_before
    assert experiment.metadata["source_file"] == str(path.resolve())
    # signals in SI with increasing axes
    assert experiment.signals
    for signal in experiment.signals.values():
        assert signal.x_axis in ("volume", "time")
        assert np.all(np.diff(signal.x) >= 0), f"{signal.name}: axis not increasing"
        assert np.all(np.isfinite(signal.values))
        assert signal.unit in INTERNAL_UNITS, f"{signal.name}: unexpected unit {signal.unit}"


# --- ÄKTA UNICORN ---------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def akta_csv():
    return read_chromatogram(AKTA_CSV)


def test_akta_files_are_detected_automatically(akta_csv):
    assert akta_csv.metadata["importer"] == "akta_csv"
    assert read_chromatogram(AKTA_XLSX).metadata["importer"] == "akta_xlsx"
    assert not get_importer("akta_csv").can_read(AKTA_XLSX)
    assert not get_importer("akta_xlsx").can_read(AKTA_CSV)
    assert not get_importer("akta_csv").can_read(GENERIC_CSV)


def test_akta_signals_and_units(akta_csv):
    assert set(akta_csv.signals) == {"UV1_280nm", "Cond", "Conc B", "pH"}  # "Chrom.1_" removed
    uv = akta_csv.signal("uv")
    assert uv.name == "UV1_280nm" and uv.original_unit == "mAU"
    assert uv.x[1] == pytest.approx(0.025e-6)  # 0.025 mL in m³
    # values_in() reverses the conversion: mAU → AU → mAU
    assert uv.values_in("mAU") == pytest.approx(uv.values * 1000.0)
    cond = akta_csv.signal("conductivity")
    assert cond.unit == "S/m"
    assert cond.values_in("mS/cm") == pytest.approx(cond.values * 10.0)  # 1 mS/cm = 0.1 S/m
    assert akta_csv.signal("concentration_b").unit == "1"  # % → fraction
    assert akta_csv.signal("ph").values.mean() == pytest.approx(5.5, abs=0.01)


def test_akta_fractions_injection_and_log(akta_csv):
    labels = [fraction.label for fraction in akta_csv.fractions]
    assert labels[:3] == ["A1", "A2", "A3"] and len(labels) == 25 and "Waste" not in labels
    first = akta_csv.fraction("A1")
    assert (first.start_m3, first.end_m3) == pytest.approx((85e-6, 87e-6))
    assert akta_csv.fraction("C1").end_m3 == pytest.approx(135e-6)  # closed by the Waste mark
    assert akta_csv.injections_m3 == pytest.approx([19.64e-6], rel=1e-3)
    assert [text for _, text in akta_csv.log][:2] == ["Method Run 10/10/2026", "Block Equilibrate"]


def test_akta_csv_and_xlsx_give_the_same_experiment(akta_csv):
    xlsx = read_chromatogram(AKTA_XLSX)
    assert set(xlsx.signals) == set(akta_csv.signals)
    for name, signal in akta_csv.signals.items():
        assert xlsx.signals[name].x == pytest.approx(signal.x)
        assert xlsx.signals[name].values == pytest.approx(signal.values)
    assert xlsx.fractions == akta_csv.fractions


@pytest.mark.parametrize(
    ("delimiter", "decimal", "encoding"),
    [("\t", ".", "utf-16"), (",", ".", "utf-8"), (";", ",", "cp1252")],
)
def test_akta_text_variants(tmp_path, delimiter, decimal, encoding):
    rows = [
        ["Chrom.1_UV 1_280", "", "Chrom.1_Cond", ""],
        ["ml", "mAU", "ml", "mS/cm"],
        ["0.0", "1.5", "0.0", "6.4"],
        ["0.5", "2.5", "1.0", "6.6"],
        ["1.0", "3.5", "", ""],
    ]
    header, data = rows[:2], rows[2:]
    lines = [delimiter.join(row) for row in header]
    lines += [delimiter.join(cell.replace(".", decimal) for cell in row) for row in data]
    text = "\n".join(lines)
    path = tmp_path / "variant.csv"
    path.write_text(text, encoding=encoding)
    experiment = read_chromatogram(path)
    uv = experiment.signal("uv")
    assert uv.name == "UV 1_280"
    assert uv.x == pytest.approx([0.0, 0.5e-6, 1.0e-6])
    assert uv.values == pytest.approx([1.5e-3, 2.5e-3, 3.5e-3])
    assert len(experiment.signal("Cond").x) == 2  # curves have their own lengths


def test_text_without_unicorn_header_is_rejected():
    rows = read_text_rows(GENERIC_CSV)
    with pytest.raises(ChromSimError, match="no UNICORN header.*generic_csv"):
        parse_unicorn_table(rows, "plain", {})


def test_unknown_file_lists_importers(tmp_path):
    path = tmp_path / "notes.txt"
    path.write_text("hello")
    with pytest.raises(ChromSimError, match="No importer recognises notes.txt.*akta_csv"):
        read_chromatogram(path)
    with pytest.raises(ChromSimError, match="not found"):
        read_chromatogram(tmp_path / "missing.csv")


# --- generic CSV ----------------------------------------------------------------------------------


def test_generic_csv_with_mapping():
    experiment = read_chromatogram(GENERIC_CSV, importer="generic_csv", mapping=GENERIC_MAPPING)
    uv = experiment.signal("UV280")
    assert uv.kind == "uv" and uv.unit == "AU"
    assert uv.values.max() == pytest.approx(0.120)
    assert experiment.signal("Cond").kind == "conductivity"
    assert uv.x[-1] == pytest.approx(4e-6)


def test_generic_csv_mapping_from_yaml(tmp_path):
    mapping_file = tmp_path / "mapping.yaml"
    mapping_file.write_text(
        "x_column: Volume (mL)\nx_unit: mL\nsignals:\n  - {column: UV 280 (mAU), unit: mAU}\n"
    )
    experiment = read_chromatogram(GENERIC_CSV, importer="generic_csv", mapping=mapping_file)
    assert experiment.signal("UV 280 (mAU)").values.max() == pytest.approx(0.120)


def test_generic_csv_reports_missing_columns():
    mapping = dict(GENERIC_MAPPING, x_column="Volume [ml]")
    with pytest.raises(ChromSimError, match="'Volume \\[ml\\]' not found.*'Volume \\(mL\\)'"):
        read_chromatogram(GENERIC_CSV, importer="generic_csv", mapping=mapping)


def test_generic_csv_needs_mapping():
    with pytest.raises(ChromSimError, match="needs mapping"):
        read_chromatogram(GENERIC_CSV, importer="generic_csv")


def test_unknown_importer_options_are_rejected():
    with pytest.raises(ChromSimError, match="does not accept option\\(s\\) sheet; accepted: none"):
        read_chromatogram(AKTA_CSV, sheet="x")


def test_duplicate_importer_name_is_rejected():
    from chromsim.data.importers.registry import register_importer
    from chromsim.errors import PluginRegistrationError

    class AnotherAkta:
        name = "akta_csv"
        description = "duplicate"
        file_extensions = (".csv",)

        def can_read(self, path):
            return False

        def read(self, path, **options):
            raise NotImplementedError

    before = list_importers()
    with pytest.raises(PluginRegistrationError, match="already registered"):
        register_importer(AnotherAkta)
    assert list_importers() == before  # the failed registration changed nothing


# --- generic tables: headerless text, Excel, concentration signals -------------------------------


def test_generic_csv_headerless_whitespace_by_position():
    mapping = {
        "header": False,
        "delimiter": "whitespace",
        "x_column": 0,
        "x_unit": "mL",
        "signals": [{"column": 1, "name": "UV", "unit": "mAU"}],
    }
    headerless = read_chromatogram(
        DATA / "generic" / "headerless.dat", "generic_csv", mapping=mapping
    )
    with_header = read_chromatogram(GENERIC_CSV, "generic_csv", mapping=GENERIC_MAPPING)
    assert np.allclose(headerless.signals["UV"].x, with_header.signals["UV280"].x)
    assert np.allclose(headerless.signals["UV"].values, with_header.signals["UV280"].values)


def test_generic_csv_position_out_of_range():
    mapping = {
        "header": False,
        "delimiter": "whitespace",
        "x_column": 0,
        "x_unit": "mL",
        "signals": [{"column": 5, "unit": "mAU"}],
    }
    with pytest.raises(ChromSimError, match="position 5 does not exist"):
        read_chromatogram(DATA / "generic" / "headerless.dat", "generic_csv", mapping=mapping)


def test_generic_xlsx_matches_csv():
    from_xlsx = read_chromatogram(
        DATA / "generic" / "plain_export.xlsx", "generic_xlsx", mapping=GENERIC_MAPPING, sheet="Run"
    )
    from_csv = read_chromatogram(GENERIC_CSV, "generic_csv", mapping=GENERIC_MAPPING)
    for name in ("UV280", "Cond"):
        assert np.allclose(from_xlsx.signals[name].values, from_csv.signals[name].values)
    with pytest.raises(ChromSimError, match="cannot read sheet"):
        read_chromatogram(
            DATA / "generic" / "plain_export.xlsx",
            "generic_xlsx",
            mapping=GENERIC_MAPPING,
            sheet="nope",
        )


def test_concentration_signals_are_converted_and_classified(tmp_path):
    table = tmp_path / "conc.csv"
    table.write_text("time [s],Salt [mM],A [g/L]\n0,50,0.0\n10,60,1.5\n")
    mapping = {
        "x_column": "time [s]",
        "x_unit": "s",
        "signals": [
            {"column": "Salt [mM]", "name": "salt", "unit": "mM"},
            {"column": "A [g/L]", "name": "A", "unit": "g/L"},
        ],
    }
    experiment = read_chromatogram(table, "generic_csv", mapping=mapping)
    salt, protein = experiment.signals["salt"], experiment.signals["A"]
    assert (salt.kind, salt.unit) == ("salt", "mol/m**3")
    assert (protein.kind, protein.unit) == ("concentration", "kg/m**3")
    assert experiment.signals["salt"].x_axis == "time"


def test_unicorn5_curve_prefix_is_removed():
    rows = [
        ["Manual run 3:10_UV", "", "Manual run 3:10_Cond", ""],
        ["ml", "mAU", "ml", "mS/cm"],
        ["0.0", "1.0", "0.0", "5.0"],
        ["0.1", "2.0", "0.1", "5.1"],
    ]
    experiment = parse_unicorn_table(rows, "run", {})
    assert set(experiment.signals) == {"UV", "Cond"}
