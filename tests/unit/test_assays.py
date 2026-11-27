"""Fraction assays read through mapping profiles (SPEC §7.2b): two differently formatted files."""

import hashlib
from pathlib import Path

import numpy as np
import pytest
from pydantic import ValidationError

from chromsim.data.assays import AssayMappingProfile, load_mapping_profile, read_assay
from chromsim.errors import ChromSimError

ROOT = Path(__file__).parents[2]
TOTAL_PROTEIN = ROOT / "tests/data/assays/total_protein_a280.csv"
SEC = ROOT / "tests/data/assays/sec_hplc.xlsx"
TOTAL_PROTEIN_PROFILE = ROOT / "examples/assay_profiles/total_protein_a280.yaml"
SEC_PROFILE = ROOT / "examples/assay_profiles/sec_hplc_long.yaml"


def test_wide_csv_with_semicolons_decimal_commas_and_volumes():
    assay = read_assay(TOTAL_PROTEIN, TOTAL_PROTEIN_PROFILE)
    assert assay.analytes == ["total_protein"]
    assert assay.fractions[:3] == ["A1", "A2", "A3"] and len(assay.fractions) == 25
    a3 = assay.results.set_index("fraction").loc["A3"]
    assert a3["value"] == pytest.approx(0.068)  # "0,068" read with decimal comma
    assert a3["sd"] == pytest.approx(0.006)
    assert a3["unit"] == "mg/mL" and a3["loq"] == "0.05 mg/mL"
    assert assay.fraction_volumes["A1"] == pytest.approx((85e-6, 87e-6))  # "85,00" mL → m³
    assert assay.readouts == {"total_protein": ["mAb", "HMW"]}


def test_below_loq_is_flagged_not_zero():
    results = read_assay(TOTAL_PROTEIN, TOTAL_PROTEIN_PROFILE).results.set_index("fraction")
    assert results.loc["A1", "below_loq"]
    assert np.isnan(results.loc["A1", "value"])
    assert not results.loc["A3", "below_loq"]


def test_long_xlsx_with_renamed_analytes_and_no_volumes():
    assay = read_assay(SEC, SEC_PROFILE)
    assert assay.analytes == ["monomer", "hmw"]
    assert assay.fraction_volumes == {}
    assert assay.fractions[0] == "Load" and assay.fractions[-1] == "A99"  # kept, reported later
    load = assay.results[assay.results["fraction"] == "Load"].set_index("analyte")
    assert load.loc["monomer", "value"] == pytest.approx(1.9)
    assert load.loc["hmw", "value"] == pytest.approx(0.1)
    assert assay.results["below_loq"].any()  # "< 0.005" entries
    assert assay.metadata["profile"] == "SEC-HPLC (long, Excel)"


def test_assay_files_are_not_modified():
    for path, profile in ((TOTAL_PROTEIN, TOTAL_PROTEIN_PROFILE), (SEC, SEC_PROFILE)):
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assay = read_assay(path, profile)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest == assay.metadata["sha256"]


def test_missing_columns_are_reported_with_available_ones(tmp_path):
    profile = load_mapping_profile(TOTAL_PROTEIN_PROFILE).model_copy(update={"skip_rows": 0})
    with pytest.raises(ChromSimError, match="Check skip_rows"):
        read_assay(TOTAL_PROTEIN, profile)


def test_unreadable_value_names_row_and_fraction(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("Fraction,Protein\nA1,0.5\nA2,lost\n")
    profile = AssayMappingProfile(
        name="t",
        layout="wide",
        fraction_column="Fraction",
        analytes=[{"column": "Protein", "analyte": "p", "unit": "g/L"}],
        readouts={"p": ["mAb"]},
    )
    with pytest.raises(ChromSimError, match="row 3.*'lost'.*fraction 'A2'"):
        read_assay(path, profile)


def test_non_concentration_units_are_rejected(tmp_path):
    path = tmp_path / "percent.csv"
    path.write_text("Fraction,Monomer\nA1,98.5\n")
    profile = AssayMappingProfile(
        name="t",
        layout="wide",
        fraction_column="Fraction",
        analytes=[{"column": "Monomer", "analyte": "monomer", "unit": "%"}],
        readouts={"monomer": ["mAb"]},
    )
    with pytest.raises(ChromSimError, match="'%', which is not a concentration"):
        read_assay(path, profile)


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        ({"layout": "wide", "analytes": []}, "needs `analytes`"),
        ({"layout": "long", "analyte_column": None}, "needs analyte_column"),
        ({"layout": "long", "unit_column": None, "unit": None}, "unit_column` or a fixed `unit`"),
        ({"start_volume_column": "Start"}, "both start_volume_column and end_volume_column"),
    ],
)
def test_profile_validation(changes, message):
    data = {
        "name": "t",
        "layout": "long",
        "fraction_column": "F",
        "analyte_column": "A",
        "value_column": "V",
        "unit": "g/L",
        "analytes": [{"column": "X", "analyte": "x", "unit": "g/L"}],
        "readouts": {},
    }
    data.update(changes)
    with pytest.raises(ValidationError, match=message):
        AssayMappingProfile.model_validate(data)


def test_load_mapping_profile_errors(tmp_path):
    with pytest.raises(ChromSimError, match="not found"):
        load_mapping_profile(tmp_path / "none.yaml")
    bad = tmp_path / "bad.yaml"
    bad.write_text("name: x\nlayout: sideways\n")
    with pytest.raises(ChromSimError, match="Invalid mapping profile"):
        load_mapping_profile(bad)


@pytest.mark.parametrize("end_cell", ["", "n/a"])
def test_unreadable_fraction_volume_is_reported(tmp_path, end_cell):
    path = tmp_path / "volumes.csv"
    path.write_text(f"Fraction,Start,End,Protein\nA1,1.0,2.0,0.5\nA2,2.0,{end_cell},0.4\n")
    profile = AssayMappingProfile(
        name="t",
        layout="wide",
        fraction_column="Fraction",
        start_volume_column="Start",
        end_volume_column="End",
        analytes=[{"column": "Protein", "analyte": "p", "unit": "g/L"}],
        readouts={"p": ["mAb"]},
    )
    with pytest.raises(ChromSimError, match="row 3.*volume.*fraction 'A2'"):
        read_assay(path, profile)
