"""Native UNICORN 3–5 result files (.res).

The fixture is our own synthetic run 01 re-encoded as .res by tests/data/generate_fixtures.py,
so it must agree with the text export of the same run up to the .res resolution (0.01 mL,
0.001 mAU, 0.001 mS/cm). Real third-party .res files are checked too when the case-study data
have been fetched (examples/case_studies/fetch_data.py); otherwise those tests are skipped.
"""

from pathlib import Path

import numpy as np
import pytest

from chromsim.data.importers import read_chromatogram
from chromsim.data.importers.unicorn_res import curve_name, read_directory
from chromsim.errors import ChromSimError

DATA = Path(__file__).parents[1] / "data"
RES = DATA / "unicorn_res" / "cex_mab_hmw_run01.res"
CSV = DATA / "akta" / "cex_mab_hmw_run01.csv"
CASE_STUDY_DATA = Path(__file__).parents[2] / "examples" / "case_studies" / "data"
ML = 1e-6


@pytest.fixture(scope="module")
def both():
    return read_chromatogram(RES), read_chromatogram(CSV)


def test_res_is_detected_and_versioned(both):
    res, _ = both
    assert res.metadata["importer"] == "unicorn_res"
    assert res.metadata["format_version"] == "UNICORN 3.10"


def test_res_curve_names_lose_the_result_prefix():
    assert curve_name("2009Jun16no001:1_UV") == "UV"
    assert curve_name("KOKO038032 1660 S200:10_UV1_280nm") == "UV1_280nm"
    assert curve_name("Methods") == "Methods"


def test_res_signals_match_the_text_export(both):
    res, csv = both
    for kind, tolerance in (("uv", 1.5e-6), ("conductivity", 1.5e-4)):  # AU, S/m
        a, b = res.signal(kind), csv.signal(kind)
        assert a.unit == b.unit
        assert np.allclose(a.x, b.x, atol=0.006 * ML)
        assert np.allclose(a.values, b.values, atol=tolerance)


def test_res_marks_match_the_text_export(both):
    res, csv = both
    assert [f.label for f in res.fractions] == [f.label for f in csv.fractions]
    assert np.allclose(
        [f.start_m3 for f in res.fractions], [f.start_m3 for f in csv.fractions], atol=1e-12
    )
    assert np.allclose(res.injections_m3, csv.injections_m3)
    assert [text for _, text in res.log] == [text for _, text in csv.log]


def test_non_res_content_is_rejected(tmp_path):
    fake = tmp_path / "fake.res"
    fake.write_bytes(b"not a unicorn file" * 50)
    with pytest.raises(ChromSimError, match="not a UNICORN .res"):
        read_directory(fake.read_bytes(), fake.name)
    with pytest.raises(ChromSimError, match="No importer recognises"):
        read_chromatogram(fake)


REAL_FILES = sorted((CASE_STUDY_DATA / "unicorn5_res").glob("*.res"))


@pytest.mark.skipif(not REAL_FILES, reason="case-study data not fetched")
@pytest.mark.parametrize("path", REAL_FILES, ids=lambda p: p.name)
def test_real_res_files(path):
    experiment = read_chromatogram(path)
    uv = experiment.signal("uv")
    assert uv.original_unit == "mAU"
    assert np.all(np.diff(uv.x) >= 0)
    assert experiment.signal("conductivity").unit == "S/m"
    assert experiment.fractions
    assert experiment.log
