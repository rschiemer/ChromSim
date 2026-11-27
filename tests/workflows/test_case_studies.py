"""The case studies on published data run end to end (examples/case_studies/).

They need the third-party data downloaded by examples/case_studies/fetch_data.py, which are not
in the repository (licenses); without them the tests are skipped. Each script asserts its own
results (porosities, Yamamoto agreement with CADET-Process, peak positions of predictions).
"""

import runpy
from pathlib import Path

import pytest

CASE_STUDIES = Path(__file__).parents[2] / "examples" / "case_studies"
DATA = CASE_STUDIES / "data"
REQUIRED_DATA = {
    "cs01_unicorn_imports.py": ("unicorn5_res", "unicorn7_csv", "unicorn5_asc", "unicorn5_xls"),
    "cs02_lysozyme_cex_calibration.py": ("fraunhofer_lysozyme",),
    "cs03_yamamoto_cross_check.py": ("cadet_workshop_yamamoto",),
}


def test_every_case_study_lists_its_data():
    scripts = {path.name for path in CASE_STUDIES.glob("cs*.py")}
    assert scripts == set(REQUIRED_DATA)


def test_fetch_manifest_is_pinned():
    manifest = runpy.run_path(str(CASE_STUDIES / "fetch_data.py"))
    for dataset, _repository, commit, path, sha256 in manifest["FILES"]:
        assert dataset in manifest["LICENSES"], f"{path}: no license recorded for {dataset}"
        assert len(commit) == 40 and len(sha256) == 64, f"{path}: not pinned"


@pytest.mark.slow
@pytest.mark.parametrize("script", sorted(REQUIRED_DATA))
def test_case_study_runs(script):
    missing = [name for name in REQUIRED_DATA[script] if not (DATA / name).is_dir()]
    if missing:
        pytest.skip(f"case-study data not fetched: {', '.join(missing)}")
    runpy.run_path(str(CASE_STUDIES / script), run_name="__main__")
