"""YAML project files: lossless round trip of the SPEC §4.3 example (Phase 1 acceptance),
schema versioning and readable errors."""

import re
from pathlib import Path

import pytest
import yaml

from chromsim.core.project import CURRENT_SCHEMA_VERSION, Project
from chromsim.errors import ChromSimError, ProjectFileError, UnitError
from chromsim.io import project_file
from chromsim.io.project_file import (
    load_project,
    migrate,
    project_from_yaml,
    project_to_yaml,
    save_project,
)
from chromsim.units import Quantity

REPO_ROOT = Path(__file__).parents[2]
EXAMPLE_PROJECT = REPO_ROOT / "examples" / "projects" / "cex_mab_lge.yaml"


def spec_example_yaml() -> str:
    """The YAML block of SPEC §4.3, read from SPEC.md so the test follows the SPEC."""
    spec = (REPO_ROOT / "SPEC.md").read_text(encoding="utf-8")
    section = spec.split("### 4.3", 1)[1]
    match = re.search(r"```yaml\n(.*?)```", section, flags=re.DOTALL)
    assert match, "SPEC §4.3 has no ```yaml block"
    return match.group(1)


def _normalise(value):
    """Compare scalars by meaning: '20 cm' as a Quantity, '1.0e2' and 100 as the float 100."""
    if isinstance(value, bool):
        return value
    if isinstance(value, int | float):
        return float(value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            pass
        try:
            return Quantity.parse(value)
        except UnitError:
            return value
    return value


def assert_contains_original(original, written, path="project"):
    """Every key and value of the original file is present, with the same meaning, in the output."""
    if isinstance(original, dict):
        assert isinstance(written, dict), path
        for key, value in original.items():
            assert key in written, f"{path}.{key} was lost"
            assert_contains_original(value, written[key], f"{path}.{key}")
    elif isinstance(original, list):
        assert isinstance(written, list) and len(written) == len(original), path
        for index, (a, b) in enumerate(zip(original, written, strict=True)):
            assert_contains_original(a, b, f"{path}[{index}]")
    else:
        assert _normalise(original) == _normalise(written), f"{path}: {original!r} != {written!r}"


# --- Phase 1 acceptance: lossless round trip ------------------------------------------------------


def test_spec_example_round_trips_losslessly():
    original_text = spec_example_yaml()
    project = project_from_yaml(original_text)
    written_text = project_to_yaml(project)
    reloaded = project_from_yaml(written_text)

    assert reloaded == project  # same Project object …
    assert project_to_yaml(reloaded) == written_text  # … the output is stable …
    # … and nothing the user wrote was lost or changed.
    assert_contains_original(yaml.safe_load(original_text), yaml.safe_load(written_text))


def test_units_are_written_as_the_user_wrote_them():
    written = project_to_yaml(project_from_yaml(spec_example_yaml()))
    for text in (
        "20 cm",
        "0.5 cm",
        "45 um",
        "1200 mM",
        "0.25 mL",
        "5 CV",
        "2 mL",
        "150 cm/h",
        "1 g/L",
    ):
        assert f"'{text}'" in written or f" {text}\n" in written, text


def test_example_project_file_matches_spec():
    assert yaml.safe_load(EXAMPLE_PROJECT.read_text(encoding="utf-8")) == yaml.safe_load(
        spec_example_yaml()
    )


def test_save_and_load_file(tmp_path):
    project = load_project(EXAMPLE_PROJECT)
    path = save_project(project, tmp_path / "copy.yaml")
    assert load_project(path) == project


def test_edited_project_round_trips(tmp_path):
    project = load_project(EXAMPLE_PROJECT)
    project.column.length = "25 cm"
    project.method[3].duration = "30 CV"
    path = save_project(project, tmp_path / "edited.yaml")
    reloaded = load_project(path)
    assert reloaded == project
    assert reloaded.column.length == Quantity.parse("25 cm")


def test_pyyaml_exponent_quirk_is_handled():
    # PyYAML reads `1.0e2` (no sign in the exponent) as the *string* "1.0e2"; it must still
    # become the number 100.
    project = project_from_yaml(spec_example_yaml())
    assert project.isotherm.parameters["mAb"]["k_kin"] == 100.0
    assert isinstance(
        yaml.safe_load(project_to_yaml(project))["isotherm"]["parameters"]["mAb"]["k_kin"], float
    )


# --- schema versioning ----------------------------------------------------------------------------


def test_current_schema_version_passes_unchanged():
    data = {"schema_version": CURRENT_SCHEMA_VERSION, "name": "x"}
    assert migrate(data, source="test") == data


def test_missing_schema_version_is_an_error():
    with pytest.raises(ProjectFileError, match="schema_version"):
        migrate({"name": "x"}, source="test")


def test_newer_schema_version_is_an_error():
    with pytest.raises(ProjectFileError, match="newer ChromSim"):
        migrate({"schema_version": CURRENT_SCHEMA_VERSION + 1}, source="test")


def test_migrations_are_applied_in_order(monkeypatch):
    def rename_key(data):
        data["title"] = data.pop("name")
        return data

    monkeypatch.setattr(project_file, "CURRENT_SCHEMA_VERSION", 2)
    monkeypatch.setattr(project_file, "MIGRATIONS", {1: rename_key})
    original = {"schema_version": 1, "name": "x"}
    assert migrate(original, source="test") == {"schema_version": 2, "title": "x"}
    assert original == {"schema_version": 1, "name": "x"}  # input is not modified


def test_project_rejects_other_schema_version_in_code():
    data = yaml.safe_load(spec_example_yaml())
    data["schema_version"] = 99
    with pytest.raises(ValueError, match="schema_version"):
        Project.model_validate(data)


# --- readable errors ------------------------------------------------------------------------------


def test_missing_file(tmp_path):
    with pytest.raises(ProjectFileError, match="not found"):
        load_project(tmp_path / "nope.yaml")


def test_invalid_yaml():
    with pytest.raises(ProjectFileError, match="not valid YAML"):
        project_from_yaml("name: [unclosed")


def test_yaml_that_is_not_a_mapping():
    with pytest.raises(ProjectFileError, match="mapping"):
        project_from_yaml("- just\n- a list\n")


def test_validation_errors_name_location_and_fix():
    text = spec_example_yaml().replace("length: 20 cm", "length: 20 mL")
    with pytest.raises(ProjectFileError) as excinfo:
        project_from_yaml(text, source="my_project.yaml")
    message = str(excinfo.value)
    assert "my_project.yaml" in message
    assert "column.length: '20 mL' is not a length; use e.g. cm or mm." in message
    assert "Value error" not in message  # pydantic's prefix is stripped


def test_project_file_error_is_a_chromsim_error():
    assert issubclass(ProjectFileError, ChromSimError)
