"""``chromsim.api`` exposes everything a notebook, the CLI or the UI needs."""

from pathlib import Path

from chromsim import api

EXAMPLE_PROJECT = Path(__file__).parents[2] / "examples" / "projects" / "cex_mab_lge.yaml"


def test_all_names_exist():
    for name in api.__all__:
        assert hasattr(api, name), name


def test_load_project_through_api():
    project = api.load_project(EXAMPLE_PROJECT)
    assert isinstance(project, api.Project)
    assert project.component_names == ["salt", "mAb"]
