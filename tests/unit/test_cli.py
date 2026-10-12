"""``chromsim doctor`` command-line behaviour, with the CADET check replaced by a fake report."""

import pytest

from chromsim import api, cli, diagnostics
from chromsim.errors import EnvironmentCheckError

VERSIONS = {
    "Python": "3.12.0",
    "ChromSim": "0.0.1",
    "CADET-Process": "0.12.0",
    "CADET-Python": "1.3.0",
    "CADET-Core": "5.1.1",
}


def _healthy_report():
    return api.DoctorReport(
        versions=dict(VERSIONS),
        platform="TestOS x86_64",
        smoke_test=api.SmokeTestResult(wall_time_s=0.1, mass_recovery=1.0),
    )


def test_doctor_prints_versions_and_exits_zero(monkeypatch, capsys):
    monkeypatch.setattr(api, "doctor", _healthy_report)
    exit_code = cli.main(["doctor"])
    output = capsys.readouterr().out
    assert exit_code == 0
    for name, version in VERSIONS.items():
        assert name in output
        assert version in output
    assert "CADET test simulation: OK" in output
    assert "All checks passed." in output


def test_doctor_reports_problem_and_exits_one(monkeypatch, capsys):
    def broken_report():
        report = _healthy_report()
        report.smoke_test = None
        report.problems.append("CADET-Core could not be found. Install it with `pixi install`.")
        return report

    monkeypatch.setattr(api, "doctor", broken_report)
    exit_code = cli.main(["doctor"])
    output = capsys.readouterr().out
    assert exit_code == 1
    assert "PROBLEM: CADET-Core could not be found" in output
    assert "FAILED" in output


def test_missing_command_is_rejected():
    with pytest.raises(SystemExit) as excinfo:
        cli.main([])
    assert excinfo.value.code == 2


def test_run_doctor_collects_smoke_test_failure(monkeypatch):
    def failing_smoke_test():
        raise EnvironmentCheckError("simulated failure; reinstall with `pixi install`")

    monkeypatch.setattr(diagnostics, "find_cadet_core_version", lambda: "5.1.1")
    monkeypatch.setattr(diagnostics, "run_cadet_smoke_test", failing_smoke_test)
    report = diagnostics.run_doctor()
    assert not report.ok
    assert report.versions["CADET-Core"] == "5.1.1"
    assert report.problems == ["simulated failure; reinstall with `pixi install`"]


def test_run_doctor_reports_missing_cadet_core(monkeypatch):
    def missing_cadet():
        raise EnvironmentCheckError("CADET-Core could not be found")

    monkeypatch.setattr(diagnostics, "find_cadet_core_version", missing_cadet)
    report = diagnostics.run_doctor()
    assert not report.ok
    assert report.smoke_test is None
    assert report.versions["CADET-Core"] == "not found"
    assert "Python" in report.versions  # versions are still reported


def test_ui_without_nicegui_explains_how_to_install(monkeypatch, capsys):
    monkeypatch.setattr(api, "ui_available", lambda: False)
    assert cli.main(["ui", "--no-browser"]) == 1
    assert "pip install 'chromsim[ui]'" in capsys.readouterr().out


def test_ui_passes_options_to_run_ui(monkeypatch):
    calls = []
    monkeypatch.setattr(api, "run_ui", lambda *args, **kwargs: calls.append((args, kwargs)))
    assert cli.main(["ui", "--workspace", "ws", "--port", "9000", "--no-browser"]) == 0
    assert calls == [(("ws",), {"port": 9000, "native": False, "show": False, "host": "127.0.0.1"})]
