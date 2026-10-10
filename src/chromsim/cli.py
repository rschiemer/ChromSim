"""Command-line interface: ``chromsim <command>``.

Like the UI, the CLI only calls functions from :mod:`chromsim.api` (SPEC §4). Phase 0 provides
``chromsim doctor``; ``simulate``, ``fit``, ``optimize`` and ``ui`` follow in later phases.
"""

import argparse
import logging

from chromsim import api


def format_doctor_report(report: api.DoctorReport) -> str:
    """Render a doctor report as plain text for the terminal.

    Parameters
    ----------
    report : DoctorReport
        Report returned by :func:`chromsim.api.doctor`.

    Returns
    -------
    str
        Multi-line, human-readable report.
    """
    lines = ["ChromSim environment check", ""]
    name_width = max(len(name) for name in report.versions)
    for name, version in report.versions.items():
        lines.append(f"  {name:<{name_width}}  {version}")
    lines.append(f"  {'Platform':<{name_width}}  {report.platform}")
    lines.append("")

    if report.smoke_test is not None:
        lines.append(
            f"CADET test simulation: OK ({report.smoke_test.wall_time_s:.2f} s wall time, "
            f"mass recovery {report.smoke_test.mass_recovery:.4%})"
        )
    else:
        lines.append("CADET test simulation: NOT RUN")

    for problem in report.problems:
        lines.append(f"PROBLEM: {problem}")

    lines.append("")
    lines.append("All checks passed." if report.ok else "Some checks FAILED (see above).")
    return "\n".join(lines)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="chromsim", description="Mechanistic chromatography modeling with CADET."
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="show log messages")
    subcommands = parser.add_subparsers(dest="command", required=True)
    subcommands.add_parser(
        "doctor", help="print package versions and run a short CADET test simulation"
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the ChromSim command line.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments without the program name; defaults to ``sys.argv[1:]``.

    Returns
    -------
    int
        Process exit code: 0 on success, 1 if a check failed.
    """
    args = _build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING)

    if args.command == "doctor":
        report = api.doctor()
        print(format_doctor_report(report))
        return 0 if report.ok else 1
    return 2  # unreachable: argparse rejects unknown commands


if __name__ == "__main__":
    raise SystemExit(main())
