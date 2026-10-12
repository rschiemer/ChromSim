"""Command-line interface: ``chromsim <command>``.

Like the UI, the CLI only calls functions from :mod:`chromsim.api` (SPEC §4). Commands:
``chromsim doctor`` (environment check) and ``chromsim ui`` (local web UI).
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
    ui_parser = subcommands.add_parser("ui", help="start the local web UI (optional extra)")
    ui_parser.add_argument("--workspace", default="chromsim-workspace", help="workspace folder")
    ui_parser.add_argument("--port", type=int, default=8080)
    ui_parser.add_argument(
        "--host",
        default="127.0.0.1",
        help="address to listen on; the default accepts only this computer (the UI has no login)",
    )
    ui_parser.add_argument("--native", action="store_true", help="open in a desktop window")
    ui_parser.add_argument("--no-browser", action="store_true", help="do not open a browser")
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
    if args.command == "ui":
        try:
            api.run_ui(
                args.workspace,
                port=args.port,
                native=args.native,
                show=not args.no_browser,
                host=args.host,
            )
        except api.UIUnavailableError as error:
            print(error)
            return 1
        return 0
    return 2  # unreachable: argparse rejects unknown commands


if __name__ == "__main__":
    raise SystemExit(main())
