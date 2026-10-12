"""The ChromSim web app: one root function, sub-pages for the workflow (SPEC §8).

Nothing is created at import time: :func:`create_root` builds the shared job runner for a
workspace folder and returns the root function NiceGUI calls once per browser session. Each
session gets its own :class:`Workspace` (same folder, shared runner), so there is no module-level
state (SPEC §8.1).
"""

from pathlib import Path

from nicegui import app, ui

from chromsim import api
from chromsim.ui import components
from chromsim.ui.pages import calibrate, data, optimize, project, results, simulate


def create_root(workspace_dir: str | Path, runner=None, max_workers: int = 2):
    """Root function for ``ui.run`` serving the workspace at ``workspace_dir``.

    Without ``runner`` a :class:`LocalProcessRunner` with ``max_workers`` processes is created and
    shut down with the app; a runner passed in is owned (and shut down) by the caller.
    """
    storage = api.FileStorage(workspace_dir)
    if runner is None:
        runner = api.LocalProcessRunner(storage, max_workers=max_workers)
        app.on_shutdown(runner.shutdown)

    def root():
        workspace = api.Workspace(storage.root, runner=runner)
        components.header(workspace)
        with ui.column().classes("w-full max-w-6xl mx-auto p-4 gap-4"):
            ui.sub_pages(
                {
                    "/": lambda: project.page(workspace),
                    "/data": lambda: data.page(workspace),
                    "/simulate": lambda: simulate.page(workspace),
                    "/calibrate": lambda: calibrate.page(workspace),
                    "/optimize": lambda: optimize.page(workspace),
                    "/results": lambda: results.page(workspace),
                }
            )

    return root


def run(
    workspace_dir: str | Path = "chromsim-workspace",
    port: int = 8080,
    native: bool = False,
    show: bool = True,
    host: str = "127.0.0.1",
):
    """Start the app (blocks). Use ``chromsim ui`` or ``api.run_ui`` rather than calling this."""
    ui.run(
        create_root(workspace_dir),
        title="ChromSim",
        host=host,
        port=port,
        native=native,
        show=show,
        reload=False,
        show_welcome_message=True,
    )
