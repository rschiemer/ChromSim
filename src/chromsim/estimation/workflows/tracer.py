"""Column and system parameters from tracer pulses (SPEC §7.3, use case U3).

Three kinds of pulse are analysed by their moments (on the volume axis):

- **system** (no column): first moment → system volume V_sys; variance → mixer volume;
- **non-pore-penetrating** tracer (e.g. large dextran): retention volume → ε_c; variance →
  axial dispersion;
- **pore-penetrating** tracer (e.g. acetone, or a salt pulse at high salt): retention volume → ε_t.

Equations
---------
Column retention volume, after removing the system and half the injection volume:

.. math:: V_R = \\mu_1 - V_{sys} - V_{inj}/2, \\qquad
          \\varepsilon_c = V_R^{np}/V_c, \\quad \\varepsilon_t = V_R^{p}/V_c, \\quad
          \\varepsilon_p = \\frac{\\varepsilon_t - \\varepsilon_c}{1 - \\varepsilon_c}

System as a stirred tank (mixer) in series with plug flow (tubing): σ²_sys = V_mix², so
V_mix = σ_sys and V_tub = V_sys − V_mix.

Axial dispersion from the column variance of the non-penetrating tracer,
σ²_col = σ² − σ²_sys − V_inj²/12, with the closed–closed dispersion model

.. math:: \\frac{\\sigma^2_{col}}{V_R^2} = \\frac{2}{Pe} - \\frac{2}{Pe^2}(1 - e^{-Pe})

.. math:: D_{ax} = \\frac{u_{int} L}{Pe}, \\qquad u_{int} = \\frac{Q}{A\\,\\varepsilon_c}

References: Levenspiel, *Chemical Reaction Engineering*, 3rd ed., Wiley 1999, ch. 13 (dispersion
model); Schmidt-Traub et al., *Preparative Chromatography*, 2nd ed., Wiley-VCH 2012, ch. 6.
"""

import logging
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from scipy.optimize import brentq

from chromsim import units
from chromsim.analysis.peaks import (
    PeakMoments,
    asymmetry_factor,
    hetp,
    peak_moments,
    plate_number,
)
from chromsim.core.column import Column
from chromsim.core.project import Project
from chromsim.data.experiment import Signal
from chromsim.errors import ChromSimError
from chromsim.units import Quantity

logger = logging.getLogger(__name__)

TracerKind = Literal["system", "non_penetrating", "pore_penetrating"]

# Below this share of the measured variance the column variance is a small difference of two
# larger, noisy numbers, and D_ax from moments is flagged as unreliable.
MIN_COLUMN_VARIANCE_SHARE = 0.5


@dataclass
class TracerRun:
    """One tracer pulse.

    Attributes
    ----------
    name : str
        Label for reports.
    kind : {"system", "non_penetrating", "pore_penetrating"}
        Run without column, with a pore-excluded tracer, or with a pore-penetrating tracer.
    signal : Signal
        Baseline-free detector signal on the volume axis (m³).
    injection_start_m3 : float
        Pumped volume at which the injection started, m³.
    injection_volume_m3 : float
        Injected volume, m³ (rectangular pulse).
    window : (float, float), optional
        Volume range (m³) used for the moments; default the whole signal.
    """

    name: str
    kind: TracerKind
    signal: Signal
    injection_start_m3: float
    injection_volume_m3: float
    window: tuple[float, float] | None = None


@dataclass
class TracerPeak:
    """Moments of one tracer peak (volumes from the injection start, m³)."""

    name: str
    kind: str
    first_moment_m3: float
    variance_m6: float
    retention_volume_m3: float | None
    plate_number: float | None
    hetp_m: float | None
    asymmetry: float | None


@dataclass
class TracerAnalysis:
    """Parameters derived from tracer pulses; ``None`` where the needed pulse was not given.

    Attributes
    ----------
    system_volume_m3, mixer_volume_m3, tubing_volume_m3 : float or None
        System dead volume and its split into mixer (CSTR) and tubing (plug flow), m³.
    col_porosity, total_porosity, par_porosity : float or None
        ε_c, ε_t, ε_p (dimensionless).
    peclet : float or None
        Column Péclet number of the non-penetrating tracer.
    axial_dispersion_m2_per_s : float or None
        D_ax, m²/s.
    peaks : list of TracerPeak
        Per-run moments, HETP and asymmetry.
    notes : list of str
        Warnings about results the data do not determine well.
    """

    system_volume_m3: float | None = None
    mixer_volume_m3: float | None = None
    tubing_volume_m3: float | None = None
    col_porosity: float | None = None
    total_porosity: float | None = None
    par_porosity: float | None = None
    peclet: float | None = None
    axial_dispersion_m2_per_s: float | None = None
    peaks: list[TracerPeak] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def apply_to(self, project: Project) -> Project:
        """A copy of ``project`` with every parameter this analysis determined.

        The system volume is measured at the detector, so it includes the project's
        ``uv_delay_volume``; that part is subtracted before the rest is split into mixer and
        tubing.
        """
        updated = project.model_copy(deep=True)
        if self.col_porosity is not None:
            updated.resin.col_porosity = self.col_porosity
        if self.par_porosity is not None:
            updated.resin.par_porosity = self.par_porosity
        if self.axial_dispersion_m2_per_s is not None:
            updated.column.axial_dispersion = Quantity(
                magnitude=self.axial_dispersion_m2_per_s, unit="m^2/s"
            )
        if self.mixer_volume_m3 is not None and self.tubing_volume_m3 is not None:
            tubing = self.tubing_volume_m3 - project.system.uv_delay_volume_m3
            if tubing < 0:
                measured_ml = units.from_si(self.system_volume_m3 or 0.0, "mL")
                raise ChromSimError(
                    f"Measured system volume {measured_ml:.3f} mL is smaller than mixer + UV "
                    "delay volume; check system.uv_delay_volume."
                )
            updated.system.mixer_volume = Quantity(
                magnitude=units.from_si(self.mixer_volume_m3, "mL"), unit="mL"
            )
            updated.system.tubing_volume = Quantity(
                magnitude=units.from_si(tubing, "mL"), unit="mL"
            )
        return Project.model_validate(updated.model_dump())


def variance_ratio(peclet: float) -> float:
    """σ²/μ² of the closed–closed axial dispersion model for Péclet number Pe."""
    pe = peclet
    return 2.0 / pe - 2.0 / pe**2 * (1.0 - np.exp(-pe))


def peclet_from_variance_ratio(ratio: float) -> float:
    """Invert :func:`variance_ratio` (Pe from σ²/μ²)."""
    if not 0.0 < ratio < 1.0:
        raise ChromSimError(
            f"Column variance ratio σ²/V_R² = {ratio:.3g} is outside (0, 1); check the system "
            "correction and the baseline of the non-penetrating tracer."
        )
    return brentq(lambda pe: variance_ratio(pe) - ratio, 1e-6, 1e7)


def _moments(run: TracerRun):
    x = run.signal.x - run.injection_start_m3
    window = None
    if run.window is not None:
        window = (run.window[0] - run.injection_start_m3, run.window[1] - run.injection_start_m3)
    moments = peak_moments(x, run.signal.values, window)
    return moments.corrected_for_injection(run.injection_volume_m3)


def _asymmetry_or_none(run: TracerRun) -> float | None:
    """Asymmetry factor, or None (with a warning) if the recorded peak is too short for it."""
    try:
        return asymmetry_factor(run.signal.x, run.signal.values, window=run.window)
    except ChromSimError as error:
        logger.warning("Tracer run '%s': no asymmetry factor (%s)", run.name, error)
        return None


def analyse_tracers(
    runs: list[TracerRun],
    column: Column,
    flow_rate_m3_per_s: float,
    system_volume_m3: float | None = None,
) -> TracerAnalysis:
    """Determine system volume, porosities and axial dispersion from tracer pulses.

    Parameters
    ----------
    runs : list of TracerRun
        Any combination of system, non-penetrating and pore-penetrating pulses (several of a
        kind are averaged).
    column : Column
        Column geometry (length, diameter).
    flow_rate_m3_per_s : float
        Flow rate of the column runs, m³/s (needed for D_ax).
    system_volume_m3 : float, optional
        Known system volume, used when no system run is given.

    Returns
    -------
    TracerAnalysis
        Derived parameters and per-peak moments.

    Raises
    ------
    ChromSimError
        If column runs are given but the system volume is unknown, or a peak cannot be analysed.
    """
    analysis = TracerAnalysis()
    by_kind: dict[str, list] = {"system": [], "non_penetrating": [], "pore_penetrating": []}
    for run in runs:
        moments = _moments(run)
        by_kind[run.kind].append((run, moments))

    system_variance = 0.0
    if by_kind["system"]:
        analysis.system_volume_m3 = float(np.mean([m.mean for _, m in by_kind["system"]]))
        system_variance = float(np.mean([m.variance for _, m in by_kind["system"]]))
        mixer = float(np.sqrt(max(system_variance, 0.0)))
        if mixer > analysis.system_volume_m3:
            logger.warning(
                "System variance larger than a single stirred tank allows; tubing set to 0."
            )
            mixer = analysis.system_volume_m3
        analysis.mixer_volume_m3 = mixer
        analysis.tubing_volume_m3 = analysis.system_volume_m3 - mixer
        system_volume = analysis.system_volume_m3
    elif system_volume_m3 is not None:
        system_volume = system_volume_m3
    elif by_kind["non_penetrating"] or by_kind["pore_penetrating"]:
        raise ChromSimError(
            "Column tracer runs need the system volume: add a run without column (kind='system') "
            "or pass system_volume_m3."
        )
    else:
        system_volume = 0.0

    column_volume = column.volume_m3
    for kind, attribute in (
        ("non_penetrating", "col_porosity"),
        ("pore_penetrating", "total_porosity"),
    ):
        if by_kind[kind]:
            retention = float(np.mean([m.mean for _, m in by_kind[kind]])) - system_volume
            setattr(analysis, attribute, retention / column_volume)

    if analysis.col_porosity is not None and analysis.total_porosity is not None:
        eps_c, eps_t = analysis.col_porosity, analysis.total_porosity
        analysis.par_porosity = (eps_t - eps_c) / (1.0 - eps_c)

    if analysis.col_porosity is not None:
        retention = analysis.col_porosity * column_volume
        measured_variance = float(np.mean([m.variance for _, m in by_kind["non_penetrating"]]))
        column_variance = measured_variance - system_variance
        share = column_variance / measured_variance
        if share < MIN_COLUMN_VARIANCE_SHARE:
            analysis.notes.append(
                f"The column contributes only {max(share, 0.0):.0%} of the measured peak variance; "
                "the rest comes from the system (mixer, tubing) and the injection. D_ax from "
                "moments is therefore unreliable: refine it with an inverse fit of the "
                "non-penetrating tracer pulse."
            )
            logger.warning(analysis.notes[-1])
        try:
            analysis.peclet = peclet_from_variance_ratio(column_variance / retention**2)
        except ChromSimError as error:
            analysis.notes.append(f"No axial dispersion from moments: {error}")
            logger.warning(analysis.notes[-1])
        else:
            u_int = flow_rate_m3_per_s / (column.cross_section_m2 * analysis.col_porosity)
            analysis.axial_dispersion_m2_per_s = u_int * column.length_m / analysis.peclet

    for kind, runs_of_kind in by_kind.items():
        for run, moments in runs_of_kind:
            column_moments = None
            if kind != "system":
                column_moments = PeakMoments(
                    moments.area, moments.mean - system_volume, moments.variance - system_variance
                )
            analysis.peaks.append(
                TracerPeak(
                    name=run.name,
                    kind=kind,
                    first_moment_m3=moments.mean,
                    variance_m6=moments.variance,
                    retention_volume_m3=None if column_moments is None else column_moments.mean,
                    plate_number=None if column_moments is None else plate_number(column_moments),
                    hetp_m=None
                    if column_moments is None
                    else hetp(column_moments, column.length_m),
                    asymmetry=_asymmetry_or_none(run),
                )
            )
    return analysis


def tracer_run_from_file(
    path, kind: TracerKind, injection_volume_ml: float = 0.1, name: str | None = None
) -> TracerRun:
    """A tracer run from a chromatogram file, with automatic UV baseline and injection mark.

    The baseline is a straight line through the mean signal before the injection and over the
    last 10 % of the run (both protein-free for a single pulse). The injection starts at the
    file's first injection mark (0 if there is none).

    Parameters
    ----------
    path : str or Path
        Chromatogram file (any importer).
    kind : {"system", "non_penetrating", "pore_penetrating"}
    injection_volume_ml : float
        Injected volume, mL.
    name : str, optional
        Label; default the file name.
    """
    from chromsim.data.importers import read_chromatogram
    from chromsim.data.preprocessing import subtract_baseline

    experiment = read_chromatogram(path)
    uv = experiment.signal("uv")
    injection = experiment.injections_m3[0] if experiment.injections_m3 else 0.0
    end = float(uv.x.max())
    windows = [(end - 0.1 * (end - injection), end)]
    if injection > uv.x.min() + 1e-12:
        windows.insert(0, (float(uv.x.min()), injection))
    corrected = subtract_baseline(uv, *windows)
    volume = units.to_si(Quantity(magnitude=injection_volume_ml, unit="mL"), "volume")
    return TracerRun(name or experiment.name, kind, corrected, injection, volume)


def tracer_summary(analysis: TracerAnalysis) -> list[dict]:
    """Results of a tracer analysis in lab units, one row per quantity (for tables and reports)."""
    rows = []

    def add(label, value, unit, convert=None):
        if value is not None:
            rows.append(
                {"quantity": label, "value": convert(value) if convert else value, "unit": unit}
            )

    def ml(value):
        return units.from_si(value, "mL")

    add("system volume (to detector)", analysis.system_volume_m3, "mL", ml)
    add("mixer volume", analysis.mixer_volume_m3, "mL", ml)
    add("tubing volume", analysis.tubing_volume_m3, "mL", ml)
    add("column porosity ε_c", analysis.col_porosity, "-")
    add("total porosity ε_t", analysis.total_porosity, "-")
    add("particle porosity ε_p", analysis.par_porosity, "-")
    add("axial dispersion D_ax", analysis.axial_dispersion_m2_per_s, "m²/s")
    add("Péclet number", analysis.peclet, "-")
    return rows
