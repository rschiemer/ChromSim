"""What an inverse fit compares: targets (measured data) and their residual functions.

A fit minimises the sum of squared residuals r over all targets of all experiments. Each target
returns its residuals already divided by √n (n its number of points), so a long UV trace and
ten assay values carry comparable weight; ``weight`` then sets their relative importance.

Residual functions (registry, extend with :func:`register_residual`):

- ``normalized`` (default for signals): r = (y_sim − y_meas) / max|y_meas|
- ``absolute``: r = y_sim − y_meas
- ``sd_weighted`` (default for assays): r = (y_sim − y_meas) / SD, SD from the assay (replicate
  standard deviation). Below-LOQ results only count when the simulation exceeds the LOQ.
"""

import logging
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

from chromsim import units
from chromsim.analysis.detectors import predict_conductivity, predict_uv
from chromsim.analysis.fractionation import compare_with_assay
from chromsim.backends.base import SimulationResult
from chromsim.core.project import Project
from chromsim.data.assays import FractionAssay
from chromsim.data.experiment import Fraction, Signal
from chromsim.data.preprocessing import ConductivityCalibration
from chromsim.errors import ChromSimError, PluginRegistrationError
from chromsim.units import Quantity

logger = logging.getLogger(__name__)

_RESIDUALS: dict[str, Callable] = {}


def register_residual(name: str):
    """Decorator adding a residual function ``f(simulated, measured, sd) -> r`` to the registry."""

    def add(function: Callable) -> Callable:
        if name in _RESIDUALS:
            raise PluginRegistrationError(f"A residual function '{name}' is already registered.")
        _RESIDUALS[name] = function
        return function

    return add


def get_residual(name: str) -> Callable:
    """The residual function registered under ``name``."""
    if name not in _RESIDUALS:
        raise ChromSimError(
            f"Unknown residual '{name}'. Available: {', '.join(sorted(_RESIDUALS))}."
        )
    return _RESIDUALS[name]


@register_residual("normalized")
def normalized_residual(simulated, measured, sd=None):
    scale = np.max(np.abs(measured))
    if scale == 0:
        raise ChromSimError("Cannot normalise residuals: the measured signal is zero everywhere.")
    return (simulated - measured) / scale


@register_residual("absolute")
def absolute_residual(simulated, measured, sd=None):
    return simulated - measured


@register_residual("sd_weighted")
def sd_weighted_residual(simulated, measured, sd):
    if sd is None or np.any(~np.isfinite(sd)) or np.any(sd <= 0):
        raise ChromSimError("sd_weighted residuals need a positive SD for every value.")
    return (simulated - measured) / sd


@dataclass
class SignalTarget:
    """A measured detector trace to fit.

    Attributes
    ----------
    signal : Signal
        Measured, preprocessed signal on the volume axis: UV in AU (baseline-corrected), or salt in
        mol/m³ (from :func:`conductivity_to_salt`), or conductivity in S/m.
    window : (float, float), optional
        Volume range to fit, m³; default the whole signal.
    weight : float
        Relative importance.
    residual : str
        Name of the residual function.
    path_length : str
        UV flow-cell path length (UV targets).
    calibration : ConductivityCalibration, optional
        Needed for conductivity targets.
    """

    signal: Signal
    window: tuple[float, float] | None = None
    weight: float = 1.0
    residual: str = "normalized"
    path_length: str = "0.2 cm"
    calibration: ConductivityCalibration | None = None

    def simulated(self, result: SimulationResult, project: Project) -> Signal:
        """The simulated counterpart of the measured signal."""
        kind = self.signal.kind
        if kind == "uv":
            return predict_uv(result, project, Quantity.parse(self.path_length))
        if kind == "salt":
            return Signal(
                "salt (simulated)",
                "salt",
                "volume",
                result.volume_m3 + project.system.conductivity_delay_volume_m3,
                result.outlet("salt"),
                "mol/m**3",
                "mM",
            )
        if kind == "conductivity":
            if self.calibration is None:
                raise ChromSimError(
                    "Conductivity targets need calibration=ConductivityCalibration(...)."
                )
            return predict_conductivity(result, project, self.calibration)
        raise ChromSimError(f"Cannot fit signals of kind '{kind}'; use uv, salt or conductivity.")

    def residuals(self, result: SimulationResult, project: Project) -> np.ndarray:
        measured = self.signal
        inside = np.ones(measured.x.size, dtype=bool)
        if self.window is not None:
            inside = (measured.x >= self.window[0]) & (measured.x <= self.window[1])
        simulated = self.simulated(result, project)
        y_sim = np.interp(measured.x[inside], simulated.x, simulated.values)
        r = get_residual(self.residual)(y_sim, measured.values[inside], None)
        return self.weight * r / np.sqrt(r.size)


@dataclass
class AssayTarget:
    """Fraction assay results to fit.

    Attributes
    ----------
    assay : FractionAssay
        Measured results; analytes are mapped to components by ``assay.readouts``.
    fractions : list of Fraction, optional
        Collected fractions (needed if the assay file has no volumes).
    analytes : list of str, optional
        Readouts to fit; default all with a readout mapping.
    weight : float
        Relative importance.
    residual : str
        Name of the residual function (default ``sd_weighted``).
    delay_volume_m3 : float
        Volume between simulated outlet and fraction collector, m³.
    """

    assay: FractionAssay
    fractions: list[Fraction] | None = None
    analytes: list[str] | None = None
    weight: float = 1.0
    residual: str = "sd_weighted"
    delay_volume_m3: float = 0.0

    def __post_init__(self):
        # Report fraction-id mismatches once, not at every fit evaluation.
        if not self.assay.fraction_volumes and self.fractions is not None:
            collected = {fraction.label for fraction in self.fractions}
            missing = [label for label in self.assay.fractions if label not in collected]
            if missing:
                logger.warning(
                    "Assay '%s': fraction(s) %s have no collected fraction and are not fitted.",
                    self.assay.name,
                    ", ".join(missing),
                )

    def residuals(self, result: SimulationResult, project: Project) -> np.ndarray:
        table = compare_with_assay(
            result, self.assay, self.fractions, self.delay_volume_m3, warn=False
        ).table
        if self.analytes is not None:
            table = table[table["analyte"].isin(self.analytes)]
        residual = get_residual(self.residual)
        values = []
        for row in table.itertuples():
            if not row.below_loq:
                values.append(
                    residual(
                        np.array([row.simulated]), np.array([row.measured]), np.array([row.sd])
                    )[0]
                )
                continue
            # Below LOQ: no residual while the simulation stays below the LOQ, a penalty above it.
            # Always one entry per row, so the residual vector keeps its length during a fit.
            loq = _loq_in_unit(self.assay, row.analyte, row.unit)
            if loq is None or row.simulated <= loq:
                values.append(0.0)
            else:
                sd = row.sd if np.isfinite(row.sd) and row.sd > 0 else 0.5 * loq
                values.append(
                    residual(np.array([row.simulated]), np.array([loq]), np.array([sd]))[0]
                )
        r = np.array(values)
        return self.weight * r / np.sqrt(max(r.size, 1))


def _loq_in_unit(assay: FractionAssay, analyte: str, unit: str) -> float | None:
    rows = assay.results[(assay.results["analyte"] == analyte) & assay.results["loq"].notna()]
    if rows.empty:
        return None
    loq = Quantity.parse(rows["loq"].iloc[0])
    return float(units.ureg.Quantity(loq.magnitude, loq.unit).to(unit).magnitude)


@dataclass
class FitExperiment:
    """One experiment in a fit: its project (method, load, …) and the data measured in it."""

    name: str
    project: Project
    targets: list = field(default_factory=list)
