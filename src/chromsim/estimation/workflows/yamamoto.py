"""SMA parameters ν and k_eq from linear-gradient elution (LGE) runs: the Yamamoto method.

Several low-load runs with different gradient slopes are needed (at least two; three or more to
judge the fit). For each run you need the gradient slope g (salt increase per pumped volume) and
the salt concentration I_R at the column outlet when the protein peak elutes.

Equations
---------
At low load the SMA isotherm is linear with distribution coefficient
K(c_s) = k_eq (Λ/c_s)^ν (bound per solid volume over pore liquid). A protein zone in a linear
gradient then obeys (Yamamoto et al. 1983; derivation for ChromSim's model in
docs/theory/yamamoto.md)

.. math:: \\int_{c_0}^{I_R} \\frac{dc_s}{K(c_s)} = g\\,(1-\\varepsilon_t)\\,V_c \\equiv GH

and for I_R ≫ c_0

.. math:: \\log_{10} GH = (\\nu + 1)\\log_{10} I_R
          - \\log_{10}\\big((\\nu+1)\\,k_{eq}\\,\\Lambda^{\\nu}\\big)

A straight-line fit of log GH against log I_R gives ν = slope − 1 and
k_eq = 10^(−intercept) / ((ν+1) Λ^ν). Concentrations are in mol/m³ (= mM).

Assumptions: low load (no steric effects), the protein accesses the same pore volume as the
salt, I_R ≫ starting salt c_0, and mass transfer does not shift the peak maximum much. With
strong mass-transfer resistance the method gives good starting values for an inverse fit.

References: Yamamoto, Nakanishi, Matsuno & Kamijubo, Biotechnol. Bioeng. 25 (1983) 1465;
Yamamoto, Nakanishi & Matsuno, *Ion-Exchange Chromatography of Proteins*, Dekker 1988;
Pedersen, Mollerup, Hansen & Jungbauer, J. Chromatogr. B 790 (2003) 161.
"""

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

from chromsim.analysis.peaks import peak_maximum
from chromsim.core.method import ResolvedStep
from chromsim.data.experiment import Signal
from chromsim.errors import ChromSimError

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class LgeRun:
    """One linear-gradient elution run.

    Attributes
    ----------
    name : str
        Label.
    gradient_slope : float
        g, increase of inlet salt per pumped volume, mol/m³ per m³.
    elution_salt : float
        I_R, outlet salt concentration at the protein peak maximum, mol/m³.
    """

    name: str
    gradient_slope: float
    elution_salt: float


def gradient_slope(step: ResolvedStep, salt_index: int = 0) -> float:
    """Gradient slope g of a method step, mol/m³ per m³ of pumped volume."""
    salt_change = step.inlet_end_mol_per_m3[salt_index] - step.inlet_start_mol_per_m3[salt_index]
    volume = step.flow_rate_m3_per_s * step.duration_s
    if salt_change <= 0:
        raise ChromSimError(f"Step '{step.name}' is not an increasing salt gradient.")
    return salt_change / volume


def elution_salt_at_peak(
    protein: Signal,
    salt: Signal,
    protein_delay_m3: float = 0.0,
    salt_delay_m3: float = 0.0,
) -> float:
    """Salt concentration at the column outlet when the protein peak maximum leaves the column.

    Parameters
    ----------
    protein : Signal
        Protein signal (e.g. baseline-corrected UV) on the volume axis.
    salt : Signal
        Salt concentration (mol/m³), e.g. from :func:`conductivity_to_salt`, on the volume axis.
    protein_delay_m3, salt_delay_m3 : float
        Detector delay volumes behind the column outlet, m³; the difference is corrected.

    Returns
    -------
    float
        I_R, mol/m³.
    """
    if salt.unit != "mol/m**3":
        raise ChromSimError(f"Salt signal is in {salt.unit}; convert it with conductivity_to_salt.")
    peak_x, _ = peak_maximum(protein.x, protein.values)
    outlet_volume = peak_x - protein_delay_m3
    return float(np.interp(outlet_volume + salt_delay_m3, salt.x, salt.values))


@dataclass
class YamamotoResult:
    """Result of the Yamamoto regression.

    Attributes
    ----------
    nu : float
        Characteristic charge ν.
    k_eq : float
        SMA equilibrium constant (dimensionless).
    slope, intercept : float
        Of log10 GH vs log10 I_R.
    r_squared : float
        Coefficient of determination of the regression.
    table : pandas.DataFrame
        Per run: ``name``, ``gradient_slope``, ``gh``, ``elution_salt``, ``log10_gh``,
        ``log10_elution_salt``, ``fitted_log10_gh``, ``residual``.
    """

    nu: float
    k_eq: float
    slope: float
    intercept: float
    r_squared: float
    table: pd.DataFrame


def yamamoto_fit(
    runs: list[LgeRun],
    column_volume_m3: float,
    total_porosity: float,
    ionic_capacity_mol_per_m3: float,
) -> YamamotoResult:
    """Estimate SMA ν and k_eq from LGE runs.

    Parameters
    ----------
    runs : list of LgeRun
        At least two runs with different gradient slopes.
    column_volume_m3 : float
        Packed-bed volume V_c, m³.
    total_porosity : float
        ε_t (from tracer analysis).
    ionic_capacity_mol_per_m3 : float
        Λ, mol/m³.

    Returns
    -------
    YamamotoResult
    """
    if len({run.gradient_slope for run in runs}) < 2:
        raise ChromSimError(
            "The Yamamoto method needs runs with at least two different gradient slopes."
        )
    if len(runs) < 3:
        logger.warning(
            "Only %d LGE runs: the regression cannot be checked for curvature.", len(runs)
        )

    gh = np.array([run.gradient_slope for run in runs]) * (1.0 - total_porosity) * column_volume_m3
    elution_salt = np.array([run.elution_salt for run in runs])
    log_gh, log_ir = np.log10(gh), np.log10(elution_salt)
    slope, intercept = np.polyfit(log_ir, log_gh, 1)
    fitted = slope * log_ir + intercept
    total_sum = np.sum((log_gh - log_gh.mean()) ** 2)
    r_squared = 1.0 - np.sum((log_gh - fitted) ** 2) / total_sum if total_sum > 0 else 1.0

    nu = slope - 1.0
    if nu <= 0:
        raise ChromSimError(
            f"Regression slope {slope:.2f} gives ν = {nu:.2f} ≤ 0; elution salt should increase "
            "with gradient slope. Check the peak positions and salt values."
        )
    k_eq = 10.0 ** (-intercept) / ((nu + 1.0) * ionic_capacity_mol_per_m3**nu)
    table = pd.DataFrame(
        {
            "name": [run.name for run in runs],
            "gradient_slope": [run.gradient_slope for run in runs],
            "gh": gh,
            "elution_salt": elution_salt,
            "log10_gh": log_gh,
            "log10_elution_salt": log_ir,
            "fitted_log10_gh": fitted,
            "residual": log_gh - fitted,
        }
    )
    return YamamotoResult(
        nu=float(nu),
        k_eq=float(k_eq),
        slope=float(slope),
        intercept=float(intercept),
        r_squared=float(r_squared),
        table=table,
    )


def lge_run_from_file(
    path, project, gradient_cv: float, calibration, gradient_step: int | None = None
) -> LgeRun:
    """An LGE run from a chromatogram file of ``project``'s method with a given gradient length.

    Parameters
    ----------
    path : str or Path
        Chromatogram with UV and conductivity.
    project : Project
        Column, method, detector delay volumes (the gradient step's length is replaced).
    gradient_cv : float
        Length of the gradient in this run, CV.
    calibration : ConductivityCalibration
        Salt ↔ conductivity.
    gradient_step : int, optional
        Index of the gradient step; default the first step with ``inlet_end``.
    """
    from chromsim.data.importers import read_chromatogram
    from chromsim.data.preprocessing import conductivity_to_salt, subtract_baseline

    if gradient_step is None:
        gradient_step = next((i for i, step in enumerate(project.method) if step.inlet_end), None)
        if gradient_step is None:
            raise ChromSimError("The project's method has no gradient step (no inlet_end).")
    method = project.model_copy(deep=True)
    method.method[gradient_step].duration = f"{gradient_cv} CV"
    slope = gradient_slope(method.resolved_method()[gradient_step])

    experiment = read_chromatogram(path)
    uv = experiment.signal("uv")
    start = float(uv.x.min())
    first_step_end = (
        method.resolved_method()[0].end_s * method.resolved_method()[0].flow_rate_m3_per_s
    )
    uv = subtract_baseline(uv, (start, max(start + 1e-9, first_step_end)))
    salt = conductivity_to_salt(experiment.signal("conductivity"), calibration)
    elution_salt = elution_salt_at_peak(
        uv, salt, project.system.uv_delay_volume_m3, project.system.conductivity_delay_volume_m3
    )
    return LgeRun(f"{experiment.name} ({gradient_cv:g} CV)", slope, elution_salt)
