"""Parameter uncertainty from the Jacobian at the optimum (SPEC §7.3, use case U5).

Linearised covariance in the scaled space, with the residual variance estimated from the fit:

.. math:: s^2 = \\frac{\\sum r_i^2}{m - p}, \\qquad C_x = s^2 \\,(J^T J)^{-1}

transformed to parameter units with D = diag(dθ/dx): C_θ = D C_x D. 95 % confidence intervals use
Student's t with m − p degrees of freedom. For log-scaled parameters the interval is computed in
log space (asymmetric in θ, always positive).

The correlation matrix R_ij = C_ij / √(C_ii C_jj) shows which parameters the data cannot separate:
|R_ij| > 0.95 is flagged. Parameters at a bound are flagged too, because the linearisation does
not hold there.

Reference: Bates & Watts, *Nonlinear Regression Analysis and Its Applications*, Wiley 1988, ch. 2;
Seber & Wild, *Nonlinear Regression*, Wiley 1989, ch. 5.
"""

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from chromsim.estimation.workflows.inverse_fit import FitResult

HIGH_CORRELATION = 0.95
# A parameter whose scaled value is within this distance of 0 or 1 is reported as "at bound".
AT_BOUND_TOLERANCE = 1e-3


@dataclass
class UncertaintyReport:
    """Confidence intervals and correlations of a fit.

    Attributes
    ----------
    table : pandas.DataFrame
        Per parameter: ``parameter``, ``value``, ``std_error``, ``ci95_low``, ``ci95_high``,
        ``at_bound``.
    correlation : pandas.DataFrame
        Correlation matrix (labels as index and columns).
    high_correlations : list of (str, str, float)
        Pairs with |correlation| > 0.95.
    identifiable : bool
        ``False`` if JᵀJ is singular (some parameter has no influence on the residuals).
    degrees_of_freedom : int
        m − p.
    """

    table: pd.DataFrame
    correlation: pd.DataFrame
    high_correlations: list[tuple[str, str, float]] = field(default_factory=list)
    identifiable: bool = True
    degrees_of_freedom: int = 0


def parameter_uncertainty(fit: FitResult, confidence: float = 0.95) -> UncertaintyReport:
    """Linearised confidence intervals and correlation matrix of a fit.

    Parameters
    ----------
    fit : FitResult
        Result of :func:`inverse_fit`.
    confidence : float
        Confidence level, default 0.95.

    Returns
    -------
    UncertaintyReport
    """
    jacobian, residuals = fit.jacobian, fit.residuals
    m, p = jacobian.shape
    dof = max(m - p, 1)
    s2 = float(residuals @ residuals) / dof
    jtj = jacobian.T @ jacobian
    identifiable = np.linalg.matrix_rank(jtj) == p
    covariance_x = s2 * np.linalg.pinv(jtj)
    t_value = float(student_t.ppf(0.5 + confidence / 2, dof))

    labels = fit.x_labels
    owners = fit.x_parameters()
    rows = []
    for i, (label, parameter, x) in enumerate(zip(labels, owners, fit.x, strict=True)):
        sd_x = float(np.sqrt(max(covariance_x[i, i], 0.0)))
        value = parameter.from_scaled(x)
        if parameter.scale == "log":
            log_range = np.log10(parameter.upper) - np.log10(parameter.lower)
            sd_log = sd_x * log_range
            low, high = value * 10 ** (-t_value * sd_log), value * 10 ** (t_value * sd_log)
        else:
            low = value - t_value * sd_x * (parameter.upper - parameter.lower)
            high = value + t_value * sd_x * (parameter.upper - parameter.lower)
        rows.append(
            {
                "parameter": label,
                "value": value,
                "std_error": abs(parameter.derivative(x)) * sd_x,
                "ci95_low": low,
                "ci95_high": high,
                "at_bound": bool(x < AT_BOUND_TOLERANCE or x > 1 - AT_BOUND_TOLERANCE),
            }
        )

    sd = np.sqrt(np.clip(np.diag(covariance_x), 1e-300, None))
    correlation = covariance_x / np.outer(sd, sd)
    np.fill_diagonal(correlation, 1.0)
    high = [
        (labels[i], labels[j], float(correlation[i, j]))
        for i in range(p)
        for j in range(i + 1, p)
        if abs(correlation[i, j]) > HIGH_CORRELATION
    ]
    return UncertaintyReport(
        table=pd.DataFrame(rows),
        correlation=pd.DataFrame(correlation, index=labels, columns=labels),
        high_correlations=high,
        identifiable=bool(identifiable),
        degrees_of_freedom=dof,
    )
