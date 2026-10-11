"""Fit parameters (paths, scaling), residual functions and the uncertainty report."""

import numpy as np
import pytest
from project_examples import make_project
from pydantic import ValidationError

from chromsim.errors import ChromSimError, PluginRegistrationError
from chromsim.estimation.objectives import get_residual, register_residual
from chromsim.estimation.parameters import FitParameter, get_value, set_values
from chromsim.estimation.uncertainty import parameter_uncertainty
from chromsim.estimation.workflows.inverse_fit import FitResult

# --- parameter paths ------------------------------------------------------------------------------


def test_get_values_of_numbers_and_quantities():
    project = make_project()
    assert get_value(project, "isotherm.parameters.mAb.k_eq") == pytest.approx(0.01)
    assert get_value(project, "resin.col_porosity") == pytest.approx(0.36)
    assert get_value(project, "column.length", unit="mm") == pytest.approx(200.0)
    assert get_value(project, "method.1.inlet.mAb", unit="mg/mL") == pytest.approx(1.0)


def test_quantity_needs_unit_and_paths_are_checked():
    project = make_project()
    with pytest.raises(ChromSimError, match='give its unit, e.g. unit="cm"'):
        get_value(project, "column.length")
    with pytest.raises(
        ChromSimError, match="'isotherm.parameters.IgG' not found. Available at that level: mAb"
    ):
        get_value(project, "isotherm.parameters.IgG.k_eq")
    with pytest.raises(
        ChromSimError, match="'colporosity' not found. Available: name, particle_radius"
    ):
        get_value(project, "resin.colporosity")


def test_set_values_returns_validated_copy():
    project = make_project()
    updated = set_values(
        project,
        [
            ("isotherm.parameters.mAb.k_eq", 0.02, None),
            ("column.length", 25.0, "cm"),
            ("method.1.inlet.mAb", 2.0, "g/L"),
        ],
    )
    assert updated.isotherm.parameters["mAb"]["k_eq"] == pytest.approx(0.02)
    assert str(updated.column.length) == "25 cm"
    assert project.isotherm.parameters["mAb"]["k_eq"] == pytest.approx(0.01)  # original untouched
    with pytest.raises(ValidationError, match="nu = 99"):
        set_values(project, [("isotherm.parameters.mAb.nu", 99.0, None)])  # bounds still checked


def test_scaling_round_trip_and_derivative():
    lin = FitParameter(path="resin.col_porosity", lower=0.2, upper=0.6)
    log = FitParameter(path="isotherm.parameters.mAb.k_eq", lower=1e-4, upper=1e-1, scale="log")
    assert lin.to_scaled(0.3) == pytest.approx(0.25)
    assert log.to_scaled(1e-2) == pytest.approx(2 / 3)
    for parameter, value in ((lin, 0.47), (log, 3.3e-3)):
        x = parameter.to_scaled(value)
        assert parameter.from_scaled(x) == pytest.approx(value)
        h = 1e-6
        numeric = (parameter.from_scaled(x + h) - parameter.from_scaled(x - h)) / (2 * h)
        assert parameter.derivative(x) == pytest.approx(numeric, rel=1e-6)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"lower": 1.0, "upper": 1.0}, "must be below"),
        ({"lower": 0.0, "upper": 1.0, "scale": "log"}, "positive lower bound"),
        ({"lower": 0.0, "upper": 1.0, "initial": 2.0}, "outside the bounds"),
    ],
)
def test_fit_parameter_validation(kwargs, message):
    with pytest.raises(ValidationError, match=message):
        FitParameter(path="x", **kwargs)


# --- residuals ------------------------------------------------------------------------------------


def test_residual_functions():
    simulated, measured, sd = np.array([1.0, 3.0]), np.array([2.0, 4.0]), np.array([0.5, 0.5])
    assert get_residual("normalized")(simulated, measured) == pytest.approx([-0.25, -0.25])
    assert get_residual("absolute")(simulated, measured) == pytest.approx([-1.0, -1.0])
    assert get_residual("sd_weighted")(simulated, measured, sd) == pytest.approx([-2.0, -2.0])
    with pytest.raises(ChromSimError, match="positive SD"):
        get_residual("sd_weighted")(simulated, measured, np.array([0.5, np.nan]))
    with pytest.raises(ChromSimError, match="Unknown residual"):
        get_residual("huber")
    with pytest.raises(PluginRegistrationError):
        register_residual("absolute")(lambda s, m, sd=None: s - m)


# --- uncertainty: exact for a linear model --------------------------------------------------------


def linear_fit_result(t, a=2.0, b=0.5, noise=0.1, seed=1):
    """FitResult of y = a + b t fitted by least squares, built by hand (no simulation)."""
    rng = np.random.default_rng(seed)
    y = a + b * t + rng.normal(0.0, noise, t.size)
    design = np.column_stack([np.ones_like(t), t])
    theta, *_ = np.linalg.lstsq(design, y, rcond=None)
    parameters = [
        FitParameter(path="a", lower=-10, upper=10, name="a"),
        FitParameter(path="b", lower=-5, upper=5, name="b"),
    ]
    x = np.array([p.to_scaled(v) for p, v in zip(parameters, theta, strict=True)])
    jacobian = design * np.array([20.0, 10.0])  # dr/dx = dr/dθ · dθ/dx (ranges 20 and 10)
    residuals = design @ theta - y
    return (
        FitResult(
            parameters=parameters,
            experiment_names=["e"],
            values={"a": theta[0], "b": theta[1]},
            x=x,
            cost=0.5 * residuals @ residuals,
            residuals=residuals,
            jacobian=jacobian,
            projects={},
            results={},
            n_evaluations=1,
            success=True,
            message="",
        ),
        design,
        residuals,
    )


def test_linear_model_matches_ordinary_least_squares():
    t = np.linspace(0.0, 10.0, 50)
    fit, design, residuals = linear_fit_result(t)
    report = parameter_uncertainty(fit)
    s2 = residuals @ residuals / (50 - 2)
    expected_se = np.sqrt(np.diag(s2 * np.linalg.inv(design.T @ design)))
    assert report.table["std_error"].to_numpy() == pytest.approx(expected_se, rel=1e-9)
    assert report.degrees_of_freedom == 48 and report.identifiable
    # 95 % CI = value ± t(0.975, 48) · SE, t = 2.0106
    half_width = report.table["ci95_high"] - report.table["value"]
    assert half_width.to_numpy() == pytest.approx(2.010635 * expected_se, rel=1e-5)


def test_correlated_parameters_are_flagged():
    fit, *_ = linear_fit_result(np.linspace(100.0, 101.0, 30))  # intercept and slope inseparable
    report = parameter_uncertainty(fit)
    assert report.high_correlations and report.high_correlations[0][:2] == ("a", "b")
    assert report.correlation.loc["a", "b"] < -0.95


def test_log_parameters_have_positive_asymmetric_intervals_and_bounds_are_flagged():
    parameter = FitParameter(path="k", lower=1e-3, upper=1.0, scale="log", name="k")
    x = np.array([parameter.to_scaled(0.05)])
    jacobian = np.array([[0.2], [0.1], [0.3]])
    residuals = np.array([0.05, -0.04, 0.02])
    fit = FitResult(
        parameters=[parameter],
        experiment_names=["e"],
        values={"k": 0.05},
        x=x,
        cost=0.0,
        residuals=residuals,
        jacobian=jacobian,
        projects={},
        results={},
        n_evaluations=1,
        success=True,
        message="",
    )
    row = parameter_uncertainty(fit).table.iloc[0]
    assert 0 < row["ci95_low"] < 0.05 < row["ci95_high"]
    assert row["ci95_high"] / 0.05 == pytest.approx(0.05 / row["ci95_low"])  # symmetric in log
    assert not row["at_bound"]
    fit.x = np.array([0.9999])
    assert parameter_uncertainty(fit).table.iloc[0]["at_bound"]


def test_get_value_in_column_volumes():
    project = make_project()
    assert get_value(project, "method.3.duration", unit="CV") == pytest.approx(20.0)
    with pytest.raises(ChromSimError, match='unit="CV"'):
        get_value(project, "method.3.duration", unit="mL")
