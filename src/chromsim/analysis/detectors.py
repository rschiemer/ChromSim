"""Predicted detector signals from a simulation (SPEC §7.1).

The simulated outlet is converted to what the instruments would record, so it can be overlaid
on measured chromatograms:

- **UV**, by Lambert–Beer summed over components with an extinction coefficient:

  .. math:: A = l \\sum_i \\varepsilon_i \\, M_i \\, c_i

  with path length l (m), mass extinction coefficient ε_i (m²/kg), molar mass M_i (kg/mol) and
  concentration c_i (mol/m³). Reference: Harris, *Quantitative Chemical Analysis*, 9th ed. 2015.
- **Conductivity**, from the salt concentration through a measured calibration curve.

Each detector sits behind its own delay volume (``System.uv_delay_volume``,
``System.conductivity_delay_volume``): the signal at pumped volume V is the outlet at
V − V_delay, i.e. the predicted curve is shifted right by V_delay.
"""

import numpy as np

from chromsim import units
from chromsim.backends.base import SimulationResult
from chromsim.core.project import Project
from chromsim.data.experiment import Signal
from chromsim.data.preprocessing import ConductivityCalibration
from chromsim.errors import ChromSimError
from chromsim.units import Quantity


def predict_uv(
    result: SimulationResult, project: Project, path_length: Quantity | str = "0.2 cm"
) -> Signal:
    """Predicted UV absorbance at the UV detector.

    Parameters
    ----------
    result : SimulationResult
        Simulation of ``project``.
    project : Project
        Supplies extinction coefficients, molecular weights and the UV delay volume.
    path_length : Quantity or str
        Flow-cell path length, default 0.2 cm.

    Returns
    -------
    Signal
        Absorbance in AU on the volume axis (displayed in mAU), kind ``"uv"``.
    """
    length = units.to_si(Quantity.parse(path_length), "length")
    absorbance = np.zeros_like(result.time_s)
    contributing = []
    for index, component in enumerate(project.components):
        epsilon = component.extinction_coeff_m2_per_kg
        molar_mass = component.molar_mass_kg_per_mol
        if epsilon is None:
            continue
        if molar_mass is None:
            raise ChromSimError(
                f"Component '{component.name}' has an extinction coefficient but no mw_kDa; "
                "both are needed to predict UV."
            )
        absorbance += length * epsilon * molar_mass * result.outlet_mol_per_m3[:, index]
        contributing.append(component.name)
    if not contributing:
        raise ChromSimError("No component has extinction_coeff_L_per_g_cm; cannot predict UV.")
    return Signal(
        name="UV (simulated)",
        kind="uv",
        x_axis="volume",
        x=result.volume_m3 + project.system.uv_delay_volume_m3,
        values=absorbance,
        unit="AU",
        original_unit="mAU",
    )


def predict_conductivity(
    result: SimulationResult,
    project: Project,
    calibration: ConductivityCalibration,
    salt_component: str = "salt",
) -> Signal:
    """Predicted conductivity at the conductivity detector.

    Parameters
    ----------
    result : SimulationResult
        Simulation of ``project``.
    project : Project
        Supplies the conductivity delay volume.
    calibration : ConductivityCalibration
        Salt ↔ conductivity calibration of the buffer system.
    salt_component : str
        Name of the salt component, default ``"salt"``.

    Returns
    -------
    Signal
        Conductivity in S/m on the volume axis (displayed in mS/cm), kind ``"conductivity"``.
    """
    if salt_component not in result.component_names:
        raise ChromSimError(
            f"No component '{salt_component}' to predict conductivity from; components are "
            f"{', '.join(result.component_names)}. Pass salt_component=..."
        )
    return Signal(
        name="Conductivity (simulated)",
        kind="conductivity",
        x_axis="volume",
        x=result.volume_m3 + project.system.conductivity_delay_volume_m3,
        values=calibration.conductivity_from_salt(result.outlet(salt_component)),
        unit="S/m",
        original_unit="mS/cm",
    )
