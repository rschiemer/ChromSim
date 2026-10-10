"""Chromatography system (instrument) dead volumes and detector delays.

The simulation backend places these in the flow sheet as (SPEC §7.1)::

    inlet → mixer (CSTR) → tubing, first half (DPFR) → column → tubing, second half (DPFR) → outlet

Zero volumes are left out. Tubing diameter and axial dispersion only shape the small tubing
contribution to peak broadening; the defaults are placeholders to be replaced by a system tracer
calibration (use case U3).
"""

from pydantic import BaseModel

from chromsim import units
from chromsim.core import MODEL_CONFIG
from chromsim.units import Diffusivity, Length, Quantity, Volume

_ZERO_VOLUME = Quantity(magnitude=0.0, unit="mL")
# Typical inner diameter of ÄKTA system tubing (0.25–0.75 mm).
_DEFAULT_TUBING_DIAMETER = Quantity(magnitude=0.5, unit="mm")
# Placeholder until calibrated. At ~0.5 mL/min in 0.5 mm tubing this gives nearly plug flow;
# the mixer dominates system broadening.
_DEFAULT_TUBING_AXIAL_DISPERSION = Quantity(magnitude=1e-6, unit="m^2/s")


class System(BaseModel):
    """Extra-column volumes of the chromatography system.

    Attributes
    ----------
    tubing_volume : Quantity
        Total tubing volume, e.g. ``"0.25 mL"``. Default 0 mL.
    mixer_volume : Quantity
        Mixer volume, e.g. ``"0.2 mL"``. Default 0 mL.
    uv_delay_volume : Quantity
        Volume between column outlet and UV detector. Default 0 mL.
    conductivity_delay_volume : Quantity
        Volume between column outlet and conductivity detector. Default 0 mL.
    tubing_diameter : Quantity
        Inner diameter of the tubing, sets the tubing length L_t = V_t / (π d_t²/4).
        Default 0.5 mm.
    tubing_axial_dispersion : Quantity
        Axial dispersion coefficient in the tubing. Default 1e-6 m²/s (placeholder; calibrate).
    """

    model_config = MODEL_CONFIG

    tubing_volume: Volume = _ZERO_VOLUME
    mixer_volume: Volume = _ZERO_VOLUME
    uv_delay_volume: Volume = _ZERO_VOLUME
    conductivity_delay_volume: Volume = _ZERO_VOLUME
    tubing_diameter: Length = _DEFAULT_TUBING_DIAMETER
    tubing_axial_dispersion: Diffusivity = _DEFAULT_TUBING_AXIAL_DISPERSION

    @property
    def tubing_volume_m3(self) -> float:
        """Tubing volume, m³."""
        return units.to_si(self.tubing_volume, "volume")

    @property
    def mixer_volume_m3(self) -> float:
        """Mixer volume, m³."""
        return units.to_si(self.mixer_volume, "volume")

    @property
    def uv_delay_volume_m3(self) -> float:
        """Delay volume to the UV detector, m³."""
        return units.to_si(self.uv_delay_volume, "volume")

    @property
    def conductivity_delay_volume_m3(self) -> float:
        """Delay volume to the conductivity detector, m³."""
        return units.to_si(self.conductivity_delay_volume, "volume")

    @property
    def tubing_diameter_m(self) -> float:
        """Tubing inner diameter, m."""
        return units.to_si(self.tubing_diameter, "length")

    @property
    def tubing_axial_dispersion_m2_per_s(self) -> float:
        """Tubing axial dispersion coefficient, m²/s."""
        return units.to_si(self.tubing_axial_dispersion, "diffusivity")
