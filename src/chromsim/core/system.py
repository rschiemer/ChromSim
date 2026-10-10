"""Chromatography system (instrument) dead volumes and detector delays.

How these volumes are placed in the flow sheet (tubing as plug flow with dispersion, mixer as a
stirred tank, before and after the column) is decided by the simulation backend (SPEC §7.1).
"""

from pydantic import BaseModel

from chromsim import units
from chromsim.core import MODEL_CONFIG
from chromsim.units import Quantity, Volume

_ZERO_VOLUME = Quantity(magnitude=0.0, unit="mL")


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
    """

    model_config = MODEL_CONFIG

    tubing_volume: Volume = _ZERO_VOLUME
    mixer_volume: Volume = _ZERO_VOLUME
    uv_delay_volume: Volume = _ZERO_VOLUME
    conductivity_delay_volume: Volume = _ZERO_VOLUME

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
