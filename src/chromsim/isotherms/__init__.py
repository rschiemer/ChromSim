"""Isotherm descriptors and the isotherm plugin registry (SPEC §5).

Importing this package registers the built-in isotherms.
"""

# Importing a module registers its isotherms (Linear, Langmuir, SMA, MPM, GIEX, HIC, BiLangmuir,
# Colloidal, and the Python isotherms MixedModeLangmuir and CPA).
from chromsim.isotherms import (  # noqa: F401
    bi_langmuir,
    cadet_native,
    colloidal,
    generalized_ion_exchange,
    hydrophobic_interaction,
    mobile_phase_modulator,
)
from chromsim.isotherms.custom import (  # noqa: F401
    colloidal_particle_adsorption,
    mixed_mode_langmuir,
)
