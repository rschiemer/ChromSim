"""Isotherm descriptors and the isotherm plugin registry (SPEC §5).

Importing this package registers the built-in isotherms.
"""

# Importing a module registers its isotherms (Linear, Langmuir, SMA, MPM, GIEX, HIC, BiLangmuir,
# and the Python example MixedModeLangmuir).
from chromsim.isotherms import (  # noqa: F401
    bi_langmuir,
    cadet_native,
    generalized_ion_exchange,
    hydrophobic_interaction,
    mobile_phase_modulator,
)
from chromsim.isotherms.custom import mixed_mode_langmuir  # noqa: F401
