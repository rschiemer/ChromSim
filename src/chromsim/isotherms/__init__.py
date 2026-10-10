"""Isotherm descriptors and the isotherm plugin registry (SPEC §5).

Importing this package registers the built-in isotherms.
"""

from chromsim.isotherms import cadet_native  # noqa: F401  (registers Linear, Langmuir, SMA)
