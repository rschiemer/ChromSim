"""Reference backend: a readable numpy/scipy solver for LRM and LRMP columns (SPEC §5b).

One file per operator (convection, dispersion, film transfer, binding) and one that assembles and
integrates them (``solver.py``). It exists for isotherms CADET does not have (Python isotherms),
for teaching, and to cross-check CADET (≤ 1 % NRMSE, tests/backends/test_reference_vs_cadet.py).
"""

from chromsim.backends.reference.solver import ReferenceBackend

__all__ = ["ReferenceBackend"]
