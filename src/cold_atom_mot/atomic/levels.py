"""Compatibility exports for the unified atomic basis."""
from .species import AtomicBasis, DipoleTransition, HyperfineState, build_atomic_basis

Rb87D2Basis = AtomicBasis
def build_rb87_d2_basis(atom=None):
    return build_atomic_basis("87Rb", "D2")
