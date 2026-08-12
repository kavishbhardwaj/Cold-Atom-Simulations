"""Exact hyperfine plus vector-Zeeman Hamiltonians for alkali manifolds.

Operators are constructed in the stable uncoupled ``|m_I,m_J>`` basis.  A
Clebsch--Gordan transformation exposes the coupled ``|F,m_F>`` ordering used by
the rate-equation and optical-Bloch models without choosing a field-dependent
quantization axis.
"""
from functools import lru_cache

import numpy as np
from scipy.constants import h, hbar
from sympy import S
from sympy.physics.wigner import clebsch_gordan

from .species import MU_B


def angular_momentum_operators(j: float):
    """Return dimensionless Jx, Jy, Jz in ascending-m basis."""
    m = np.arange(-j, j + 1, dtype=float)
    raising = np.zeros((len(m), len(m)), complex)
    for index, value in enumerate(m[:-1]):
        raising[index + 1, index] = np.sqrt(j * (j + 1) - value * (value + 1))
    lowering = raising.T
    return (raising + lowering) / 2, (raising - lowering) / (2j), np.diag(m)


@lru_cache(maxsize=None)
def manifold_operators(nuclear_spin: float, electronic_j: float):
    """Return I and J Cartesian operators on ``|m_I,m_J>``."""
    ix, iy, iz = angular_momentum_operators(nuclear_spin)
    jx, jy, jz = angular_momentum_operators(electronic_j)
    eye_i, eye_j = np.eye(ix.shape[0]), np.eye(jx.shape[0])
    i_ops = tuple(np.kron(op, eye_j) for op in (ix, iy, iz))
    j_ops = tuple(np.kron(eye_i, op) for op in (jx, jy, jz))
    return i_ops, j_ops


@lru_cache(maxsize=None)
def coupled_transformation(nuclear_spin: float, electronic_j: float):
    """Columns transform ascending ``|F,m_F>`` states to the uncoupled basis."""
    mi = np.arange(-nuclear_spin, nuclear_spin + 1)
    mj = np.arange(-electronic_j, electronic_j + 1)
    coupled = [(f, m) for f in range(int(abs(nuclear_spin-electronic_j)),
                                     int(nuclear_spin+electronic_j)+1)
               for m in range(-f, f+1)]
    transform = np.zeros((len(mi)*len(mj), len(coupled)), complex)
    for column, (f, m) in enumerate(coupled):
        for ii, m_i in enumerate(mi):
            for jj, m_j in enumerate(mj):
                if np.isclose(m_i + m_j, m):
                    transform[ii*len(mj)+jj, column] = float(clebsch_gordan(
                        S(nuclear_spin), S(electronic_j), S(f), S(m_i), S(m_j), S(m)))
    return transform, tuple(coupled)


def hyperfine_zeeman_hamiltonian(line, manifold: str, magnetic_field_t,
                                 *, basis: str = "coupled", angular=True):
    """Return H_hfs+H_Z for one fine-structure manifold.

    The magnetic term is ``mu_B (g_J J + g_I I).B``.  Hyperfine constants are
    specified in Hz.  The quadrupole expression is included only when both I
    and J admit a rank-2 moment.  Results are angular frequencies by default.
    """
    if manifold not in ("ground", "excited"):
        raise ValueError("manifold must be 'ground' or 'excited'")
    i = line.nuclear_spin
    j = line.ground_j if manifold == "ground" else line.excited_j
    a = line.species.ground_hyperfine_a_hz if manifold == "ground" else line.excited_hyperfine_a_hz
    b = 0.0 if manifold == "ground" else line.excited_hyperfine_b_hz
    g_j = line.ground_g_j if manifold == "ground" else line.excited_g_j
    i_ops, j_ops = manifold_operators(i, j)
    identity = np.eye(i_ops[0].shape[0], dtype=complex)
    i_dot_j = sum(ii @ jj for ii, jj in zip(i_ops, j_ops))
    matrix_hz = a * i_dot_j
    if b and i >= 1 and j >= 1:
        numerator = (3 * (i_dot_j @ i_dot_j) + 1.5 * i_dot_j -
                     i*(i+1)*j*(j+1)*identity)
        matrix_hz += b * numerator / (2*i*(2*i-1)*j*(2*j-1))
    magnetic_field_t = np.asarray(magnetic_field_t, float)
    if magnetic_field_t.shape != (3,):
        raise ValueError("magnetic_field_t must be a three-vector")
    zeeman_joule = MU_B * sum(component * (g_j*jop + line.species.nuclear_g_factor*iop)
                              for component, iop, jop in zip(magnetic_field_t, i_ops, j_ops))
    matrix = h * matrix_hz + zeeman_joule
    if basis == "coupled":
        transform, _ = coupled_transformation(i, j)
        matrix = transform.conj().T @ matrix @ transform
    elif basis != "uncoupled":
        raise ValueError("basis must be 'coupled' or 'uncoupled'")
    matrix = (matrix + matrix.conj().T) / 2
    return matrix / hbar if angular else matrix


def linear_zeeman_energies(line, manifold: str, magnetic_field_t):
    """Weak-field coupled-basis energies in angular-frequency units."""
    magnitude = np.linalg.norm(magnetic_field_t)
    states = [(f, m) for f in (line.ground_f if manifold == "ground" else line.excited_f)
              for m in range(-f, f+1)]
    return np.array([2*np.pi*line.hyperfine_energy_hz(manifold, f) +
                     MU_B/hbar*line.lande_gf(manifold, f)*m*magnitude
                     for f, m in states])
