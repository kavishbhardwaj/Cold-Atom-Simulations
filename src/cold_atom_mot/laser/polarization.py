"""Propagation-relative polarization and local spherical-basis utilities."""

import numpy as np


def unit(vector: np.ndarray) -> np.ndarray:
    vector = np.asarray(vector, dtype=float)
    norm = np.linalg.norm(vector)
    if norm == 0:
        raise ValueError("direction must be non-zero")
    return vector / norm


def transverse_basis(k_hat: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return a right-handed transverse basis (e1,e2,k_hat)."""
    k_hat = unit(k_hat)
    reference = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(reference, k_hat)) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])
    e1 = unit(np.cross(reference, k_hat))
    e2 = np.cross(k_hat, e1)
    return e1, e2


def circular_polarization(k_hat: np.ndarray, helicity: int) -> np.ndarray:
    """Complex unit polarization; helicity is defined relative to propagation."""
    if helicity not in (-1, 1):
        raise ValueError("helicity must be -1 or +1")
    e1, e2 = transverse_basis(k_hat)
    return (e1 + 1j * helicity * e2) / np.sqrt(2.0)


def spherical_fractions(polarization: np.ndarray, quantization_axis: np.ndarray) -> dict[int, float]:
    """Return local q=-1,0,+1 intensity fractions in a spherical basis."""
    qaxis = unit(quantization_axis)
    e1, e2 = transverse_basis(qaxis)
    basis = {
        +1: -(e1 + 1j * e2) / np.sqrt(2.0),
        0: qaxis.astype(complex),
        -1: (e1 - 1j * e2) / np.sqrt(2.0),
    }
    epsilon = np.asarray(polarization, dtype=complex)
    epsilon /= np.linalg.norm(epsilon)
    fractions = {q: float(abs(np.vdot(component, epsilon)) ** 2) for q, component in basis.items()}
    total = sum(fractions.values())
    return {q: value / total for q, value in fractions.items()}
