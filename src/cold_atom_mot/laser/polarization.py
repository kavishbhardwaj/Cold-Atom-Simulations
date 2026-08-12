"""Propagation-relative polarization and local spherical-basis utilities."""

import numpy as np
from dataclasses import dataclass


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


def jones_rotation(angle):
    """Real rotation in a beam's transverse Jones basis."""
    return np.array([[np.cos(angle), -np.sin(angle)],
                     [np.sin(angle), np.cos(angle)]], complex)


@dataclass(frozen=True)
class JonesElement:
    """Lossless retarder, or an ideal linear polarizer, in radians."""
    kind: str
    angle: float = 0.0
    retardance: float | None = None
    retardance_error: float = 0.0

    def matrix(self):
        rotation = jones_rotation(self.angle)
        if self.kind == "polarizer":
            native = np.diag([1, 0]).astype(complex)
        else:
            nominal = {"quarter_wave": np.pi/2, "half_wave": np.pi,
                       "retarder": self.retardance}.get(self.kind)
            if nominal is None:
                raise ValueError("Jones element must be polarizer, quarter_wave, half_wave, or retarder")
            native = np.diag([np.exp(-.5j*(nominal+self.retardance_error)),
                              np.exp(.5j*(nominal+self.retardance_error))])
        return rotation @ native @ rotation.T


def propagate_jones(initial, elements=()):
    """Propagate and normalize a two-component Jones vector."""
    vector = np.asarray(initial, complex)
    if vector.shape != (2,):
        raise ValueError("Jones vector must have two components")
    for element in elements:
        vector = element.matrix() @ vector
    norm = np.linalg.norm(vector)
    if norm == 0:
        raise ValueError("optical train extinguishes the input")
    return vector/norm
