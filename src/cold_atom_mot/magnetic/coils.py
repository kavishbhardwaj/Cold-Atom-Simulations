"""Segmented-wire Biot-Savart circular coils with arbitrary pose.

A loop is represented by straight current elements evaluated at their
midpoints.  This converges as O(N^-2) away from the wire; convergence is tested
and documented.  Points on a conductor are outside the model domain.
"""

from dataclasses import dataclass
import numpy as np
from scipy.constants import mu_0
from scipy.optimize import least_squares
from ..laser.polarization import transverse_basis, unit


@dataclass(frozen=True)
class CircularCoil:
    centre: np.ndarray
    normal: np.ndarray
    radius: float
    current: float
    turns: int = 1
    segments: int = 256

    def __post_init__(self) -> None:
        if self.radius <= 0 or self.turns <= 0 or self.segments < 16:
            raise ValueError("radius/turns must be positive and segments >= 16")
        object.__setattr__(self, "centre", np.asarray(self.centre, dtype=float))
        object.__setattr__(self, "normal", unit(self.normal))

    def wire_elements(self) -> tuple[np.ndarray, np.ndarray]:
        e1, e2 = transverse_basis(self.normal)
        edges = np.linspace(0.0, 2 * np.pi, self.segments + 1)
        wire = self.centre + self.radius * (np.cos(edges)[:, None] * e1 + np.sin(edges)[:, None] * e2)
        dl = np.diff(wire, axis=0)
        midpoint = 0.5 * (wire[:-1] + wire[1:])
        return midpoint, dl

    def field(self, positions: np.ndarray, time: float = 0.0) -> np.ndarray:
        del time
        points = np.asarray(positions, dtype=float)
        original_shape = points.shape
        flat = points.reshape(-1, 3)
        midpoint, dl = self.wire_elements()
        displacement = flat[:, None, :] - midpoint[None, :, :]
        distance = np.linalg.norm(displacement, axis=-1)
        if np.any(distance == 0):
            raise ValueError("Biot-Savart field is singular on the wire")
        integrand = np.cross(dl[None, :, :], displacement) / distance[..., None] ** 3
        result = mu_0 / (4 * np.pi) * self.current * self.turns * np.sum(integrand, axis=1)
        return result.reshape(original_shape)


@dataclass(frozen=True)
class AntiHelmholtzPair:
    first: CircularCoil
    second: CircularCoil

    @classmethod
    def symmetric(cls, radius: float, separation: float, current: float, turns: int, segments: int = 256, tilt_y: float = 0.0, lateral_offset: float = 0.0, current_imbalance: float = 0.0):
        """Build coaxial opposed-current coils; tilt_y applies to the second coil."""
        if separation <= 0 or abs(current_imbalance) >= 1:
            raise ValueError("separation must be positive and |imbalance| < 1")
        normal1 = np.array([0.0, 0.0, 1.0])
        normal2 = np.array([np.sin(tilt_y), 0.0, np.cos(tilt_y)])
        first = CircularCoil(np.array([0.0, 0.0, -separation / 2]), normal1, radius, current * (1 + current_imbalance), turns, segments)
        second = CircularCoil(np.array([lateral_offset, 0.0, separation / 2]), normal2, radius, -current * (1 - current_imbalance), turns, segments)
        return cls(first, second)

    def field(self, positions: np.ndarray, time: float = 0.0) -> np.ndarray:
        return self.first.field(positions, time) + self.second.field(positions, time)

    def jacobian(self, position: np.ndarray, step: float = 1e-5) -> np.ndarray:
        point = np.asarray(position, dtype=float)
        jac = np.empty((3, 3))
        for axis in range(3):
            offset = np.zeros(3)
            offset[axis] = step
            jac[:, axis] = (self.field(point + offset) - self.field(point - offset)) / (2 * step)
        return jac

    def field_zero(self, guess: np.ndarray | None = None) -> np.ndarray:
        initial = np.zeros(3) if guess is None else np.asarray(guess, dtype=float)
        result = least_squares(lambda point: self.field(point), initial, xtol=1e-13, ftol=1e-13, gtol=1e-13)
        return result.x
