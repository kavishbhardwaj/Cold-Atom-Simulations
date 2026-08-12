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

    @property
    def is_time_independent(self) -> bool:
        return True


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

    @property
    def is_time_independent(self) -> bool:
        return self.first.is_time_independent and self.second.is_time_independent


@dataclass(frozen=True)
class HelmholtzPair:
    """Two approximately co-directed circular bias coils with independent geometry."""
    first: CircularCoil
    second: CircularCoil

    @classmethod
    def imperfect(cls, axis, radius, separation, current, turns, *, segments=256,
                   separation_error=0.0, radius_mismatch=0.0,
                   current_imbalance=0.0, turns_mismatch=0,
                   lateral_displacement=(0, 0, 0), tilt=(0, 0, 0)):
        axis = unit(axis); separation += separation_error
        lateral = np.asarray(lateral_displacement, float)
        tilt = np.asarray(tilt, float)
        second_axis = unit(axis+tilt)
        return cls(
            CircularCoil(-axis*separation/2, axis, radius*(1-radius_mismatch/2),
                         current*(1-current_imbalance/2), turns, segments),
            CircularCoil(axis*separation/2+lateral, second_axis,
                         radius*(1+radius_mismatch/2), current*(1+current_imbalance/2),
                         turns+turns_mismatch, segments))

    def field(self, positions, time=0.0):
        return self.first.field(positions, time)+self.second.field(positions, time)

    def jacobian(self, position=np.zeros(3), step=1e-5):
        point=np.asarray(position,float); output=np.empty((3,3))
        for axis in range(3):
            shift=np.zeros(3);shift[axis]=step
            output[:,axis]=(self.field(point+shift)-self.field(point-shift))/(2*step)
        return output

    def curvature(self, position=np.zeros(3), step=1e-4):
        """Return d2 B_i / dx_j dx_k as a 3x3x3 tensor."""
        point=np.asarray(position,float); output=np.empty((3,3,3))
        for axis in range(3):
            shift=np.zeros(3);shift[axis]=step
            output[:,:,axis]=(self.jacobian(point+shift,step/5)-self.jacobian(point-shift,step/5))/(2*step)
        return output

    @property
    def is_time_independent(self): return True


@dataclass(frozen=True)
class ThreeAxisBiasCoils:
    """Physical x/y/z pairs plus an empirical B=M I+B_offset calibration."""
    pairs: tuple
    calibration_matrix: np.ndarray
    offset: np.ndarray

    def __post_init__(self):
        if len(self.pairs)!=3: raise ValueError("three bias-coil pairs are required")
        matrix=np.asarray(self.calibration_matrix,float); offset=np.asarray(self.offset,float)
        if matrix.shape!=(3,3) or offset.shape!=(3,): raise ValueError("calibration must be 3x3 with a three-vector offset")
        object.__setattr__(self,"calibration_matrix",matrix);object.__setattr__(self,"offset",offset)

    def compensation_currents(self, background, target=np.zeros(3)):
        """Least-squares currents minimizing M I+B_offset+B_background-target."""
        return np.linalg.lstsq(self.calibration_matrix,
            np.asarray(target)-self.offset-np.asarray(background),rcond=None)[0]

    def calibrated_field(self, currents, background=np.zeros(3)):
        return self.calibration_matrix@np.asarray(currents)+self.offset+np.asarray(background)

    def physical_field(self, positions, currents):
        currents=np.asarray(currents,float)
        total=np.zeros_like(np.asarray(positions,float))
        for pair,current in zip(self.pairs,currents):
            reference=pair.first.current
            total += pair.field(positions)*current/reference
        return total
