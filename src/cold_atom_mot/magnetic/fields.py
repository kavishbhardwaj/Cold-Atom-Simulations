"""Ideal quadrupole and composable residual magnetic fields, in tesla."""

from dataclasses import dataclass, field
import numpy as np


@dataclass(frozen=True)
class IdealQuadrupole:
    """Maxwell-consistent linear field B=G r with trace(G)=0."""
    radial_gradient: float
    centre: np.ndarray = field(default_factory=lambda: np.zeros(3))
    rotation: np.ndarray = field(default_factory=lambda: np.eye(3))

    def __post_init__(self) -> None:
        if self.radial_gradient <= 0:
            raise ValueError("radial_gradient must be positive")

    @property
    def gradient(self) -> np.ndarray:
        local = np.diag([self.radial_gradient, self.radial_gradient, -2 * self.radial_gradient])
        rotation = np.asarray(self.rotation, dtype=float)
        return rotation @ local @ rotation.T

    def field(self, positions: np.ndarray, time: float = 0.0) -> np.ndarray:
        del time
        return (np.asarray(positions) - np.asarray(self.centre)) @ self.gradient.T

    def jacobian(self, position: np.ndarray | None = None) -> np.ndarray:
        del position
        return self.gradient.copy()


@dataclass(frozen=True)
class ResidualField:
    uniform: np.ndarray = field(default_factory=lambda: np.zeros(3))
    gradient: np.ndarray = field(default_factory=lambda: np.zeros((3, 3)))
    ac_amplitude: np.ndarray = field(default_factory=lambda: np.zeros(3))
    ac_frequency: float = 0.0
    ac_phase: float = 0.0

    def field(self, positions: np.ndarray, time: float = 0.0) -> np.ndarray:
        points = np.asarray(positions, dtype=float)
        oscillation = self.ac_amplitude * np.sin(2 * np.pi * self.ac_frequency * time + self.ac_phase)
        return points @ np.asarray(self.gradient).T + np.asarray(self.uniform) + oscillation


@dataclass(frozen=True)
class CompositeField:
    components: tuple

    def field(self, positions: np.ndarray, time: float = 0.0) -> np.ndarray:
        return sum((component.field(positions, time) for component in self.components), start=np.zeros_like(np.asarray(positions, dtype=float)))

    def jacobian(self, position: np.ndarray, step: float = 1e-5) -> np.ndarray:
        point = np.asarray(position, dtype=float)
        result = np.empty((3, 3))
        for axis in range(3):
            offset = np.zeros(3)
            offset[axis] = step
            result[:, axis] = (self.field(point + offset) - self.field(point - offset)) / (2 * step)
        return result
