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

    @property
    def is_time_independent(self) -> bool:
        return True


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

    @property
    def is_time_independent(self) -> bool:
        return self.ac_frequency == 0.0 or not np.any(self.ac_amplitude)


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

    @property
    def is_time_independent(self) -> bool:
        return all(getattr(component, "is_time_independent", False)
                   for component in self.components)


@dataclass(frozen=True)
class HarmonicResidualField:
    """Uniform/gradient background plus configurable AC harmonics."""
    uniform: np.ndarray
    gradient: np.ndarray = field(default_factory=lambda: np.zeros((3,3)))
    # tuples (vector amplitude [T], frequency [Hz], phase [rad])
    harmonics: tuple = ()

    def field(self, positions, time=0.0):
        result=np.asarray(positions,float)@np.asarray(self.gradient,float).T+np.asarray(self.uniform,float)
        for amplitude,frequency,phase in self.harmonics:
            result=result+np.asarray(amplitude,float)*np.sin(2*np.pi*frequency*time+phase)
        return result

    @property
    def is_time_independent(self): return len(self.harmonics)==0


@dataclass(frozen=True)
class SwitchingTransientField:
    """L/R current and multiple eddy exponentials, or a measured waveform."""
    switch_time: float
    dc_field: np.ndarray
    coil_amplitude: np.ndarray
    coil_time_constant: float
    eddy_components: tuple = ()  # (amplitude vector, tau)
    waveform_time: np.ndarray | None = None
    waveform_field: np.ndarray | None = None

    def __post_init__(self):
        if self.coil_time_constant<=0 or any(tau<=0 for _,tau in self.eddy_components):
            raise ValueError("transient time constants must be positive")
        if (self.waveform_time is None)!=(self.waveform_field is None):
            raise ValueError("waveform time and field must be supplied together")

    def field(self, positions, time=0.0):
        del positions
        elapsed=max(0,time-self.switch_time); result=np.asarray(self.dc_field,float)
        if self.waveform_time is not None:
            return result+np.array([np.interp(elapsed,self.waveform_time,np.asarray(self.waveform_field)[:,i]) for i in range(3)])
        result=result+np.asarray(self.coil_amplitude)*np.exp(-elapsed/self.coil_time_constant)
        for amplitude,tau in self.eddy_components: result=result+np.asarray(amplitude)*np.exp(-elapsed/tau)
        return result

    @property
    def is_time_independent(self): return False
