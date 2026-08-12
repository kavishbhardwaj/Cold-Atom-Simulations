"""Time-sequenced controls for cold-atom experimental cycles.

The sequence engine is deliberately solver-agnostic: force and internal-state
solvers consume immutable snapshots at the requested laboratory time. All
times are seconds, fields tesla, gradients T/m, powers watts and detunings rad/s.
"""
from dataclasses import dataclass, replace, field
from functools import cached_property
import numpy as np


@dataclass(frozen=True)
class Ramp:
    """Interpolate between endpoints using a step, linear, or smoothstep law."""
    start: object
    stop: object | None = None
    kind: str = "step"

    def __post_init__(self):
        if self.kind not in ("step", "linear", "smooth"):
            raise ValueError("ramp kind must be step, linear, or smooth")
        if self.stop is None:
            object.__setattr__(self, "stop", self.start)

    def value(self, fraction):
        fraction = float(np.clip(fraction, 0, 1))
        if fraction == 0:
            return np.asarray(self.start).copy() if np.asarray(self.start).ndim else self.start
        if fraction == 1:
            return np.asarray(self.stop).copy() if np.asarray(self.stop).ndim else self.stop
        if self.kind == "step":
            fraction = 0.0 if fraction < 1 else 1.0
        elif self.kind == "smooth":
            fraction = fraction*fraction*(3-2*fraction)
        start, stop = np.asarray(self.start), np.asarray(self.stop)
        value = start+(stop-start)*fraction
        return float(value) if value.ndim == 0 else value


@dataclass(frozen=True)
class Stage:
    name: str
    duration: float
    cooling_detuning: Ramp
    cooling_power: Ramp
    repump_power: Ramp
    quadrupole_gradient: Ramp
    bias_field: Ramp = field(default_factory=lambda: Ramp(np.zeros(3)))
    cooling_frequency_offset: Ramp = field(default_factory=lambda: Ramp(0.0))
    repump_frequency_offset: Ramp = field(default_factory=lambda: Ramp(0.0))
    cooling_coherence_groups: tuple | None = None
    repump_coherence_groups: tuple | None = None
    cooling_phases: tuple | None = None
    polarization_purity: Ramp = field(default_factory=lambda: Ramp(1.0))

    def __post_init__(self):
        if self.duration < 0:
            raise ValueError("stage duration must be non-negative")
        for values in (self.cooling_coherence_groups, self.repump_coherence_groups,
                       self.cooling_phases):
            if values is not None and len(values) != 6:
                raise ValueError("beam controls require six entries")


@dataclass(frozen=True)
class ControlSnapshot:
    time: float
    stage: str
    fraction: float
    cooling_detuning: float
    cooling_power: float
    repump_power: float
    quadrupole_gradient: float
    bias_field: np.ndarray
    cooling_frequency_offset: float
    repump_frequency_offset: float
    polarization_purity: float
    cooling_coherence_groups: tuple | None
    repump_coherence_groups: tuple | None
    cooling_phases: tuple | None


@dataclass(frozen=True)
class EddyCurrentResponse:
    """Finite coil switch-off plus eddy, DC, and mains-frequency residuals."""
    switch_time: float
    initial_gradient: float
    coil_time_constant: float
    eddy_amplitude: np.ndarray
    eddy_time_constant: float
    residual_dc: np.ndarray
    ac_amplitude: np.ndarray = None
    ac_frequency: float = 0.0
    ac_phase: float = 0.0

    def __post_init__(self):
        if self.coil_time_constant <= 0 or self.eddy_time_constant <= 0:
            raise ValueError("decay time constants must be positive")
        for name in ("eddy_amplitude", "residual_dc"):
            object.__setattr__(self, name, np.asarray(getattr(self, name), float))
        amplitude = np.zeros(3) if self.ac_amplitude is None else self.ac_amplitude
        object.__setattr__(self, "ac_amplitude", np.asarray(amplitude, float))

    def gradient(self, time):
        elapsed = max(0.0, time-self.switch_time)
        return self.initial_gradient if time < self.switch_time else self.initial_gradient*np.exp(-elapsed/self.coil_time_constant)

    def field(self, position, time=0.0):
        elapsed = max(0.0, time-self.switch_time)
        eddy = self.eddy_amplitude*(1.0 if time < self.switch_time else np.exp(-elapsed/self.eddy_time_constant))
        ac = self.ac_amplitude*np.sin(2*np.pi*self.ac_frequency*time+self.ac_phase)
        gradient = np.diag([self.gradient(time), self.gradient(time), -2*self.gradient(time)])
        return gradient@np.asarray(position, float)+self.residual_dc+eddy+ac

    @property
    def is_time_independent(self):
        return False


@dataclass(frozen=True)
class ExperimentalSequence:
    stages: tuple[Stage, ...]
    magnetic_response: EddyCurrentResponse | None = None

    def __post_init__(self):
        if not self.stages:
            raise ValueError("sequence requires at least one stage")

    @cached_property
    def boundaries(self):
        return np.cumsum([0.0]+[stage.duration for stage in self.stages])

    @property
    def duration(self):
        return float(self.boundaries[-1])

    def locate(self, time):
        if time < 0 or time > self.duration:
            raise ValueError("time lies outside sequence")
        index = min(np.searchsorted(self.boundaries[1:], time, side="right"), len(self.stages)-1)
        stage = self.stages[index]; start = self.boundaries[index]
        fraction = 1.0 if stage.duration == 0 else (time-start)/stage.duration
        return stage, float(np.clip(fraction, 0, 1))

    def snapshot(self, time):
        stage, fraction = self.locate(time)
        gradient = stage.quadrupole_gradient.value(fraction)
        bias = stage.bias_field.value(fraction)
        if self.magnetic_response is not None and time >= self.magnetic_response.switch_time:
            gradient = self.magnetic_response.gradient(time)
            bias = bias+self.magnetic_response.field(np.zeros(3), time)
        return ControlSnapshot(time, stage.name, fraction,
            stage.cooling_detuning.value(fraction), stage.cooling_power.value(fraction),
            stage.repump_power.value(fraction), gradient, np.asarray(bias),
            stage.cooling_frequency_offset.value(fraction),
            stage.repump_frequency_offset.value(fraction),
            stage.polarization_purity.value(fraction), stage.cooling_coherence_groups,
            stage.repump_coherence_groups, stage.cooling_phases)

    def apply_beams(self, beams, time):
        """Return beam copies carrying the controls at laboratory ``time``."""
        snapshot = self.snapshot(time); output = []
        for index, beam in enumerate(beams):
            cooling = index < 6
            groups = snapshot.cooling_coherence_groups if cooling else snapshot.repump_coherence_groups
            phase = (snapshot.cooling_phases[index] if cooling and snapshot.cooling_phases else beam.phase)
            output.append(replace(beam,
                power=snapshot.cooling_power if cooling else snapshot.repump_power,
                detuning=snapshot.cooling_detuning if cooling else beam.detuning,
                frequency_offset=(snapshot.cooling_frequency_offset if cooling else snapshot.repump_frequency_offset),
                coherence_group=(groups[index % 6] if groups else beam.coherence_group),
                phase=phase, polarization_purity=snapshot.polarization_purity))
        return output

    def apply_beam_families(self, families, time):
        """Apply controls to cooling/repump BeamFamily objects for rate/OBE solvers."""
        beams = self.apply_beams([family.beam for family in families], time)
        return [replace(family, beam=beam) for family, beam in zip(families, beams)]


class SequencedForce:
    """Adapter allowing deterministic trajectories to consume a sequence."""
    def __init__(self, base_model, sequence):
        self.base_model, self.sequence, self.atom = base_model, sequence, base_model.atom

    def _model(self, time):
        model = replace(self.base_model)
        model.beams = self.sequence.apply_beams(self.base_model.beams, time)
        if self.sequence.magnetic_response is not None:
            model.magnetic_field = self.sequence.magnetic_response
        return model

    def force(self, position, velocity, time=0.0):
        return self._model(time).force(position, velocity, time)

    def scattering_rates(self, position, velocity, time=0.0):
        return self._model(time).scattering_rates(position, velocity, time)

    @property
    def beams(self):
        return self.base_model.beams

    @property
    def gravity(self):
        return self.base_model.gravity
