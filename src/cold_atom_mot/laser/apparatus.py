"""Configurable physical six-beam apparatus and retroreflection helpers."""
from dataclasses import dataclass, replace
import numpy as np
from .beam import GaussianBeam
from .polarization import propagate_jones


@dataclass(frozen=True)
class Retroreflection:
    power_transmission: float = 1.0
    phase_offset: float = 0.0
    mirror_tilt: tuple = (0.0, 0.0, 0.0)
    double_pass_elements: tuple = ()
    coherence_group: str | None = None

    def reflected(self, incident: GaussianBeam, label=None):
        if not 0 <= self.power_transmission <= 1:
            raise ValueError("retroreflection transmission must be in [0,1]")
        direction = -incident.direction+np.asarray(self.mirror_tilt, float)
        direction /= np.linalg.norm(direction)
        jones = incident.jones_vector
        if jones is not None:
            # Reflection reverses propagation handedness; the same transverse
            # lab field is represented by conjugating the second basis axis.
            jones = propagate_jones([jones[0], -jones[1]], self.double_pass_elements)
        return replace(incident, direction=direction,
            power=incident.power*self.power_transmission,
            phase=incident.phase+self.phase_offset, jones_vector=jones,
            optical_elements=(), coherence_group=self.coherence_group,
            label=label or f"{incident.label}-retro")


@dataclass(frozen=True)
class SixBeamApparatus:
    """Exactly six independently specified physical beams."""
    beams: tuple[GaussianBeam, ...]
    topology: str = "six_independent"

    def __post_init__(self):
        if len(self.beams) != 6:
            raise ValueError("a six-beam apparatus requires exactly six beams")
        if self.topology not in ("six_independent", "three_retroreflected", "three_pairs"):
            raise ValueError("unsupported apparatus topology")

    def local_polarizations(self, quantization_axis):
        from .polarization import spherical_fractions
        return tuple(spherical_fractions(beam.polarization, quantization_axis) for beam in self.beams)

    @classmethod
    def three_retroreflected(cls, inputs, retroreflectors):
        if len(inputs) != 3 or len(retroreflectors) != 3:
            raise ValueError("three inputs and retroreflectors are required")
        beams=[]
        for beam, retro in zip(inputs, retroreflectors):
            beams.extend((beam,retro.reflected(beam)))
        return cls(tuple(beams),"three_retroreflected")
