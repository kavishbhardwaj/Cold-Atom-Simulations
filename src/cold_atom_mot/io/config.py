"""Validated YAML configuration and model construction."""

from pathlib import Path
import numpy as np
import yaml
from ..atomic.rb87 import Rb87D2
from ..laser.beam import six_beam_mot
from ..magnetic.fields import CompositeField, IdealQuadrupole, ResidualField
from ..physics.force import EffectiveMOTForce


def load_config(path: str | Path) -> dict:
    with Path(path).open(encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    validate_config(config)
    return config


def validate_config(config: dict) -> None:
    """Reject missing or manifestly unphysical Phase-1 inputs."""
    required = ("laser", "magnetic_field", "gravity", "simulation", "monte_carlo")
    if any(section not in config for section in required):
        raise ValueError(f"configuration requires sections: {required}")
    laser = config["laser"]
    if laser["power_per_beam_w"] < 0 or laser["waist_m"] <= 0:
        raise ValueError("laser power must be non-negative and waist positive")
    if config["simulation"]["duration_s"] <= 0 or config["monte_carlo"]["time_step_s"] <= 0:
        raise ValueError("simulation times must be positive")


def build_effective_model(config: dict) -> EffectiveMOTForce:
    atom = Rb87D2()
    laser = config["laser"]
    beams = six_beam_mot(
        laser["power_per_beam_w"], laser["waist_m"],
        laser["detuning_gamma"] * atom.gamma, atom.wavelength,
    )
    magnetic = config["magnetic_field"]
    quadrupole = IdealQuadrupole(magnetic["radial_gradient_t_per_m"])
    residual = ResidualField(
        uniform=np.asarray(magnetic.get("uniform_stray_t", [0, 0, 0]), dtype=float),
        gradient=np.asarray(magnetic.get("stray_gradient_t_per_m", np.zeros((3, 3))), dtype=float),
    )
    return EffectiveMOTForce(atom, beams, CompositeField((quadrupole, residual)), np.asarray(config["gravity"]["vector_m_per_s2"], dtype=float))
