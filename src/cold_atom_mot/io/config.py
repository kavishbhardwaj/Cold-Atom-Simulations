"""Validated YAML configuration and model construction."""

from pathlib import Path
import numpy as np
import yaml
from ..atomic.rb87 import Rb87D2
from ..laser.beam import six_beam_mot
from ..magnetic.fields import CompositeField, IdealQuadrupole, ResidualField
from ..physics.force import EffectiveMOTForce
from ..physics.rate_equation import BeamFamily, MultilevelRateEquationMOT
from ..atomic.levels import build_rb87_d2_basis


def load_config(path: str | Path, *, validate: bool = True) -> dict:
    with Path(path).open(encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    if validate:
        validate_config(config)
    return config


def validate_config(config: dict) -> None:
    """Reject missing or manifestly unphysical inputs for each fidelity level."""
    if config.get("model") == "level_c_reduced_two_level_obe":
        required = ("atom", "obe", "output")
        if any(section not in config for section in required):
            raise ValueError(f"Level-C configuration requires sections: {required}")
        obe = config["obe"]
        if obe["saturation"] < 0 or obe["duration_lifetimes"] <= 0:
            raise ValueError("OBE saturation must be non-negative and duration positive")
        if obe["rtol"] <= 0 or obe["atol"] <= 0 or obe["max_step_lifetimes"] <= 0:
            raise ValueError("OBE tolerances and maximum step must be positive")
        return
    required = ("laser", "magnetic_field", "gravity", "simulation", "monte_carlo")
    if any(section not in config for section in required):
        raise ValueError(f"configuration requires sections: {required}")
    laser = config["laser"]
    if laser["power_per_beam_w"] < 0 or laser["waist_m"] <= 0:
        raise ValueError("laser power must be non-negative and waist positive")
    if config.get("model") == "level_b_multilevel_rate_equation":
        repump = config.get("repump")
        if not repump or repump["power_per_beam_w"] < 0 or repump["waist_m"] <= 0:
            raise ValueError("Level-B repump power must be non-negative and waist positive")
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


def build_multilevel_model(config: dict) -> MultilevelRateEquationMOT:
    """Build the Level-B cooling+repump population model."""
    atom = Rb87D2()
    magnetic = config["magnetic_field"]
    field = CompositeField((
        IdealQuadrupole(magnetic["radial_gradient_t_per_m"]),
        ResidualField(
            uniform=np.asarray(magnetic.get("uniform_stray_t", [0, 0, 0]), dtype=float),
            gradient=np.asarray(magnetic.get("stray_gradient_t_per_m", np.zeros((3, 3))), dtype=float),
        ),
    ))
    families = []
    for section, ground_f, target_f in (("laser", 2, 3), ("repump", 1, 2)):
        parameters = config[section]
        beams = six_beam_mot(
            parameters["power_per_beam_w"],
            parameters["waist_m"],
            parameters["detuning_gamma"] * atom.gamma,
            atom.wavelength,
        )
        families.extend(BeamFamily(beam, ground_f, target_f, section) for beam in beams)
    return MultilevelRateEquationMOT(
        atom,
        build_rb87_d2_basis(atom),
        families,
        field,
        np.asarray(config["gravity"]["vector_m_per_s2"], dtype=float),
    )
