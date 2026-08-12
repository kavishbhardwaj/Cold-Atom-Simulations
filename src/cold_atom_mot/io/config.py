"""Validated YAML configuration and model construction."""

from pathlib import Path
import numpy as np
import yaml
from ..atomic.species import get_atomic_line, build_atomic_basis
from ..laser.beam import six_beam_mot
from ..magnetic.fields import CompositeField, IdealQuadrupole, ResidualField
from ..physics.force import EffectiveMOTForce
from ..physics.rate_equation import BeamFamily, MultilevelRateEquationMOT
from ..physics.subdoppler import coherent_six_beam_field, PolarizationGradientModel
from ..vacuum import VaporState, rubidium_vapor_pressure_pa


def load_config(path: str | Path, *, validate: bool = True) -> dict:
    with Path(path).open(encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    if validate:
        validate_config(config)
    return config


def validate_config(config: dict) -> None:
    """Reject missing or manifestly unphysical inputs for each fidelity level."""
    if config.get("model") == "configurable_vapor_loading":
        for section in ("mot_config", "vacuum", "capture", "loading", "output"):
            if section not in config:
                raise ValueError(f"vapor-loading configuration requires {section}")
        vacuum, capture, loading = config["vacuum"], config["capture"], config["loading"]
        temperatures=(vacuum["rb_reservoir_temperature_k"],vacuum["vapor_temperature_k"],
                      vacuum["background_temperature_k"])
        if min(temperatures) <= 0 or vacuum["background_gas_pressure_pa"] < 0:
            raise ValueError("all temperatures must be positive and background pressure non-negative")
        if vacuum.get("rb_partial_pressure_pa") is not None and vacuum["rb_partial_pressure_pa"] < 0:
            raise ValueError("Rb partial pressure must be non-negative")
        if not np.isclose(sum(vacuum["isotope_fractions"].values()), 1.0):
            raise ValueError("isotope fractions must sum to one")
        edges = np.asarray(capture["speed_bin_edges_m_s"])
        criterion = capture["criterion"]
        if (capture["surface_radius_m"] <= 0 or capture["atoms_per_bin"] <= 0 or
                capture["max_step_s"] <= 0 or edges[0] != 0 or np.any(np.diff(edges) <= 0)):
            raise ValueError("capture surface, sampling, step, and speed edges must be physical")
        if min(criterion.values()) <= 0:
            raise ValueError("capture criterion values must be positive")
        if (not 0 < capture["tail_relative_loading_tolerance"] < 1 or
                capture["maximum_speed_m_s"] <= edges[-1]):
            raise ValueError("tail tolerance must be in (0,1) and maximum speed above initial edges")
        if min(loading["background_one_body_loss_s"], loading["hot_rb_one_body_loss_s"],
               loading["two_body_loss_m3_s"]) < 0:
            raise ValueError("loss coefficients must be non-negative")
        collision = loading.get("background_collision_model")
        components = loading.get("background_gas_components")
        if collision is not None and components:
            raise ValueError("background collision model and component list would double count")
        if (collision is not None or components) and loading["background_one_body_loss_s"] > 0:
            raise ValueError(
                "background_one_body_loss_s must be zero when kinetic background-gas inputs are present"
            )
        if collision is not None and min(collision["particle_mass_kg"],
                                         collision["effective_loss_cross_section_m2"]) <= 0:
            raise ValueError("kinetic background collision inputs must be positive")
        for component in components or []:
            if min(component["partial_pressure_pa"],component["particle_mass_kg"],
                   component["effective_loss_cross_section_m2"]) <= 0:
                raise ValueError("background component inputs must be positive")
        return
    if config.get("model") == "polarization_gradient":
        for section in ("atom", "laser", "magnetic_field", "simulation", "output"):
            if section not in config:
                raise ValueError(f"polarization-gradient configuration requires {section}")
        laser, simulation = config["laser"], config["simulation"]
        if laser["saturation_per_beam"] < 0 or laser["detuning_gamma"] == 0:
            raise ValueError("PGC saturation must be non-negative and detuning non-zero")
        if len(laser["phases_rad"]) != 6 or simulation["periods"] <= simulation["discard_periods"]:
            raise ValueError("PGC requires six phases and periods greater than discarded periods")
        if simulation["steps_per_period"] < 8 or simulation["velocity_m_per_s"] == 0:
            raise ValueError("PGC resolution must be >=8 and probe velocity non-zero")
        return
    if config.get("model") == "two_level_obe":
        required = ("atom", "obe", "output")
        if any(section not in config for section in required):
            raise ValueError(f"coherent-model configuration requires sections: {required}")
        obe = config["obe"]
        if obe["saturation"] < 0 or obe["duration_lifetimes"] <= 0:
            raise ValueError("OBE saturation must be non-negative and duration positive")
        if obe["rtol"] <= 0 or obe["atol"] <= 0 or obe["max_step_lifetimes"] <= 0:
            raise ValueError("OBE tolerances and maximum step must be positive")
        if obe.get("pure_dephasing_gamma", 0.0) < 0:
            raise ValueError("OBE pure dephasing must be non-negative")
        return
    required = ("laser", "magnetic_field", "gravity", "simulation", "monte_carlo")
    if any(section not in config for section in required):
        raise ValueError(f"configuration requires sections: {required}")
    laser = config["laser"]
    if laser["power_per_beam_w"] < 0 or laser["waist_m"] <= 0:
        raise ValueError("laser power must be non-negative and waist positive")
    if config.get("model") == "multilevel_rate_equation":
        repump = config.get("repump")
        if not repump or repump["power_per_beam_w"] < 0 or repump["waist_m"] <= 0:
            raise ValueError("rate-equation repump power must be non-negative and waist positive")
    if config["simulation"]["duration_s"] <= 0 or config["monte_carlo"]["time_step_s"] <= 0:
        raise ValueError("simulation times must be positive")


def build_effective_model(config: dict) -> EffectiveMOTForce:
    atom_config = config.get("atom", {"isotope": "87Rb", "line": "D2"})
    atom = get_atomic_line(atom_config.get("isotope", atom_config.get("species", "87Rb")), atom_config.get("line", "D2"))
    laser = config["laser"]
    beams = six_beam_mot(
        laser["power_per_beam_w"], laser["waist_m"],
        laser["detuning_gamma"] * atom.gamma_rad_s, atom.wavelength_m,
    )
    magnetic = config["magnetic_field"]
    quadrupole = IdealQuadrupole(magnetic["radial_gradient_t_per_m"])
    residual = ResidualField(
        uniform=np.asarray(magnetic.get("uniform_stray_t", [0, 0, 0]), dtype=float),
        gradient=np.asarray(magnetic.get("stray_gradient_t_per_m", np.zeros((3, 3))), dtype=float),
    )
    return EffectiveMOTForce(atom, beams, CompositeField((quadrupole, residual)), np.asarray(config["gravity"]["vector_m_per_s2"], dtype=float))


def build_multilevel_model(config: dict) -> MultilevelRateEquationMOT:
    """Build the rate-equation cooling+repump population model."""
    atom_config = config["atom"]
    atom = get_atomic_line(atom_config.get("isotope", "87Rb"), atom_config.get("line", "D2"))
    if not atom.rate_equation_mot or atom.cooling_transition is None:
        raise ValueError(f"rate-equation MOT is not supported for {atom.isotope} {atom.line}")
    for section, transition in (("laser", atom.cooling_transition), ("repump", atom.repump_transition)):
        expected=f"f{transition[0]}_to_fprime{transition[1]}"
        configured=config[section].get("target_transition",expected).lower()
        if configured != expected:
            raise ValueError(f"{atom.isotope} {atom.line} {section} requires {expected}, not {configured}")
    magnetic = config["magnetic_field"]
    field = CompositeField((
        IdealQuadrupole(magnetic["radial_gradient_t_per_m"]),
        ResidualField(
            uniform=np.asarray(magnetic.get("uniform_stray_t", [0, 0, 0]), dtype=float),
            gradient=np.asarray(magnetic.get("stray_gradient_t_per_m", np.zeros((3, 3))), dtype=float),
        ),
    ))
    families = []
    for section, transition in (("laser", atom.cooling_transition), ("repump", atom.repump_transition)):
        ground_f, target_f = transition
        parameters = config[section]
        beams = six_beam_mot(
            parameters["power_per_beam_w"],
            parameters["waist_m"],
            parameters["detuning_gamma"] * atom.gamma_rad_s,
            atom.wavelength_m,
        )
        families.extend(BeamFamily(beam, ground_f, target_f, section) for beam in beams)
    return MultilevelRateEquationMOT(
        atom,
        build_atomic_basis(atom.isotope, atom.line),
        families,
        field,
        np.asarray(config["gravity"]["vector_m_per_s2"], dtype=float),
    )


def build_subdoppler_model(config: dict) -> PolarizationGradientModel:
    """Build the explicitly phase-coherent polarization-gradient F=2 -> F'=3 model."""
    atom_config=config["atom"]; atom=get_atomic_line(atom_config["isotope"],atom_config["line"]); laser = config["laser"]
    groups = laser.get("coherence_groups", ["all"] * 6)
    beams = coherent_six_beam_field(atom.wave_number_rad_m, laser["saturation_per_beam"], laser["phases_rad"], groups)
    ground_f, excited_f = atom.cooling_transition
    return PolarizationGradientModel(build_atomic_basis(atom.isotope, atom.line), ground_f, excited_f,
                                     laser["detuning_gamma"] * atom.gamma_rad_s, beams,
                                     magnetic_field_t=config["magnetic_field"]["uniform_t"])


def build_vapor_state(config: dict) -> VaporState:
    """Build thermodynamic Rb state without conflating background pressure."""
    vacuum = config["vacuum"]
    reservoir_temperature = vacuum["rb_reservoir_temperature_k"]
    rb_pressure = vacuum.get("rb_partial_pressure_pa")
    if rb_pressure is None:
        rb_pressure = rubidium_vapor_pressure_pa(
            reservoir_temperature,
            allow_extrapolation=vacuum.get("allow_vapor_pressure_extrapolation", False),
        )
    return VaporState(
        reservoir_temperature,
        vacuum["vapor_temperature_k"],
        vacuum["background_temperature_k"],
        rb_pressure,
        vacuum["background_gas_pressure_pa"],
        vacuum["isotope_fractions"],
    )
