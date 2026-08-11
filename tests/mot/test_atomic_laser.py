import numpy as np
import pytest
from scipy.integrate import quad
from cold_atom_mot.atomic.rb87 import Rb87D2
from cold_atom_mot.laser.beam import GaussianBeam, six_beam_mot
from cold_atom_mot.laser.polarization import circular_polarization, spherical_fractions


def test_atomic_strengths_and_recoil_are_physical():
    atom = Rb87D2(); atom.validate()
    assert atom.recoil_velocity > 0
    assert 0 < atom.recoil_temperature < atom.doppler_temperature


def test_gaussian_power_integral():
    beam = GaussianBeam(np.array([1, 0, 0]), np.zeros(3), 0.012, 0.008, -1.0, 780e-9, 1)
    integral = quad(lambda radius: beam.peak_intensity * np.exp(-2 * radius**2 / beam.waist**2) * 2 * np.pi * radius, 0, np.inf)[0]
    assert integral == pytest.approx(beam.power, rel=1e-8)


def test_propagation_relative_helicity_and_local_decomposition():
    plus = circular_polarization(np.array([0, 0, 1]), 1)
    minus_direction = circular_polarization(np.array([0, 0, -1]), 1)
    assert np.allclose(np.vdot(plus, minus_direction), 0.0)
    fractions = spherical_fractions(plus, np.array([0, 0, 1]))
    assert sum(fractions.values()) == pytest.approx(1.0)
    assert max(fractions.values()) == pytest.approx(1.0)


def test_six_beams_are_independent_and_opposite_momenta():
    beams = six_beam_mot(0.01, 0.008, -1.0, 780e-9)
    assert len(beams) == 6 and len({beam.label for beam in beams}) == 6
    assert np.allclose(sum((beam.k_vector for beam in beams), start=np.zeros(3)), 0)


def test_counterpropagating_pair_has_opposite_local_sigma_components():
    beams = six_beam_mot(0.01, 0.008, -1.0, 780e-9)
    minus_x = spherical_fractions(beams[0].polarization, np.array([1, 0, 0]))
    plus_x = spherical_fractions(beams[1].polarization, np.array([1, 0, 0]))
    assert minus_x[-1] == pytest.approx(1.0)
    assert plus_x[+1] == pytest.approx(1.0)
