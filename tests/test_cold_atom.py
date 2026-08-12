import math

import numpy as np

from cold_atom_mot.foundations import (
    ballistic_trajectory,
    gaussian_beam_waist,
    gravitational_sag,
    optical_dipole_potential,
    radial_trap_frequency,
    rayleigh_range,
    thermal_velocity_sigma,
)


def test_thermal_velocity_zero_temperature():
    assert thermal_velocity_sigma(0.0, 1.0) == 0.0


def test_ballistic_trajectory():
    t = np.array([0.0, 1.0, 2.0])
    x = ballistic_trajectory(t, initial_position=1.0, initial_velocity=2.0, acceleration=-10.0)
    np.testing.assert_allclose(x, np.array([1.0, -2.0, -15.0]))


def test_gaussian_waist_at_focus():
    waist = 50e-6
    wavelength = 1064e-9
    assert math.isclose(gaussian_beam_waist(0.0, waist, wavelength), waist)
    assert rayleigh_range(waist, wavelength) > 0


def test_dipole_potential_at_center_equals_depth():
    depth = 1e-28
    u = optical_dipole_potential(
        0.0,
        0.0,
        trap_depth=depth,
        waist=50e-6,
        wavelength=1064e-9,
    )
    assert math.isclose(u, -depth)


def test_trap_frequency_and_sag_are_positive():
    omega = radial_trap_frequency(1e-28, 1.44e-25, 50e-6)
    assert omega > 0
    assert gravitational_sag(9.81, omega) > 0
