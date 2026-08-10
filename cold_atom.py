"""Compact numerical models for idealized cold-atom systems."""

from __future__ import annotations

from math import exp, pi, sqrt
from typing import Iterable

import numpy as np

BOLTZMANN = 1.380649e-23  # J/K


def thermal_velocity_sigma(temperature: float, mass: float) -> float:
    """Return 1D thermal velocity standard deviation sqrt(k_B T / m)."""

    if temperature < 0:
        raise ValueError("temperature must be non-negative")
    if mass <= 0:
        raise ValueError("mass must be positive")
    return sqrt(BOLTZMANN * temperature / mass)


def ballistic_trajectory(
    time: Iterable[float],
    *,
    initial_position: float = 0.0,
    initial_velocity: float = 0.0,
    acceleration: float = -9.80665,
) -> np.ndarray:
    """Return x(t) for constant acceleration."""

    t = np.asarray(list(time), dtype=float)
    return initial_position + initial_velocity * t + 0.5 * acceleration * t**2


def rayleigh_range(waist: float, wavelength: float) -> float:
    """Return Gaussian-beam Rayleigh range z_R = pi w0^2 / lambda."""

    if waist <= 0 or wavelength <= 0:
        raise ValueError("waist and wavelength must be positive")
    return pi * waist**2 / wavelength


def gaussian_beam_waist(z: float, waist: float, wavelength: float) -> float:
    """Return Gaussian beam radius w(z)."""

    z_r = rayleigh_range(waist, wavelength)
    return waist * sqrt(1.0 + (z / z_r) ** 2)


def optical_dipole_potential(
    radius: float,
    z: float,
    *,
    trap_depth: float,
    waist: float,
    wavelength: float,
) -> float:
    """Ideal attractive Gaussian dipole potential in joules.

    ``trap_depth`` is supplied as a positive magnitude; the returned potential
    is negative, with -trap_depth at the beam focus.
    """

    if trap_depth <= 0:
        raise ValueError("trap_depth must be positive")
    w_z = gaussian_beam_waist(z, waist, wavelength)
    intensity_factor = (waist / w_z) ** 2 * exp(-2.0 * radius**2 / w_z**2)
    return -trap_depth * intensity_factor


def radial_trap_frequency(trap_depth: float, mass: float, waist: float) -> float:
    """Return radial angular trap frequency near a Gaussian trap minimum."""

    if trap_depth <= 0 or mass <= 0 or waist <= 0:
        raise ValueError("trap_depth, mass, and waist must be positive")
    return sqrt(4.0 * trap_depth / (mass * waist**2))


def axial_trap_frequency(
    trap_depth: float,
    mass: float,
    waist: float,
    wavelength: float,
) -> float:
    """Return axial angular trap frequency near the focus."""

    if trap_depth <= 0 or mass <= 0:
        raise ValueError("trap_depth and mass must be positive")
    z_r = rayleigh_range(waist, wavelength)
    return sqrt(2.0 * trap_depth / (mass * z_r**2))


def gravitational_sag(gravity: float, trap_angular_frequency: float) -> float:
    """Return equilibrium displacement g/omega^2 in a harmonic trap."""

    if trap_angular_frequency <= 0:
        raise ValueError("trap_angular_frequency must be positive")
    return gravity / trap_angular_frequency**2
