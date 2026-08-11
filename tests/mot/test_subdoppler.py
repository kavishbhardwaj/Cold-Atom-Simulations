import numpy as np
import pytest
from cold_atom_mot.atomic.rb87 import Rb87D2
from cold_atom_mot.physics.subdoppler import coherent_six_beam_field, PolarizationGradientModel


def model(**kwargs):
    atom = Rb87D2()
    beams = coherent_six_beam_field(atom.wave_number, 0.04, [0, 0, 0, np.pi/2, 0, np.pi/4])
    return PolarizationGradientModel(atom.gamma, -3 * atom.gamma, atom.wave_number, beams, **kwargs)


def test_phase_resolved_field_is_periodic_and_polarization_normalized():
    pgc = model(); period = 2 * np.pi / pgc.wave_number
    np.testing.assert_allclose(pgc.electric_field([0, 0, 0]), pgc.electric_field([period, 0, 0]), atol=1e-12)
    _, fractions = pgc.polarization_components([0.137 * period, 0, 0])
    assert sum(fractions.values()) == pytest.approx(1.0)


def test_optical_pumping_is_probability_conserving_and_physical():
    pgc = model(); generator = pgc.pumping_generator([1e-8, 2e-8, 0])
    np.testing.assert_allclose(generator.sum(axis=0), 0, atol=1e-9)
    off_diagonal = generator.copy(); np.fill_diagonal(off_diagonal, 0)
    assert np.all(off_diagonal >= 0)
    p = pgc.stationary_populations([1e-8, 2e-8, 0])
    assert p.sum() == pytest.approx(1); assert np.all(p >= 0)


def test_light_shifts_and_forces_are_finite_and_state_resolved():
    pgc = model()
    assert pgc.light_shifts([0, 0, 0]).shape == (5,)
    force = pgc.state_forces([1e-8, 0, 0])
    assert force.shape == (5, 3); assert np.isfinite(force).all()
    assert np.ptp(pgc.light_shifts([1e-8, 0, 0])) > 0


def test_residual_axial_field_gives_correct_zeeman_spacing():
    atom = Rb87D2(); zero = model(); biased = model(magnetic_field_t=[0, 0, 2e-5])
    difference = biased.light_shifts([0, 0, 0]) - zero.light_shifts([0, 0, 0])
    # Direct ground Zeeman spacing dominates; the small nonuniform correction
    # is the physically retained excited-minus-ground transition shift.
    np.testing.assert_allclose(np.diff(difference), np.diff(difference)[0], rtol=5e-4)
    assert np.diff(difference)[0] > 0


def test_velocity_force_and_resolution_refinement_converge():
    pgc = model()
    coarse = pgc.moving_average_force(0.03, periods=8, discard=4, steps_per_period=24)
    fine = pgc.moving_average_force(0.03, periods=8, discard=4, steps_per_period=48)
    assert coarse == pytest.approx(fine, rel=0.08, abs=2e-24)
