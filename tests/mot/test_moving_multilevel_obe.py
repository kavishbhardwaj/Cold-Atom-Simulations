import numpy as np
import pytest

from cold_atom_mot.io.config import build_multilevel_model, load_config
from cold_atom_mot.physics.multilevel_obe import MultilevelOBE


def solver():
    rate = build_multilevel_model(load_config("configs/rb87_d2_multilevel.yaml"))
    return MultilevelOBE(rate.basis, rate.beam_families[:1], rate.magnetic_field)


def test_rb87_basis_is_complete_24_state_d2_basis():
    obe = solver()
    assert obe.basis.state_count == 24
    assert {state.F for state in obe.basis.ground} == {1, 2}
    assert {state.F for state in obe.basis.excited} == {0, 1, 2, 3}


def test_trajectory_phase_has_per_beam_doppler_sign():
    obe = solver(); family = obe.beam_families[0]
    velocity = 0.17*family.beam.direction
    dt = 1e-10
    c0 = obe._beam_coupling(family, np.zeros(3), velocity, 0)
    c1 = obe._beam_coupling(family, np.zeros(3), velocity, dt)
    ng = len(obe.basis.ground)
    index = np.argwhere(np.abs(c0[ng:, :ng]) > 0)[0] + np.array([ng, 0])
    observed = np.angle(c1[tuple(index)]/c0[tuple(index)])/dt
    expected_sign = np.sign(np.dot(family.beam.k_vector, velocity)-obe._laser_offset(family))
    assert np.sign(observed) == expected_sign


def test_velocity_changes_single_beam_force():
    obe = solver(); direction = obe.beam_families[0].beam.direction
    zero = np.dot(obe.force(np.zeros(3), np.zeros(3)), direction)
    shifted = np.dot(obe.force(np.zeros(3), 4.0*direction), direction)
    assert abs(shifted-zero) > 1e-3*max(abs(zero), abs(shifted))


def test_density_matrix_is_hermitian_normalized_and_positive_within_tolerance():
    rho = solver().steady_state(np.zeros(3), np.array([0.2, 0, 0]))
    np.testing.assert_allclose(rho, rho.conj().T, atol=1e-11)
    assert np.trace(rho) == pytest.approx(1)
    assert np.linalg.eigvalsh(rho).min() > -2e-9


def test_off_resonant_excited_hyperfine_couplings_are_present():
    obe = solver(); family = obe.beam_families[0]
    coupling = obe._beam_coupling(family, np.zeros(3), np.zeros(3), 0)
    ng = len(obe.basis.ground)
    coupled_f = {state.F for i, state in enumerate(obe.basis.excited)
                 if np.any(np.abs(coupling[ng+i, :]) > 0)}
    assert family.target_excited_f in coupled_f
    assert coupled_f - {family.target_excited_f}
