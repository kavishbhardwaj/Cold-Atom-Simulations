import numpy as np
import pytest
from dataclasses import replace

from cold_atom_mot.io.config import build_multilevel_model, load_config
from cold_atom_mot.physics.multilevel_obe import MultilevelOBE


def solver():
    rate = build_multilevel_model(load_config("configs/rb87_d2_multilevel.yaml"))
    # A single circular beam has an exactly degenerate dark kernel.  Request a
    # tiny, explicit CP mixing rate for tests that need a unique stationary
    # representative; production defaults remain exactly zero.
    return MultilevelOBE(rate.basis, rate.beam_families[:1], rate.magnetic_field,
                         ground_relaxation_rate=1e-4*rate.atom.gamma_rad_s)


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


def full_solver(*, phase_samples=2):
    rate = build_multilevel_model(load_config("configs/rb87_d2_multilevel.yaml"))
    return MultilevelOBE(rate.basis, rate.beam_families, rate.magnetic_field,
                         phase_samples=phase_samples)


def test_six_cooling_plus_six_repump_block_frame_and_short_evolution():
    obe = full_solver()
    assert len(obe.beam_families) == 12
    carriers = obe._manifold_carriers()
    removed_beat = abs(carriers[2]-carriers[1])
    assert removed_beat/(2*np.pi) > 6e9
    np.testing.assert_allclose(obe.retained_beat_frequencies(), 0, atol=1e-6)
    times = np.linspace(0, 3/obe.basis.line.gamma_rad_s, 9)
    density, diagnostics = obe.evolve(np.zeros(3), np.zeros(3), times,
                                      max_step=0.2/obe.basis.line.gamma_rad_s,
                                      return_diagnostics=True)
    np.testing.assert_allclose(np.trace(density, axis1=1, axis2=2), 1, atol=2e-8)
    np.testing.assert_allclose(density, density.swapaxes(1, 2).conj(), atol=2e-9)
    assert min(np.linalg.eigvalsh(r).min() for r in density) > -2e-7
    assert diagnostics["max_step_s"] <= 0.21/obe.basis.line.gamma_rad_s
    stationary = obe.steady_state(np.zeros(3), np.zeros(3))
    assert np.trace(stationary) == pytest.approx(1)
    assert np.linalg.eigvalsh(stationary).min() > -2e-9


def test_incoherent_group_phase_is_invariant_but_coherent_pair_is_phase_sensitive():
    base = solver(); first = base.beam_families[0]
    opposite = replace(first, beam=replace(first.beam, direction=-first.beam.direction,
                                            helicity=-first.beam.helicity))
    independent_a = MultilevelOBE(base.basis, [first, opposite], base.magnetic_field, phase_samples=4)
    independent_b = MultilevelOBE(base.basis,
        [first, replace(opposite, beam=replace(opposite.beam, phase=1.234))],
        base.magnetic_field, phase_samples=4)
    point = np.array([0, 0, 0.13*first.beam.wavelength])
    # Absolute phases of singleton incoherent groups are removed before the
    # identical deterministic phase-cycling ensemble is formed.
    for sample in range(4):
        a = independent_a._phase_realization(sample)
        b = independent_b._phase_realization(sample)
        np.testing.assert_allclose(a.hamiltonian(point), b.hamiltonian(point), atol=1e-8)
    coherent_a = MultilevelOBE(base.basis,
        [replace(first, beam=replace(first.beam, coherence_group="pair")),
         replace(opposite, beam=replace(opposite.beam, coherence_group="pair"))], base.magnetic_field)
    coherent_b = MultilevelOBE(base.basis,
        [coherent_a.beam_families[0], replace(coherent_a.beam_families[1],
          beam=replace(coherent_a.beam_families[1].beam, phase=np.pi/2))], base.magnetic_field)
    assert np.linalg.norm(coherent_a.hamiltonian(point)-coherent_b.hamiltonian(point)) > 1e3


def test_analytic_interaction_gradient_matches_finite_difference():
    obe = solver(); point = np.array([1.2e-4, -2.1e-4, 0.7e-4])
    analytic = obe.force_operators(point)
    for scale in (3e-4, 1e-4, 3e-5):
        finite = obe.force_operators(point, finite_difference=True,
                                     step=scale*obe.basis.line.wavelength_m)
        np.testing.assert_allclose(finite, analytic, rtol=2e-5, atol=2e-32)


def test_ground_manifold_rwa_bound_is_small_and_regularizer_is_explicit():
    obe = full_solver()
    assert max(obe.cross_ground_rwa_bound(f) for f in obe.beam_families) < 2e-5
    assert obe.ground_relaxation_rate == 0
    with pytest.raises(ValueError, match="non-negative"):
        replace(obe, ground_relaxation_rate=-1)
    tiny = replace(obe, ground_relaxation_rate=1e-10*obe.basis.line.gamma_rad_s)
    rho_zero = obe.steady_state(np.zeros(3), np.zeros(3))
    rho_tiny = tiny.steady_state(np.zeros(3), np.zeros(3))
    np.testing.assert_allclose(rho_zero, rho_tiny, atol=2e-7)
