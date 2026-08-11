import numpy as np
import pytest
from cold_atom_mot.physics.optical_bloch import TwoLevelOBE


def test_steady_state_matches_analytic_two_level_limit():
    for saturation in (0.01, 0.2, 1.0, 10.0):
        for detuning_gamma in (-3.0, -0.5, 0.0, 1.0):
            model = TwoLevelOBE.from_saturation(2.0, detuning_gamma * 2.0, saturation)
            rho = model.steady_state()
            assert rho[1, 1].real == pytest.approx(model.analytic_excited_population(), rel=2e-12)
            assert np.trace(rho) == pytest.approx(1.0)
            np.testing.assert_allclose(rho, rho.conj().T, atol=1e-13)
            assert np.linalg.eigvalsh(rho).min() > -1e-12


def test_time_evolution_preserves_trace_and_approaches_steady_state():
    model = TwoLevelOBE.from_saturation(1.0, -0.5, 2.0)
    time, density = model.evolve(np.array([[1, 0], [0, 0]], complex), 30.0, max_step=0.1)
    np.testing.assert_allclose(np.trace(density, axis1=1, axis2=2), 1.0, atol=2e-9)
    np.testing.assert_allclose(density[-1], model.steady_state(), atol=2e-6)
    assert time[-1] == pytest.approx(30.0)


def test_zero_drive_has_ground_steady_state_and_zero_force():
    model = TwoLevelOBE.from_saturation(1.0, 0.0, 0.0)
    np.testing.assert_allclose(model.steady_state(), [[1, 0], [0, 0]], atol=1e-14)
    np.testing.assert_allclose(model.scattering_force([1, 0, 0]), 0.0)


def test_tolerance_refinement_converges():
    model = TwoLevelOBE.from_saturation(1.0, -1.0, 1.5)
    initial = np.array([[1, 0], [0, 0]], complex)
    _, coarse = model.evolve(initial, 8.0, rtol=1e-6, atol=1e-8, max_step=0.2)
    _, fine = model.evolve(initial, 8.0, rtol=1e-9, atol=1e-11, max_step=0.05)
    np.testing.assert_allclose(coarse[-1], fine[-1], atol=2e-6)
