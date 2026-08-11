import numpy as np
import pytest
from cold_atom_mot.io.config import build_effective_model, load_config
from cold_atom_mot.solvers.deterministic import integrate_trajectory
from cold_atom_mot.solvers.monte_carlo import isotropic_directions, simulate_photon_events


@pytest.fixture
def model():
    config = load_config('configs/rb87_standard_mot.yaml')
    config['gravity']['vector_m_per_s2'] = [0, 0, 0]
    return build_effective_model(config)


def test_force_restores_and_damps_near_origin(model):
    for axis in range(3):
        position = np.zeros(3); position[axis] = 1e-4
        velocity = np.zeros(3); velocity[axis] = 0.02
        assert model.force(position, np.zeros(3))[axis] * position[axis] < 0
        assert model.force(np.zeros(3), velocity)[axis] * velocity[axis] < 0
    damping, restoring = model.linear_coefficients()
    assert np.all(damping > 0) and np.all(restoring > 0)


def test_force_symmetry(model):
    position = np.array([1e-4, -2e-4, 3e-4]); velocity = np.array([.01, -.02, .03])
    assert np.allclose(model.force(position, velocity), -model.force(-position, -velocity), rtol=1e-11, atol=1e-27)


def test_deterministic_trajectory_is_three_dimensional(model):
    trajectory = integrate_trajectory(model, [2e-4, 1e-4, -1e-4], [0, 0, 0], 2e-4, max_step=1e-5)
    assert trajectory.position.shape[1] == 3
    assert np.linalg.norm(trajectory.position[-1]) < np.linalg.norm(trajectory.position[0])


def test_isotropic_recoil_statistics():
    directions = isotropic_directions(np.random.default_rng(4), 100_000)
    assert np.linalg.norm(directions.mean(axis=0)) < 0.01
    assert np.allclose(np.var(directions, axis=0), 1 / 3, atol=0.01)


def test_monte_carlo_seed_and_recoil(model):
    kwargs = dict(force_model=model, position=np.zeros((8, 3)), velocity=np.zeros((8, 3)), duration=2e-6, time_step=5e-9, seed=42, store_every=40)
    first = simulate_photon_events(**kwargs); second = simulate_photon_events(**kwargs)
    assert np.array_equal(first.position, second.position)
    assert np.array_equal(first.velocity, second.velocity)
    assert first.scattering_events > 0


def test_timestep_refinement_statistical_mean(model):
    position = np.zeros((300, 3)); velocity = np.tile([0.05, 0, 0], (300, 1))
    coarse = simulate_photon_events(model, position, velocity, 1e-6, 5e-9, seed=9, store_every=200)
    fine = simulate_photon_events(model, position, velocity, 1e-6, 2.5e-9, seed=9, store_every=400)
    assert abs(coarse.velocity[-1, :, 0].mean() - fine.velocity[-1, :, 0].mean()) < 0.003
