import copy
import numpy as np
import pytest
from cold_atom_mot.atomic.levels import build_rb87_d2_basis
from cold_atom_mot.io.config import build_multilevel_model, load_config


def test_full_rb87_basis_and_selection_rules():
    basis = build_rb87_d2_basis()
    assert len(basis.ground) == 8
    assert len(basis.excited) == 16
    assert basis.state_count == 24
    assert all(t.q in (-1, 0, 1) for t in basis.transitions)
    assert all(basis.excited[t.excited_index].m - basis.ground[t.ground_index].m == t.q for t in basis.transitions)
    stretched = [t for t in basis.transitions if basis.ground[t.ground_index].F == 2 and basis.ground[t.ground_index].m == 2 and basis.excited[t.excited_index].F == 3 and basis.excited[t.excited_index].m == 3]
    assert len(stretched) == 1 and stretched[0].strength == pytest.approx(1.0)


def test_spontaneous_branching_normalizes_per_excited_state():
    branching = build_rb87_d2_basis().spontaneous_branching
    np.testing.assert_allclose(branching.sum(axis=1), 1.0, atol=1e-14)
    assert np.all(branching >= 0)


def test_generator_conserves_probability_and_steady_state_is_physical():
    model = build_multilevel_model(load_config('configs/rb87_d2_multilevel.yaml'))
    generator = model.generator(np.zeros(3), np.zeros(3))
    np.testing.assert_allclose(generator.sum(axis=0), 0.0, atol=2e-8)
    population = model.steady_state(np.zeros(3), np.zeros(3))
    assert population.sum() == pytest.approx(1.0)
    assert np.all(population >= 0)
    assert model.manifold_populations(population)['ground_F2'] > 0.5


def test_repump_recovers_f2_population():
    config = load_config('configs/rb87_d2_multilevel.yaml')
    with_repump = build_multilevel_model(config).steady_state(np.zeros(3), np.zeros(3))
    weak = copy.deepcopy(config)
    weak['repump']['power_per_beam_w'] = 1e-9
    weak_model = build_multilevel_model(weak)
    without_repump = weak_model.steady_state(np.zeros(3), np.zeros(3))
    full_model = build_multilevel_model(config)
    assert full_model.manifold_populations(with_repump)['ground_F2'] > weak_model.manifold_populations(without_repump)['ground_F2']


def test_multilevel_force_restores_and_damps():
    config = load_config('configs/rb87_d2_multilevel.yaml')
    config['gravity']['vector_m_per_s2'] = [0, 0, 0]
    model = build_multilevel_model(config)
    assert model.force(np.array([1e-4, 0, 0]), np.zeros(3))[0] < 0
    assert model.force(np.zeros(3), np.array([0.02, 0, 0]))[0] < 0
