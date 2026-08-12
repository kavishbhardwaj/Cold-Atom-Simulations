import numpy as np
import pytest
from scipy.constants import hbar

from cold_atom_mot.atomic.species import MU_B, get_atomic_line
from cold_atom_mot.atomic.zeeman import (
    coupled_transformation,
    hyperfine_zeeman_hamiltonian,
    linear_zeeman_energies,
)
from cold_atom_mot.io.config import build_multilevel_model, load_config
from cold_atom_mot.physics.multilevel_obe import MultilevelOBE


@pytest.mark.parametrize("isotope", ["87Rb", "85Rb"])
@pytest.mark.parametrize("manifold", ["ground", "excited"])
def test_zero_field_exactly_recovers_hyperfine_energies(isotope, manifold):
    line = get_atomic_line(isotope, "D2")
    matrix = hyperfine_zeeman_hamiltonian(line, manifold, np.zeros(3))
    expected = []
    manifolds = line.ground_f if manifold == "ground" else line.excited_f
    for f in manifolds:
        expected.extend([2*np.pi*line.hyperfine_energy_hz(manifold, f)] * (2*f+1))
    np.testing.assert_allclose(matrix, np.diag(expected), rtol=2e-13, atol=2e-5)


def test_vector_hamiltonian_is_hermitian_and_direction_independent_spectrum():
    line = get_atomic_line("87Rb", "D2")
    magnitude = 3.7e-4
    directions = ([magnitude, 0, 0], [0, magnitude, 0],
                  np.asarray([1, 2, -3])*magnitude/np.sqrt(14))
    spectra = []
    for field in directions:
        matrix = hyperfine_zeeman_hamiltonian(line, "excited", field)
        np.testing.assert_allclose(matrix, matrix.conj().T, atol=1e-7)
        spectra.append(np.linalg.eigvalsh(matrix))
    np.testing.assert_allclose(spectra[0], spectra[1], rtol=2e-12, atol=1e-4)
    np.testing.assert_allclose(spectra[0], spectra[2], rtol=2e-12, atol=1e-4)


def test_weak_field_diagonal_slopes_recover_lande_gf():
    line = get_atomic_line("87Rb", "D2")
    step = 1e-6
    plus = hyperfine_zeeman_hamiltonian(line, "ground", [0, 0, step])
    minus = hyperfine_zeeman_hamiltonian(line, "ground", [0, 0, -step])
    slope = np.diag((plus-minus)/(2*step)).real
    expected = [MU_B/hbar*line.lande_gf("ground", f)*m
                for f in line.ground_f for m in range(-f, f+1)]
    np.testing.assert_allclose(slope, expected, rtol=2e-10, atol=1e-2)


def test_continuity_at_zero_and_nonlinear_large_field_correction():
    line = get_atomic_line("87Rb", "D2")
    zero = np.linalg.eigvalsh(hyperfine_zeeman_hamiltonian(line, "ground", [0,0,0]))
    tiny = np.linalg.eigvalsh(hyperfine_zeeman_hamiltonian(line, "ground", [1e-14,-2e-14,3e-14]))
    np.testing.assert_allclose(tiny, zero, rtol=0, atol=4e-3)
    field = np.array([0, 0, 0.1])
    exact = np.linalg.eigvalsh(hyperfine_zeeman_hamiltonian(line, "ground", field))
    linear = np.sort(linear_zeeman_energies(line, "ground", field))
    assert np.max(np.abs(exact-linear))/(2*np.pi) > 1e8


def test_coupled_transformation_is_unitary():
    transform, states = coupled_transformation(1.5, 1.5)
    assert len(states) == 16
    np.testing.assert_allclose(transform.conj().T@transform, np.eye(16), atol=1e-14)


def test_multilevel_obe_with_transverse_field_preserves_trace_and_hermiticity():
    rate = build_multilevel_model(load_config("configs/rb87_d2_multilevel.yaml"))
    rate.magnetic_field = type("Uniform", (), {"field": lambda self, position, time=0: np.array([2e-5,-3e-5,1e-5])})()
    obe = MultilevelOBE(rate.basis, rate.beam_families[:1], rate.magnetic_field)
    h = obe.hamiltonian(np.zeros(3))
    np.testing.assert_allclose(h, h.conj().T, atol=1e-5)
    liouvillian = obe.liouvillian(np.zeros(3))
    n = rate.basis.state_count
    trace = np.zeros(n*n); trace[::n+1] = 1
    np.testing.assert_allclose(trace@liouvillian.toarray(), 0, atol=3e-6)
