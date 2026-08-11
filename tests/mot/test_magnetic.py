import numpy as np
import pytest
from cold_atom_mot.magnetic.coils import AntiHelmholtzPair
from cold_atom_mot.magnetic.fields import IdealQuadrupole, ResidualField


def test_ideal_quadrupole_is_divergence_free_and_symmetric():
    field = IdealQuadrupole(0.1)
    assert np.trace(field.gradient) == pytest.approx(0.0, abs=1e-15)
    assert np.allclose(field.field([1e-3, 2e-3, 3e-3]), [1e-4, 2e-4, -6e-4])


def test_symmetric_antihelmholtz_zero_and_maxwell_gradient():
    pair = AntiHelmholtzPair.symmetric(0.04, 0.04, 2.0, 20, segments=256)
    assert np.linalg.norm(pair.field(np.zeros(3))) < 1e-12
    jac = pair.jacobian(np.zeros(3), step=2e-5)
    assert np.trace(jac) == pytest.approx(0.0, abs=2e-7)
    assert jac[0, 0] == pytest.approx(jac[1, 1], rel=2e-3)
    assert jac[2, 2] == pytest.approx(-2 * jac[0, 0], rel=3e-3)


def test_segmented_coil_converges_and_tilt_moves_zero():
    coarse = AntiHelmholtzPair.symmetric(0.04, 0.04, 2.0, 20, segments=128)
    fine = AntiHelmholtzPair.symmetric(0.04, 0.04, 2.0, 20, segments=256)
    point = np.array([0.003, 0.001, 0.002])
    assert np.linalg.norm(coarse.field(point) - fine.field(point)) / np.linalg.norm(fine.field(point)) < 1e-3
    tilted = AntiHelmholtzPair.symmetric(0.04, 0.04, 2.0, 20, segments=128, tilt_y=np.deg2rad(1), lateral_offset=5e-4)
    assert np.linalg.norm(tilted.field_zero()) > 1e-6


def test_uniform_stray_field_displaces_ideal_zero():
    quadrupole = IdealQuadrupole(0.1)
    stray = ResidualField(uniform=np.array([2e-5, 0, 0]))
    expected = np.array([-2e-4, 0, 0])
    assert np.linalg.norm(quadrupole.field(expected) + stray.field(expected)) < 1e-15
