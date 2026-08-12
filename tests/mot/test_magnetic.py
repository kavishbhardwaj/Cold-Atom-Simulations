import numpy as np
import pytest
from cold_atom_mot.magnetic.coils import AntiHelmholtzPair, HelmholtzPair, ThreeAxisBiasCoils
from cold_atom_mot.magnetic.fields import IdealQuadrupole, ResidualField, HarmonicResidualField, SwitchingTransientField


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


def test_ideal_helmholtz_is_uniform_symmetric_and_linear_in_current():
    pair=HelmholtzPair.imperfect([0,0,1],.04,.04,1,20,segments=192)
    centre=pair.field(np.zeros(3)); assert centre[2]>0
    assert np.linalg.norm(pair.jacobian()) < 2e-7
    scaled=HelmholtzPair.imperfect([0,0,1],.04,.04,2,20,segments=192)
    np.testing.assert_allclose(scaled.field([.001,0,0]),2*pair.field([.001,0,0]))


def test_nonorthogonal_calibration_least_squares_compensates_background():
    pair=HelmholtzPair.imperfect([0,0,1],.04,.04,1,10,segments=64)
    matrix=np.array([[2,.1,0],[.05,1.8,.1],[0,.07,2.1]])*1e-5
    coils=ThreeAxisBiasCoils((pair,pair,pair),matrix,[1e-7,-2e-7,3e-7])
    background=np.array([2e-5,-1e-5,4e-6]); currents=coils.compensation_currents(background)
    np.testing.assert_allclose(coils.calibrated_field(currents,background),0,atol=1e-18)


def test_multiple_eddy_exponentials_and_measured_waveform():
    field=SwitchingTransientField(1,[1e-6,0,0],[4e-6,0,0],.01,
                                  ((np.array([2e-6,0,0]),.02),))
    assert field.field([0,0,0],1+.02)[0]==pytest.approx(1e-6+4e-6*np.exp(-2)+2e-6/np.e)
    measured=SwitchingTransientField(0,[0,0,0],[0,0,0],1,waveform_time=np.array([0,1]),waveform_field=np.array([[1,0,0],[3,0,0]]))
    assert measured.field([0,0,0],.5)[0]==pytest.approx(2)


def test_harmonic_background_and_rotation_covariance():
    residual=HarmonicResidualField([1e-6,0,0],harmonics=((np.array([0,2e-6,0]),50,0),))
    np.testing.assert_allclose(residual.field([0,0,0],1/200),[1e-6,2e-6,0])
    pair=HelmholtzPair.imperfect([1,0,0],.04,.04,1,10,segments=96)
    rotation=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
    rotated=HelmholtzPair.imperfect(rotation@[1,0,0],.04,.04,1,10,segments=96)
    np.testing.assert_allclose(rotated.field(rotation@[.001,0,0]),rotation@pair.field([.001,0,0]),atol=1e-12)


def test_helmholtz_maxwell_divergence_and_curvature_are_finite():
    pair=HelmholtzPair.imperfect([0,0,1],.04,.04,1,20,segments=192)
    assert np.trace(pair.jacobian([.001,.002,.003]))==pytest.approx(0,abs=2e-8)
    assert np.isfinite(pair.curvature()).all()
