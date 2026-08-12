from dataclasses import replace
import numpy as np
import pytest

from cold_atom_mot.laser.apparatus import Retroreflection, SixBeamApparatus
from cold_atom_mot.laser.beam import GaussianBeam, six_beam_mot
from cold_atom_mot.laser.polarization import JonesElement, propagate_jones


def beam(**kwargs):
    values=dict(direction=[1,0,0],origin=[0,0,0],power=.01,waist=.008,
                detuning=-1,wavelength=780e-9,helicity=1)
    values.update(kwargs); return GaussianBeam(**values)


def test_ideal_apparatus_exactly_recovers_symmetric_constructor():
    old=six_beam_mot(.01,.008,-2,780e-9)
    apparatus=SixBeamApparatus(tuple(old))
    for left,right in zip(apparatus.beams,old):
        assert left==right


def test_jones_normalization_and_ideal_qwp_circularity():
    vector=propagate_jones([1,0],[JonesElement("quarter_wave",np.pi/4)])
    assert np.linalg.norm(vector)==pytest.approx(1)
    assert abs(abs(vector[0])-1/np.sqrt(2))<1e-14
    assert abs(abs(vector[1])-1/np.sqrt(2))<1e-14
    assert abs(np.imag(vector[0]*vector[1].conjugate()))==pytest.approx(.5)


def test_qwp_angle_and_retardance_errors_are_continuous():
    angles=np.linspace(np.pi/4-.01,np.pi/4+.01,9)
    values=np.array([propagate_jones([1,0],[JonesElement("quarter_wave",a,
                          retardance_error=.02)]) for a in angles])
    assert np.max(np.linalg.norm(np.diff(values,axis=0),axis=1))<.01


def test_rotation_covariance_of_local_spherical_fractions():
    original=beam(jones_vector=np.array([1,1j])/np.sqrt(2))
    rotation=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
    rotated=replace(original,direction=rotation@original.direction)
    # Jones coefficients are intrinsic to the transported transverse basis.
    a=SixBeamApparatus(tuple([original]*6)).local_polarizations([0,0,1])[0]
    b=SixBeamApparatus(tuple([rotated]*6)).local_polarizations([0,0,1])[0]
    assert a==pytest.approx(b)


def test_retroreflection_reverses_k_and_applies_loss_and_tilt():
    incoming=beam(jones_vector=[1,0],coherence_group="pair")
    exact=Retroreflection(.8,np.pi,coherence_group="pair").reflected(incoming)
    np.testing.assert_allclose(exact.k_vector,-incoming.k_vector)
    assert exact.power==pytest.approx(.8*incoming.power)
    tilted=Retroreflection(mirror_tilt=(0,1e-3,0)).reflected(incoming)
    assert tilted.direction[1]>0 and np.dot(tilted.direction,incoming.direction)<0


def test_full_gaussian_propagation_and_elliptical_waists():
    gaussian=beam(waist=.001,waist_y=.002,propagation_mode="gaussian")
    assert gaussian.intensity([gaussian.rayleigh_range,0,0]) < gaussian.intensity([0,0,0])
    assert gaussian.intensity([0,.001,0]) != gaussian.intensity([0,0,.001])


def test_linewidth_changes_physical_scattering_rate():
    from cold_atom_mot.io.config import build_effective_model,load_config
    model=build_effective_model(load_config("configs/rb87_d2_mot.yaml"))
    narrow=model.scattering_rates(np.zeros(3),np.zeros(3)).sum()
    model.beams=[replace(b,linewidth=model.atom.gamma_rad_s) for b in model.beams]
    broad=model.scattering_rates(np.zeros(3),np.zeros(3)).sum()
    assert broad != pytest.approx(narrow)
