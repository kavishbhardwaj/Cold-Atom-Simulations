import numpy as np
import pytest
from scipy.constants import k as k_B

from cold_atom_mot.physics.collective import GaussianCloud,MultipleScatteringModel,CollectiveLoading
from cold_atom_mot.vacuum import loading_curve


def cloud(n=1e7,sigma=1e-3): return GaussianCloud(n,[sigma]*3,[100e-6]*3,1.443e-25)


def scattering(reabs=1e-13): return MultipleScatteringModel(2e-13,reabs,100,8.05e6,1e6)


def test_gaussian_normalization_two_body_integral_and_optical_depth():
    c=cloud(); assert c.peak_density_m3==(pytest.approx(1e7/((2*np.pi)**1.5*1e-9)))
    assert c.two_body_integral_m3()==pytest.approx(c.atom_number**2/(8*np.pi**1.5*1e-9))
    np.testing.assert_allclose(c.optical_depth(2e-13),1e7*2e-13/(2*np.pi*1e-6))


def test_reabsorption_and_radiation_trapping_vanish_at_low_optical_depth():
    dilute=cloud(1e-9); model=scattering()
    assert model.reabsorption_probability(dilute)<1e-15
    assert model.radiation_trapping_diffusion(dilute)<1e-60


def test_multiple_scattering_is_outward_and_shadow_opposes_transmission():
    model=scattering(); c=cloud()
    assert model.repulsive_force(1e-3,c)>0 and model.repulsive_force(-1e-3,c)<0
    force=model.attenuated_pair_force(1e-20,c)
    assert 0<force<1e-20


def test_density_limited_cloud_expands_as_cube_root_of_population():
    model=scattering(); a=model.equilibrium_sigma(1e6,1e-18,100e-6,1.443e-25)
    b=model.equilibrium_sigma(8e6,1e-18,100e-6,1.443e-25)
    assert b/a==pytest.approx(2)


def test_two_body_beta_requires_source_and_loss_is_exact_gaussian_integral():
    args=dict(loading_rate_s=1e6,background_loss_s=.1,hot_rb_loss_s=.2,
              two_body_coefficient_m3_s=1e-16,temperature_k=100e-6,
              restoring_coefficient_n_m=1e-18,atom_mass_kg=1.443e-25,fixed_sigma_m=1e-3)
    with pytest.raises(ValueError,match="source"): CollectiveLoading(**args)
    model=CollectiveLoading(**args,beta_source="user supplied")
    n=1e7; expected=1e6-.3*n-1e-16*cloud(n).two_body_integral_m3()
    assert model.derivative(n)==pytest.approx(expected)


def test_collective_low_density_limit_recovers_independent_loading_curve():
    times=np.linspace(0,5,51); model=CollectiveLoading(12,.1,.2,0,100e-6,1e-18,1.443e-25,fixed_sigma_m=1e-3)
    result=model.evolve(times)
    np.testing.assert_allclose(result["population"],loading_curve(times,12,.3),rtol=2e-7)


def test_quick_and_research_collective_evolution_converge():
    model=CollectiveLoading(1e6,.1,.05,1e-16,100e-6,1e-18,1.443e-25,
                            scattering=scattering(),beta_source="user supplied")
    times=np.linspace(0,2,51); quick=model.evolve(times,mode="quick")
    research=model.evolve(times,mode="research")
    np.testing.assert_allclose(quick["population"],research["population"],rtol=2e-6)
