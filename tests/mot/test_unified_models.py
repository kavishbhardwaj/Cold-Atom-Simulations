import numpy as np
import pytest
from scipy.linalg import eigvalsh
from cold_atom_mot.atomic.species import build_atomic_basis,get_atomic_line
from cold_atom_mot.io.config import load_config,build_multilevel_model
from cold_atom_mot.laser.beam import GaussianBeam,grouped_intensity
from cold_atom_mot.physics.multilevel_obe import MultilevelOBE
from cold_atom_mot.solvers.monte_carlo import at_least_one_event_probability
from cold_atom_mot.vacuum import rubidium_vapor_pressure_pa,number_density,loading_curve

def test_coherence_groups_control_cross_terms():
    args=dict(direction=np.array([0,0,1]),origin=np.zeros(3),power=.001,waist=.01,detuning=0,wavelength=780e-9,helicity=1)
    a=GaussianBeam(**args,phase=0,coherence_group="pair")
    b=GaussianBeam(**args,phase=np.pi,coherence_group="pair")
    independent=GaussianBeam(**args,phase=np.pi,coherence_group="other")
    assert grouped_intensity([a,b],[0,0,0]) < 1e-20
    assert grouped_intensity([a,independent],[0,0,0]) == pytest.approx(2*a.peak_intensity)

def test_exact_poisson_event_probability():
    assert at_least_one_event_probability(2.0,.5)==pytest.approx(1-np.exp(-1))
    assert at_least_one_event_probability(2.0,1e-9)==pytest.approx(2e-9)

def test_85rb_rate_model_uses_genuine_manifolds():
    model=build_multilevel_model(load_config("configs/rb85_d2_multilevel.yaml"))
    assert model.basis.state_count==36
    assert model.atom.cooling_transition==(3,4) and model.atom.repump_transition==(2,3)
    assert model.steady_state(np.zeros(3),np.zeros(3)).sum()==pytest.approx(1)

def test_sparse_multilevel_obe_operators_are_trace_preserving():
    rate=build_multilevel_model(load_config("configs/rb87_d2_multilevel.yaml"))
    obe=MultilevelOBE(rate.basis,rate.beam_families[:1],rate.magnetic_field)
    h=obe.hamiltonian(np.zeros(3)); np.testing.assert_allclose(h,h.conj().T)
    L=obe.liouvillian(np.zeros(3)); n=rate.basis.state_count
    trace=np.zeros(n*n); trace[::n+1]=1
    np.testing.assert_allclose(trace@L.toarray(),0,atol=2e-7)

def test_vapour_pressure_density_and_loading_solution():
    pressure=rubidium_vapor_pressure_pa(300); assert 1e-9 < pressure < 1e-4
    assert number_density(pressure,300)>0
    t=np.linspace(0,3,20); n=loading_curve(t,100,2)
    np.testing.assert_allclose(n,50*(1-np.exp(-2*t)))
