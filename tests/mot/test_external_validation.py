"""Optional independent-software validation; public APIs only."""
import numpy as np
import pytest

from cold_atom_mot.physics.optical_bloch import TwoLevelOBE
from cold_atom_mot.atomic.species import get_atomic_line,MU_B
from cold_atom_mot.atomic.zeeman import hyperfine_zeeman_hamiltonian
from scipy.constants import hbar


qutip = pytest.importorskip("qutip")


@pytest.mark.parametrize("saturation,detuning", [(.02,-2.0),(.4,-.5),(2.,0.)])
def test_two_level_steady_state_and_liouvillian_against_qutip(saturation,detuning):
    gamma=1.7; internal=TwoLevelOBE.from_saturation(gamma,detuning*gamma,saturation)
    omega=internal.rabi_frequency; h=qutip.Qobj([[0,np.conjugate(omega)/2],[omega/2,-internal.detuning]])
    collapse=np.sqrt(gamma)*qutip.basis(2,0)*qutip.basis(2,1).dag()
    external=qutip.steadystate(h,[collapse]).full()
    np.testing.assert_allclose(internal.steady_state(),external,rtol=2e-12,atol=2e-13)
    # Internal TwoLevelOBE uses row stacking; QuTiP uses column stacking.
    permutation=np.eye(4)[[0,2,1,3]]
    external_liouvillian=permutation@qutip.liouvillian(h,[collapse]).full()@permutation.T
    np.testing.assert_allclose(internal.liouvillian(),external_liouvillian,atol=2e-14)


def test_rabi_oscillation_and_spontaneous_decay_against_qutip():
    gamma=.3; model=TwoLevelOBE(gamma,0,1.1); initial=np.array([[1,0],[0,0]],complex)
    times=np.linspace(0,8,101); _,density=model.evolve(initial,times[-1],max_step=.02)
    # Internal adaptive output is not on t_eval; compare through an independent
    # solve at its returned times rather than interpolating populations.
    internal_time,internal=model.evolve(initial,times[-1],max_step=.02)
    h=qutip.Qobj(model.hamiltonian_over_hbar); c=np.sqrt(gamma)*qutip.basis(2,0)*qutip.basis(2,1).dag()
    external=qutip.mesolve(h,qutip.Qobj(initial),internal_time,c_ops=[c],
                           options={"rtol":1e-11,"atol":1e-13}).states
    np.testing.assert_allclose(internal,np.array([rho.full() for rho in external]),atol=2e-8)
    decay=TwoLevelOBE(1.3,0,0); t,rho=decay.evolve([[0,0],[0,1]],4,max_step=.02)
    np.testing.assert_allclose(rho[:,1,1].real,np.exp(-1.3*t),atol=3e-8)


def test_pylcp_two_beam_force_matches_internal_formula():
    pylcp=pytest.importorskip("pylcp")
    s,delta=.05,-2.; params=[dict(kvec=np.array([1.,0,0]),pol=1,s=s,delta=delta),
                            dict(kvec=np.array([-1.,0,0]),pol=1,s=s,delta=delta)]
    beams=pylcp.laserBeams(params,beam_type=pylcp.infinitePlaneWaveBeam)
    external=pylcp.heuristiceq(beams,pylcp.constantMagneticField([0,0,0]),gamma=1,k=1,mass=1)
    velocities=np.linspace(-.3,.3,13); pylcp_force=[]; internal=[]
    for velocity in velocities:
        pylcp_force.append(external.force(np.zeros(3),[velocity,0,0],0)[0][0])
        plus=.5*s/(1+2*s+(2*(delta-velocity))**2)
        minus=.5*s/(1+2*s+(2*(delta+velocity))**2)
        internal.append(plus-minus)
    np.testing.assert_allclose(internal,pylcp_force,rtol=2e-14,atol=2e-16)


def test_rb87_vector_zeeman_spectrum_against_pylcp():
    pylcp=pytest.importorskip("pylcp"); line=get_atomic_line("87Rb","D2")
    # PyLCP's public hyperfine helper defines gI with the opposite sign to this
    # repository's Steck Hamiltonian convention; make that convention mapping explicit.
    h0,mu=pylcp.hamiltonians.hyperfine_coupled(.5,1.5,line.ground_g_j,
        -line.species.nuclear_g_factor,line.species.ground_hyperfine_a_hz,
        muB=MU_B/hbar/(2*np.pi))
    for field in (0,1e-6,1e-4):
        external=np.linalg.eigvalsh(h0+mu[1]*field)
        internal=np.linalg.eigvalsh(hyperfine_zeeman_hamiltonian(line,"ground",[0,0,field]))/(2*np.pi)
        np.testing.assert_allclose(internal,external,rtol=0,atol=.7)
