"""Sparse multilevel optical-Bloch operators for a generated atomic basis."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import hbar, physical_constants
from scipy.sparse import csr_matrix, eye, kron
from scipy.sparse.linalg import spsolve
from ..laser.polarization import spherical_fractions

MU_B=physical_constants["Bohr magneton"][0]


@dataclass
class MultilevelOBE:
    """Hyperfine/Zeeman density-matrix model with explicit laser couplings.

    Beam families use the same basis and saturation convention as the population
    solver. The rotating-frame construction assumes one frequency per addressed
    ground manifold; arbitrary multi-frequency Floquet physics is out of scope.
    """
    basis: object
    beam_families: list
    magnetic_field: object

    def hamiltonian(self, position, velocity=(0,0,0), time=0.0):
        n=len(self.basis.ground)+len(self.basis.excited); ng=len(self.basis.ground)
        h=np.zeros((n,n),complex); b=np.asarray(self.magnetic_field.field(position,time)); bmag=np.linalg.norm(b)
        for i,state in enumerate(self.basis.ground): h[i,i]=2*np.pi*state.frequency_offset_hz+MU_B/hbar*state.g_factor*state.m*bmag
        for i,state in enumerate(self.basis.excited): h[ng+i,ng+i]=2*np.pi*state.frequency_offset_hz+MU_B/hbar*state.g_factor*state.m*bmag
        for family in self.beam_families:
            axis=b/bmag if bmag>1e-12 else np.array([0.,0.,1.]); fractions=spherical_fractions(family.beam.polarization,axis)
            s=float(family.beam.intensity(position)/self.basis.line.saturation_intensity_w_m2)
            target=next(e.frequency_offset_hz for e in self.basis.excited if e.F==family.target_excited_f)
            for transition in self.basis.transitions:
                g=self.basis.ground[transition.ground_index]; e=self.basis.excited[transition.excited_index]
                if g.F != family.ground_f: continue
                omega=self.basis.line.gamma_rad_s*np.sqrt(s*transition.strength*fractions[transition.q]/2)
                optical_phase=np.dot(family.beam.k_vector,np.asarray(position))+family.beam.phase
                coupling=.5*omega*np.exp(1j*optical_phase)
                h[ng+transition.excited_index,transition.ground_index]+=coupling
                h[transition.ground_index,ng+transition.excited_index]+=np.conjugate(coupling)
                h[ng+transition.excited_index,ng+transition.excited_index] += -(family.beam.detuning+2*np.pi*target)/max(1,len(self.beam_families))
        return h

    def collapse_operators(self):
        n=self.basis.state_count; ng=len(self.basis.ground); operators=[]
        for ei,row in enumerate(self.basis.spontaneous_branching):
            for gi,branch in enumerate(row):
                if branch>0:
                    op=np.zeros((n,n),complex); op[gi,ng+ei]=np.sqrt(self.basis.line.gamma_rad_s*branch); operators.append(csr_matrix(op))
        return operators

    def liouvillian(self, position, velocity=(0,0,0), time=0.0):
        h=csr_matrix(self.hamiltonian(position,velocity,time)); n=h.shape[0]; ident=eye(n,format="csr",dtype=complex)
        L=-1j*(kron(ident,h)-kron(h.T,ident))
        for c in self.collapse_operators():
            cd_c=c.getH()@c
            L += kron(c.conjugate(),c)-.5*kron(ident,cd_c)-.5*kron(cd_c.T,ident)
        return L.tocsr()

    def steady_state(self, position=(0,0,0), velocity=(0,0,0)):
        n=self.basis.state_count; matrix=self.liouvillian(position,velocity).tolil(); rhs=np.zeros(n*n,complex)
        trace_row=np.zeros(n*n,complex); trace_row[::n+1]=1
        matrix[-1,:]=trace_row; rhs[-1]=1
        rho=spsolve(matrix.tocsr(),rhs).reshape((n,n),order="F")
        return (rho+rho.conj().T)/2
