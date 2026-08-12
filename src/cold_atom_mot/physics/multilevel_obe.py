"""Sparse, moving-atom multilevel optical Bloch equations.

Hamiltonians in this module are angular-frequency Hamiltonians (``H / hbar``).
The excited manifold is transformed at one *reference* optical frequency.  A
laser whose frequency differs from that reference consequently remains an
explicitly time-dependent term.  This avoids assigning an atom -- which has
only one excited-state phase -- a different rotating frame for every laser.
"""
from dataclasses import dataclass, field
import numpy as np
from scipy.constants import hbar
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix, eye, kron
from scipy.sparse.linalg import spsolve
from ..laser.polarization import spherical_fractions
from ..atomic.zeeman import hyperfine_zeeman_hamiltonian


@dataclass
class MultilevelOBE:
    """Hyperfine/Zeeman OBE with phase-resolved travelling-wave lasers.

    ``position`` is the position at ``time=0`` and the optical phase is
    ``k.(position + velocity*time) - omega_offset*time + phase``.  Thus every
    beam has its own Doppler shift ``k.velocity`` and frequency/AOM offset.
    ``mode='quick'`` uses fewer periods/samples than ``mode='research'`` when
    averaging a genuinely time-dependent, multi-frequency configuration.
    """
    basis: object
    beam_families: list
    magnetic_field: object
    mode: str = "quick"
    _collapses: list = field(default=None, init=False, repr=False)

    def __post_init__(self):
        if self.mode not in ("quick", "research"):
            raise ValueError("mode must be 'quick' or 'research'")

    @property
    def _ground_reference_hz(self):
        # The cycling-transition ground state is a convenient but arbitrary
        # global energy origin; physics is invariant under this scalar shift.
        cooling_f = max(self.basis.line.ground_f)
        return self.basis.line.hyperfine_energy_hz("ground", cooling_f)

    def _laser_offset(self, family):
        """Laser angular frequency relative to the cycling reference."""
        target = self.basis.line.hyperfine_energy_hz("excited", family.target_excited_f)
        excited_ref = self.basis.line.hyperfine_energy_hz("excited", max(self.basis.line.excited_f))
        ground = self.basis.line.hyperfine_energy_hz("ground", family.ground_f)
        return 2*np.pi*((target-excited_ref) - (ground-self._ground_reference_hz)) + \
            family.beam.detuning + family.beam.frequency_offset

    def _bare_hamiltonian(self, position, time):
        ng = len(self.basis.ground); n = self.basis.state_count
        b = np.asarray(self.magnetic_field.field(position, time), float)
        h = np.zeros((n, n), complex)
        h[:ng, :ng] = hyperfine_zeeman_hamiltonian(self.basis.line, "ground", b)
        h[:ng, :ng] -= 2*np.pi*self._ground_reference_hz*np.eye(ng)
        h[ng:, ng:] = hyperfine_zeeman_hamiltonian(self.basis.line, "excited", b)
        h[ng:, ng:] -= 2*np.pi*self.basis.line.hyperfine_energy_hz(
            "excited", max(self.basis.line.excited_f))*np.eye(n-ng)
        return h

    def _beam_coupling(self, family, position, velocity, time):
        """Return one beam's RWA coupling matrix, including its trajectory phase."""
        ng = len(self.basis.ground); matrix = np.zeros((self.basis.state_count,)*2, complex)
        beam = family.beam
        r = np.asarray(position, float) + np.asarray(velocity, float)*time
        fractions = spherical_fractions(beam.polarization, np.array([0., 0., 1.]))
        saturation = float(beam.intensity(r)/self.basis.line.saturation_intensity_w_m2)
        phase = np.dot(beam.k_vector, r-beam.origin) - self._laser_offset(family)*time + beam.phase
        for transition in self.basis.transitions:
            ground = self.basis.ground[transition.ground_index]
            if ground.F != family.ground_f:
                continue
            omega = self.basis.line.gamma_rad_s*np.sqrt(
                saturation*transition.strength*fractions[transition.q]/2)
            value = .5*omega*np.exp(1j*phase)
            ei, gi = ng+transition.excited_index, transition.ground_index
            matrix[ei, gi] += value
            matrix[gi, ei] += value.conjugate()
        return matrix

    def hamiltonian(self, position, velocity=(0, 0, 0), time=0.0):
        """Return ``H/hbar`` at time along ``r(t)=position+velocity*time``."""
        r = np.asarray(position, float) + np.asarray(velocity, float)*time
        h = self._bare_hamiltonian(r, time)
        for family in self.beam_families:
            h += self._beam_coupling(family, position, velocity, time)
        return h

    def collapse_operators(self):
        if self._collapses is None:
            ng = len(self.basis.ground); operators = []
            for ei, row in enumerate(self.basis.spontaneous_branching):
                for gi, branch in enumerate(row):
                    if branch > 0:
                        op = csr_matrix(([np.sqrt(self.basis.line.gamma_rad_s*branch)],
                                         ([gi], [ng+ei])), shape=(self.basis.state_count,)*2)
                        operators.append(op)
            self._collapses = operators
        return self._collapses

    def liouvillian(self, position, velocity=(0, 0, 0), time=0.0):
        h = csr_matrix(self.hamiltonian(position, velocity, time)); n = h.shape[0]
        ident = eye(n, format="csr", dtype=complex)
        operator = -1j*(kron(ident, h)-kron(h.T, ident))
        for collapse in self.collapse_operators():
            cd_c = collapse.getH()@collapse
            operator += kron(collapse.conjugate(), collapse)-.5*kron(ident, cd_c)-.5*kron(cd_c.T, ident)
        return operator.tocsr()

    def steady_state(self, position=(0, 0, 0), velocity=(0, 0, 0), time=0.0):
        """Instantaneous stationary state (appropriate for a stationary Hamiltonian).

        Multi-frequency/multi-direction moving configurations should use
        :meth:`evolve` and time-average :meth:`per_beam_force` instead.
        """
        effective = np.asarray([self._laser_offset(f)-np.dot(f.beam.k_vector, velocity)
                                for f in self.beam_families])
        if len(effective) and not np.allclose(effective, effective[0], rtol=0,
                                              atol=1e-7*self.basis.line.gamma_rad_s):
            raise ValueError("no stationary rotating frame: use evolve() for unequal laser/Doppler frequencies")
        n = self.basis.state_count
        # When all fields share a trajectory-frame frequency, one additional
        # excited-manifold rotation makes the Hamiltonian stationary.  Its
        # diagonal -delta is precisely the per-beam Doppler detuning.
        h = self.hamiltonian(position, velocity, time)
        if len(effective):
            ng = len(self.basis.ground)
            h[ng:, ng:] -= effective[0]*np.eye(n-ng)
        ident = eye(n, format="csr", dtype=complex); hs = csr_matrix(h)
        matrix = (-1j*(kron(ident, hs)-kron(hs.T, ident))).tolil()
        for collapse in self.collapse_operators():
            cd_c = collapse.getH()@collapse
            matrix += (kron(collapse.conjugate(), collapse)-.5*kron(ident, cd_c)-
                       .5*kron(cd_c.T, ident)).tolil()
        # Infinitesimal isotropic ground relaxation selects a unique member of
        # an exactly degenerate dark-state kernel (a common single-beam test
        # case).  At 1e-10 Gamma it is far below reported force precision and
        # represents the zero-relaxation limit, while avoiding an arbitrary,
        # possibly non-positive sparse-nullspace solution.
        ng = len(self.basis.ground); relaxation = 1e-10*self.basis.line.gamma_rad_s
        for source in range(ng):
            matrix[source+n*source, source+n*source] -= relaxation
            for target in range(ng):
                matrix[target+n*target, source+n*source] += relaxation/ng
        rhs = np.zeros(n*n, complex)
        trace_row = np.zeros(n*n, complex); trace_row[::n+1] = 1
        matrix[-1, :] = trace_row; rhs[-1] = 1
        rho = spsolve(matrix.tocsr(), rhs).reshape((n, n), order="F")
        rho = (rho+rho.conj().T)/2
        trace = np.trace(rho).real
        if not np.isfinite(rho).all() or abs(trace) < 1e-12:
            # A single polarization can have an exactly degenerate dark-state
            # kernel.  Select the physical state reached from an unpolarized
            # ground ensemble rather than an arbitrary null-space vector.
            duration = (20 if self.mode == "quick" else 60)/self.basis.line.gamma_rad_s
            return self.evolve(position, velocity, np.linspace(0, duration, 81))[-1]
        rho /= trace
        eigenvalues, vectors = np.linalg.eigh(rho)
        if eigenvalues.min() < -1e-10:
            # Singular dark manifolds also make row-replaced null solves
            # numerically indefinite.  Projection selects a positive density
            # matrix while preserving its supported stationary subspace.
            eigenvalues = np.maximum(eigenvalues, 0)
            rho = (vectors*eigenvalues)@vectors.conj().T
            rho /= np.trace(rho)
        return rho

    def evolve(self, position, velocity, times, rho0=None, *, rtol=None, atol=None):
        """Integrate the explicitly time-dependent master equation."""
        times = np.asarray(times, float)
        if times.ndim != 1 or len(times) < 2 or np.any(np.diff(times) <= 0):
            raise ValueError("times must be a strictly increasing 1D array")
        n = self.basis.state_count
        if rho0 is None:
            rho0 = np.zeros((n, n), complex)
            rho0[np.arange(len(self.basis.ground)), np.arange(len(self.basis.ground))] = 1/len(self.basis.ground)
        tolerances = (2e-5, 2e-8) if self.mode == "quick" else (2e-7, 2e-10)
        solution = solve_ivp(lambda t, y: self.liouvillian(position, velocity, t)@y,
                             (times[0], times[-1]), np.asarray(rho0).reshape(-1, order="F"),
                             t_eval=times, rtol=rtol or tolerances[0], atol=atol or tolerances[1],
                             method="DOP853")
        if not solution.success:
            raise RuntimeError(solution.message)
        density = solution.y.T.reshape((-1, n, n), order="F")
        return np.asarray([(r+r.conj().T)/(2*np.trace(r).real) for r in density])

    def force_operators(self, position, velocity=(0, 0, 0), time=0.0):
        """Per-beam optical force operators ``-gradient(H_int)`` in newtons.

        The derivative includes both travelling-wave phase and the Gaussian
        envelope.  A symmetric spatial derivative is used so displaced beams
        retain their dipole-force contribution as well as radiation pressure.
        """
        position = np.asarray(position, float); step = max(self.basis.line.wavelength_m*1e-4, 1e-11)
        output = np.empty((len(self.beam_families), 3, self.basis.state_count, self.basis.state_count), complex)
        for bi, family in enumerate(self.beam_families):
            for axis in range(3):
                shift = np.zeros(3); shift[axis] = step
                derivative = (self._beam_coupling(family, position+shift, velocity, time) -
                              self._beam_coupling(family, position-shift, velocity, time))/(2*step)
                output[bi, axis] = -hbar*derivative
        return output

    def per_beam_force(self, position, velocity, rho=None, time=0.0):
        """Expectation of each beam's rigorously differentiated interaction."""
        rho = self.steady_state(position, velocity, time) if rho is None else np.asarray(rho)
        operators = self.force_operators(position, velocity, time)
        return np.real(np.einsum("ba,ikab->ik", rho, operators))

    def force(self, position, velocity, rho=None, time=0.0):
        """Return ``(Fx,Fy,Fz)``.

        Supplying ``rho`` gives the instantaneous force.  With no state, a
        stationary common-frequency problem is solved directly.  Otherwise
        the explicit moving-atom master equation is integrated and the final
        quarter of the trajectory is averaged (QUICK: 12 lifetimes, RESEARCH:
        40 lifetimes).  This is intentionally a point solver, not a 3-D-grid
        routine.
        """
        if rho is not None:
            return self.per_beam_force(position, velocity, rho, time).sum(axis=0)
        try:
            state = self.steady_state(position, velocity, time)
            return self.per_beam_force(position, velocity, state, time).sum(axis=0)
        except ValueError:
            lifetimes = 12 if self.mode == "quick" else 40
            samples = 49 if self.mode == "quick" else 241
            times = time + np.linspace(0, lifetimes/self.basis.line.gamma_rad_s, samples)
            density = self.evolve(position, velocity, times)
            start = 3*samples//4
            values = [self.per_beam_force(position, velocity, density[i], times[i]).sum(axis=0)
                      for i in range(start, samples)]
            return np.mean(values, axis=0)
