"""Sparse, moving-atom multilevel optical Bloch equations.

Hamiltonians in this module are angular-frequency Hamiltonians (``H / hbar``).
The excited manifold uses one optical reference and each ground hyperfine block
uses its own carrier. Residual same-manifold frequency/Doppler differences stay
explicitly time dependent, without numerically following the removed GHz
cooling--repump beat.
"""
from dataclasses import dataclass, field, replace
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
    ground_relaxation_rate: float = 0.0
    phase_samples: int = 4
    phase_tolerance: float = 0.01
    _phase_resolved: bool = field(default=False, repr=False)
    _collapses: list = field(default=None, init=False, repr=False)
    _dissipator: object = field(default=None, init=False, repr=False)
    last_force_convergence: dict | None = field(default=None, init=False, repr=False)
    last_phase_diagnostics: dict | None = field(default=None, init=False, repr=False)

    def __post_init__(self):
        if self.mode not in ("quick", "research"):
            raise ValueError("mode must be 'quick' or 'research'")
        if self.ground_relaxation_rate < 0:
            raise ValueError("ground_relaxation_rate must be non-negative")
        if self.phase_samples < 1:
            raise ValueError("phase_samples must be positive")
        if self.phase_tolerance <= 0:
            raise ValueError("phase_tolerance must be positive")

    def _coherence_groups(self):
        groups = {}
        for index, family in enumerate(self.beam_families):
            key = family.beam.coherence_group
            # None means a mutually incoherent singleton, not one common group.
            key = ("independent", index) if key is None else ("coherent", key)
            groups.setdefault(key, []).append(index)
        return list(groups.values())

    def _phase_realization(self, sample, total_samples=None):
        """Return a solver realization for deterministic optical phase cycling."""
        groups = self._coherence_groups(); families = list(self.beam_families)
        total_samples = total_samples or self.phase_samples
        for group_index, members in enumerate(groups):
            reference = families[members[0]].beam.phase
            common = 2*np.pi*((sample*(group_index+1)) % total_samples)/total_samples
            for member in members:
                family = families[member]
                # Remove the arbitrary absolute phase, preserve within-group
                # phase differences, and cycle each incoherent group.
                beam = replace(family.beam, phase=family.beam.phase-reference+common)
                families[member] = replace(family, beam=beam)
        return replace(self, beam_families=families, _phase_resolved=True)

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
        # Independent F-block optical rotations make cross-F magnetic elements
        # oscillate at the 6.835-GHz ground splitting.  Apply the corresponding
        # secular approximation while retaining the complete vector Hamiltonian
        # inside each F block.  Its omitted mixing amplitude is O(mu_B B/Delta_hfs).
        for i, left in enumerate(self.basis.ground):
            for j, right in enumerate(self.basis.ground):
                if left.F != right.F:
                    h[i, j] = 0
        h[:ng, :ng] -= 2*np.pi*self._ground_reference_hz*np.eye(ng)
        h[ng:, ng:] = hyperfine_zeeman_hamiltonian(self.basis.line, "excited", b)
        h[ng:, ng:] -= 2*np.pi*self.basis.line.hyperfine_energy_hz(
            "excited", max(self.basis.line.excited_f))*np.eye(n-ng)
        return h

    def _manifold_carriers(self, velocity=(0, 0, 0)):
        """Carrier angular frequency chosen for each addressed ground manifold."""
        carriers = {}; velocity = np.asarray(velocity, float)
        for family in self.beam_families:
            carriers.setdefault(family.ground_f,
                                self._laser_offset(family)-np.dot(family.beam.k_vector, velocity))
        return carriers

    def retained_beat_frequencies(self, velocity=(0, 0, 0)):
        """Residual beat frequencies after the ground-manifold block rotation.

        Cooling/repump's GHz separation is removed exactly by independently
        rotating F=2/F=1.  Only multiple frequencies addressing the *same* F,
        and their independent Doppler shifts, remain explicitly time dependent.
        """
        carriers = self._manifold_carriers(velocity); velocity = np.asarray(velocity, float)
        return np.asarray([self._laser_offset(f)-carriers[f.ground_f]-np.dot(f.beam.k_vector, velocity)
                           for f in self.beam_families])

    def cross_ground_rwa_diagnostics(self, family):
        """Actual-transition bound for the discarded other-ground-F drive.

        Detunings use every generated dipole-allowed transition originating in
        the discarded ground manifold and the real configured laser frequency.
        """
        other = [f for f in self.basis.line.ground_f if f != family.ground_f]
        if not other:
            return {"omega_max_rad_s": 0.0, "delta_min_rad_s": np.inf,
                    "amplitude_ratio": 0.0, "population_bound": 0.0}
        laser = self._laser_offset(family)
        excited_ref = self.basis.line.hyperfine_energy_hz("excited", max(self.basis.line.excited_f))
        detunings = []
        for transition in self.basis.transitions:
            ground = self.basis.ground[transition.ground_index]
            if ground.F not in other:
                continue
            excited = self.basis.excited[transition.excited_index]
            transition_offset = 2*np.pi*((excited.frequency_offset_hz) -
                (ground.frequency_offset_hz))
            detunings.append(abs(laser-transition_offset))
        delta = min(detunings)
        saturation = family.beam.peak_intensity/self.basis.line.saturation_intensity_w_m2
        omega = self.basis.line.gamma_rad_s*np.sqrt(saturation/2)
        ratio = float(omega/delta)
        return {"omega_max_rad_s": float(omega), "delta_min_rad_s": float(delta),
                "amplitude_ratio": ratio, "population_bound": ratio**2}

    def cross_ground_rwa_bound(self, family):
        """Backward-compatible population bound for the ground-manifold RWA."""
        return self.cross_ground_rwa_diagnostics(family)["population_bound"]

    def is_hamiltonian_time_independent(self, velocity=(0, 0, 0), *, tolerance=None):
        """Determine stationarity from declared physical sources, never samples."""
        tolerance = tolerance or 1e-9*self.basis.line.gamma_rad_s
        optical_static = np.all(np.abs(self.retained_beat_frequencies(velocity)) <= tolerance)
        magnetic_static = getattr(self.magnetic_field, "is_time_independent", False)
        return bool(optical_static and magnetic_static)

    def _beam_coupling(self, family, position, velocity, time, phase_shift=0.0):
        """Return one beam's RWA coupling matrix, including its trajectory phase."""
        ng = len(self.basis.ground); matrix = np.zeros((self.basis.state_count,)*2, complex)
        beam = family.beam
        r = np.asarray(position, float) + np.asarray(velocity, float)*time
        fractions = spherical_fractions(beam.polarization, np.array([0., 0., 1.]))
        saturation = float(beam.intensity(r)/self.basis.line.saturation_intensity_w_m2)
        carrier = self._manifold_carriers(velocity)[family.ground_f]
        residual = self._laser_offset(family)-carrier
        phase = np.dot(beam.k_vector, r-beam.origin) - residual*time + beam.phase + phase_shift
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
        for index, state in enumerate(self.basis.ground):
            h[index, index] += self._manifold_carriers(velocity).get(state.F, 0.0)
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

    def _ground_relaxation_operators(self):
        """Optional explicitly requested CP isotropic ground mixing."""
        if not self.ground_relaxation_rate:
            return []
        ng = len(self.basis.ground); rate = self.ground_relaxation_rate/ng
        return [csr_matrix(([np.sqrt(rate)], ([target], [source])),
                           shape=(self.basis.state_count,)*2)
                for source in range(ng) for target in range(ng)]

    def liouvillian(self, position, velocity=(0, 0, 0), time=0.0):
        h = csr_matrix(self.hamiltonian(position, velocity, time)); n = h.shape[0]
        ident = eye(n, format="csr", dtype=complex)
        operator = -1j*(kron(ident, h)-kron(h.T, ident))
        if self._dissipator is None:
            dissipator = csr_matrix((n*n, n*n), dtype=complex)
            for collapse in self.collapse_operators()+self._ground_relaxation_operators():
                cd_c = collapse.getH()@collapse
                dissipator += kron(collapse.conjugate(), collapse)-.5*kron(ident, cd_c)-.5*kron(cd_c.T, ident)
            self._dissipator = dissipator.tocsr()
        operator += self._dissipator
        return operator.tocsr()

    def steady_state_realization(self, position=(0, 0, 0), velocity=(0, 0, 0), time=0.0):
        """Stationary state of one explicitly phase-resolved realization.

        Multi-frequency/multi-direction moving configurations should use
        :meth:`evolve` and time-average :meth:`per_beam_force` instead.
        """
        effective = self.retained_beat_frequencies(velocity)
        if len(effective) and not np.allclose(effective, 0, rtol=0,
                                              atol=1e-9*self.basis.line.gamma_rad_s):
            raise ValueError("no stationary rotating frame: use evolve() for unequal laser/Doppler frequencies")
        if not getattr(self.magnetic_field, "is_time_independent", False):
            raise ValueError("magnetic field is not explicitly time independent: use evolve()")
        n = self.basis.state_count
        h = self.hamiltonian(position, velocity, time)
        ident = eye(n, format="csr", dtype=complex); hs = csr_matrix(h)
        matrix = (-1j*(kron(ident, hs)-kron(hs.T, ident))).tolil()
        for collapse in self.collapse_operators()+self._ground_relaxation_operators():
            cd_c = collapse.getH()@collapse
            matrix += (kron(collapse.conjugate(), collapse)-.5*kron(ident, cd_c)-
                       .5*kron(cd_c.T, ident)).tolil()
        rhs = np.zeros(n*n, complex)
        trace_row = np.zeros(n*n, complex); trace_row[::n+1] = 1
        matrix[-1, :] = trace_row; rhs[-1] = 1
        rho = spsolve(matrix.tocsr(), rhs).reshape((n, n), order="F")
        rho = (rho+rho.conj().T)/2
        trace = np.trace(rho).real
        if not np.isfinite(rho).all() or abs(trace) < 1e-12:
            return self.long_time_state(position, velocity)
        rho /= trace
        if np.linalg.eigvalsh(rho).min() < -1e-8:
            return self.long_time_state(position, velocity)
        return rho

    def _phase_average(self, observable, *, base_samples=None):
        """Refine deterministic phase ensembles and return value plus diagnostics."""
        groups = self._coherence_groups()
        if len(groups) <= 1 or self._phase_resolved:
            return observable(self), {"coherence_groups": len(groups), "phase_samples": 1,
                                      "relative_change": 0.0, "converged": True}
        n = max(base_samples or self.phase_samples, len(groups)+1)
        previous = np.mean([observable(self._phase_realization(i, n)) for i in range(n)], axis=0)
        refinements = 2 if self.mode == "research" else 1
        change = np.inf
        for _ in range(refinements):
            n *= 2
            current = np.mean([observable(self._phase_realization(i, n)) for i in range(n)], axis=0)
            scale = max(np.linalg.norm(current), np.linalg.norm(previous), 1e-30)
            change = float(np.linalg.norm(current-previous)/scale)
            previous = current
        diagnostics = {"coherence_groups": len(groups), "phase_samples": n,
                       "relative_change": change, "converged": change <= self.phase_tolerance}
        if self.mode == "research" and not diagnostics["converged"]:
            raise RuntimeError(f"incoherent phase average not converged: relative change {change:.3g}")
        return previous, diagnostics

    def phase_averaged_steady_state(self, position=(0, 0, 0), velocity=(0, 0, 0), time=0.0):
        """Physical ensemble density matrix for mutually incoherent groups."""
        state, diagnostics = self._phase_average(
            lambda realization: realization.steady_state_realization(position, velocity, time))
        self.last_phase_diagnostics = diagnostics
        return state

    def steady_state(self, position=(0, 0, 0), velocity=(0, 0, 0), time=0.0):
        """Return a phase average, refusing ambiguity for incoherent apparatus."""
        if len(self._coherence_groups()) > 1 and not self._phase_resolved:
            return self.phase_averaged_steady_state(position, velocity, time)
        return self.steady_state_realization(position, velocity, time)

    def evolve(self, position, velocity, times, rho0=None, *, rtol=None, atol=None,
               max_step=np.inf, return_diagnostics=False):
        """Integrate the explicitly time-dependent master equation."""
        times = np.asarray(times, float)
        if times.ndim != 1 or len(times) < 2 or np.any(np.diff(times) <= 0):
            raise ValueError("times must be a strictly increasing 1D array")
        n = self.basis.state_count
        if rho0 is None:
            rho0 = np.zeros((n, n), complex)
            rho0[np.arange(len(self.basis.ground)), np.arange(len(self.basis.ground))] = 1/len(self.basis.ground)
        tolerances = (2e-5, 2e-8) if self.mode == "quick" else (2e-7, 2e-10)
        # Cache only when all known physical sources explicitly declare
        # stationarity. Equal endpoints do not establish this for periodic H(t).
        stationary = self.is_hamiltonian_time_independent(velocity)
        fixed = self.liouvillian(position, velocity, times[0]) if stationary else None
        rhs = ((lambda _t, y: fixed@y) if fixed is not None else
               (lambda t, y: self.liouvillian(position, velocity, t)@y))
        solution = solve_ivp(rhs,
                             (times[0], times[-1]), np.asarray(rho0).reshape(-1, order="F"),
                             t_eval=None if return_diagnostics else times,
                             dense_output=return_diagnostics,
                             rtol=rtol or tolerances[0], atol=atol or tolerances[1],
                             max_step=max_step, method="DOP853")
        if not solution.success:
            raise RuntimeError(solution.message)
        evaluated = solution.sol(times).T if return_diagnostics else solution.y.T
        density = evaluated.reshape((-1, n, n), order="F")
        density = np.asarray([(r+r.conj().T)/(2*np.trace(r).real) for r in density])
        if return_diagnostics:
            steps = np.diff(solution.t)
            diagnostics = {"internal_steps": len(steps), "min_step_s": float(steps.min()),
                           "max_step_s": float(steps.max()), "median_step_s": float(np.median(steps)),
                           "retained_max_beat_rad_s": float(np.max(np.abs(self.retained_beat_frequencies(velocity)), initial=0)),
                           "fixed_liouvillian": stationary}
            return density, diagnostics
        return density

    def long_time_state(self, position, velocity, *, lifetimes=80, samples=161, rho0=None):
        """Physical asymptotic state selected by an initial density matrix."""
        times = np.linspace(0, lifetimes/self.basis.line.gamma_rad_s, samples)
        return self.evolve(position, velocity, times, rho0)[-1]

    def force_operators(self, position, velocity=(0, 0, 0), time=0.0, *, finite_difference=False, step=None):
        """Per-beam optical force operators ``-gradient(H_int)`` in newtons.

        The derivative includes both travelling-wave phase and the Gaussian
        envelope.  A symmetric spatial derivative is used so displaced beams
        retain their dipole-force contribution as well as radiation pressure.
        """
        position = np.asarray(position, float)
        step = step or max(self.basis.line.wavelength_m*1e-4, 1e-11)
        output = np.empty((len(self.beam_families), 3, self.basis.state_count, self.basis.state_count), complex)
        for bi, family in enumerate(self.beam_families):
            for axis in range(3):
                if finite_difference:
                    shift = np.zeros(3); shift[axis] = step
                    derivative = (self._beam_coupling(family, position+shift, velocity, time) -
                                  self._beam_coupling(family, position-shift, velocity, time))/(2*step)
                else:
                    coupling = self._beam_coupling(family, position, velocity, time)
                    r = position+np.asarray(velocity, float)*time-family.beam.origin
                    transverse = r-np.dot(r, family.beam.direction)*family.beam.direction
                    factor = 1j*family.beam.k_vector[axis]-2*transverse[axis]/family.beam.waist**2
                    derivative = np.zeros_like(coupling)
                    ng = len(self.basis.ground)
                    derivative[ng:, :ng] = coupling[ng:, :ng]*factor
                    derivative[:ng, ng:] = derivative[ng:, :ng].conj().T
                output[bi, axis] = -hbar*derivative
        return output

    def per_beam_force(self, position, velocity, rho=None, time=0.0):
        """Expectation of each beam's rigorously differentiated interaction."""
        rho = self.steady_state(position, velocity, time) if rho is None else np.asarray(rho)
        operators = self.force_operators(position, velocity, time)
        return np.real(np.einsum("ba,ikab->ik", rho, operators))

    def _window_average(self, values, *, lifetimes=None, samples=None):
        """Average the last window and enforce the RESEARCH convergence policy."""
        values = np.asarray(values); split = len(values)//2
        early, late = values[:split].mean(axis=0), values[split:].mean(axis=0)
        scale = max(np.linalg.norm(late), np.max(np.linalg.norm(values, axis=1)), 1e-30)
        metric = float(np.linalg.norm(late-early)/scale)
        self.last_force_convergence = {"relative_window_change": metric,
                                       "lifetimes": lifetimes, "samples": samples,
                                       "averaging_fraction": 0.25}
        if self.mode == "research" and metric > 0.01:
            raise RuntimeError(f"time-averaged force not converged: relative window change {metric:.3g}")
        return late

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
        groups = self._coherence_groups()
        if len(groups) > 1 and not self._phase_resolved:
            value, diagnostics = self._phase_average(
                lambda realization: realization.force(position, velocity, time=time))
            self.last_phase_diagnostics = diagnostics
            return value
        try:
            state = self.steady_state(position, velocity, time)
            return self.per_beam_force(position, velocity, state, time).sum(axis=0)
        except ValueError:
            lifetimes = 12 if self.mode == "quick" else 40
            samples = 49 if self.mode == "quick" else 241
            times = time + np.linspace(0, lifetimes/self.basis.line.gamma_rad_s, samples)
            density = self.evolve(position, velocity, times)
            middle = samples//2
            values = np.asarray([self.per_beam_force(position, velocity, density[i], times[i]).sum(axis=0)
                                 for i in range(middle, samples)])
            return self._window_average(values, lifetimes=lifetimes, samples=samples)
