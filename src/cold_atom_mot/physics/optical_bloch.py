"""Two-level OBE reduced two-level optical Bloch equations.

This backend resolves one explicitly selected |g> <-> |e> transition.  It is a
coherence-capable validation and teaching backend, not a 24-state six-beam OBE.
"""

from dataclasses import dataclass
import numpy as np
from scipy.constants import hbar
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class TwoLevelOBE:
    gamma: float
    detuning: float
    rabi_frequency: complex
    dephasing_rate: float = 0.0

    def __post_init__(self) -> None:
        if self.gamma <= 0:
            raise ValueError("gamma must be positive")
        if self.dephasing_rate < 0:
            raise ValueError("dephasing_rate must be non-negative")

    @classmethod
    def from_saturation(
        cls,
        gamma: float,
        detuning: float,
        saturation: float,
        *,
        dephasing_rate: float = 0.0,
    ) -> "TwoLevelOBE":
        if saturation < 0:
            raise ValueError("saturation must be non-negative")
        # s = 2 |Omega|^2 / Gamma^2 for the stated two-level convention.
        return cls(
            gamma,
            detuning,
            gamma * np.sqrt(saturation / 2),
            dephasing_rate,
        )

    @property
    def hamiltonian_over_hbar(self) -> np.ndarray:
        omega = self.rabi_frequency
        return np.array([[0.0, np.conjugate(omega) / 2], [omega / 2, -self.detuning]], dtype=complex)

    def derivative(self, density: np.ndarray) -> np.ndarray:
        """Return dρ/dt=-i[H/ℏ,ρ]+D[sqrt(Γ)|g><e|]ρ."""
        rho = np.asarray(density, dtype=complex).reshape(2, 2)
        hamiltonian = self.hamiltonian_over_hbar
        coherent = -1j * (hamiltonian @ rho - rho @ hamiltonian)
        collapse = np.array([[0.0, np.sqrt(self.gamma)], [0.0, 0.0]], dtype=complex)
        cdag_c = collapse.conj().T @ collapse
        dissipative = collapse @ rho @ collapse.conj().T - 0.5 * (cdag_c @ rho + rho @ cdag_c)
        if self.dephasing_rate:
            # C_phi=sqrt(gamma_phi/2)*sigma_z makes rho_ge decay at gamma_phi.
            sigma_z = np.diag([1.0, -1.0]).astype(complex)
            collapse_phi = np.sqrt(self.dephasing_rate / 2) * sigma_z
            dissipative += (
                collapse_phi @ rho @ collapse_phi.conj().T
                - 0.5 * (collapse_phi.conj().T @ collapse_phi @ rho + rho @ collapse_phi.conj().T @ collapse_phi)
            )
        return coherent + dissipative

    def liouvillian(self) -> np.ndarray:
        """Construct the 4×4 linear superoperator in row-major vectorization."""
        basis = []
        for index in range(4):
            matrix = np.zeros((2, 2), dtype=complex)
            matrix.flat[index] = 1.0
            basis.append(self.derivative(matrix).ravel())
        return np.column_stack(basis)

    def steady_state(self) -> np.ndarray:
        matrix = self.liouvillian().copy()
        rhs = np.zeros(4, dtype=complex)
        matrix[-1] = np.array([1, 0, 0, 1], dtype=complex)
        rhs[-1] = 1
        rho = np.linalg.solve(matrix, rhs).reshape(2, 2)
        return 0.5 * (rho + rho.conj().T)

    def evolve(self, initial_density: np.ndarray, duration: float, *, rtol: float = 1e-9, atol: float = 1e-11, max_step: float | None = None):
        if duration <= 0:
            raise ValueError("duration must be positive")
        initial = np.asarray(initial_density, dtype=complex).reshape(4)

        def rhs(_time, vector):
            rho = vector[:4] + 1j * vector[4:]
            derivative = self.derivative(rho.reshape(2, 2)).ravel()
            return np.concatenate([derivative.real, derivative.imag])

        packed = np.concatenate([initial.real, initial.imag])
        options = {} if max_step is None else {"max_step": max_step}
        solution = solve_ivp(rhs, (0, duration), packed, rtol=rtol, atol=atol, **options)
        density = (solution.y[:4] + 1j * solution.y[4:]).T.reshape(-1, 2, 2)
        return solution.t, density

    def analytic_excited_population(self) -> float:
        omega_squared = abs(self.rabi_frequency) ** 2
        transverse_decay = self.gamma / 2 + self.dephasing_rate
        return (
            omega_squared * transverse_decay
            / (
                2 * self.gamma * (self.detuning**2 + transverse_decay**2)
                + 2 * omega_squared * transverse_decay
            )
        )

    def scattering_force(self, wave_vector: np.ndarray) -> np.ndarray:
        """Single travelling-wave force ℏ k Γ ρee at steady state."""
        return hbar * np.asarray(wave_vector, dtype=float) * self.gamma * self.steady_state()[1, 1].real
