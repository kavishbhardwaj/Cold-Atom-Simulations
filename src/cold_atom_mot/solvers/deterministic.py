"""Adaptive deterministic mean-force trajectory integration."""

from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Trajectory:
    time: np.ndarray
    position: np.ndarray
    velocity: np.ndarray


def integrate_trajectory(force_model, position, velocity, duration, *, max_step=1e-5,
                         rtol=1e-7, atol=1e-10, sample_step=None) -> Trajectory:
    """Integrate dr/dt=v and m dv/dt=F using adaptive RK45."""
    if duration <= 0 or max_step <= 0:
        raise ValueError("duration and max_step must be positive")
    initial = np.concatenate([np.asarray(position, dtype=float), np.asarray(velocity, dtype=float)])

    def derivative(time, state):
        return np.concatenate([state[3:], force_model.force(state[:3], state[3:], time=time) / force_model.atom.mass])

    if sample_step is not None and sample_step <= 0:
        raise ValueError("sample_step must be positive")
    t_eval = None
    if sample_step is not None:
        t_eval = np.r_[np.arange(0, duration, sample_step), duration]
    solution = solve_ivp(derivative, (0.0, duration), initial, max_step=max_step,
                         rtol=rtol, atol=atol, t_eval=t_eval)
    return Trajectory(solution.t, solution.y[:3].T, solution.y[3:].T)
