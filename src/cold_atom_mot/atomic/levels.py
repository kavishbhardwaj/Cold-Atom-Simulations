"""Explicit 87Rb D2 hyperfine/Zeeman basis and electric-dipole graph."""

from dataclasses import dataclass
import numpy as np
from .angular_momentum import clebsch_gordan
from .rb87 import Rb87D2


@dataclass(frozen=True)
class HyperfineState:
    manifold: str
    F: int
    m: int
    g_factor: float
    frequency_offset_hz: float


@dataclass(frozen=True)
class DipoleTransition:
    ground_index: int
    excited_index: int
    q: int
    strength: float


@dataclass(frozen=True)
class Rb87D2Basis:
    ground: tuple[HyperfineState, ...]
    excited: tuple[HyperfineState, ...]
    transitions: tuple[DipoleTransition, ...]
    spontaneous_branching: np.ndarray

    @property
    def state_count(self) -> int:
        return len(self.ground) + len(self.excited)


def build_rb87_d2_basis(atom: Rb87D2 | None = None) -> Rb87D2Basis:
    """Construct all 24 Zeeman states and allowed q=-1,0,+1 couplings.

    Hyperfine strengths are normalized to the Steck line-strength convention,
    then globally scaled so the stretched F=2,m=2 -> F'=3,m'=3 cycling
    transition has unit strength, matching the configured saturation intensity.
    Spontaneous branching is obtained from the same squared matrix elements and
    normalized separately for every excited Zeeman state.
    """
    atom = Rb87D2() if atom is None else atom
    ground = []
    for F, g_factor, offset in ((1, atom.ground_g_f[0], -atom.ground_hyperfine_splitting_hz), (2, atom.ground_g_f[1], 0.0)):
        ground.extend(HyperfineState("ground", F, m, g_factor, offset) for m in range(-F, F + 1))
    excited = []
    for Fp, (g_factor, offset) in enumerate(zip(atom.excited_g_f, atom.excited_hyperfine_offsets_hz)):
        excited.extend(HyperfineState("excited", Fp, m, g_factor, offset) for m in range(-Fp, Fp + 1))

    line_strength = {
        (1, 0): atom.f1_strengths_to_fprime_0_1_2[0],
        (1, 1): atom.f1_strengths_to_fprime_0_1_2[1],
        (1, 2): atom.f1_strengths_to_fprime_0_1_2[2],
        (2, 1): atom.f2_strengths_to_fprime_1_2_3[0],
        (2, 2): atom.f2_strengths_to_fprime_1_2_3[1],
        (2, 3): atom.f2_strengths_to_fprime_1_2_3[2],
    }
    raw = []
    for gi, g in enumerate(ground):
        for ei, e in enumerate(excited):
            if (g.F, e.F) not in line_strength:
                continue
            q = e.m - g.m
            if q not in (-1, 0, 1):
                continue
            cg2 = clebsch_gordan(g.F, g.m, 1, q, e.F, e.m) ** 2
            if cg2 > 1e-15:
                raw.append([gi, ei, q, line_strength[(g.F, e.F)] * cg2])

    # Correct each hyperfine line so its m-averaged total equals Steck's S_FF'.
    for F, Fp in line_strength:
        indices = [i for i, (gi, ei, _, _) in enumerate(raw) if ground[gi].F == F and excited[ei].F == Fp]
        current_average = sum(raw[i][3] for i in indices) / (2 * F + 1)
        correction = line_strength[(F, Fp)] / current_average
        for i in indices:
            raw[i][3] *= correction

    stretched = next(value for gi, ei, q, value in raw if ground[gi].F == 2 and ground[gi].m == 2 and excited[ei].F == 3 and excited[ei].m == 3 and q == 1)
    transitions = tuple(DipoleTransition(gi, ei, q, value / stretched) for gi, ei, q, value in raw)

    branching = np.zeros((len(excited), len(ground)))
    for transition in transitions:
        branching[transition.excited_index, transition.ground_index] += transition.strength
    branching /= branching.sum(axis=1, keepdims=True)
    return Rb87D2Basis(tuple(ground), tuple(excited), transitions, branching)
