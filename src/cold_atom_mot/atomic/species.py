"""Single-source rubidium isotope, fine-structure and hyperfine data.

Constants follow Steck's 85Rb/87Rb D-line data, revision 2.3.2 (2021).
Angular frequencies use rad/s; tabulated hyperfine constants use Hz.
"""
from dataclasses import dataclass
from functools import lru_cache
from types import MappingProxyType
import numpy as np
from scipy.constants import hbar, k as boltzmann, physical_constants
from sympy import S
from sympy.physics.wigner import clebsch_gordan, wigner_6j

MU_B = physical_constants["Bohr magneton"][0]


@dataclass(frozen=True)
class AtomicSpecies:
    isotope: str
    mass_kg: float
    abundance: float
    nuclear_spin: float
    nuclear_g_factor: float
    ground_hyperfine_a_hz: float
    source: str


@dataclass(frozen=True)
class FineStructureLine:
    species: AtomicSpecies
    line: str
    wavelength_m: float
    lifetime_s: float
    ground_j: float
    excited_j: float
    ground_g_j: float
    excited_g_j: float
    excited_hyperfine_a_hz: float
    excited_hyperfine_b_hz: float
    saturation_intensity_w_m2: float
    cooling_transition: tuple[int, int] | None
    repump_transition: tuple[int, int] | None
    rate_equation_mot: bool
    multilevel_obe: bool = False
    gray_molasses: bool = False
    @property
    def model_support(self):
        capabilities = ["atomic data", "two-level OBE benchmark"]
        if self.rate_equation_mot: capabilities += ["effective MOT", "multilevel rate-equation MOT"]
        return "; ".join(capabilities)

    @property
    def isotope(self): return self.species.isotope
    @property
    def mass(self): return self.mass_kg
    @property
    def gamma(self): return self.gamma_rad_s
    @property
    def wave_number(self): return self.wave_number_rad_m
    @property
    def wavelength(self): return self.wavelength_m
    @property
    def saturation_intensity(self): return self.saturation_intensity_w_m2
    @property
    def mass_kg(self): return self.species.mass_kg
    @property
    def natural_abundance(self): return self.species.abundance
    @property
    def nuclear_spin(self): return self.species.nuclear_spin
    @property
    def gamma_rad_s(self): return 1 / self.lifetime_s
    @property
    def wave_number_rad_m(self): return 2 * np.pi / self.wavelength_m
    @property
    def recoil_velocity_m_s(self): return hbar * self.wave_number_rad_m / self.mass_kg
    @property
    def recoil_temperature_k(self): return (hbar * self.wave_number_rad_m) ** 2 / (2 * self.mass_kg * boltzmann)
    @property
    def doppler_temperature_k(self): return hbar * self.gamma_rad_s / (2 * boltzmann)
    @property
    def ground_f(self): return tuple(range(int(abs(self.nuclear_spin-self.ground_j)), int(self.nuclear_spin+self.ground_j)+1))
    @property
    def excited_f(self): return tuple(range(int(abs(self.nuclear_spin-self.excited_j)), int(self.nuclear_spin+self.excited_j)+1))

    def hyperfine_energy_hz(self, manifold: str, f: int) -> float:
        j = self.ground_j if manifold == "ground" else self.excited_j
        a = self.species.ground_hyperfine_a_hz if manifold == "ground" else self.excited_hyperfine_a_hz
        b = 0.0 if manifold == "ground" else self.excited_hyperfine_b_hz
        i = self.nuclear_spin; K = f*(f+1)-i*(i+1)-j*(j+1)
        energy = 0.5*a*K
        if b and i >= 1 and j >= 1:
            energy += b * (0.75*K*(K+1)-i*(i+1)*j*(j+1)) / (2*i*(2*i-1)*j*(2*j-1))
        return energy

    def lande_gf(self, manifold: str, f: int) -> float:
        if f == 0: return 0.0
        j = self.ground_j if manifold == "ground" else self.excited_j
        gj = self.ground_g_j if manifold == "ground" else self.excited_g_j
        i = self.nuclear_spin; denominator = 2*f*(f+1)
        return (gj*(f*(f+1)+j*(j+1)-i*(i+1)) +
                self.species.nuclear_g_factor*(f*(f+1)+i*(i+1)-j*(j+1))) / denominator


@dataclass(frozen=True)
class HyperfineState:
    manifold: str; F: int; m: int; g_factor: float; frequency_offset_hz: float


@dataclass(frozen=True)
class DipoleTransition:
    ground_index: int; excited_index: int; q: int; strength: float


@dataclass(frozen=True)
class AtomicBasis:
    line: FineStructureLine
    ground: tuple[HyperfineState, ...]
    excited: tuple[HyperfineState, ...]
    transitions: tuple[DipoleTransition, ...]
    spontaneous_branching: np.ndarray
    @property
    def state_count(self): return len(self.ground)+len(self.excited)


@lru_cache(maxsize=None)
def hyperfine_reduced_strength(i, jg, fg, je, fe):
    """Relative |<Fe||d||Fg>|² including the Wigner-6j factor."""
    if abs(fe-fg) > 1 or (fg == 0 and fe == 0): return 0.0
    six = float(wigner_6j(S(je), S(fe), S(i), S(fg), S(jg), S(1)))
    return (2*fe+1)*(2*jg+1)*six*six


@lru_cache(maxsize=None)
def zeeman_strength(i, jg, fg, mg, je, fe, me):
    q = me-mg
    if q not in (-1, 0, 1): return 0.0
    cg = float(clebsch_gordan(S(fg), S(1), S(fe), S(mg), S(q), S(me)))
    return hyperfine_reduced_strength(i, jg, fg, je, fe)*cg*cg


@lru_cache(maxsize=None)
def build_atomic_basis(isotope: str, line_name: str) -> AtomicBasis:
    line = get_atomic_line(isotope, line_name)
    ground_zero = line.hyperfine_energy_hz("ground", max(line.ground_f))
    excited_zero = line.hyperfine_energy_hz("excited", max(line.excited_f))
    ground = tuple(HyperfineState("ground", f, m, line.lande_gf("ground", f),
                                 line.hyperfine_energy_hz("ground", f)-ground_zero)
                   for f in line.ground_f for m in range(-f, f+1))
    excited = tuple(HyperfineState("excited", f, m, line.lande_gf("excited", f),
                                   line.hyperfine_energy_hz("excited", f)-excited_zero)
                    for f in line.excited_f for m in range(-f, f+1))
    raw = []
    for gi,g in enumerate(ground):
        for ei,e in enumerate(excited):
            value = zeeman_strength(line.nuclear_spin,line.ground_j,g.F,g.m,line.excited_j,e.F,e.m)
            if value > 1e-15: raw.append((gi,ei,e.m-g.m,value))
    maximum = max(value for *_,value in raw)
    transitions = tuple(DipoleTransition(gi,ei,q,value/maximum) for gi,ei,q,value in raw)
    branching = np.zeros((len(excited),len(ground)))
    for t in transitions: branching[t.excited_index,t.ground_index] += t.strength
    branching /= branching.sum(axis=1,keepdims=True)
    return AtomicBasis(line,ground,excited,transitions,branching)


S87=AtomicSpecies("87Rb",1.44316060e-25,0.2783,1.5,-0.0009951414,3.417341305452145e9,"Steck Rb87 rev 2.3.2")
S85=AtomicSpecies("85Rb",1.409993199e-25,0.7217,2.5,-0.0002936400,1.0119108130e9,"Steck Rb85 rev 2.3.2")
ATOMIC_LINES=MappingProxyType({
 ("87Rb","D2"):FineStructureLine(S87,"D2",780.241209686e-9,26.2348e-9,.5,1.5,2.00233113,1.334,84.7185e6,12.4965e6,16.69,(2,3),(1,2),True),
 ("87Rb","D1"):FineStructureLine(S87,"D1",794.978851156e-9,27.679e-9,.5,.5,2.00233113,.666,408.328e6,0,35.77,None,None,False),
 ("85Rb","D2"):FineStructureLine(S85,"D2",780.241368271e-9,26.2348e-9,.5,1.5,2.00233113,1.334,25.0020e6,25.790e6,16.67,(3,4),(2,3),True),
 ("85Rb","D1"):FineStructureLine(S85,"D1",794.9788509e-9,27.679e-9,.5,.5,2.00233113,.666,120.527e6,0,35.75,None,None,False),
})


def get_atomic_line(isotope, line):
    token = isotope.lower()
    normalized = {"rb87": "87Rb", "87rb": "87Rb", "rb85": "85Rb", "85rb": "85Rb"}.get(token, isotope)
    try: return ATOMIC_LINES[(normalized,line.upper())]
    except KeyError as error:
        choices = ", ".join(f"{i} {l}" for i,l in ATOMIC_LINES)
        raise ValueError(f"unknown isotope/line; choose one of: {choices}") from error
