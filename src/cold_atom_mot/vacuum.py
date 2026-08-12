"""Rubidium vapour, thermal surface flux, and configurable MOT loading/loss."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import k as k_B, pi
from scipy.special import gammainc, gammaincinv


@dataclass(frozen=True)
class VaporState:
    """Thermodynamic Rb vapour state, kept separate from non-Rb background gas."""

    temperature_k: float
    rb_partial_pressure_pa: float
    background_gas_pressure_pa: float
    isotope_fractions: dict[str, float]

    def __post_init__(self):
        if self.temperature_k <= 0:
            raise ValueError("temperature must be positive")
        if self.rb_partial_pressure_pa < 0 or self.background_gas_pressure_pa < 0:
            raise ValueError("partial and background pressures must be non-negative")
        if not np.isclose(sum(self.isotope_fractions.values()), 1.0):
            raise ValueError("isotope fractions must sum to one")
        if any(value < 0 for value in self.isotope_fractions.values()):
            raise ValueError("isotope fractions must be non-negative")

    @property
    def rb_number_density_m3(self):
        return number_density(self.rb_partial_pressure_pa, self.temperature_k)

    def isotope_number_density_m3(self, isotope):
        return self.rb_number_density_m3 * self.isotope_fractions[isotope]

def rubidium_vapor_pressure_pa(temperature_k: float) -> float:
    """Natural-Rb vapour pressure using the Alcock-Itkin-Horrigan fit.

    Solid (T<312.46 K): log10(P/Pa)=7.738-4215/T.
    Liquid: log10(P/Pa)=7.193-4040/T. Valid only over the source's fit range.
    """
    if temperature_k <= 0: raise ValueError("temperature must be positive")
    return 10**((7.738-4215/temperature_k) if temperature_k < 312.46 else (7.193-4040/temperature_k))

def number_density(pressure_pa, temperature_k):
    if pressure_pa < 0 or temperature_k <= 0: raise ValueError("pressure must be non-negative and temperature positive")
    return pressure_pa/(k_B*temperature_k)


def one_sided_thermal_flux_m2_s(number_density_m3, temperature_k, mass_kg):
    """Equilibrium flux crossing a plane from one side: n*sqrt(kT/(2*pi*m))."""
    if number_density_m3 < 0 or temperature_k <= 0 or mass_kg <= 0:
        raise ValueError("density must be non-negative; temperature and mass positive")
    return number_density_m3 * np.sqrt(k_B * temperature_k / (2 * pi * mass_kg))


def background_collision_loss_rate_s(
    pressure_pa,
    temperature_k,
    trapped_atom_mass_kg,
    background_particle_mass_kg,
    effective_loss_cross_section_m2,
):
    """Kinetic one-body loss n_bg*sigma_loss*<v_relative>.

    ``sigma_loss`` is an experiment/model input, not a package default. It must
    include the probability that a collision ejects a trapped atom.
    """
    if min(pressure_pa, effective_loss_cross_section_m2) < 0:
        raise ValueError("pressure and loss cross section must be non-negative")
    if min(temperature_k, trapped_atom_mass_kg, background_particle_mass_kg) <= 0:
        raise ValueError("temperature and masses must be positive")
    reduced_mass = trapped_atom_mass_kg * background_particle_mass_kg / (
        trapped_atom_mass_kg + background_particle_mass_kg
    )
    mean_relative_speed = np.sqrt(8 * k_B * temperature_k / (pi * reduced_mass))
    return number_density(pressure_pa, temperature_k) * effective_loss_cross_section_m2 * mean_relative_speed


def sample_flux_speeds(temperature_k, mass_kg, count, rng):
    """Sample the effusive/surface-flux speed law p(v)=2*a²*v³*exp(-a*v²).

    This is the distribution incident on a capture surface. It is weighted by
    one extra factor of speed relative to the bulk Maxwell distribution.
    """
    if temperature_k <= 0 or mass_kg <= 0 or count < 0:
        raise ValueError("temperature/mass must be positive and count non-negative")
    # y=a*v² follows Gamma(shape=2, scale=1), a=m/(2*kT).
    y = rng.gamma(shape=2.0, scale=1.0, size=count)
    return np.sqrt(2 * k_B * temperature_k * y / mass_kg)


def flux_speed_cdf(speed_m_s, temperature_k, mass_kg):
    """CDF of the surface-flux speed distribution."""
    if temperature_k <= 0 or mass_kg <= 0:
        raise ValueError("temperature and mass must be positive")
    y = mass_kg * np.asarray(speed_m_s)**2 / (2 * k_B * temperature_k)
    return gammainc(2.0, y)


def sample_flux_speeds_between(low, high, temperature_k, mass_kg, count, rng):
    """Inverse-CDF samples conditional on low <= v < high."""
    if low < 0 or high <= low:
        raise ValueError("speed interval must satisfy 0 <= low < high")
    cdf_low, cdf_high = flux_speed_cdf([low, high], temperature_k, mass_kg)
    probability = cdf_high - cdf_low
    if probability <= 0:
        raise ValueError("speed interval has negligible floating-point probability")
    y = gammaincinv(2.0, rng.uniform(cdf_low, cdf_high, count))
    return np.sqrt(2 * k_B * temperature_k * y / mass_kg), float(probability)


def sample_spherical_inward_flux(radius_m, temperature_k, mass_kg, count, *, seed):
    """Sample positions and cosine-law inward velocities on a sphere."""
    if radius_m <= 0:
        raise ValueError("capture-surface radius must be positive")
    rng = np.random.default_rng(seed)
    normals = rng.normal(size=(count, 3))
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    positions = radius_m * normals
    # Lambert/cosine incidence: mu=cos(theta) has p(mu)=2*mu.
    mu = np.sqrt(rng.random(count))
    azimuth = 2 * np.pi * rng.random(count)
    reference = np.tile([0.0, 0.0, 1.0], (count, 1))
    near = np.abs(normals[:, 2]) > 0.9
    reference[near] = [0.0, 1.0, 0.0]
    tangent_1 = np.cross(normals, reference)
    tangent_1 /= np.linalg.norm(tangent_1, axis=1)[:, None]
    tangent_2 = np.cross(normals, tangent_1)
    directions = (-mu[:, None] * normals + np.sqrt(1 - mu**2)[:, None] *
                  (np.cos(azimuth)[:, None] * tangent_1 + np.sin(azimuth)[:, None] * tangent_2))
    speeds = sample_flux_speeds(temperature_k, mass_kg, count, rng)
    return positions, speeds[:, None] * directions, speeds

def loading_curve(time_s, loading_rate_s, one_body_loss_s, *, two_body_coefficient=0.0, effective_volume_m3=None):
    """Integrate dN/dt=R-gamma*N-(beta/V)*N²; beta is never invented."""
    from scipy.integrate import solve_ivp
    if min(loading_rate_s,one_body_loss_s,two_body_coefficient)<0: raise ValueError("rates must be non-negative")
    if two_body_coefficient and (effective_volume_m3 is None or effective_volume_m3<=0): raise ValueError("two-body loss requires positive effective volume")
    t=np.asarray(time_s,float)
    if t.ndim != 1 or len(t) == 0 or np.any(np.diff(t) < 0) or t[0] < 0:
        raise ValueError("time must be a non-empty, non-decreasing 1D array")
    if loading_rate_s == 0:
        return np.zeros_like(t)
    if one_body_loss_s == 0 and two_body_coefficient == 0:
        return loading_rate_s * t
    if two_body_coefficient==0 and one_body_loss_s>0: return loading_rate_s/one_body_loss_s*(1-np.exp(-one_body_loss_s*t))
    coefficient=0 if not two_body_coefficient else two_body_coefficient/effective_volume_m3
    return solve_ivp(lambda _,n: loading_rate_s-one_body_loss_s*n-coefficient*n*n,(0,float(t[-1])),[0.],t_eval=t).y[0]
