"""Rubidium vapour, thermal surface flux, and configurable MOT loading/loss."""
from dataclasses import dataclass
import numpy as np
from scipy.constants import k as k_B, pi
from scipy.special import gammainc, gammaincinv


@dataclass(frozen=True)
class VaporState:
    """Thermodynamic Rb vapour state, kept separate from non-Rb background gas."""

    rb_reservoir_temperature_k: float
    vapor_temperature_k: float
    background_temperature_k: float
    rb_partial_pressure_pa: float
    background_gas_pressure_pa: float
    isotope_fractions: dict[str, float]

    def __post_init__(self):
        if min(self.rb_reservoir_temperature_k, self.vapor_temperature_k,
               self.background_temperature_k) <= 0:
            raise ValueError("reservoir, vapor, and background temperatures must be positive")
        if self.rb_partial_pressure_pa < 0 or self.background_gas_pressure_pa < 0:
            raise ValueError("partial and background pressures must be non-negative")
        if not np.isclose(sum(self.isotope_fractions.values()), 1.0):
            raise ValueError("isotope fractions must sum to one")
        if any(value < 0 for value in self.isotope_fractions.values()):
            raise ValueError("isotope fractions must be non-negative")

    @property
    def rb_number_density_m3(self):
        return number_density(self.rb_partial_pressure_pa, self.vapor_temperature_k)

    def isotope_number_density_m3(self, isotope):
        return self.rb_number_density_m3 * self.isotope_fractions[isotope]

def rubidium_vapor_pressure_pa(temperature_k: float, *, allow_extrapolation=False) -> float:
    """Natural-Rb vapour pressure using the Alcock-Itkin-Horrigan fit.

    Solid, 298.15 <= T < 312.46 K: log10(P/Pa)=7.738-4215/T.
    Liquid, 312.46 <= T <= 550 K: log10(P/Pa)=7.193-4040/T.
    Pressure is returned directly in pascals; these are the source coefficients
    after expressing pressure in Pa (not torr or atmospheres).
    """
    if not 298.15 <= temperature_k <= 550.0:
        if not allow_extrapolation:
            raise ValueError("Alcock Rb vapour-pressure fit requires 298.15 <= T <= 550 K")
        import warnings
        warnings.warn("extrapolating Alcock Rb vapour-pressure fit outside 298.15-550 K", RuntimeWarning)
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


def wilson_interval(successes, trials, confidence=0.95):
    """Two-sided Wilson score interval, including finite bounds at k=0,n."""
    from scipy.stats import norm
    if trials <= 0 or not 0 <= successes <= trials:
        raise ValueError("require trials>0 and 0<=successes<=trials")
    z = norm.ppf(0.5 + confidence / 2)
    p = successes / trials
    denominator = 1 + z*z/trials
    centre = (p + z*z/(2*trials))/denominator
    half = z*np.sqrt(p*(1-p)/trials + z*z/(4*trials*trials))/denominator
    low=max(0.0,centre-half)
    if successes == 0: low=0.0
    high=min(1.0,centre+half)
    if successes == trials: high=1.0
    return low,high


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


def flux_speed_pdf(speed_m_s, temperature_k, mass_kg):
    speed=np.asarray(speed_m_s,float)
    a=mass_kg/(2*k_B*temperature_k)
    return 2*a*a*speed**3*np.exp(-a*speed**2)


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

def gaussian_two_body_effective_volume(sigma_x_m, sigma_y_m, sigma_z_m):
    """V2=N²/integral(n²dV)=8*pi^(3/2)*sigma_x*sigma_y*sigma_z."""
    if min(sigma_x_m, sigma_y_m, sigma_z_m) <= 0:
        raise ValueError("Gaussian cloud widths must be positive")
    return 8 * pi**1.5 * sigma_x_m * sigma_y_m * sigma_z_m


def steady_state_population(loading_rate_s, one_body_loss_s, *,
                            two_body_coefficient=0.0, effective_volume_m3=None):
    if min(loading_rate_s, one_body_loss_s, two_body_coefficient) < 0:
        raise ValueError("rates must be non-negative")
    if loading_rate_s == 0:
        return 0.0
    if two_body_coefficient == 0:
        return np.inf if one_body_loss_s == 0 and loading_rate_s > 0 else (
            0.0 if loading_rate_s == 0 else loading_rate_s/one_body_loss_s)
    if effective_volume_m3 is None or effective_volume_m3 <= 0:
        raise ValueError("two-body loss requires positive effective volume")
    coefficient = two_body_coefficient/effective_volume_m3
    discriminant=np.sqrt(one_body_loss_s**2+4*coefficient*loading_rate_s)
    return 2*loading_rate_s/(one_body_loss_s+discriminant)


def loading_curve(time_s, loading_rate_s, one_body_loss_s, *, initial_population=0.0,
                  two_body_coefficient=0.0, effective_volume_m3=None):
    """Integrate dN/dt=R-gamma*N-(beta/V)*N²; beta is never invented."""
    from scipy.integrate import solve_ivp
    if min(loading_rate_s,one_body_loss_s,two_body_coefficient,initial_population)<0: raise ValueError("rates and N0 must be non-negative")
    if two_body_coefficient and (effective_volume_m3 is None or effective_volume_m3<=0): raise ValueError("two-body loss requires positive effective volume")
    t=np.asarray(time_s,float)
    if t.ndim != 1 or len(t) == 0 or np.any(np.diff(t) < 0) or t[0] < 0:
        raise ValueError("time must be a non-empty, non-decreasing 1D array")
    if t[-1] == 0:
        return np.full_like(t, initial_population)
    if loading_rate_s == 0:
        if two_body_coefficient == 0:
            return initial_population*np.exp(-one_body_loss_s*t)
    if one_body_loss_s == 0 and two_body_coefficient == 0:
        return initial_population + loading_rate_s * t
    if two_body_coefficient==0 and one_body_loss_s>0:
        steady=loading_rate_s/one_body_loss_s
        return steady+(initial_population-steady)*np.exp(-one_body_loss_s*t)
    coefficient=0 if not two_body_coefficient else two_body_coefficient/effective_volume_m3
    return solve_ivp(lambda _,n: loading_rate_s-one_body_loss_s*n-coefficient*n*n,
                     (0,float(t[-1])),[initial_population],t_eval=t,
                     rtol=1e-9,atol=1e-10).y[0]
