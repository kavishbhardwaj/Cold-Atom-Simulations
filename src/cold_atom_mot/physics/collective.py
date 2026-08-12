"""Optional Gaussian-cloud mean-field physics for optically thick MOTs.

The validated single-atom force remains unchanged.  This module adds a
continuum approximation; it never represents N atoms one by one.
"""
from dataclasses import dataclass
import numpy as np
from scipy.constants import c, hbar, k as k_B, pi
from scipy.integrate import solve_ivp
from scipy.special import erf


@dataclass(frozen=True)
class GaussianCloud:
    atom_number: float
    sigma_m: np.ndarray
    temperature_k: np.ndarray
    atom_mass_kg: float

    def __post_init__(self):
        sigma=np.broadcast_to(np.asarray(self.sigma_m,float),(3,)).copy()
        temperature=np.broadcast_to(np.asarray(self.temperature_k,float),(3,)).copy()
        if self.atom_number<0 or np.any(sigma<=0) or np.any(temperature<0) or self.atom_mass_kg<=0:
            raise ValueError("N and temperature must be non-negative; widths and mass positive")
        object.__setattr__(self,"sigma_m",sigma);object.__setattr__(self,"temperature_k",temperature)

    @property
    def peak_density_m3(self):
        return self.atom_number/((2*pi)**1.5*np.prod(self.sigma_m))

    def density(self, position):
        position=np.asarray(position,float)
        return self.peak_density_m3*np.exp(-.5*np.sum((position/self.sigma_m)**2,axis=-1))

    @property
    def rms_velocity_m_s(self): return np.sqrt(k_B*self.temperature_k/self.atom_mass_kg)

    def central_column_density_m2(self, axis):
        transverse=np.delete(self.sigma_m,axis)
        return self.atom_number/(2*pi*np.prod(transverse))

    def optical_depth(self, cross_section_m2):
        """Central Beer-Lambert OD along x/y/z."""
        return np.array([cross_section_m2*self.central_column_density_m2(i) for i in range(3)])

    def two_body_integral_m3(self):
        """Exact integral n(r)^2 d3r for the Gaussian profile."""
        return self.atom_number**2/(8*pi**1.5*np.prod(self.sigma_m))


@dataclass(frozen=True)
class MultipleScatteringModel:
    """Walker/Sesko-style Coulomb mean field plus Beer-Lambert shadowing.

    ``laser_cross_section_m2`` and ``reabsorption_cross_section_m2`` must be
    calculated or supplied/cited by the caller.  No arbitrary repulsion
    constant is accepted.  The model assumes an isotropic spherical cloud and
    single-frequency effective cross sections; it is not radiative transport.
    """
    laser_cross_section_m2: float
    reabsorption_cross_section_m2: float
    total_laser_intensity_w_m2: float
    wave_number_rad_m: float
    scattering_rate_s: float

    def __post_init__(self):
        if min(self.laser_cross_section_m2,self.reabsorption_cross_section_m2,
               self.total_laser_intensity_w_m2,self.wave_number_rad_m,self.scattering_rate_s)<0:
            raise ValueError("collective optical parameters must be non-negative")

    @property
    def coulomb_coefficient_n_m2(self):
        """Q=sigma_L sigma_R I/(4 pi c), giving F=Q N_enclosed/r^2."""
        return self.laser_cross_section_m2*self.reabsorption_cross_section_m2*self.total_laser_intensity_w_m2/(4*pi*c)

    def reabsorption_probability(self, cloud):
        return 1-np.exp(-np.mean(cloud.optical_depth(self.reabsorption_cross_section_m2)))

    def radiation_trapping_diffusion(self, cloud):
        """Isotropic extra recoil Dpp per Cartesian axis (single-reabsorption)."""
        return (hbar*self.wave_number_rad_m)**2*self.scattering_rate_s*self.reabsorption_probability(cloud)/3

    def enclosed_fraction(self, radius, sigma):
        x=np.asarray(radius)/(np.sqrt(2)*sigma)
        return erf(x)-np.sqrt(2/pi)*(np.asarray(radius)/sigma)*np.exp(-x*x)

    def repulsive_force(self, radius, cloud):
        radius=np.asarray(radius,float); sigma=float(np.mean(cloud.sigma_m))
        enclosed=cloud.atom_number*self.enclosed_fraction(np.abs(radius),sigma)
        return np.sign(radius)*self.coulomb_coefficient_n_m2*enclosed/np.maximum(radius*radius,1e-30)

    def attenuated_pair_force(self, incident_force_n, cloud, axis=0):
        """Maximum on-axis shadow imbalance from Beer-Lambert attenuation."""
        transmission=np.exp(-cloud.optical_depth(self.laser_cross_section_m2)[axis])
        return incident_force_n*(1-transmission)

    def density_limited_value(self, restoring_coefficient_n_m):
        """Uniform-sphere force-balance density 3*kappa/(4*pi*Q)."""
        if restoring_coefficient_n_m<=0 or self.coulomb_coefficient_n_m2<=0:
            return np.inf
        return 3*restoring_coefficient_n_m/(4*pi*self.coulomb_coefficient_n_m2)

    def equilibrium_sigma(self, atom_number, restoring_coefficient_n_m,
                          temperature_k, atom_mass_kg):
        """Larger of thermal harmonic width and uniform density-limit RMS size."""
        del atom_mass_kg
        thermal=np.sqrt(k_B*temperature_k/restoring_coefficient_n_m)
        density=self.density_limited_value(restoring_coefficient_n_m)
        radius=(3*atom_number/(4*pi*density))**(1/3) if np.isfinite(density) and atom_number else 0
        return max(thermal,radius/np.sqrt(5))


@dataclass(frozen=True)
class CollectiveLoading:
    loading_rate_s: float
    background_loss_s: float
    hot_rb_loss_s: float
    two_body_coefficient_m3_s: float
    temperature_k: float
    restoring_coefficient_n_m: float
    atom_mass_kg: float
    scattering: MultipleScatteringModel | None = None
    fixed_sigma_m: float | None = None
    beta_source: str | None = None

    def __post_init__(self):
        if min(self.loading_rate_s,self.background_loss_s,self.hot_rb_loss_s,
               self.two_body_coefficient_m3_s,self.temperature_k,self.restoring_coefficient_n_m)<0:
            raise ValueError("loading/loss/temperature/restoring inputs must be non-negative")
        if self.two_body_coefficient_m3_s and not self.beta_source:
            raise ValueError("two-body beta requires an experimental/literature source label")

    def cloud(self, population):
        if self.fixed_sigma_m is not None: sigma=self.fixed_sigma_m
        elif self.scattering is not None:
            sigma=self.scattering.equilibrium_sigma(population,self.restoring_coefficient_n_m,
                                                     self.temperature_k,self.atom_mass_kg)
        else: sigma=np.sqrt(k_B*self.temperature_k/self.restoring_coefficient_n_m)
        return GaussianCloud(population,[sigma]*3,[self.temperature_k]*3,self.atom_mass_kg)

    def derivative(self, population):
        cloud=self.cloud(max(float(population),0)); two_body=self.two_body_coefficient_m3_s*cloud.two_body_integral_m3()
        return self.loading_rate_s-(self.background_loss_s+self.hot_rb_loss_s)*population-two_body

    def evolve(self, times, initial_population=0, mode="quick"):
        times=np.asarray(times,float)
        if mode not in ("quick","research") or np.any(np.diff(times)<0) or times[0]<0:
            raise ValueError("mode must be quick/research and times non-decreasing")
        if times[-1]==0: populations=np.full_like(times,initial_population)
        else:
            tol=(1e-7,1e-9) if mode=="quick" else (1e-10,1e-12)
            populations=solve_ivp(lambda _,n:self.derivative(n[0]),(0,times[-1]),[initial_population],t_eval=times,rtol=tol[0],atol=tol[1]).y[0]
        clouds=[self.cloud(n) for n in populations]
        return {"time_s":times,"population":populations,
                "sigma_m":np.array([c.sigma_m for c in clouds]),
                "peak_density_m3":np.array([c.peak_density_m3 for c in clouds])}
