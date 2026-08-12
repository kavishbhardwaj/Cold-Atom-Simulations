"""Backward-compatible 87Rb D2 view of the unified atomic registry."""
from dataclasses import dataclass
from .species import get_atomic_line

@dataclass(frozen=True)
class Rb87D2:
    """Deprecated compatibility adapter; new code uses ``FineStructureLine``."""
    @property
    def _line(self): return get_atomic_line("87Rb","D2")
    @property
    def mass(self): return self._line.mass_kg
    @property
    def wavelength(self): return self._line.wavelength_m
    @property
    def lifetime(self): return self._line.lifetime_s
    @property
    def saturation_intensity(self): return self._line.saturation_intensity_w_m2
    @property
    def gamma(self): return self._line.gamma_rad_s
    @property
    def wave_number(self): return self._line.wave_number_rad_m
    @property
    def recoil_velocity(self): return self._line.recoil_velocity_m_s
    @property
    def recoil_temperature(self): return self._line.recoil_temperature_k
    @property
    def doppler_temperature(self): return self._line.doppler_temperature_k
    ground_hyperfine_splitting_hz=6.83468261090429e9
    excited_hyperfine_offsets_hz=(-495.815e6,-423.600e6,-266.650e6,0.0)
    ground_g_f=(-0.501827,0.499836)
    excited_g_f=(0.0,0.666,0.666,0.667)
    f2_strengths_to_fprime_1_2_3=(0.05,0.25,0.70)
    f1_strengths_to_fprime_0_1_2=(1/6,5/12,5/12)
    def validate(self): return None
