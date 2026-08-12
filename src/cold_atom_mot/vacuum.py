"""Rubidium vapour density and configurable MOT loading/loss models."""
import numpy as np
from scipy.constants import k as k_B

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

def loading_curve(time_s, loading_rate_s, one_body_loss_s, *, two_body_coefficient=0.0, effective_volume_m3=None):
    """Integrate dN/dt=R-gamma*N-(beta/V)*N²; beta is never invented."""
    from scipy.integrate import solve_ivp
    if min(loading_rate_s,one_body_loss_s,two_body_coefficient)<0: raise ValueError("rates must be non-negative")
    if two_body_coefficient and (effective_volume_m3 is None or effective_volume_m3<=0): raise ValueError("two-body loss requires positive effective volume")
    t=np.asarray(time_s,float)
    if two_body_coefficient==0 and one_body_loss_s>0: return loading_rate_s/one_body_loss_s*(1-np.exp(-one_body_loss_s*t))
    coefficient=0 if not two_body_coefficient else two_body_coefficient/effective_volume_m3
    return solve_ivp(lambda _,n: loading_rate_s-one_body_loss_s*n-coefficient*n*n,(0,float(t[-1])),[0.],t_eval=t).y[0]
