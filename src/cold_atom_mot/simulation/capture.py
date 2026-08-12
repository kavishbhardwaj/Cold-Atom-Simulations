"""Trajectory-derived capture metrics with an explicit numerical criterion."""
from dataclasses import dataclass
import numpy as np
from ..solvers.deterministic import integrate_trajectory

@dataclass(frozen=True)
class CaptureCriterion:
    radius_m: float
    speed_m_s: float
    duration_s: float
    dwell_s: float

def evaluate_capture(force_model, positions, velocities, criterion, *, max_step):
    captured=[]; capture_time=[]
    for r,v in zip(np.asarray(positions),np.asarray(velocities)):
        tr=integrate_trajectory(force_model,r,v,criterion.duration_s,max_step=max_step)
        inside=(np.linalg.norm(tr.position,axis=1)<=criterion.radius_m)&(np.linalg.norm(tr.velocity,axis=1)<=criterion.speed_m_s)
        required=max(1,int(np.ceil(criterion.dwell_s/np.median(np.diff(tr.time)))))
        hits=np.convolve(inside.astype(int),np.ones(required,dtype=int),mode="valid") if len(inside)>=required else []
        ok=bool(len(hits) and np.any(hits==required)); captured.append(ok)
        capture_time.append(float(tr.time[np.argmax(hits==required)]) if ok else np.nan)
    return np.asarray(captured),np.asarray(capture_time)
