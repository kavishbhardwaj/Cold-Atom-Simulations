# Parameter studies

`generate_parameter_studies.py` records every held-fixed constraint in NPZ
metadata. The damping study compares numerical `beta=-dF/dv|0` with the analytic
two-beam derivative, while also plotting scattering. At low saturation beta
rises; saturation and power broadening flatten the velocity slope so beta turns
over even as scattering approaches its ceiling. **More scattering is not more
damping.** Derivative-step refinement is tested.

The waist study holds 10 mW per beam, the YAML detuning and gradient fixed. It
reports peak intensity and a trajectory-derived speed-grid acceptance proxy
under an explicit radius/speed/dwell criterion. Narrow beams supply stronger
local intensity but less interaction distance; wide beams cover more volume but
have lower intensity. This is not an atom-number optimum: vapour flux and loss
calibration are separate.
