# Time-sequenced experimental cycles

`simulation.sequence` represents a laboratory timeline rather than applying one
laser and magnetic configuration forever. A stage carries cooling detuning and
power, repump power, gradient, bias, AOM offsets, polarization purity, phases,
and coherence groups. Each scalar/vector control is a step, linear ramp, or
smoothstep. Zero-duration stages are skipped deterministically at boundaries.

The illustrative 87Rb D2 configuration in
`configs/rb87_d2_reference_sequence.yaml` contains vapor-MOT loading, compressed
MOT, gradient switch-off, field settling, PGC ramp, molasses hold, and release / 
time of flight. Values are representative simulation inputs, **not a literature-
validated optimum**. In particular, the mains-field amplitude defaults to zero:
the simulator does not invent unmeasured laboratory noise.

## Magnetic response

`EddyCurrentResponse` supplies a finite coil response
`G(t)=G0 exp[-(t-t_off)/tau_coil]` and residual field
`B=B_DC+B_eddy exp[-(t-t_off)/tau_eddy]+B_AC sin(2 pi f t+phi)` after switch-off.
Before switch-off it returns the initial gradient and eddy amplitude. Position
dependence uses the Maxwell-consistent `diag(G,G,-2G)` quadrupole. The response
declares itself time dependent, preventing stationary-Liouvillian caching.

`SequencedForce` lets deterministic and photon-event trajectory solvers consume
the time-varying effective force. `snapshot(t)` and `apply_beams(...)` also let
rate-equation and OBE callers construct internal-state solvers at each requested
time without mutating the base configuration.

## Interpretation of generated timing scans

`generate_sequence_results.py` plots the full control timeline and an effective-
model cooling proxy versus field-settling delay, eddy time constant, residual
field, and molasses duration. The proxy uses local linear damping and the
effective model's recoil-event diffusion; it is not a coherent PGC temperature.
It demonstrates the timing/field connection and establishes regression data,
but a quantitative experimental prediction still needs measured coil transfer
functions, beam extinction/rise traces, density-dependent losses, and a
converged multilevel OBE force-noise calculation.
