# Laboratory magnetic apparatus

The magnetic package now contains physical circular Helmholtz pairs for three
independent bias axes as well as the segmented anti-Helmholtz MOT pair. Each
bias pair permits unequal radii, separation error, current and turns imbalance,
lateral displacement, and vector tilt. Numerical derivatives provide the field
Jacobian and rank-3 curvature tensor away from conductors.

`ThreeAxisBiasCoils` separates the Biot-Savart geometry from an empirical
calibration `B=M I+B_offset`. `compensation_currents` solves the full, possibly
nonorthogonal 3x3 system by least squares, so measured cross-axis response is
retained. The residual must be checked with `calibrated_field`; singular or
poorly conditioned calibrations cannot guarantee cancellation.

Backgrounds support a uniform Earth/component field, linear stray gradient and
any number of sinusoidal harmonics (including 50/60 Hz). `SwitchingTransientField`
adds an L/R exponential and multiple vector eddy exponentials, or interpolates
a user-supplied measured three-axis waveform. Amplitudes are never invented by
the model.

The reference generator uses a clearly synthetic nonorthogonal calibration and
configured background. It produces Bx/By/Bz/|B| maps, a Maxwell divergence
check, compensation currents, anti-Helmholtz zero displacement, and a table of
controlled imperfections. MOT-centre shifts use the effective linear gradient;
the PGC column intentionally does not report temperature. Static vector-OBE
coherence is useful as an internal-state diagnostic, but a “10% PGC degradation”
or allowed current imbalance requires converged moving-OBE friction and matched
force-noise diffusion. The current code does not scientifically support that
tolerance, and the table says so rather than substituting a scalar detuning.
