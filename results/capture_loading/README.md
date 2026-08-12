# Vapour capture and loading results

Generated with `python examples/generate_capture_loading_results.py`. The model
combines a sourced natural-Rb vapour-pressure fit, isotope-resolved equilibrium
surface flux, deterministic MOT capture trajectories, and explicitly calibrated
loss scenarios. It is not an experimentally validated atom-number prediction.

The stored response map resolves capture versus incident speed and impact
parameter. Reservoir temperature controls equilibrium Rb pressure, vapour
temperature controls the incident flux distribution, and background-gas
temperature is reserved for kinetic collision-loss estimates. Wilson intervals
and an adaptive high-speed tail bound retain finite uncertainty when no capture
is observed. See the [visual gallery](../README.md) and
[`capture_loading_reference.npz`](capture_loading_reference.npz).
