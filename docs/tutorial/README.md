# Tutorial: From Rubidium Atoms to a MOT Simulation

This documentation is written to teach the physics **as a calculation**, not to list APIs.

> **Start here:** [**Continuous equation-to-result walkthrough**](continuous_walkthrough.md)
>
> It reads from one $^{87}$Rb atom to a MOT, defines every symbol before using it, derives the governing equations, explains why each model/tool was chosen, states the approximations, shows the corresponding calculation result or pedagogical plot, and then explains why the next model is needed.

Two companions make the tutorial easier to audit:

- [**Notation and complete governing-equation inventory**](00_notation_and_equation_inventory.md) — symbol glossary plus the core physical equations that directly determine repository results. It distinguishes easy-to-confuse symbols such as $A_{\rm hfs}$ versus the rate-generator matrix, nuclear spin $I$ versus optical intensity, and damping $\beta_v$ versus two-body loss $\beta_2$.
- [**Equation visual atlas**](equation_visual_atlas.md) — pedagogical plots for the hyperfine equation, Gaussian beam, Doppler damping, MOT restoring force, Lindblad decay, Rabi oscillations, OBE steady state, thermal flux and loading, alongside the actual generated repository results.

> **Scope / provenance.** This repository is independent after-hours work developed from personal scientific interest and kept as a reproducible record and backup. Laboratory control, acquisition, and other lab codes are not kept here.

## Deeper reference chapters

After the continuous walkthrough, these chapters provide more detail without interrupting the main narrative:

1. [Part I — Physical system and apparatus](01_physical_system_and_apparatus.md) — $^{87}$Rb D2 choice, hyperfine/Zeeman basis, fields, coils, six physical beams, polarization and coherence.
2. [Part II — MOT force models, OBEs, and sub-Doppler physics](02_mot_force_models_and_obes.md) — effective MOT force, multilevel rate equations, two-level and 24-state OBEs, PGC and residual magnetic fields.
3. [Part III — Motion, loading, and collective physics](03_motion_loading_and_collective_physics.md) — RK45 motion, recoil, time sequence, vapour flux, capture, loading/loss and collective effects.
4. [Part IV — Cross-cutting interpretation, approximations, and validation](04_results_validation_and_scope.md) — complete pipeline, approximation boundaries and claim language.
5. [Part V — Digital twin, learning, and reproduction](05_digital_twin_learning_and_reproduction.md) — experiment-specific inputs, clean-checkout reproduction, equation-to-code map and reproducibility checklist.

## The modelling principle

The repository intentionally does **not** use one monolithic “exact MOT model.” It uses a hierarchy:

**atomic constants and basis** → **six beams + magnetic field** → **effective scattering force** → **multilevel rate equations** → **coherent OBEs** → **polarization-gradient physics** → **classical motion/recoil** → **vapour capture/loading** → **time sequence** → **collective-cloud effects**.

Each layer exists because it answers a different physical question at a different computational cost. A result is interpreted only at the fidelity of the model that generated it.

## What “all equations” means here

The [equation inventory](00_notation_and_equation_inventory.md) includes the **core governing physical equations that directly determine the scientific outputs**, together with numerical/statistical definitions that materially affect those outputs (for example phase averaging, capture statistics, loading/loss and sequence transients).

It does not promote trivial array reshaping, plotting transforms or ordinary software bookkeeping into “physics equations.” If an implementation equation is scientifically consequential, it belongs in the inventory and is explained in the walkthrough or a reference chapter.

The quick result-only companion remains the [scientific results gallery](../../results/README.md), and the current external-validation status is in [validation.md](../validation.md).
