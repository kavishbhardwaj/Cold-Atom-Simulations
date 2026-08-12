# Tutorial: From Rubidium Atoms to a MOT Simulation

This is the main technical tutorial for the repository. It explains the code as a **continuous physics calculation**, not as an API reference: what is calculated, which equations are solved, why each model was chosen, which approximations were accepted, what each result means, and how to reproduce the calculations.

> **Scope / provenance.** This repository is independent after-hours work developed from personal scientific interest and kept as a reproducible record and backup. Laboratory control, acquisition, and other lab codes are not kept here.

## How to read it

1. [Part I — Physical system and apparatus](01_physical_system_and_apparatus.md)  
   87Rb D2 choice, hyperfine/Zeeman basis, magnetic field, coils, six physical beams, polarization and coherence.
2. [Part II — MOT force models, OBEs, and sub-Doppler physics](02_mot_force_models_and_obes.md)  
   Effective MOT force, multilevel rate equations, two-level and 24-state OBEs, polarization-gradient cooling, residual magnetic fields.
3. [Part III — Motion, loading, and collective physics](03_motion_loading_and_collective_physics.md)  
   Newton/RK45 trajectories, photon recoil, time sequences, vapour flux, capture, loading/loss, multiple scattering and collective MOT effects.
4. [Part IV — Results, validation, approximations, and scope](04_results_validation_and_scope.md)  
   What the main figures mean, the complete simulation flow, every major approximation, validation philosophy, and what is not yet claimed.
5. [Part V — From framework to digital twin: learning and reproduction](05_digital_twin_learning_and_reproduction.md)  
   Experimental inputs needed for a real digital twin, a student learning path, repository map, clean-checkout reproduction, equation-to-code map, and reproducibility checklist.

## The central modelling decision

The repository intentionally does **not** use one monolithic “exact MOT model.” It uses a hierarchy:

**atomic constants and basis** → **six beams + magnetic field** → **effective scattering force** → **multilevel rate equations** → **coherent OBEs** → **polarization-gradient physics** → **classical motion/recoil** → **vapour capture/loading** → **time sequence** → **collective-cloud effects**.

Each layer exists because it answers a different question at a different computational cost. The tutorial repeatedly states the validity boundary before interpreting a result.

## What a reader should be able to reproduce

After completing the tutorial, a reader should be able to:

- reconstruct the 87Rb D2 24-state hyperfine/Zeeman basis and dipole graph;
- derive and reproduce the effective MOT scattering force;
- solve the multilevel population equations with cooling and repump;
- reproduce the two-level OBE/Lindblad benchmark and understand the 24-state moving OBE;
- understand the reduced Sisyphus model and why it cannot yet provide a fully defensible `T(B)`;
- propagate deterministic/stochastic trajectories;
- sample the correct one-sided thermal flux and calculate trajectory-derived loading;
- reproduce the committed validation checks against QuTiP and PyLCP;
- identify which quantities are calculated, sourced, user supplied, or require experimental calibration.

The quick visual companion is the [scientific results gallery](../../results/README.md); current validation status is in [validation.md](../validation.md).
