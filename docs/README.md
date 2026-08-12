# Documentation map

Use this page as a short guide to the repository. The code is organized as a hierarchy of models rather than one monolithic “truth model”; each layer trades speed for physical detail.

## Start here

- **[Continuous equation-to-result tutorial](tutorial/continuous_walkthrough.md)** — the recommended first read. It defines symbols before use and follows the actual scientific workflow: question → equation → modelling/tool decision → approximation → result/plot → interpretation → next model.
- **[Notation and governing-equation inventory](tutorial/00_notation_and_equation_inventory.md)** — glossary and completeness audit of the core physical equations that determine repository results.
- **[Equation visual atlas](tutorial/equation_visual_atlas.md)** — pedagogical equation plots plus the corresponding generated simulation results.
- [Tutorial index and deeper reference chapters](tutorial/README.md) — chaptered extensions of the main walkthrough.
- [Model hierarchy](model_hierarchy.md) — which solver to use, and what each one neglects.
- [Validation](validation.md) — what is independently verified, internally tested, or not yet validated.
- [Scientific results](../results/README.md) — figures, numerical data, provenance, and fidelity notes.

The tutorial is intentionally **not** an expanded README. The primary walkthrough follows the calculation continuously from atomic constants through hyperfine structure, laser/magnetic apparatus, MOT force, rate equations, Lindblad/OBE physics, PGC, trajectories, vapour loading, experimental timing, collective effects and validation.

## Atomic and cooling physics

- [Atomic systems](atomic_systems.md) — 85Rb/87Rb, D1/D2, hyperfine structure and supported transitions.
- [Cooling physics](cooling_physics.md) — Doppler force, optical pumping, recoil, and polarization-gradient cooling.
- [Numerical methods](numerical_methods.md) — integration, convergence, sampling, and solver choices.

## Laboratory apparatus

- [Six-beam apparatus](six_beam_apparatus.md) — six physical beams, Gaussian propagation, Jones optics, QWP errors, retroreflection, and coherence groups.
- [Magnetic apparatus](magnetic_apparatus.md) — MOT/bias coils, calibration matrices, stray fields, 50/60-Hz components, switch-off, and eddy currents.
- [Experimental sequences](experimental_sequences.md) — time-dependent MOT load, compression, field settling, PGC/molasses, and TOF stages.

## Loading and many-atom extensions

- [Collective MOT physics](collective_mot.md) — optional Gaussian density, two-body loss, shadowing, multiple scattering, and radiation-trapping proxy.
- Vapour/capture/loading physics is derived continuously in the [main tutorial](tutorial/continuous_walkthrough.md) and documented further through the model hierarchy, configuration files, and the [capture/loading results](../results/README.md#vapour-capture-loading-and-loss).

## How to interpret claims

Three words matter throughout the repository:

1. **Implemented** — code exists and is exercised by internal tests.
2. **Externally verified** — a matched result agrees with an independent implementation such as QuTiP or PyLCP.
3. **Experimentally calibrated** — measured apparatus inputs and quantitative agreement with a specified experiment support the prediction.

The repository deliberately keeps these categories separate. In particular, the 24-state moving 87Rb D2 OBE is implemented and internally tested, while its full multilevel force and quantitative PGC temperature are still external-validation targets.
