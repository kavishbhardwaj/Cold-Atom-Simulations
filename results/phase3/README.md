# Phase-3 reduced optical-Bloch results

Generated with `python examples/generate_phase3_results.py` using package 0.3.1.
The OBE figures use the effective 87Rb D2 stretched transition. The conditioned
engineering scans use the documented Level-A six-beam force so that the effect
of held-fixed power/intensity is isolated transparently.

This is not a full 24-state, six-beam coherent OBE. `phase3_reference.npz` stores
all plotted arrays, maximum analytical comparison error, model fidelity, basis,
collapse operator and units. Every PNG has a same-figure SVG counterpart.

Phase 3 also compares registered 85Rb/87Rb D1/D2 atomic scales without claiming
that all four systems have full MOT solver support, and adds a configurable
pure-dephasing Lindblad channel for homogeneous coherence broadening.
