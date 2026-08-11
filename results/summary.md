# Representative simulation results

These deterministic results use rubidium-87 atoms in an ideal 1064 nm Gaussian
optical dipole trap. They can be regenerated with `python results/generate_results.py`.

| Quantity | Value |
| --- | ---: |
| Atom temperature | 20.0 µK |
| Trap depth | 1.0 mK × k_B |
| Beam waist | 50.0 µm |
| 1D thermal velocity σ | 43.74 mm/s |
| Radial trap frequency ω_r / 2π | 1969.1 Hz |
| Axial trap frequency ω_z / 2π | 9.4 Hz |
| Radial gravitational sag | 0.064 µm |
| Cloud σ after 30 ms time of flight | 1.312 mm |

The calculations use the idealized models in `cold_atom.py`; they omit atom-atom
interactions, photon scattering, trap anharmonicity, and technical noise.
