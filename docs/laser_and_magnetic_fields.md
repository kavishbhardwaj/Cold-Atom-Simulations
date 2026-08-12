# Laser and magnetic fields

Gaussian beams record direction, origin, power, waist, wavelength, detuning,
frequency offset, phase, helicity, polarization purity, linewidth and coherence
group. The Rayleigh range and waist evolution are available; the MOT force uses
the collimated local-intensity approximation. Fields interfere only within the
same coherence group. `None` groups are mutually incoherent.

The effective Zeeman force projects circular photon angular momentum onto the
actual local magnetic-field direction; it no longer classifies a beam by its
nearest Cartesian axis. The rate model uses the local B direction as a
quantization axis, with a documented fixed-axis convention near zero. The sparse
OBE currently uses |B| diagonal shifts rather than a full transverse vector
Zeeman Hamiltonian, so arbitrary transverse coherent predictions remain a
limitation.

Ideal rotated quadrupoles, uniform/gradient/AC residual fields and independently
positioned/tilted segmented circular coils remain available.
