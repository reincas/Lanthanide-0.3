# Lanthanide-0.3

This is the Python software which I developed and used for my PhD work and some time after finishing my doctoral
dissertation in 2001 [1].
The software provides the calculation of all matrix elements of angular spherical tensor operators for the electron
configurations of lanthanide ions.
Included is also code to fit radial integrals to measured energy level spectra,
as well as the code for Judd-Ofelt fits,
and the calculation of all radiative electric and magnetic dipole transitions inside the energy level spectrum
of lanthanide ions.

## Important Note

This software is outdated and intended for reference only.
If you really want to execute it, you can use a virtual machine running SuSE Linux 9.3.
This is how I restarted it in 2025 to compare it with the results of its successor package [Lanthanide](https://github.com/reincas/Lanthanide)
which you should use instead.
Important: Place a symbolic link to the folder `matrix` in each of your working directories
to avoid the repeated calculation of all required matrices,

## Documentation

You find a copy of the [thesis](docs/Dissertation.pdf) [1] and [corrections](docs/errata5.pdf) in the folder `docs`.
There are also some more documents included which I prepared in that time for presentations or as memos:
An [introduction](docs/terms.pdf) to states and energy levels of lanthanide ions (in German),
an [introduction](docs/tensors.pdf) to tensor operators for $l^N$ configurations,
a [collection](docs/wigner.pdf) of useful properties of Wigner 3j- and 6j-symbols,
and a [summary](docs/hamiltonian.pdf) of the Hamiltonians for $l^N$ configurations.
For the generation of nice term schemata with transitions I developed the Python package `termdia`.
It was required for the generation of `docs/terms.pdf` and you find it in the folder `termdia-1.3` together
with a [manual](termdia-1.3/Anleitung.pdf) (in German).

## Reference

[1] Reinhard Caspary: "Applied Rare-Earth Spectroscopy for Fiber Laser Optimization", doctoral dissertation at
Technische Universität Braunschweig, published with Shaker, Aachen, 2002