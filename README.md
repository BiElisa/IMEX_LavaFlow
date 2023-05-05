# Depth-averaged lava flow model

DOI: https://zenodo.org/badge/latestdoi/453471156

Shallow water model for lava flow with vertical profiles of velocity and temperature and temperature-dependent viscosity. 

To compile:

> ./configure

To compile the code with OpenMP add the following flag in src/Makefile:
1) with gfortran: -fopenmp
2) with intel: -qopenmp

> make

> make install

The executable is copied in the bin folder.

Several examples can be found in the EXAMPLES folder.
