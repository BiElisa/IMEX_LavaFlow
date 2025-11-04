# Depth-averaged lava flow model 🌋

[![DOI](https://zenodo.org/badge/453471156.svg)](https://doi.org/10.5281/zenodo.7900929)

IMEX_LavaFlow is a Fortran 90 numerical simulation code designed to model the propagation of lava flows on topographic surfaces. The model is based on the numerical solution of a two-dimensional system of Shallow Water Equations (SWE), integrated with source and relaxation terms to capture the effects of flow resistance and lava rheology. The model includes Additional characteristics are the implementation of vertical profiles of velocity and temperature and temperature-dependent viscosity. 

#### Numerical method
The code uses a Runge-Kutta Explicit-IMplicit (IMEX) scheme for time integration:

* Hyperbolic Part (Explicit): Solved explicitly with a Finite Volume method (semidiscrete, central scheme), ensuring stability and precision in flow transport.

* Stiff Terms (Implicit): Treated implicitly, essential for managing the numerical rigidity introduced by nonlinear terms such as rheology.

* Implicit Resolution: The nonlinear system is solved using the Newton-Raphson algorithm, with the Jacobian matrix computed numerically using the extremely precise Complex Step Derivative technique (implemented in complexify.f90).

### Prerequisites
| Category  | Requirement  | Details  | Installation  |
|-----------|--------------|----------|--------------|
| Language  | Fortran 90/95/2003  |    |  |
| Compiler        | GNU Fortran (gfortran) | Required.     | [Fortran-Lang](https://fortran-lang.org/),  [GNU Fortran docs](https://gcc.gnu.org/fortran/)|
| Parallelization | OpenMP  | Used for computational acceleration.  | Included with gfortran
| Libraries       | LAPACK / BLAS  | Essential for high-performance linear algebra; used to solve the implicit system. | [LAPACK site](https://www.netlib.org/lapack/)
| Build System    | Autotools   | The project uses `Makefile.am` and `Makefile.in` (Autoconf/Automake).  | [GNU Autotools](https://www.gnu.org/software/automake/)

Install all required components on Debian/Ubuntu-based systems with:
```bash
$ sudo apt install gfortran
$ sudo apt install liblapack-dev libblas-dev
$ sudo apt install make
$ sudo apt install automake autoconf libtool
```

### Installation and Compiling
Download locally the repository (if needed):
```bash
$ git clone https://github.com/BiElisa/IMEX_LavaFlow.git
$ cd IMEX_LavaFlow
```
Configure (you might need to make the file executable):
```bash
$ chmod +x configure
$ ./configure
```

Compile:
```bash
$ make
$ make install
```

To compile the code with OpenMP add the following flag in src/Makefile:
1) with gfortran: `-fopenmp`
2) with intel: `-qopenmp`

The executable `IMEX_LavaFlow` is copied in the `bin` folder.

Examples are located in the EXAMPLES folder.

## Run an example
To run an example you need to have in the same folder an input file `IMEX_LavaFlow.inp`, an asc file for the topography (for example `topography_dem.asc` contained into the folder `EXAMPLE_DEM_FOR_STUDENTS`), and a copy or link of the executable file contained into the `bin` folder.
Execute the simulation by launching:
```bash
$ ./IMEX_LavaFlow
```
> *Note:*
> If your simulation requires a 