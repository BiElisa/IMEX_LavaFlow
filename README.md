# Depth-averaged lava flow model 🌋

[![DOI](https://zenodo.org/badge/453471156.svg)](https://doi.org/10.5281/zenodo.7900929)

IMEX_LavaFlow is a Fortran 90 numerical simulation code designed to model the propagation of lava flows on topographic surfaces. The model is based on the numerical solution of a two-dimensional system of Shallow Water Equations (SWE), integrated with source and relaxation terms to capture the effects of flow resistance and lava rheology. The model includes Additional characteristics are the implementation of vertical profiles of velocity and temperature and temperature-dependent viscosity. 

#### Numerical method
The code uses a Runge-Kutta Explicit-IMplicit (IMEX) scheme for time integration:

* Hyperbolic Part (Explicit): Solved explicitly with a Finite Volume method (semidiscrete, central scheme), ensuring stability and precision in flow transport.

* Stiff Terms (Implicit): Treated implicitly, essential for managing the numerical rigidity introduced by nonlinear terms such as rheology.

* Implicit Resolution: The nonlinear system is solved using the Newton-Raphson algorithm, with the Jacobian matrix computed numerically using the extremely precise Complex Step Derivative technique (implemented in complexify.f90).

### Bibliographic references
<a id="ref1"></a>
[1] Biagioli, E., de' Michieli Vitturi, M. & Di Benedetto, F. (2021). *Modified shallow water model for viscous fluids and positivity preserving numerical approximation*. Applied Matheematical Modelling, 94(1-2), 482–505.  
   DOI: [10.1016/j.apm.2020.12.036](https://doi.org/10.1016/j.apm.2020.12.036)

<a id="ref2"></a>
[2] Biagioli, E., de' Michieli Vitturi, M., Di Benedetto, F., & Polacci, M. (2023). *Benchmarking a new 2.5D shallow water model for lava flows*. Journal of Volcanology and Geothermal Research, 444(7), 107935.  
   DOI: [110.1016/j.jvolgeores.2023.107935](https://doi.org/10.1016/j.jvolgeores.2023.107935)

## Prerequisites
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

## Installation and Compiling
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

## Read and modify an input file

The input file `IMEX_LavaFlow.inp` defines all the physical, numerical and geometrical parameters required to execute a numerical simulation.

The file is composed of blocks, Fortran NAMELISTS, each introduced by `&BLOCK_NAME` and terminated by `/`. Each block controls a specific aspect of the simulation.

### RUN_PARAMETERS

| Variable  |  Meaning  |
|-----------|-----------|
| `RUN_NAME`  | Run name, used for output files.  |
| `RESTART`   | If `T`, resume from a saved state; if `F`, start a new run. |
| `T_START`  | Initial simulation time (s). |
| `T_END`  | Final simulation time (s). |
| `DT_OUTPUT` | Time interval between saving results (s). |
| `OUTPUT_CONS_FLAG`  | Save *conservative variables* such as `h, hU, hT` in a file `.q_2d`; these variables are used for numerical purpose, see [[1]](#ref1). |
| `OUTPUT_ESRI_FLAG`  | Export results in ESRI ASCII format (in a file`.asc`). |
| `OUTPUT_PHYS_FLAG`  | Save *physical variables* such as `h, U, T` (thickness, depth-averaged velocity, depth-averaged temperature) in a file `.p_2d`. |
| `OUTPUT_RUNOUT_FLAG` | If `T`, save the maximum distances traveled by the flow. |
| `VERBOSE_LEVEL` | On-screen message level (0 = silent). |

### NEWRUN_PARAMETERS

| Variable  |  Meaning  |
|-----------|-----------|
| `X0`, `Y0`  | UTM coordinates (m) of the southwest corner of the grid. |
| `COMP_CELLS_X`, `COMP_CELLS_Y` | Number of cells in the X and Y direction. |
| `CELL_SIZE` | Cell width (m). |
| `RHEOLOGY_FLAG` | Activates rheological laws (e.g. non-Newtonian viscosity).  |
| `ENERGY_FLAG` | Activate the energy equation (set `F` for lava flows). |
| `LIQUID_FLAG` | `T` lava is considered a liquid |
| `RADIAL_SOURCE_FLAG` | Activate radial sources from the surface (set `F` for lava flows). |
| `BOTTOM_RADIAL_SOURCE_FLAG` | Activate radial sources from the bottom (set `T` for lava flows). |
| `COLLAPSING_VOLUME_FLAG` | Activate initial volume collapse (set `F` for lava flows).|
| `VELOCITY_PROFILE_FLAG`, `TEMPERATURE_PROFILE_FLAG` | Activate velocity and temperature vertical profiles. |

### BOUNDARY CONDITIONS

There are 4 identical blocks for the borders West, East, South and North, `&WEST_BOUNDARY_CONDITIONS`, etc. For each direction, we define boundary conditions for

| Variable  |  Meaning  |
|-----------|-----------|
| `H_BC%*` | Boundary condition for `h`, the thickness (m) |
| `HU_BC%*`, `HV_BC%*` | Boundary conditions for `HU` and `HV`, the momentum along the two directions (m^2/s) |
| `T_BC%*` | Boundary condition for `T`, the temperature (K) |

For each variable, the `FLAG` must be set (for example `H_BC%FLAG`). It can assume the values:
* `1` → Dirichlet condition (imposed value),
* `0` → Neumann condition (zero gradient),
* `-1` → condition not set.

### TEMPERATURE_PARAMETERS

| Variable  |  Meaning  |
|-----------|-----------|
| `EMISSIVITY` | Surface emissivity (0–1).  |
| `ATM_HEAT_TRANSF_COEFF`  | Heat exchange coefficient with the atmosphere (W m⁻² K⁻¹).  |
| `EXP_AREA_FRACT` | Fraction of area exposed to cooling. |
| `C_P`  | Specific heat (J kg⁻¹ K⁻¹). |
| `ENNE`  | Relative thickness of the thermal layer (dimensionless units). |
| `T_ENV`, `T_GROUND` | Ambient and soil temperatures (K). |
| `THERMAL_CONDUCTIVITY_FLUID`, `THERMAL_CONDUCTIVITY_SOIL` | Thermal conductivity of fluid and soil (W m⁻¹ K⁻¹). |
| `EMME`  | Relative thickness of the conductive layer in the soil. |
| `RHO_SOIL`, `C_P_SOIL`, `T_SOIL`  | Density, specific heat and deep soil temperature. |

### NUMERIC_PARAMETERS

| Variable  |  Meaning  |
|-----------|-----------|
| `SOLVER_SCHEME` | Numerical scheme for the iperbolic part of the equations: `"KT"`, `"LxF"`, `"GFORCE"`, `"UP"` (recommended `KT`). |
| `DT0`, `MAX_DT` | Initial and maximum time-step. |
| `CFL` | Courant number (for numerical stability). For 1D simulations: <0.5; for 2D simulations: <0.25.  |
| `LIMITER` | Select the *limiter* for the flux reconstruction (hyperbolic terms): 0 – 3.                  |
| `THETA`, `RECONSTR_COEFF` | Coefficients for the Total Variation Diminishing (TVD)/WENO schemes.                               |
| `N_RK` | Number of steps for the IM-EX Runge-Kutta scheme. |

