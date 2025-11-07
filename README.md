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
| `HU_BC%*`, `HV_BC%*` | Boundary conditions for `HU` and `HV`, the momentum along the two directions (m²s⁻¹) |
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

### EXPL_TERMS_PARAMETERS

| Variable  |  Meaning  |
|-----------|-----------|
| `GRAV`    | Gravitational acceleration (m s⁻²). |

### RHEOLOGY_PARAMETERS

| Variable  |  Meaning  |
|-----------|-----------|
| `RHEOLOGY_MODEL`  | Model type. 3 = Temperature-dependent viscosity (Eq. 6 from Costa & Macedonio, 2005). |
| `VISC_PAR`  | Exponential parameter of temperature-dependent model. Corresponds to `b` (K⁻¹) in Eq. 6 from Costa & Macedonio, 2005. |
| `NU_REF` | Reference kinematic viscosity (m²s⁻¹) - Eq. 6 from Costa & Macedonio, 2005. Could be set equal to the viscosity at the vent. |
| `T_REF` | Reference temperature (K) - Eq. 6 from Costa & Macedonio, 2005. Could be set equal to the temperature at the vent. `T_REF` and  `NU_REF` should be chosen consistently, i.e., `NU_REF` should correspond to the viscosity at `T_REF`. |
| `TAU0` | Yield stress for Bingham fluids (Pa). Set equal to 0 to have Netwonian fluid.|

### RADIAL_SOURCE_PARAMETERS

| Variable  |  Meaning  |
|-----------|-----------|
| `X_SOURCE`, `Y_SOURCE` | Source Coordinate (m). |
| `R_SOURCE` | Radius of the source (m). |
| `VEL_SOURCE` | Injection velocity from the source (ms⁻¹). `VEL_SOURCE` is  the ratio between the effusion rate (which is a volume flow rate) (m³s⁻¹) and the source area (m²), that is `π R_SOURCE²`|
| `T_SOURCE` | Temperature of the emitted material (K). |
| `TIME_PARAM` | Source time parameters (see below for more details). |

The input parameter `TIME_PARAM` is a 4-element array defining the temporal behavior of the source. The four entries are:
| Index | Name | Description    |
| ----- | ---- | -------------- |
| TIME_PARAM(1) | Period (`T_period`) | Total period after which the source repeats. |
| TIME_PARAM(2) | Duration (`T_active`) | Duration within the period during which the source is active. |
| TIME_PARAM(3) | Ramp (`T_ramp`) | Time to ramp up from 0 to full strength and ramp down from full strength to 0. |
| TIME_PARAM(4) | Phase shift (`T_shift`) | Time offset applied to shift the source activation within the period. |



A key concept linking these parameters and the simulation is the variable `t_coeff` used in the code, which serves as a scaling factor for time-based adjustments of the flow emission.


Below are three example cases showing how `t_coeff` varies over time for different `TIME_PARAM` settings.

#### Case 1 – Continuous emission
**TIME_PARAM:** `[50.0, 50.0, 0.0, 50.0]`  
<img 
   width="800" 
   height="400" 
   alt="t_coeff_Case_1" 
   src="https://github.com/user-attachments/assets/2cbe5d67-d1fa-47dc-a322-7946f93fc6e0" 
/>

#### Case 2 – Shorter active duration with ramp
**TIME_PARAM:** `[50.0, 30.0, 10.0, 0.0]`  
<img 
   width="800" 
   height="400" 
   alt="t_coeff_Case_2" 
   src="https://github.com/user-attachments/assets/daa4fa83-cf04-4e7f-a459-c757d9c35dec" 
/>

#### Case 3 – Phase-shifted, partial ramp
**TIME_PARAM:** `[60.0, 40.0, 5.0, 10.0]`  
<img 
   width="800" 
   height="400" 
   alt="t_coeff_Case_3" 
   src="https://github.com/user-attachments/assets/be669d84-1c49-4722-9dde-1b95982e891e" 
/>

#### Case 4 – Sawtooth source
**TIME_PARAM:** `[50.0, 50.0, 50.0, 0.0]`  
This profile exhibits a gradual ramp-up to full strength, followed by a reset to zero, repeating periodically.  
<img 
   width="800" 
   height="400" 
   alt="t_coeff_Case_4" 
   src="https://github.com/user-attachments/assets/3c0dcf4e-772c-45a3-887a-88e87cec20aa" 
/>

#### Case 5 – Symmetric Sawtooth source
**TIME_PARAM:** `[50.0, 50.0, 25.0, 0.0]`  
<img 
   width="800" 
   height="400" 
   alt="t_coeff_Case_5" 
   src="https://github.com/user-attachments/assets/9b673666-697d-4bc9-beec-c0968e4313b1" 
/>

#### Case 6 – Square waves
**TIME_PARAM:** `[50.0, 30.0, 0.0, 10.0]`  
<img 
   width="800" 
   height="400" 
   alt="t_coeff_Case_6" 
   src="https://github.com/user-attachments/assets/debab8f9-43d0-4e88-862d-7059968c9241" 
/>



### GAS_TRANSPORT_PARAMETERS 
No need to touch these parameters for lava flow simulations.

| Variable  |  Meaning  |
|-----------|-----------|
| `SP_HEAT_A`      | Specific heat (J kg⁻¹ K⁻¹). |
| `SP_GAS_CONST_A` | Gas specific constant (J kg⁻¹ K⁻¹). |
| `KIN_VISC_A`     | Kinematic viscosity (m²s⁻¹). |
| `PRES`           | Ambient pressure (Pa). |
| `T_AMBIENT`      | Ambient temperature (K). |


### LIQUID_TRANSPORT_PARAMETERS 
| Variable  |  Meaning  |
|-----------|-----------|
| `SP_HEAT_L`  | Specific heat of the liquid (J kg⁻¹ K⁻¹). |
| `RHO_L`      | Density (kg m⁻³).                         |
| `KIN_VISC_L` | Kinematic viscosity (m² s⁻¹).             |

### Notes
* All numeric variables can be written in E-notation format (e.g. `1.0E-3`).
* Logical values ​​must be `T` (true) or `F` (false).
