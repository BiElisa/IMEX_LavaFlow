
AUTOMAKE_OPTIONS = foreign

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
bin_PROGRAMS = IMEX_LavaFlow

bindir=@bindir@

# list of sources for the binary
IMEX_LavaFlow_SOURCES = \
                parameters_2d.f90 \
                complexify.f90 \
                geometry_2d.f90 \
                constitutive_2d.f90 \
                solver_2d.f90 \
                init_2d.f90 \
                inpout_2d.f90 \
                IMEX_LavaFlow.f90

FCFLAGS         = -fopenmp -Ofast -funroll-all-loops
FCDEBUGFLAGS    = -Wall -g
# Use these flags to DEBUG
#FCFLAGS = -fopenmp -g -O0 -fcheck=all -fbacktrace -Wall -Wextra
#FCFLAGS = -fopenmp -g -O0 -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow

LIBS            = -lm -llapack

FC              = gfortran
MPFC            = gfortran
LINKER          = gfortran

clean-local:
	-rm -f *.mod *.o



