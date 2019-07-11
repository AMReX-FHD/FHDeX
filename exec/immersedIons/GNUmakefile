# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex/

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 3
DSMC      = FALSE


TINY_PROFILE  = FALSE
USE_PARTICLES = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include ../../src_rng/Make.package
VPATH_LOCATIONS   += ../../src_rng
INCLUDE_LOCATIONS += ../../src_rng

include ../../src_rng/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_rng/src_F90
INCLUDE_LOCATIONS += ../../src_rng/src_F90

include ../../src_hydro/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_hydro/src_F90
INCLUDE_LOCATIONS += ../../src_hydro/src_F90

include ../../src_hydro/Make.package
VPATH_LOCATIONS   += ../../src_hydro/
INCLUDE_LOCATIONS += ../../src_hydro/

include ../../src_analysis/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_analysis/src_F90
INCLUDE_LOCATIONS += ../../src_analysis/src_F90

#include ../../src_analysis/Make.package
#VPATH_LOCATIONS   += ../../src_analysis/
#INCLUDE_LOCATIONS += ../../src_analysis/

#include $(AMREX_HOME)/Src/test/FParallelMG_RNG.mak
#VPATH_LOCATIONS   += $(AMREX_HOME)/Src/test/
#INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/test/

include ../../src_gmres/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_gmres/src_F90
INCLUDE_LOCATIONS += ../../src_gmres/src_F90

include ../../src_gmres/Make.package
VPATH_LOCATIONS   += ../../src_gmres/
INCLUDE_LOCATIONS += ../../src_gmres/

include ../../src_common/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_common/src_F90
INCLUDE_LOCATIONS += ../../src_common/src_F90

include ../../src_common/Make.package
VPATH_LOCATIONS   += ../../src_common/
INCLUDE_LOCATIONS += ../../src_common/

include ../../src_particles/Make.package
VPATH_LOCATIONS   += ../../src_particles/
INCLUDE_LOCATIONS += ../../src_particles/

include ../../src_particles/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_particles/src_F90
INCLUDE_LOCATIONS += ../../src_particles/src_F90

include ../../src_surfaces/Make.package
VPATH_LOCATIONS   += ../../src_surfaces
INCLUDE_LOCATIONS += ../../src_surfaces

include ../../src_electrostatic/Make.package
VPATH_LOCATIONS   += ../../src_electrostatic
INCLUDE_LOCATIONS += ../../src_electrostatic

include ./Make_precheck.package

include ../../src_surfaces/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_surfaces/src_F90
INCLUDE_LOCATIONS += ../../src_surfaces/src_F90

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package
include $(AMREX_HOME)/Src/Extern/SWFFT/Make.package
INCLUDE_LOCATIONS	+= $(AMREX_HOME)/Src/Extern/SWFFT
VPATH_LOCATIONS		+= $(AMREX_HOME)/Src/Extern/SWFFT

LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3_omp -lfftw3

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

ifeq ($(DSMC), TRUE)
  DEFINES += -DDSMC
endif

#DEFINES	+= -DDSMC=$(DSMC)				

