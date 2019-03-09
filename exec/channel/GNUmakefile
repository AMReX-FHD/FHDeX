# AMREX_HOME defines the directory in which we will find all the AMReX code. It
# can be defined in tools/Path.local (we have to check if it exists first though)
ifneq ("$(wildcard ../../tools/Path.local)","")
  $(info Loading tools/Path.local...)
  include        ../../tools/Path.local
endif

# ... alternatively if you set AMREX_HOME as an environment variable, this line
# will be ignored => this line acts as a default
AMREX_HOME ?= ../../../amrex
FBOXLIB_HOME ?= ../../../FBoxLib

DEBUG     = FALSE
PROFILE   = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 3
DSMC      = TRUE

USE_PARTICLES = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

# ifneq ("$(wildcard $(AMREX_HOME)/Src/F_BaseLib)", "")
#  $(info Loading AMReX F_BaseLib)
#  include $(AMREX_HOME)/Src/F_BaseLib/FParallelMG2.mak
#  VPATH_LOCATIONS   += $(AMREX_HOME)/Src/F_BaseLib/
#  INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/F_BaseLib/
# else
#   $(info Loading FBoxLib)
#   include $(FBOXLIB_HOME)/Src/BaseLib/FParallelMG_RNG.mak
#   VPATH_LOCATIONS   += $(FBOXLIB_HOME)/Src/BaseLib/
#   INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/Src/BaseLib/
# endif

include ../../src_hydro/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_hydro/src_F90
INCLUDE_LOCATIONS += ../../src_hydro/src_F90

include ../../src_hydro/Make.package
VPATH_LOCATIONS   += ../../src_hydro/
INCLUDE_LOCATIONS += ../../src_hydro/

include ../../src_analysis/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_analysis/src_F90
INCLUDE_LOCATIONS += ../../src_analysis/src_F90

include ../../src_analysis/Make.package
VPATH_LOCATIONS   += ../../src_analysis/
INCLUDE_LOCATIONS += ../../src_analysis/

include ../../src_rng/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_rng/src_F90
INCLUDE_LOCATIONS += ../../src_rng/src_F90

include ../../src_rng/Make.package
VPATH_LOCATIONS   += ../../src_rng/
INCLUDE_LOCATIONS += ../../src_rng/

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

#include ../../src_immersed-boundary/src_F90/Make.package
#VPATH_LOCATIONS   += ../../src_immersed-boundary/src_F90
#INCLUDE_LOCATIONS += ../../src_immersed-boundary/src_F90

# include ../../src_immersed-boundary/Make.package
# VPATH_LOCATIONS   += ../../src_immersed-boundary
# INCLUDE_LOCATIONS += ../../src_immersed-boundary


include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Src/Extern/SWFFT/Make.package
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/SWFFT
VPATH_LOCATIONS	  += $(AMREX_HOME)/Src/Extern/SWFFT


# check if running on macOS (as there are some subtle differences to linux).
# UNAME_S = $(shell uname -s)
# ifeq ($(UNAME_S), Darwin)
#     # on macOS (if libfftw3 has been installed via homebrew, there are no
#     # separated omp and mpi libraries)
#     LIBRARIES += -L$(FFTW_DIR) -lfftw3
# else
    LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3_omp -lfftw3
# endif

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
