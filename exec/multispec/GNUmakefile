# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex/

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
USE_CUDA  = FALSE
COMP      = gnu
DIM       = 2
DSMC      = FALSE

TINY_PROFILE = TRUE
USE_PARTICLES = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include ../../src_hydro/Make.package
VPATH_LOCATIONS   += ../../src_hydro/
INCLUDE_LOCATIONS += ../../src_hydro/

include ../../src_multispec/Make.package
VPATH_LOCATIONS   += ../../src_multispec/
INCLUDE_LOCATIONS += ../../src_multispec/

include ../../src_analysis/Make.package
VPATH_LOCATIONS   += ../../src_analysis/
INCLUDE_LOCATIONS += ../../src_analysis/

include ../../src_rng/Make.package
VPATH_LOCATIONS   += ../../src_rng/
INCLUDE_LOCATIONS += ../../src_rng/

include ../../src_gmres/Make.package
VPATH_LOCATIONS   += ../../src_gmres/
INCLUDE_LOCATIONS += ../../src_gmres/


include ../../src_common/Make.package
VPATH_LOCATIONS   += ../../src_common/
INCLUDE_LOCATIONS += ../../src_common/

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

ifeq ($(findstring cgpu, $(HOST)), cgpu)
  CXXFLAGS += $(FFTW)
endif

ifeq ($(USE_CUDA),TRUE)
  LIBRARIES += -lcufft
else
  LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3
endif
