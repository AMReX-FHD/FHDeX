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
USE_PARTICLES = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include ../../src_rng/Make.package
VPATH_LOCATIONS   += ../../src_rng/
INCLUDE_LOCATIONS += ../../src_rng/

include ../../src_common/Make.package
VPATH_LOCATIONS   += ../../src_common/
INCLUDE_LOCATIONS += ../../src_common/

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

ifeq ($(findstring cgpu, $(HOST)), cgpu)
  CXXFLAGS += $(FFTW)
endif
