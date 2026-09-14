# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME ?= ../../../amrex

DEBUG        = FALSE
USE_MPI      = TRUE
USE_OMP      = FALSE
COMP         = gnu
DIM          = 2
USE_FFT      = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include ../../src_rng/Make.package
VPATH_LOCATIONS   += ../../src_rng/
INCLUDE_LOCATIONS += ../../src_rng/

include ../../src_common/Make.package
VPATH_LOCATIONS   += ../../src_common/
INCLUDE_LOCATIONS += ../../src_common/

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

MAX_SPEC  = 2
MAXSPECIES := $(strip $(MAX_SPEC))
DEFINES += -DMAX_SPECIES=$(MAXSPECIES)
