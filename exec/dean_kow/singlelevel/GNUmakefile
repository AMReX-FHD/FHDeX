# AMREX_HOME defines the directory in which we will find all the AMReX code.
AMREX_HOME ?= ../../../amrex
USE_MPI   = TRUE
USE_OMP   = FALSE
COMP      = gcc
DIM       = 2
MAX_SPEC      = 8

USE_CUDA = FALSE
USE_HIP  = FALSE
USE_SYCL = FALSE

DEBUG     = TRUE

USE_PARTICLES = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ../../../src_deankow/singlelevel/Make.package
VPATH_LOCATIONS  += ../../../src_deankow/singlelevel
INCLUDE_LOCATIONS += ../../../src_deankow/singlelevel

include ../../../src_rng/Make.package
VPATH_LOCATIONS   += ../../../src_rng/
INCLUDE_LOCATIONS += ../../../src_rng/

include ../../../src_common/Make.package
VPATH_LOCATIONS   += ../../../src_common/
INCLUDE_LOCATIONS += ../../../src_common/

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

MAXSPECIES := $(strip $(MAX_SPEC))
DEFINES += -DMAX_SPECIES=$(MAXSPECIES)
