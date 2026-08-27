FHDEX_ROOT ?= ../../..
AMREX_HOME ?= $(FHDEX_ROOT)/../amrex

DEBUG         = FALSE
USE_MPI       = TRUE
USE_OMP       = FALSE
USE_CUDA      = FALSE
USE_HIP       = FALSE
COMP          = gnu
DIM           = 3
TINY_PROFILE  = FALSE
PRECISION     = DOUBLE
BL_NO_FORT    = TRUE
USE_FFT       = TRUE
MAX_SPEC      = 2

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include $(FHDEX_ROOT)/src_analysis/Make.package
VPATH_LOCATIONS   += $(FHDEX_ROOT)/src_analysis
INCLUDE_LOCATIONS += $(FHDEX_ROOT)/src_analysis

include $(FHDEX_ROOT)/src_common/Make.package
VPATH_LOCATIONS   += $(FHDEX_ROOT)/src_common
INCLUDE_LOCATIONS += $(FHDEX_ROOT)/src_common

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

ifeq ($(findstring cgpu, $(HOST)), cgpu)
  CXXFLAGS += $(FFTW)
endif

ifeq ($(USE_CUDA),TRUE)
  LIBRARIES += -lcufft
else
  LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3
endif

MAXSPECIES := $(strip $(MAX_SPEC))
DEFINES += -DMAX_SPECIES=$(MAXSPECIES)
