# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex/

DEBUG         = FALSE
USE_MPI       = TRUE
USE_OMP       = FALSE
USE_CUDA      = FALSE
USE_HIP       = FALSE
COMP          = gnu
DIM           = 3
TINY_PROFILE  = FALSE
MAX_SPEC      = 8
MAX_REAC      = 5

USE_PARTICLES = FALSE
DO_TURB       = FALSE

USE_HEFFTE_FFTW   = FALSE
USE_HEFFTE_CUFFT  = FALSE
USE_HEFFTE_ROCFFT = FALSE
USE_DISTRIBUTED_FFT = TRUE

ifeq ($(USE_HEFFTE_FFTW),TRUE)
  HEFFTE_HOME ?= ../../../heffte/
else ifeq ($(USE_HEFFTE_CUFFT),TRUE)
  HEFFTE_HOME ?= ../../../heffte-org/build_aware/
else ifeq ($(USE_HEFFTE_ROCFFT),TRUE)
  HEFFTE_HOME ?= ../../../heffte-org/build_noaware/
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include ../../src_compressible_stag/Make.package
VPATH_LOCATIONS   += ../../src_compressible_stag/
INCLUDE_LOCATIONS += ../../src_compressible_stag/

include ../../src_compressible/Make.package
VPATH_LOCATIONS   += ../../src_compressible/
INCLUDE_LOCATIONS += ../../src_compressible/

include ../../src_chemistry/Make.package
VPATH_LOCATIONS   += ../../src_chemistry/
INCLUDE_LOCATIONS += ../../src_chemistry/

include ../../src_MFsurfchem/Make.package
VPATH_LOCATIONS   += ../../src_MFsurfchem/
INCLUDE_LOCATIONS += ../../src_MFsurfchem/

include ../../src_rng/Make.package
VPATH_LOCATIONS   += ../../src_rng/
INCLUDE_LOCATIONS += ../../src_rng/

include ../../src_common/Make.package
VPATH_LOCATIONS   += ../../src_common/
INCLUDE_LOCATIONS += ../../src_common/

#ifeq ($(USE_HEFFTE_FFTW),TRUE)
#  include $(HEFFTE_HOME)/src/Make.package
#else ifeq ($(USE_HEFFTE_CUFFT),TRUE)
#  include $(HEFFTE_HOME)/src/Make.package
#else ifeq ($(USE_HEFFTE_ROCFFT),TRUE)
#  include $(HEFFTE_HOME)/src/Make.package
#endif

ifeq ($(USE_HEFFTE_FFTW),TRUE)
  DEFINES += -DHEFFTE_FFTW
  LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3 -lfftw3f
else ifeq ($(USE_HEFFTE_CUFFT),TRUE)
  DEFINES += -DHEFFTE_CUFFT
else ifeq ($(USE_HEFFTE_ROCFFT),TRUE)
  DEFINES += -DHEFFTE_ROCFFT
endif

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/FFT/Make.package

include ../../src_analysis/Make.package
VPATH_LOCATIONS   += ../../src_analysis/
INCLUDE_LOCATIONS += ../../src_analysis/
 
include $(AMREX_HOME)/Tools/GNUMake/Make.rules

ifeq ($(USE_CUDA),TRUE)
  LIBRARIES += -lcufft
else ifeq ($(USE_HIP),TRUE)
  # Use rocFFT.  ROC_PATH is defined in amrex
  INCLUDE_LOCATIONS += $(ROC_PATH)/rocfft/include
  LIBRARY_LOCATIONS += $(ROC_PATH)/rocfft/lib
  LIBRARIES += -L$(ROC_PATH)/rocfft/lib -lrocfft
else
  LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3 -lfftw3f
endif

ifeq ($(DO_TURB), TRUE)
  DEFINES += -DTURB
endif

MAXSPECIES := $(strip $(MAX_SPEC))
DEFINES += -DMAX_SPECIES=$(MAXSPECIES)

MAXREACTION := $(strip $(MAX_REAC))
DEFINES += -DMAX_REACTION=$(MAXREACTION)
