AMREX_HOME ?= ../../../../amrex/

DEBUG         = FALSE
USE_MPI       = TRUE
USE_OMP       = FALSE
USE_CUDA      = FALSE
USE_HIP       = FALSE
COMP          = gnu
DIM           = 3
TINY_PROFILE  = FALSE

USE_HEFFTE_FFTW   = FALSE
USE_HEFFTE_CUFFT  = FALSE
USE_HEFFTE_ROCFFT = FALSE

ifeq ($(USE_HEFFTE_FFTW),TRUE)
  HEFFTE_HOME ?= ../../../../heffte/
else ifeq ($(USE_HEFFTE_CUFFT),TRUE)
  HEFFTE_HOME ?= ../../../../heffte/
else ifeq ($(USE_HEFFTE_ROCFFT),TRUE)
  HEFFTE_HOME ?= ../../../../heffte/
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

ifeq ($(USE_HEFFTE_FFTW),TRUE)
  include $(HEFFTE_HOME)/src/Make.package
else ifeq ($(USE_HEFFTE_CUFFT),TRUE)
  include $(HEFFTE_HOME)/src/Make.package
else ifeq ($(USE_HEFFTE_ROCFFT),TRUE)
  include $(HEFFTE_HOME)/src/Make.package
endif

include ./Make.package

include $(AMREX_HOME)/Src/Base/Make.package

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

ifeq ($(USE_HEFFTE_FFTW),TRUE)
  DEFINES += -DHEFFTE_FFTW
  LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3 -lfftw3f
else ifeq ($(USE_HEFFTE_CUFFT),TRUE)
  DEFINES += -DHEFFTE_CUFFT
else ifeq ($(USE_HEFFTE_ROCFFT),TRUE)
  DEFINES += -DHEFFTE_ROCFFT
endif
