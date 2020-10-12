# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex/
MUI_HOME ?= MUI/

DEBUG         = FALSE
USE_MPI       = TRUE
USE_OMP       = FALSE
COMP          = gnu
DIM           = 3
TINY_PROFILE  = FALSE

USE_PARTICLES = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include ../../src_compressible/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_compressible/src_F90
INCLUDE_LOCATIONS += ../../src_compressible/src_F90

include ../../src_compressible/Make.package
VPATH_LOCATIONS   += ../../src_compressible/
INCLUDE_LOCATIONS += ../../src_compressible/

include ../../src_rng/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_rng/src_F90
INCLUDE_LOCATIONS += ../../src_rng/src_F90

include ../../src_rng/Make.package
VPATH_LOCATIONS   += ../../src_rng/
INCLUDE_LOCATIONS += ../../src_rng/

include ../../src_common/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_common/src_F90
INCLUDE_LOCATIONS += ../../src_common/src_F90

include ../../src_common/Make.package
VPATH_LOCATIONS   += ../../src_common/
INCLUDE_LOCATIONS += ../../src_common/

VPATH_LOCATIONS   += $(MUI_HOME)
INCLUDE_LOCATIONS += $(MUI_HOME)

include $(AMREX_HOME)/Src/Base/Make.package

# structure factor and SWFFT stuff; does not work with CUDA build
ifeq ($(USE_CUDA),FALSE)
 include ../../src_analysis/src_F90/Make.package
 VPATH_LOCATIONS   += ../../src_analysis/src_F90
 INCLUDE_LOCATIONS += ../../src_analysis/src_F90

 include ../../src_analysis/Make.package
 VPATH_LOCATIONS   += ../../src_analysis/
 INCLUDE_LOCATIONS += ../../src_analysis/

 include $(AMREX_HOME)/Src/Extern/SWFFT/Make.package
 INCLUDE_LOCATIONS    += $(AMREX_HOME)/Src/Extern/SWFFT
 VPATH_LOCATIONS      += $(AMREX_HOME)/Src/Extern/SWFFT
endif

LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3

include $(AMREX_HOME)/Tools/GNUMake/Make.rules