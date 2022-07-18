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
USE_MUI       = TRUE

USE_PARTICLES = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs
VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include ../../src_compressible/Make.package
VPATH_LOCATIONS   += ../../src_compressible/
INCLUDE_LOCATIONS += ../../src_compressible/

include ../../src_chemistry/Make.package
VPATH_LOCATIONS   += ../../src_chemistry/
INCLUDE_LOCATIONS += ../../src_chemistry/

include ../../src_MFsurfchem/Make.package
VPATH_LOCATIONS   += ../../src_MFsurfchem/
INCLUDE_LOCATIONS += ../../src_MFsurfchem/

include ../../src_surfchem_mui/Make.package
VPATH_LOCATIONS   += ../../src_surfchem_mui/
INCLUDE_LOCATIONS += ../../src_surfchem_mui/

include ../../src_rng/Make.package
VPATH_LOCATIONS   += ../../src_rng/
INCLUDE_LOCATIONS += ../../src_rng/


include ../../src_common/Make.package
VPATH_LOCATIONS   += ../../src_common/
INCLUDE_LOCATIONS += ../../src_common/

VPATH_LOCATIONS   += $(MUI_HOME)
INCLUDE_LOCATIONS += $(MUI_HOME)

include $(AMREX_HOME)/Src/Base/Make.package

include ../../src_analysis/Make.package
VPATH_LOCATIONS   += ../../src_analysis/
INCLUDE_LOCATIONS += ../../src_analysis/

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

ifeq ($(findstring cgpu, $(HOST)), cgpu)
  CXXFLAGS += $(FFTW)
endif

ifeq ($(USE_CUDA),TRUE)
  LIBRARIES += -lcufft
else
  LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3
endif

ifeq ($(USE_MUI), TRUE)
  DEFINES += -DMUI
endif
