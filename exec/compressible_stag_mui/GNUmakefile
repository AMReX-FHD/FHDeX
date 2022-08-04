# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex/
MUI_HOME ?= ../../src_mui/

BL_NO_FORT    = TRUE
DEBUG         = FALSE
USE_MPI       = TRUE
USE_OMP       = FALSE
USE_CUDA      = FALSE
COMP          = gnu
DIM           = 3
TINY_PROFILE  = FALSE
USE_MUI       = TRUE
USE_SLURM     = FALSE
USE_AMREX_MPMD = TRUE

USE_PARTICLES = FALSE

ifeq ($(MAKECMDGOALS),spk)
  EBASE = spk
else ifeq ($(USE_AMREX_MPMD),TRUE)
  EBASE = fhd
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

ifeq ($(USE_AMREX_MPMD),TRUE)
  USE_MUI := FALSE
endif

ifeq ($(MAKECMDGOALS),spk)

  include SPPARKS_MUI/Make.package
  VPATH_LOCATIONS   += SPPARKS_MUI
  INCLUDE_LOCATIONS += SPPARKS_MUI

else

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

  ifeq ($(USE_MUI),TRUE)
    include ../../src_surfchem_mui/Make.package
    VPATH_LOCATIONS   += ../../src_surfchem_mui/
    INCLUDE_LOCATIONS += ../../src_surfchem_mui/
  else
    include ../../src_surfchem_amrex/Make.package
    VPATH_LOCATIONS   += ../../src_surfchem_amrex/
    INCLUDE_LOCATIONS += ../../src_surfchem_amrex/
  endif

  include ../../src_rng/Make.package
  VPATH_LOCATIONS   += ../../src_rng/
  INCLUDE_LOCATIONS += ../../src_rng/

  include ../../src_common/Make.package
  VPATH_LOCATIONS   += ../../src_common/
  INCLUDE_LOCATIONS += ../../src_common/

  ifeq ($(USE_MUI),TRUE)
    VPATH_LOCATIONS   += $(MUI_HOME)
    INCLUDE_LOCATIONS += $(MUI_HOME)
  endif

  include ../../src_analysis/Make.package
  VPATH_LOCATIONS   += ../../src_analysis/
  INCLUDE_LOCATIONS += ../../src_analysis/

endif

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

ifeq ($(MAKECMDGOALS),spk)

vpath main_spk.cpp $(srcTempDir)

$(srcTempDir)/main_spk.cpp: SPPARKS_MUI/main.cpp
	@if [ ! -d $(srcTempDir) ]; then mkdir -p $(srcTempDir); fi
	@$(SHELL) -ec "ln -fs ../../../$< $@"

spk: $(executable)

else
  ifeq ($(findstring cgpu, $(HOST)), cgpu)
    CXXFLAGS += $(FFTW)
  endif

  ifeq ($(USE_CUDA),TRUE)
    LIBRARIES += -lcufft
  else
    LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3
  endif
endif

ifeq ($(USE_MUI), TRUE)
  DEFINES += -DMUI
endif

ifeq ($(USE_AMREX_MPMD),TRUE)
  DEFINES += -DUSE_AMREX_MPMD
endif

ifeq ($(USE_SLURM), TRUE)
  DEFINES += -DSLURM
endif
