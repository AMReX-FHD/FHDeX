# AMREX_HOME defines the directory in which we will find all the AMReX code. It
# can be defined in tools/Path.local (we have to check if it exists first though)
ifneq ("$(wildcard ../../../tools/Path.local)","")
  $(info Loading tools/Path.local...)
  include        ../../../tools/Path.local
endif

# ... alternatively if you set AMREX_HOME as an environment variable, this line
# will be ignored => this line acts as a default
AMREX_HOME ?= ../../../../amrex
FBOXLIB_HOME ?= ../../../../FBoxLib

DEBUG     = FALSE
PROFILE   = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 3
DSMC      = TRUE

USE_PARTICLES = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include ../src_models/Make.package
VPATH_LOCATIONS   += ../src_models
INCLUDE_LOCATIONS += ../src_models

include ../../../src_hydro/Make.package
VPATH_LOCATIONS   += ../../../src_hydro/
INCLUDE_LOCATIONS += ../../../src_hydro/

include ../../../src_rng/Make.package
VPATH_LOCATIONS   += ../../../src_rng/
INCLUDE_LOCATIONS += ../../../src_rng/

include ../../../src_gmres/Make.package
VPATH_LOCATIONS   += ../../../src_gmres/
INCLUDE_LOCATIONS += ../../../src_gmres/

include ../../../src_common/src_F90/Make.package
VPATH_LOCATIONS   += ../../../src_common/src_F90
INCLUDE_LOCATIONS += ../../../src_common/src_F90

include ../../../src_common/Make.package
VPATH_LOCATIONS   += ../../../src_common/
INCLUDE_LOCATIONS += ../../../src_common/

include ../../../src_immersed-boundary/src_F90/Make.package
VPATH_LOCATIONS   += ../../../src_immersed-boundary/src_F90
INCLUDE_LOCATIONS += ../../../src_immersed-boundary/src_F90

include ../../../src_immersed-boundary/Make.package
VPATH_LOCATIONS   += ../../../src_immersed-boundary
INCLUDE_LOCATIONS += ../../../src_immersed-boundary

include ../../../src_particles/MakeKernels.package
VPATH_LOCATIONS   += ../../../src_particles
INCLUDE_LOCATIONS += ../../../src_particles


include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

# include $(AMREX_HOME)/Src/Extern/SWFFT/Make.package
# INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/SWFFT
# VPATH_LOCATIONS	  += $(AMREX_HOME)/Src/Extern/SWFFT

# LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
