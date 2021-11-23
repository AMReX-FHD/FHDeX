# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../../programs/amrex/
FHDeX ?= ../../../
FIONS = $(FHDeX)/exec/tests/RFDanalysis

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
USE_CUDA  = FALSE
COMP      = gnu
DIM       = 3
DSMC      = FALSE

TINY_PROFILE  = FALSE
USE_PARTICLES = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

VPATH_LOCATIONS   += $(FIONS)
INCLUDE_LOCATIONS += $(FIONS)

include $(FHDeX)/src_rng/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_rng
INCLUDE_LOCATIONS += $(FHDeX)/src_rng


include $(FHDeX)/src_hydro/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_hydro/
INCLUDE_LOCATIONS += $(FHDeX)/src_hydro/

include $(FHDeX)/src_analysis/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_analysis/
INCLUDE_LOCATIONS += $(FHDeX)/src_analysis/


include $(FHDeX)/src_gmres/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_gmres/
INCLUDE_LOCATIONS += $(FHDeX)/src_gmres/

include $(FHDeX)/src_common/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_common/
INCLUDE_LOCATIONS += $(FHDeX)/src_common/


include $(FHDeX)/src_particles/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_particles/
INCLUDE_LOCATIONS += $(FHDeX)/src_particles/

include $(FHDeX)/src_geometry/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_geometry
INCLUDE_LOCATIONS += $(FHDeX)/src_geometry

include $(FHDeX)/src_electrostatic/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_electrostatic
INCLUDE_LOCATIONS += $(FHDeX)/src_electrostatic

include $(FHDeX)/src_immersed-boundary/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDeX)src_immersed-boundary/src_F90
INCLUDE_LOCATIONS += $(FHDeX)src_immersed-boundary/src_F90

include $(FHDeX)/src_immersed-boundary/MakeBase.package
VPATH_LOCATIONS   += $(FHDeX)src_immersed-boundary
INCLUDE_LOCATIONS += $(FHDeX)src_immersed-boundary

include $(FIONS)/Make.package
include $(FIONS)/Make_precheck.package

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

ifeq ($(DSMC), TRUE)
  DEFINES += -DDSMC
endif

#DEFINES	+= -DDSMC=$(DSMC)

ifeq ($(findstring cgpu, $(HOST)), cgpu)
  CXXFLAGS += $(FFTW)
endif

ifeq ($(USE_CUDA),TRUE)
  LIBRARIES += -lcufft
else
  LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3
endif
