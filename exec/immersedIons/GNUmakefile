# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex/
FHDeX ?= ../../
FMEMBRANE = $(FHDeX)/exec/fluctuatingMembrane
HYDROLIB_HOME := ../../../LowMachFHD/HydroGrid

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 3
DSMC      = FALSE

TINY_PROFILE  = FALSE
USE_PARTICLES = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

VPATH_LOCATIONS   += $(FMEMBRANE)
INCLUDE_LOCATIONS += $(FMEMBRANE)

include $(FMEMBRANE)/Make.package

include $(FHDeX)/src_rng/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_rng
INCLUDE_LOCATIONS += $(FHDeX)/src_rng

include $(FHDeX)/src_rng/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_rng/src_F90
INCLUDE_LOCATIONS += $(FHDeX)/src_rng/src_F90

include $(FHDeX)/src_hydro/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_hydro/
INCLUDE_LOCATIONS += $(FHDeX)/src_hydro/

include $(FHDeX)/src_hydro/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_hydro/src_F90
INCLUDE_LOCATIONS += $(FHDeX)/src_hydro/src_F90

include $(FHDeX)/src_analysis/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_analysis/
INCLUDE_LOCATIONS += $(FHDeX)/src_analysis/

include $(FHDeX)/src_analysis/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_analysis/src_F90
INCLUDE_LOCATIONS += $(FHDeX)/src_analysis/src_F90

#include $(AMREX_HOME)/Src/test/FParallelMG_RNG.mak
#VPATH_LOCATIONS   += $(AMREX_HOME)/Src/test/
#INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/test/

include $(FHDeX)/src_gmres/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_gmres/
INCLUDE_LOCATIONS += $(FHDeX)/src_gmres/

include $(FHDeX)/src_gmres/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_gmres/src_F90
INCLUDE_LOCATIONS += $(FHDeX)/src_gmres/src_F90

include $(FHDeX)/src_common/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_common/
INCLUDE_LOCATIONS += $(FHDeX)/src_common/

include $(FHDeX)/src_common/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_common/src_F90
INCLUDE_LOCATIONS += $(FHDeX)/src_common/src_F90

include $(FHDeX)/src_particles/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_particles/
INCLUDE_LOCATIONS += $(FHDeX)/src_particles/

include $(FHDeX)/src_particles/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_particles/src_F90
INCLUDE_LOCATIONS += $(FHDeX)/src_particles/src_F90

include $(FHDeX)/src_surfaces/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_surfaces
INCLUDE_LOCATIONS += $(FHDeX)/src_surfaces

include $(FHDeX)/src_surfaces/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_surfaces/src_F90
INCLUDE_LOCATIONS += $(FHDeX)/src_surfaces/src_F90

include $(FHDeX)/src_electrostatic/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_electrostatic
INCLUDE_LOCATIONS += $(FHDeX)/src_electrostatic

include $(FHDeX)/src_immersed-boundary/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_immersed-boundary/src_F90
INCLUDE_LOCATIONS += $(FHDeX)/src_immersed-boundary/src_F90

include $(FHDeX)/src_immersed-boundary/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_immersed-boundary
INCLUDE_LOCATIONS += $(FHDeX)/src_immersed-boundary

include $(HYDROLIB_HOME)/GPackage_noRNG.mak
VPATH_LOCATIONS += $(HYDROLIB_HOME)

include $(FMEMBRANE)/Make_precheck.package
include $(FMEMBRANE)/Make_debug.package

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package
include $(AMREX_HOME)/Src/Extern/SWFFT/Make.package

INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/Extern/SWFFT
VPATH_LOCATIONS	  += $(AMREX_HOME)/Src/Extern/SWFFT

LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3_omp -lfftw3

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

ifeq ($(DSMC), TRUE)
  DEFINES += -DDSMC
endif

#DEFINES	+= -DDSMC=$(DSMC)				

