# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../../amrex/
FHDEX_HOME ?= ../../..

DEBUG     = FALSE
PROFILE   = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 2

USE_PARTICLES = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include $(FHDEX_HOME)/src_hydro/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_hydro/src_F90
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_hydro/src_F90

include $(FHDEX_HOME)/src_hydro/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_hydro/
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_hydro/

include $(FHDEX_HOME)/src_multispec/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_multispec/src_F90
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_multispec/src_F90

include $(FHDEX_HOME)/src_multispec/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_multispec/
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_multispec/

include $(FHDEX_HOME)/src_analysis/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_analysis/src_F90
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_analysis/src_F90

include $(FHDEX_HOME)/src_analysis/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_analysis/
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_analysis/

include $(FHDEX_HOME)/src_rng/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_rng/src_F90
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_rng/src_F90

include $(FHDEX_HOME)/src_rng/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_rng/
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_rng/

include $(FHDEX_HOME)/src_gmres/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_gmres/src_F90
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_gmres/src_F90

include $(FHDEX_HOME)/src_gmres/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_gmres/
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_gmres/

include $(FHDEX_HOME)/src_common/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_common/src_F90
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_common/src_F90

include $(FHDEX_HOME)/src_common/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_common/
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_common/

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

include $(AMREX_HOME)/Src/Extern/SWFFT/Make.package
INCLUDE_LOCATIONS	+= $(AMREX_HOME)/Src/Extern/SWFFT
VPATH_LOCATIONS		+= $(AMREX_HOME)/Src/Extern/SWFFT

LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
