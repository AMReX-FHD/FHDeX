# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex/
FBOXLIB_HOME ?= ../../../../../programs/FBoxLib/

DEBUG     = FASLE
USE_MPI   = TRUE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 3
DSMC      = TRUE

USE_PARTICLES = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include $(AMREX_HOME)/Src/F_BaseLib/FParallelMG2.mak
VPATH_LOCATIONS   += $(AMREX_HOME)/Src/F_BaseLib/
INCLUDE_LOCATIONS += $(AMREX_HOME)/Src/F_BaseLib/

#include $(FBOXLIB_HOME)/Src/BaseLib/FParallelMG_RNG.mak
#VPATH_LOCATIONS   += $(FBOXLIB_HOME)/Src/BaseLib/
#INCLUDE_LOCATIONS += $(FBOXLIB_HOME)/Src/BaseLib/


include ../../src_compressible/Make.package
VPATH_LOCATIONS   += ../../src_compressible/
INCLUDE_LOCATIONS += ../../src_compressible/

include ../../src_compressible/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_compressible/src_F90
INCLUDE_LOCATIONS += ../../src_compressible/src_F90

#include ../../src_analysis/src_F90/Make.package
#VPATH_LOCATIONS   += ../../src_analysis/src_F90
#INCLUDE_LOCATIONS += ../../src_analysis/src_F90

#include ../../src_analysis/Make.package
#VPATH_LOCATIONS   += ../../src_analysis/
#INCLUDE_LOCATIONS += ../../src_analysis/

include ../../src_rng/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_rng/src_F90
INCLUDE_LOCATIONS += ../../src_rng/src_F90

include ../../src_rng/Make.package
VPATH_LOCATIONS   += ../../src_rng/
INCLUDE_LOCATIONS += ../../src_rng/

include ../../src_gmres/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_gmres/src_F90
INCLUDE_LOCATIONS += ../../src_gmres/src_F90

include ../../src_gmres/Make.package
VPATH_LOCATIONS   += ../../src_gmres/
INCLUDE_LOCATIONS += ../../src_gmres/

include ../../src_common/src_F90/Make.package
VPATH_LOCATIONS   += ../../src_common/src_F90
INCLUDE_LOCATIONS += ../../src_common/src_F90

include ../../src_common/Make.package
VPATH_LOCATIONS   += ../../src_common/
INCLUDE_LOCATIONS += ../../src_common/

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/C_CellMG/Make.package
include $(AMREX_HOME)/Src/LinearSolvers/MLMG/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package

 include $(AMREX_HOME)/Src/Extern/SWFFT/Make.package
 INCLUDE_LOCATIONS	+= $(AMREX_HOME)/Src/Extern/SWFFT
 VPATH_LOCATIONS		+= $(AMREX_HOME)/Src/Extern/SWFFT

 LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3_omp -lfftw3

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
