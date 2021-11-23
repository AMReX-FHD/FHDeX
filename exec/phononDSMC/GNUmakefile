# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex/
FHDeX ?= ../../
FDSMC = $(FHDeX)/exec/phononDSMC

DEBUG     = FALSE
PROFILE   = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
USE_CUDA  = FALSE
COMP      = gnu
DIM       = 3

TINY_PROFILE  = FALSE
USE_PARTICLES = TRUE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

VPATH_LOCATIONS   += $(FDSMC)
INCLUDE_LOCATIONS += $(FDSMC)

include $(FHDeX)/src_rng/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_rng
INCLUDE_LOCATIONS += $(FHDeX)/src_rng

include $(FHDeX)/src_analysis/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_analysis/
INCLUDE_LOCATIONS += $(FHDeX)/src_analysis/

include $(FHDeX)/src_common/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_common/
INCLUDE_LOCATIONS += $(FHDeX)/src_common/


include $(FHDeX)/src_particles/MakeDSMC.package
VPATH_LOCATIONS   += $(FHDeX)/src_particles/
INCLUDE_LOCATIONS += $(FHDeX)/src_particles/

include $(FHDeX)/src_geometry/Make.package
VPATH_LOCATIONS   += $(FHDeX)/src_geometry
INCLUDE_LOCATIONS += $(FHDeX)/src_geometry


include $(FDSMC)/Make.package

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package
include $(AMREX_HOME)/Src/Particle/Make.package


#ifeq ($(NERSC_HOST),) # check if NERSC_HOST is empty => hack to test if on cori
#    LIBRARIES += -L$(FFTW_DIR) -lfftw3
#    ifeq ($(USE_OMP), TRUE)
#	LIBRARIES += -lfftw3_omp
#    endif
#    ifeq ($(USE_MPI), TRUE)
#	LIBRARIES += -lfftw3_mpi
#    endif
#else
#    $(info NERSC_HOST not empty => relying on compiler wrapper for fftw3)
#endif

LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3_omp -lfftw3 -lgomp
#LIBRARIES += -L$(FFTW_DIR) -lfftw3

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

#DEFINES	+= -DDSMC=$(DSMC)

ifeq ($(findstring cgpu, $(HOST)), cgpu)
CXXFLAGS += $(FFTW)
endif
