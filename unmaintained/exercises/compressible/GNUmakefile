# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../../amrex/
FHDEX_HOME ?= ../../..

DEBUG     = FALSE
USE_MPI   = TRUE
USE_OMP   = FALSE
COMP      = gnu
DIM       = 3
USE_PARTICLES = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include $(FHDEX_HOME)/src_compressible/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_compressible/src_F90
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_compressible/src_F90

include $(FHDEX_HOME)/src_compressible/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_compressible/
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_compressible/


include $(FHDEX_HOME)/src_analysis/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_analysis/
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_analysis/

include $(FHDEX_HOME)/src_rng/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_rng/src_F90
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_rng/src_F90

include $(FHDEX_HOME)/src_rng/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_rng/
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_rng/

include $(FHDEX_HOME)/src_common/src_F90/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_common/src_F90
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_common/src_F90

include $(FHDEX_HOME)/src_common/Make.package
VPATH_LOCATIONS   += $(FHDEX_HOME)/src_common/
INCLUDE_LOCATIONS += $(FHDEX_HOME)/src_common/

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package

 include $(AMREX_HOME)/Src/Extern/SWFFT/Make.package
 INCLUDE_LOCATIONS	+= $(AMREX_HOME)/Src/Extern/SWFFT
 VPATH_LOCATIONS	+= $(AMREX_HOME)/Src/Extern/SWFFT

 LIBRARIES += -L$(FFTW_DIR) -lfftw3_mpi -lfftw3

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
