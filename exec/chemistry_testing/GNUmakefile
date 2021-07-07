# AMREX_HOME defines the directory in which we will find all the AMReX code.
# On this machine I created an extra directory so I needed to go one level up
# compared to Dr. Nonaka's code

AMREX_HOME ?= ../../../amrex

DEBUG        = FALSE
USE_MPI      = TRUE 
USE_OMP      = FALSE
COMP         = gnu
DIM          = 2

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
VPATH_LOCATIONS  += .
INCLUDE_LOCATIONS += .

include ../../src_chemistry/Make.package
VPATH_LOCATIONS   += ../../src_chemistry/
INCLUDE_LOCATIONS += ../../src_chemistry/

include $(AMREX_HOME)/Src/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
