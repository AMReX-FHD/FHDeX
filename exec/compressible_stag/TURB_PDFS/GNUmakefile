AMREX_HOME ?= ../../../../amrex/

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gcc

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = FALSE
USE_CUDA  = TRUE

TINY_PROFILE = FALSE

###################################################

#EBASE     = main
#EBASE     = main_decomp
EBASE     = main_multisteps

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

vpath %.c   : . $(vpathdir)
vpath %.h   : . $(vpathdir)
vpath %.cpp : . $(vpathdir)
vpath %.H   : . $(vpathdir)
vpath %.F   : . $(vpathdir)
vpath %.f   : . $(vpathdir)
vpath %.f90 : . $(vpathdir)

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
