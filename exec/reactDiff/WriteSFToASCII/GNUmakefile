
AMREX_HOME ?= ../../../../amrex/

DEBUG	= TRUE
#DEBUG	= FALSE

DIM	= 3
DIM	= 2

COMP    = gcc

PRECISION = DOUBLE

USE_MPI   = FALSE
USE_OMP   = FALSE

###################################################

EBASE     = WriteSFToASCII

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

