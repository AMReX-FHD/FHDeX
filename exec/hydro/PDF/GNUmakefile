AMREX_HOME ?= ../../..

DEBUG	= TRUE
DEBUG	= FALSE

DIM	= 3

COMP    = gcc

PRECISION = DOUBLE

USE_MPI   = TRUE
USE_OMP   = FALSE
USE_CUDA  = FALSE

TINY_PROFILE = FALSE

###################################################

#EBASE     = PDF
EBASE     = PDF_all

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
