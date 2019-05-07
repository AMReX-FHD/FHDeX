# AMREX_HOME defines the directory in which we will find all the AMReX code. It
# can be defined in tools/Path.local (we have to check if it exists first though)
ifneq ("$(wildcard ../../../tools/Path.local)","")
  $(info Loading tools/Path.local...)
  include       ../../../tools/Path.local
endif

# ... alternatively if you set AMREX_HOME as an environment variable, this line
# will be ignored => this line acts as a default
AMREX_HOME  ?= ../../../../amrex
DIM         ?= 3


include_dirs = -I.. -I$(AMREX_HOME)/Src/Base

test_kernel_6p   : test_kernel_6p.o
	gfortran -o test_kernel_6p test_kernel_6p.o

test_kernel_6p.o : test_kernel_6p.F90
	gfortran -O3 -c test_kernel_6p.F90 $(include_dirs) -DAMREX_SPACEDIM=$(DIM)

clean :
	trash *.mod *.o
