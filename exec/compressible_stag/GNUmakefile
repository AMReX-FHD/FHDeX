# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex/

DEBUG         = FALSE
USE_MPI       = TRUE
USE_OMP       = FALSE
USE_CUDA      = FALSE
CUDA_ARCH     = 70
USE_HIP       = FALSE
COMP          = gnu
DIM           = 3
TINY_PROFILE  = FALSE
MAX_SPEC      = 8
MAX_REAC      = 5
USE_FFT       = TRUE

USE_PARTICLES = FALSE
DO_TURB       = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

include ../../src_compressible_stag/Make.package
VPATH_LOCATIONS   += ../../src_compressible_stag/
INCLUDE_LOCATIONS += ../../src_compressible_stag/

include ../../src_compressible/Make.package
VPATH_LOCATIONS   += ../../src_compressible/
INCLUDE_LOCATIONS += ../../src_compressible/

include ../../src_chemistry/Make.package
VPATH_LOCATIONS   += ../../src_chemistry/
INCLUDE_LOCATIONS += ../../src_chemistry/

include ../../src_MFsurfchem/Make.package
VPATH_LOCATIONS   += ../../src_MFsurfchem/
INCLUDE_LOCATIONS += ../../src_MFsurfchem/

include ../../src_rng/Make.package
VPATH_LOCATIONS   += ../../src_rng/
INCLUDE_LOCATIONS += ../../src_rng/

include ../../src_common/Make.package
VPATH_LOCATIONS   += ../../src_common/
INCLUDE_LOCATIONS += ../../src_common/

include $(AMREX_HOME)/Src/Base/Make.package

include ../../src_analysis/Make.package
VPATH_LOCATIONS   += ../../src_analysis/
INCLUDE_LOCATIONS += ../../src_analysis/

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

ifeq ($(DO_TURB), TRUE)
  DEFINES += -DTURB
endif

MAXSPECIES := $(strip $(MAX_SPEC))
DEFINES += -DMAX_SPECIES=$(MAXSPECIES)

MAXREACTION := $(strip $(MAX_REAC))
DEFINES += -DMAX_REACTION=$(MAXREACTION)

#ifeq ($(USE_CUDA),TRUE)
#    PYTORCH_ROOT := ./test_MFsurfchem_COAr_eq_ML/libtorch_cuda
#else
#    PYTORCH_ROOT := ./test_MFsurfchem_COAr_eq_ML/libtorch_cpu
#endif

#ifeq ($(USE_CUDA),TRUE)
#    PYTORCH_ROOT := $(abspath $(CURDIR))/test_MFsurfchem_COAr_eq_ML/libtorch_cuda
#    TORCH_RPATH  := $$ORIGIN/test_MFsurfchem_COAr_eq_ML/libtorch_cuda/lib
#else
#    PYTORCH_ROOT := $(abspath $(CURDIR))/test_MFsurfchem_COAr_eq_ML/libtorch_cpu
#    TORCH_RPATH  := $$ORIGIN/test_MFsurfchem_COAr_eq_ML/libtorch_cpu/lib
#endif

#TORCH_LIBPATH = $(PYTORCH_ROOT)/lib
#
#ifeq ($(USE_CUDA),TRUE)
#    TORCH_LIBS = -ltorch -ltorch_cpu -lc10 -lc10_cuda -lcuda
#else
#    TORCH_LIBS = -ltorch -ltorch_cpu -lc10
#endif


#INCLUDE_LOCATIONS += $(PYTORCH_ROOT)/include \
#                     $(PYTORCH_ROOT)/include/torch/csrc/api/include
#LIBRARY_LOCATIONS += $(TORCH_LIBPATH)

#DEFINES += -D_GLIBCXX_USE_CXX11_ABI=1
#ifeq ($(USE_CUDA),TRUE)
#    LDFLAGS += -Xlinker "--no-as-needed,-rpath $(TORCH_RPATH) $(TORCH_LIBS)"
#else
#    LDFLAGS += -Wl,--no-as-needed,-rpath=$(TORCH_RPATH) $(TORCH_LIBS)
#endif



ifeq ($(USE_CUDA),TRUE)
    PYTORCH_ROOT := $(abspath $(CURDIR))/test_MFsurfchem_COAr_eq_ML/libtorch_cuda
else
    PYTORCH_ROOT := $(abspath $(CURDIR))/test_MFsurfchem_COAr_eq_ML/libtorch_cpu
endif

TORCH_LIBPATH = $(PYTORCH_ROOT)/lib

ifeq ($(USE_CUDA),TRUE)
    TORCH_LIBS = -ltorch -ltorch_cpu -lc10 -lc10_cuda -lcuda
else
    TORCH_LIBS = -ltorch -ltorch_cpu -lc10
endif

INCLUDE_LOCATIONS += $(PYTORCH_ROOT)/include \
                     $(PYTORCH_ROOT)/include/torch/csrc/api/include
LIBRARY_LOCATIONS += $(TORCH_LIBPATH)

DEFINES += -D_GLIBCXX_USE_CXX11_ABI=1

# FIX: Use TORCH_LIBPATH for the rpath, and ensure commas/spaces are correct!
ifeq ($(USE_CUDA),TRUE)
    LDFLAGS += -Xlinker --no-as-needed -Xlinker -rpath=$(TORCH_LIBPATH) $(TORCH_LIBS)
else
    LDFLAGS += -Wl,--no-as-needed,-rpath=$(TORCH_LIBPATH) $(TORCH_LIBS)
endif
