# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
# AMREX_HOME ?= ../../../amrex/
PELE_PHYSICS_HOME ?= ../../../PelePhysics/
AMREX_HOME        ?= $(abspath $(lastword $(PELE_PHYSICS_HOME)/Submodules/amrex))
SUNDIALS_HOME     ?= $(abspath $(lastword $(PELE_PHYSICS_HOME)/Submodules/sundials))

# FHDeX compiler flags
DEBUG         = FALSE
USE_MPI       = TRUE
USE_OMP       = FALSE
USE_CUDA      = FALSE
USE_HIP       = FALSE
USE_SUNDIALS  = TRUE
COMP          = gnu
DIM           = 3
TINY_PROFILE  = FALSE
MAX_SPEC      = 8
MAX_REAC      = 5
USE_FFT       = TRUE
USE_PARTICLES = FALSE
DO_TURB       = FALSE

# Pele Physics compiler flags
USE_PELE_PHYSICS = TRUE
PELE_CVODE_FORCE_YCORDER = FALSE
PELE_USE_MAGMA = FALSE
PELE_COMPILE_AJACOBIAN = FALSE
Eos_Model := Fuego
Chemistry_Model := LiDryer
Transport_Model := Simple

# FHDeX Definition
MAXSPECIES := $(strip $(MAX_SPEC))
DEFINES += -DMAX_SPECIES=$(MAXSPECIES)

MAXREACTION := $(strip $(MAX_REAC))
DEFINES += -DMAX_REACTION=$(MAXREACTION)

# Turbulence
ifeq ($(DO_TURB), TRUE)
  DEFINES += -DTURB
endif

# Set up the compiler locations
include $(AMREX_HOME)/Tools/GNUMake/Make.defs
ifeq ($(USE_AMREX_LIB), TRUE)
  CEXE_sources :=
endif

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

# PelePhysics
ifeq ($(USE_PELE_PHYSICS), TRUE)

  DEFINES           += -DPELEPHYSICS
  
  UTILITY_HOME = $(PELE_PHYSICS_HOME)/Source/Utility/
  VPATH_LOCATIONS += $(UTILITY_HOME)
  INCLUDE_LOCATIONS += $(UTILITY_HOME)
  include $(UTILITY_HOME)/BlackBoxFunction/Make.package
  VPATH_LOCATIONS += $(UTILITY_HOME)/BlackBoxFunction
  INCLUDE_LOCATIONS += $(UTILITY_HOME)/BlackBoxFunction
  include $(UTILITY_HOME)/Utilities/Make.package
  VPATH_LOCATIONS += $(UTILITY_HOME)/Utilities
  INCLUDE_LOCATIONS += $(UTILITY_HOME)/Utilities

  PP_SRC_HOME 	     = $(PELE_PHYSICS_HOME)/Source
  EXTERN_CORE       += $(PP_SRC_HOME)
  include $(PP_SRC_HOME)/Make.package
  VPATH_LOCATIONS   += $(PP_SRC_HOME)
  INCLUDE_LOCATIONS += $(PP_SRC_HOME)

  ifeq ($(PELE_COMPILE_AJACOBIAN), TRUE)
    DEFINES += -DPELE_COMPILE_AJACOBIAN
  endif
  ifeq ($(PELE_CVODE_FORCE_YCORDER), TRUE)
    DEFINES += -DPELE_CVODE_FORCE_YCORDER
  endif

  # EOS
  EOS_HOME = $(PELE_PHYSICS_HOME)/Source/Eos/
  ifeq ($(Eos_Model),$(filter $(Eos_Model),GammaLaw))
    DEFINES += -DUSE_GAMMALAW_EOS
  endif
  ifeq ($(Eos_Model),$(filter $(Eos_Model),Fuego))
    DEFINES += -DUSE_FUEGO_EOS
  endif
  ifeq ($(Eos_Model),$(filter $(Eos_Model),Soave-Redlich-Kwong))
    DEFINES += -DUSE_SRK_EOS
  endif
  EXTERN_CORE       += $(EOS_HOME)
  include $(EOS_HOME)/Make.package
  INCLUDE_LOCATIONS += $(EOS_HOME)
  VPATH_LOCATIONS   += $(EOS_HOME)
  ifeq ($(Eos_Model), Fuego)
    TRANSPORT_TYPE := IDEAL_GAS
  else
    ifeq ($(Eos_Model), GammaLaw)
      TRANSPORT_TYPE := IDEAL_GAS
    else
      TRANSPORT_TYPE := REAL_GAS
    endif
  endif
  
  # Reactions
  REACTIONS_HOME ?= $(PELE_PHYSICS_HOME)/Source/Reactions
  include $(REACTIONS_HOME)/Make.package
  INCLUDE_LOCATIONS += $(REACTIONS_HOME)
  VPATH_LOCATIONS   += $(REACTIONS_HOME)
  ifeq ($(Eos_Model), GammaLaw)
    ifneq ($(Chemistry_Model), Null)
      $(error Chemistry_Model definition not compatible with Eos_Model=GammaLaw)
    endif
  endif
  CHEM_HOME = $(PELE_PHYSICS_HOME)/Mechanisms/$(Chemistry_Model)
  CHEM_ALL  = $(PELE_PHYSICS_HOME)/Mechanisms
  VPATH_LOCATIONS += $(CHEM_HOME) $(CHEM_ALL)
  include $(CHEM_HOME)/Make.package
  include $(CHEM_ALL)/Make.package
	Blocs += $(CHEM_HOME) $(CHEM_ALL)
  REACTIONS_HOME = $(PELE_PHYSICS_HOME)/Source/Reactions
  EXTERN_CORE       += $(REACTIONS_HOME)
  include $(PELE_PHYSICS_HOME)/ThirdParty/Make.ThirdParty

  ifdef Chemistry_Model
    CHEM_HOME = $(PELE_PHYSICS_HOME)/Mechanisms/$(Chemistry_Model)
  endif
  include $(CHEM_HOME)/Make.package
  VPATH_LOCATIONS += $(CHEM_HOME)
  INCLUDE_LOCATIONS += $(CHEM_HOME)

  # Transport
  TRANSPORT_HOME = $(PELE_PHYSICS_HOME)/Source/Transport
  ifeq ($(Transport_Model),$(filter $(Transport_Model),Constant))
    DEFINES += -DUSE_CONSTANT_TRANSPORT
  endif
  ifeq ($(Transport_Model),$(filter $(Transport_Model),Simple))
    DEFINES += -DUSE_SIMPLE_TRANSPORT
  endif
  ifeq ($(Transport_Model),$(filter $(Transport_Model),Sutherland))
    DEFINES += -DUSE_SUTHERLAND_TRANSPORT
  endif
  EXTERN_CORE       += $(TRANSPORT_HOME)
  include $(TRANSPORT_HOME)/Make.package
  INCLUDE_LOCATIONS += $(TRANSPORT_HOME)
  VPATH_LOCATIONS   += $(TRANSPORT_HOME)
endif


TPL:
	@echo "==> Building SUNDIALS library"
	@cd $(PELE_PHYSICS_HOME)/ThirdParty; $(MAKE) sundials SUNDIALS_HOME=$(SUNDIALS_HOME) AMREX_HOME=$(AMREX_HOME) USE_SYCL=$(USE_SYCL) USE_HIP=$(USE_HIP) USE_CUDA=$(USE_CUDA) PELE_USE_KLU=$(PELE_USE_KLU) PELE_USE_MAGMA=$(PELE_USE_MAGMA) DEBUG=$(DEBUG) COMP=$(HOSTCC) NVCC=$(COMP) CUDA_ARCH=$(CUDA_ARCH) PRECISION=$(PRECISION)

TPLclean:
	@echo "==> Removing SUNDIALS library"
	@cd $(PELE_PHYSICS_HOME)/ThirdParty; $(MAKE) SUNDIALS_HOME=$(SUNDIALS_HOME) AMREX_HOME=$(AMREX_HOME) USE_SYCL=$(USE_SYCL) USE_HIP=$(USE_HIP) USE_CUDA=$(USE_CUDA) PELE_USE_KLU=$(PELE_USE_KLU) PELE_USE_MAGMA=$(PELE_USE_MAGMA) DEBUG=$(DEBUG) COMP=$(HOSTCC) NVCC=$(COMP) CUDA_ARCH=$(CUDA_ARCH) PRECISION=$(PRECISION) clean

TPLrealclean:
	@echo "==> Removing SUNDIALS library"
	@cd $(PELE_PHYSICS_HOME)/ThirdParty; $(MAKE) SUNDIALS_HOME=$(SUNDIALS_HOME) AMREX_HOME=$(AMREX_HOME) USE_SYCL=$(USE_SYCL) USE_HIP=$(USE_HIP) USE_CUDA=$(USE_CUDA) PELE_USE_KLU=$(PELE_USE_KLU) PELE_USE_MAGMA=$(PELE_USE_MAGMA) DEBUG=$(DEBUG) COMP=$(HOSTCC) NVCC=$(COMP) CUDA_ARCH=$(CUDA_ARCH) PRECISION=$(PRECISION) realclean

