# AMREX_HOME defines the directory in which we will find all the AMReX code.
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex/
PELE_PHYSICS_HOME ?= ../../../PelePhysics/

# FHDeX compiler flags
DEBUG         = FALSE
USE_MPI       = TRUE
USE_OMP       = FALSE
USE_CUDA      = FALSE
USE_HIP       = FALSE
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

VPATH_LOCATIONS   += .
INCLUDE_LOCATIONS += .

VPATH_LOCATIONS   += ../../src_compressible_stag/
INCLUDE_LOCATIONS += ../../src_compressible_stag/
Bpack += ../../src_compressible_stag/Make.package
Blocs += ../../src_compressible_stag

VPATH_LOCATIONS   += ../../src_compressible/
INCLUDE_LOCATIONS += ../../src_compressible/
Bpack += ../../src_compressible/Make.package
Bpack += ../../src_compressible

VPATH_LOCATIONS   += ../../src_chemistry/
INCLUDE_LOCATIONS += ../../src_chemistry/
Bpack += ../../src_chemistry/Make.package
Blocs += ../../src_chemistry

VPATH_LOCATIONS   += ../../src_rng/
INCLUDE_LOCATIONS += ../../src_rng/
Bpack += ../../src_rng/Make.package
Blocs += ../../src_rng

VPATH_LOCATIONS   += ../../src_common/
INCLUDE_LOCATIONS += ../../src_common/
Bpack += ../../src_common/Make.package
Blocs += ../../src_common

Bpack += $(AMREX_HOME)/Src/Base/Make.package
Blocs += $(AMREX_HOME)/Src/Base

VPATH_LOCATIONS   += ../../src_analysis/
INCLUDE_LOCATIONS += ../../src_analysis/
Bpack += ../../src_analysis/Make.package
Blocs += ../../src_analysis

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

# PelePhysics
ifeq ($(USE_PELE_PHYSICS), TRUE)

  DEFINES           += -DPELEPHYSICS
  PP_SRC_HOME 	     = $(PELE_PHYSICS_HOME)/Source
  EXTERN_CORE       += $(PP_SRC_HOME)
  INCLUDE_LOCATIONS += $(PP_SRC_HOME)
  VPATH_LOCATIONS   += $(PP_SRC_HOME)
  Bpack             += $(PP_SRC_HOME)/Make.package
  Blocs             += $(PP_SRC_HOME)

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
  INCLUDE_LOCATIONS += $(EOS_HOME)
  VPATH_LOCATIONS   += $(EOS_HOME)
  Bpack             += $(EOS_HOME)/Make.package
  Blocs             += $(EOS_HOME)

  # Reactions
  REACTIONS_HOME = $(PELE_PHYSICS_HOME)/Source/Reactions
  EXTERN_CORE       += $(REACTIONS_HOME)
  INCLUDE_LOCATIONS += $(REACTIONS_HOME)
  VPATH_LOCATIONS   += $(REACTIONS_HOME)
  Bpack += $(REACTIONS_HOME)/Make.package
  Blocs += $(REACTIONS_HOME)
  include $(PELE_PHYSICS_HOME)/ThirdParty/Make.ThirdParty

  ifdef Chemistry_Model
    CHEM_HOME = $(PELE_PHYSICS_HOME)/Mechanisms/$(Chemistry_Model)
  endif
  VPATH_LOCATIONS += $(CHEM_HOME)
  Bpack += $(CHEM_HOME)/Make.package
  Blocs += $(CHEM_HOME)

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
  INCLUDE_LOCATIONS += $(TRANSPORT_HOME)
  VPATH_LOCATIONS   += $(TRANSPORT_HOME)
  Bpack             += $(TRANSPORT_HOME)/Make.package
  Blocs             += $(TRANSPORT_HOME)
endif

$(info $${Blocs} is [${Blocs}])

include $(Bpack)

INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)
