# Flocking Dean-Kawasaki: one source directory, two builds.
#   make DIM=3   SPDE on (x, y, theta)       -> main3d.gnu.MPI.ex
#   make DIM=2   particles with theta        -> main2d.gnu.MPI.ex
PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM       ?= 3

MAX_SPEC   = 8

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE
USE_CUDA   = FALSE

ifeq ($(DIM),2)
USE_PARTICLES = TRUE
else
USE_PARTICLES = FALSE
endif

Bpack   := ./Make.package
Blocs   := .

include Make.Adv
