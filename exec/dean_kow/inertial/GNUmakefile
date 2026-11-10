# Inertial Dean-Kawasaki in (x, v): one 2D build runs the SPDE or the
# particles (inertial.run_mode).
#   make -j8 && mpirun -np 4 ./main2d.gnu.MPI.ex inputs_spde
PRECISION  = DOUBLE
PROFILE    = FALSE

DEBUG      = FALSE

DIM        = 2

MAX_SPEC   = 8

COMP	   = gnu

USE_MPI    = TRUE
USE_OMP    = FALSE
USE_CUDA   = FALSE

USE_PARTICLES = TRUE

Bpack   := ./Make.package
Blocs   := .

include Make.Adv
