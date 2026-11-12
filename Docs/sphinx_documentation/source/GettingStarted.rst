.. _chap:getting-started:

Getting Started
===============

This chapter walks through downloading FHDeX, building the incompressible
hydrodynamics application ``exec/hydro`` and running a short test problem.
All other applications are built and run the same way.

Prerequisites
-------------

* A C++17 compiler (GCC is the default; ``COMP`` selects others) and GNU Make.
* MPI. All applications default to ``USE_MPI = TRUE``.
* FFTW3 with MPI support, for applications built with ``USE_FFT = TRUE``
  (structure factors and spectral analysis). The ``GNUmakefile`` links
  ``-lfftw3_mpi -lfftw3`` from ``$(FFTW_DIR)``; CUDA builds link cuFFT
  instead.
* Optionally CUDA or HIP for GPU builds (``USE_CUDA`` / ``USE_HIP``); see
  :ref:`amrex:sec:gpu:build`.

AMReX itself does not need to be installed separately. Each FHDeX application
compiles the AMReX sources it needs directly, as described in
:ref:`amrex:sec:build:make`.

Downloading the code
--------------------

FHDeX and AMReX are cloned side by side, because the ``GNUmakefile`` in every
application looks for AMReX in ``../../../amrex`` by default:

.. code-block:: bash

   git clone https://github.com/AMReX-Codes/amrex.git
   git clone https://github.com/AMReX-FHD/FHDeX.git

which gives

.. code-block:: text

   amrex/
   FHDeX/
     exec/hydro/      <- AMREX_HOME ?= ../../../amrex/

If AMReX lives somewhere else, set the ``AMREX_HOME`` environment variable;
it overrides the default in the ``GNUmakefile``.

.. note::

   FHDeX is developed against the ``development`` branch of AMReX (the
   default branch). If the build fails inside AMReX sources, first update
   both repositories.

Building ``exec/hydro``
-----------------------

Each directory under ``exec/`` holds one application. Its ``GNUmakefile``
sets the build options and pulls in the FHDeX source libraries it needs:

.. code-block:: make

   AMREX_HOME ?= ../../../amrex/

   DEBUG     = FALSE
   USE_MPI   = TRUE
   USE_OMP   = FALSE
   USE_CUDA  = FALSE
   COMP      = gnu
   DIM       = 3
   MAX_SPEC  = 2
   MAX_ELEM  = 1     # must be MAX_SPEC*(MAX_SPEC-1)/2
   USE_FFT   = TRUE

   include $(AMREX_HOME)/Tools/GNUMake/Make.defs

   include ../../src_hydro/Make.package
   VPATH_LOCATIONS   += ../../src_hydro/
   INCLUDE_LOCATIONS += ../../src_hydro/
   # ... the same for src_multispec, src_analysis, src_rng, src_gmres, src_common

   include $(AMREX_HOME)/Src/Base/Make.package
   # ... AMReX FFT, Boundary and MLMG components
   include $(AMREX_HOME)/Tools/GNUMake/Make.rules

Any option can be overridden on the command line. To build the 2D version:

.. code-block:: bash

   cd FHDeX/exec/hydro
   make -j8 DIM=2

The build produces an executable whose name encodes the options, e.g.
``main2d.gnu.MPI.ex``. The AMReX documentation lists all standard build
options (``DEBUG``, ``COMP``, ``USE_OMP``, ``TINY_PROFILE``, ...); see
:ref:`amrex:sec:build:make`.

Two options are specific to FHDeX and are compiled in as preprocessor
constants:

``MAX_SPEC``
   Maximum number of species (``-DMAX_SPECIES``). Fixed-size per-species
   arrays such as ``molmass`` or ``rhobar`` are sized by it, so the
   ``nspecies`` runtime parameter must not exceed it.

``MAX_ELEM``
   Size of the arrays holding the off-diagonal entries of symmetric
   species-pair matrices, such as the Maxwell--Stefan diffusion coefficients
   ``Dbar`` (``-DMAX_ELEMENT``). It must equal ``MAX_SPEC*(MAX_SPEC-1)/2``;
   ``src_multispec`` aborts at start-up otherwise.

Compressible and chemistry applications also use ``MAX_REAC``
(``-DMAX_REACTION``), the maximum number of reactions.

Running a test problem
----------------------

Every application directory contains a set of ``inputs_*`` files. Pass one
as the first argument of the executable:

.. code-block:: bash

   mpiexec -n 4 ./main2d.gnu.MPI.ex inputs_regression_vortex_2d

This runs a 32 x 32 vortex problem for three time steps, split into four
16 x 16 boxes (one per MPI rank). The inputs file uses the AMReX
``ParmParse`` format (:ref:`amrex:sec:basics:parmparse`):

.. code-block:: text

   prob_lo   = 0.0 0.0       # physical lo coordinate
   prob_hi   = 3200. 3200.   # physical hi coordinate
   prob_type = 1             # 1 = vortex
   n_cells       = 32 32
   max_grid_size = 16 16
   fixed_dt  = 10
   max_step  = 3
   plot_int  = 3
   bc_vel_lo = 1 2           # -1 periodic, 1 slip, 2 no-slip
   bc_vel_hi = 1 2

Any parameter can also be overridden on the command line, which is handy for
parameter scans:

.. code-block:: bash

   mpiexec -n 4 ./main2d.gnu.MPI.ex inputs_regression_vortex_2d max_step=10 plot_int=5

Chapter :ref:`chap:params` lists every parameter with its default value.

Output
------

The hydro application writes:

``plt0000003/``
   AMReX plotfile with cell-averaged velocities and pressure, written every
   ``plot_int`` steps. The prefix is set by ``plot_base_name``.

``stagx0000003/``, ``stagy...``, ``stagz...``
   The face-centred (staggered) velocity components, when ``plot_stag = 1``.

``plt_SF_mag000000003/``, ``plt_SF_real_imag000000003/``
   Static structure factor of the fluctuating fields (magnitude, and real and
   imaginary parts), when ``struct_fact_int > 0``; see
   :ref:`sec:code:analysis`.

``chk0000003/``
   Checkpoint files every ``chk_int`` steps (prefix ``chk_base_name``).
   Restart from one with ``restart = 3``.

Plotfiles can be viewed with AMReX's ``amrvis``, VisIt, ParaView or yt; see
:ref:`amrex:sec:visit` and :ref:`amrex:sec:postprocessing`.
