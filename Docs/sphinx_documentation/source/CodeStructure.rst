.. _chap:code-structure:

Code Structure
==============

FHDeX has no single executable. It is a set of source libraries
(``src_*``) that applications (``exec/*``) combine at compile time. This
chapter describes how an application is put together and what each library
provides.

How an application is assembled
-------------------------------

Source selection through ``VPATH``
   An application's ``GNUmakefile`` includes the ``Make.package`` of each
   library it needs and appends that directory to ``VPATH_LOCATIONS`` and
   ``INCLUDE_LOCATIONS`` (see :ref:`chap:getting-started` and
   :ref:`amrex:sec:build:make`). Make searches ``VPATH_LOCATIONS`` **in
   order** and uses the first file with a matching name. The application's
   own directory (``.``) is listed first, so an application can replace any
   library file by providing a file with the same name. This is how every
   application supplies its own ``main_driver.cpp``. Library order matters
   too: ``exec/compressible_stag`` lists ``src_compressible_stag`` before
   ``src_compressible``, so the staggered versions of shared file names
   win.

Program entry
   ``src_common/main.cpp`` calls ``amrex::Initialize``, passes the inputs
   file name to ``main_driver()`` and calls ``amrex::Finalize``. The default
   ``src_common/main_driver.cpp`` just aborts with *"write a custom routine
   and put in local source directory"*. A real driver comes either from the
   application directory (e.g. ``exec/hydro/main_driver.cpp``) or from a
   solver library (e.g. ``src_compressible_stag/main_driver.cpp``).

Parameters as namespaces
   Runtime parameters live in C++ namespaces (``common``, ``gmres``,
   ``multispec``, ``compressible``, ...). Each namespace is declared in
   ``src_*/<name>_namespace.H`` and filled by ``Initialize<Name>Namespace()``
   in ``<name>_functions.cpp``, which first sets defaults and then reads the
   inputs file with ``ParmParse`` (:ref:`amrex:sec:basics:parmparse`). Many
   variables are ``AMREX_GPU_MANAGED`` so that GPU kernels can read them
   directly. A driver typically begins with

   .. code-block:: c++

      InitializeCommonNamespace();
      InitializeGmresNamespace();

   Chapter :ref:`chap:params` documents every namespace.

Data layout
   Cell-centred quantities (densities, concentrations, pressure) are
   ``MultiFab`` objects. Velocities and momenta on the staggered (MAC) grid
   are stored as ``std::array<MultiFab, AMREX_SPACEDIM>``, one face-centred
   ``MultiFab`` per direction. The index types are predefined in the
   ``common`` namespace: ``nodal_flag_x``/``_y``/``_z`` for faces,
   ``nodal_flag_xy``/``_xz``/``_yz`` for edges, ``nodal_flag`` for nodes.
   See :ref:`amrex:sec:basics:box` for AMReX index types and
   :ref:`amrex:sec:basics:multifab` for the containers.

Source libraries
----------------

Infrastructure
~~~~~~~~~~~~~~

``src_common``
   Used by every application. Contains the ``common`` parameter namespace,
   the program entry point, and grid utilities: averaging and shifting
   between cell, face, edge and node centres (``AverageFaceToCC``,
   ``AverageCCToFace``, ...), discrete divergence, gradient, Laplacian and
   curl (``ComputeDiv``, ``ComputeGrad``, ``ComputeStagLap``, ...),
   statistics (``ComputeBasicStats``, ``WriteHorizontalAverage``) and
   physical boundary conditions (``MultiFabPhysBC*``). The file
   ``src_common/README_FILLROUTINES`` explains what each ghost-cell filling
   routine does at walls, reservoirs and periodic boundaries; it complements
   :ref:`amrex:sec:basics:boundary`.

``src_rng``
   Random-number fields for the stochastic fluxes, e.g.
   ``MultiFabFillRandom`` and ``MultiFabFillRandomNormal``. These build on
   the AMReX random-number engine (:ref:`amrex:sec:basics:random`).

.. _sec:code:analysis:

``src_analysis``
   Run-time analysis. The ``StructFact`` class accumulates static structure
   factors (spatial power spectra and cross-correlations of fluctuating
   fields) with FFTs and writes them as plotfiles. ``TurbSpectra`` computes
   shell-integrated kinetic-energy and scalar spectra for turbulence runs.

``src_gmres``
   Krylov solver for the saddle-point (Stokes) system of the staggered
   discretisation: GMRES with block preconditioners (``precon_type``), a
   staggered multigrid solver for the velocity block (``StagMGSolver``) and a
   MAC projection for the pressure block. It is controlled by the ``gmres``
   namespace.

Physics solvers
~~~~~~~~~~~~~~~

``src_hydro``
   Incompressible fluctuating hydrodynamics on the staggered grid: the time
   advance, the stochastic momentum flux (``StochMomFlux``), advective
   fluxes (centred or BDS), vorticity, and conversion between momentum and
   velocity.

``src_multispec``
   Low-Mach multispecies mixtures: diffusive and stochastic mass fluxes
   (``DiffusiveMassFluxdiv``, ``StochMassFlux``), mixture properties,
   correction fluxes that keep the mass fractions summing to one, the
   reversible stress for capillary effects, and electrodiffusion of charged
   species. It is controlled by the ``multispec`` namespace.

``src_compressible``
   Compressible fluctuating Navier--Stokes with all variables cell-centred:
   fluxes, conversion between conserved and primitive variables, boundary
   conditions, membranes and statistics. It is controlled by the
   ``compressible`` namespace.

``src_compressible_stag``
   Compressible fluctuating Navier--Stokes with momentum on the staggered
   grid, advanced with a three-stage Runge--Kutta scheme
   (``RK3stepStag``). It adds particle reservoirs, membranes, turbulence
   forcing and spectral output, and provides the ``main_driver`` used by
   ``exec/compressible_stag``.

``src_reactDiff``
   Stochastic reaction--diffusion of species number densities, with several
   splitting and integration choices (``temporal_integrator``,
   ``reactDiff_diffusion_type``, ``reaction_type``).

``src_chemistry``
   Bulk chemical reactions (stoichiometry, rate constants) for the
   compressible and reaction--diffusion solvers (``chemistry`` namespace).

``src_MFsurfchem``
   Mean-field surface chemistry (adsorption and desorption at a wall),
   coupled to the compressible staggered solver (``MFsurfchem``
   namespace).

``src_electrostatic``
   Poisson solve for the electric potential and field from a charge density
   (``esSolve``), using AMReX's ``MLPoisson`` multigrid solver.

Particles and immersed boundaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``src_particles``
   AMReX particle containers for FHD particle methods
   (``FhdParticleContainer``) and DSMC (``DsmcParticleContainer``), and the
   ``dsmc`` namespace. Background on AMReX particle containers is in
   :ref:`amrex:sec:Particles:ParticleContainer`.

``src_geometry``
   Parametric planes (``paramPlane``) used as walls and surfaces by the
   particle codes.

``src_immersed-boundary``
   Immersed-boundary markers, flagella and colloids coupled to the fluid
   solver (``IBCore``, ``IBMarkerContainer``, ``IBGMRES``). Its parameters
   are read with a Fortran namelist (``src_F90/immbdy_namelist.F90``)
   rather than ``ParmParse``, so they are not part of Chapter
   :ref:`chap:params`.

Multiscale coupling
~~~~~~~~~~~~~~~~~~~

``src_mui``
   A bundled copy of MUI (Multiscale Universal Interface), a header-only
   library for coupling concurrently running solvers (see
   ``src_mui/README.md`` for the reference and licence).

``src_surfchem_mui`` and ``src_surfchem_amrex``
   Surface-chemistry coupling between FHDeX and an external particle code,
   through MUI or through AMReX's MPMD interface respectively.

Other directories
~~~~~~~~~~~~~~~~~

``src_deankow``
   Dean--Kawasaki (fluctuating particle density) solvers. Each subdirectory
   (``singlelevel``, ``multilevel``, ``ensemble``, ``helfrich_membrane``,
   ...) is a self-contained program with its own ``main.cpp``, built from
   the matching directory under ``exec/dean_kow/``.

``unmaintained/``
   Older applications kept for reference. They are not guaranteed to build
   with the current libraries.
