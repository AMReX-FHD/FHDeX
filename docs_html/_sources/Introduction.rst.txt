.. _chap:introduction:

Introduction
============

What is fluctuating hydrodynamics?
----------------------------------

At mesoscopic scales (from nanometres to micrometres) thermal fluctuations
change how fluids behave. They drive diffusion enhancement, giant
concentration fluctuations, the onset of instabilities and spontaneous
pattern formation. Fluctuating hydrodynamics models these effects by adding
stochastic fluxes to the usual conservation laws of continuum fluid
dynamics. The noise amplitude is fixed by the fluctuation--dissipation
theorem, so that the equations reproduce the correct equilibrium
fluctuations.

The simplest example is the incompressible fluctuating Stokes / Navier--Stokes
system for the velocity :math:`\mathbf{v}` and pressure :math:`\pi`,

.. math::

   \rho\,\partial_t \mathbf{v} + \nabla\pi
     = -\rho\,\nabla\cdot(\mathbf{v}\mathbf{v}^T)
       + \eta\nabla^2\mathbf{v}
       + \nabla\cdot\left(\sqrt{2\eta k_B T}\;\boldsymbol{\mathcal{W}}\right),
   \qquad \nabla\cdot\mathbf{v} = 0,

where :math:`\boldsymbol{\mathcal{W}}` is a white-noise random tensor field.
FHDeX solves this system and extends it to multispecies mixtures, compressible
flows, chemical reactions, charged species, and coupling to particles and
immersed boundaries.

The equations and numerical methods of all solver families are described in
:ref:`chap:theory`.

Solver families
---------------

FHDeX contains several groups of solvers. Each group is a set of source
libraries (``src_*``) together with one or more applications (``exec/*``):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Family
     - Description
   * - Incompressible hydrodynamics
     - Staggered-grid (MAC) fluctuating Navier--Stokes with a GMRES Stokes
       solver and multigrid preconditioner (``src_hydro``, ``src_gmres``;
       application ``exec/hydro``).
   * - Low-Mach multispecies
     - Multispecies mixtures with stochastic mass fluxes, optional
       electrodiffusion for charged species, inertial and Boussinesq time
       integrators (``src_multispec``; application ``exec/multispec``).
   * - Compressible
     - Compressible fluctuating Navier--Stokes for multispecies gas mixtures,
       in a collocated (``src_compressible``) and a staggered
       (``src_compressible_stag``) formulation, with optional chemistry and
       surface chemistry.
   * - Reaction--diffusion
     - Stochastic reaction--diffusion of species number densities, with
       deterministic, chemical-Langevin, SSA or tau-leaping reactions and
       explicit, semi-implicit or multinomial diffusion
       (``src_reactDiff``; application ``exec/reactDiff``).
   * - Particle methods
     - DSMC, Brownian / immersed ions and immersed-boundary markers built on
       AMReX particle containers (``src_particles``,
       ``src_immersed-boundary``).
   * - Multiscale coupling
     - Coupling of continuum solvers to external particle codes through the
       MUI library (``src_mui``, ``src_surfchem_mui``).

Relationship to AMReX
---------------------

FHDeX applications are AMReX applications. All grid data are stored in AMReX
``MultiFab`` objects, staggered velocities are face-centred (nodal in one
direction) ``MultiFab`` objects, parallelism uses AMReX's MPI domain
decomposition and GPU offloading, and input parameters are read with
``ParmParse``. Readers new to AMReX should first read:

* :ref:`amrex:sec:basics:box`, :ref:`amrex:sec:basics:ba` and
  :ref:`amrex:sec:basics:dm` -- how the domain is split into boxes and
  distributed over MPI ranks;
* :ref:`amrex:sec:basics:multifab` and :ref:`amrex:sec:basics:mfiter` -- the
  data containers and loops used in every FHDeX kernel;
* :ref:`amrex:sec:basics:parmparse` -- the input-file format;
* :ref:`amrex:sec:build:make` -- the GNU Make build system used by every
  ``exec/*`` directory.

Repository layout
-----------------

.. code-block:: text

   FHDeX/
     src_*/          source libraries, one per physics component (see Code Structure)
     exec/           applications: one directory per executable, each with a
                     GNUmakefile and a set of inputs_* files (see Applications)
     python/         post-processing scripts (AMReX particle output -> ParaView)
     tools/          template for a local Path.local (sets AMREX_HOME)
     unmaintained/   old applications kept for reference; not documented here
     Docs/           this guide
