.. _sec:apps:compressible-stag:

``exec/compressible_stag`` -- compressible staggered-grid FHD
=============================================================

Purpose
   Solves the compressible fluctuating Navier--Stokes equations for an
   ideal-gas mixture of up to ``MAX_SPEC`` species, with momentum on the
   staggered grid. Optional features are bulk chemistry
   (``src_chemistry``), mean-field surface chemistry (``src_MFsurfchem``),
   particle reservoirs, membranes and turbulence forcing.

Libraries
   ``src_compressible_stag`` (which provides ``main_driver``),
   ``src_compressible``, ``src_chemistry``, ``src_MFsurfchem``,
   ``src_rng``, ``src_analysis``, ``src_common``.

Build options
   ``MAX_SPEC`` (default 8) and ``MAX_REAC`` (default 5) bound the number of
   species and reactions. ``DO_TURB = TRUE`` compiles in turbulence forcing.
   The scripts ``build_perlmutter.sh`` and ``build_frontier.sh`` show GPU
   builds on NERSC Perlmutter (CUDA) and OLCF Frontier (HIP).

Algorithm
   Explicit three-stage Runge--Kutta time stepping (``RK3stepStag`` in
   ``timeStepStag.cpp``). The state is stored as conserved variables
   ``cu`` (density, momentum, total energy, species densities) and derived
   primitive variables ``prim``. The driver checks that

   * ``nvars = AMREX_SPACEDIM + 2 + nspecies``, and
   * ``nprimvars = AMREX_SPACEDIM + 3 + 2*nspecies``

   and aborts otherwise. Both must therefore be set consistently in the
   inputs file (e.g. ``nvars = 9``, ``nprimvars = 14`` for four species in
   3D), together with ``ngc = 2 2 2`` ghost cells.

Key parameters
   .. list-table::
      :header-rows: 1
      :widths: 30 70

      * - Parameter
        - Meaning
      * - ``prob_type``
        - Initial condition, e.g. 1 = constant concentration,
          2 = Rayleigh--Taylor, 3 = diffusion barrier. Further cases
          (4, 5 and 100--111) are defined in ``src_compressible_stag``.
      * - ``algorithm_type``
        - 1 = single component, 2 = multispecies.
      * - ``visc_type``
        - 1 = non-symmetric, 2 = symmetric, 3 = symmetric with bulk
          viscosity (1D/2D runs require 1 or 2).
      * - ``advection_type``
        - 1 = interpolate primitive, 2 = interpolate conserved quantities.
      * - ``nspecies``, ``molmass``, ``diameter``, ``dof``, ``hcv``,
          ``hcp``
        - Species properties (CGS units in the examples).
      * - ``k_B``, ``Runiv``, ``T_init``, ``rho0``
        - Physical constants and initial state.
      * - ``bc_vel_lo/hi``, ``bc_mass_lo/hi``, ``bc_therm_lo/hi``
        - Velocity (-1 periodic, 1 slip, 2 no-slip), mass (1 wall,
          2 concentration, 3 reservoir) and thermal (1 adiabatic,
          2 isothermal with ``t_lo/hi``) boundary conditions. Periodic
          velocity forces all boundary conditions in that direction to be
          periodic.
      * - ``transport_type``
        - Transport-coefficient model (1 = Giovangigli, the default);
          see :ref:`params:compressible`.
      * - ``plot_means``, ``plot_vars``, ``plot_covars``, ``plot_cross``,
          ``cross_cell``
        - Statistics written to plotfiles.

Example inputs
   ``inputs_regression_*`` cover periodic, slip/adiabatic, temperature
   gradient, concentration and reservoir boundaries in 1D and 3D. The
   ``test_*`` subdirectories contain chemistry and surface-chemistry test
   cases.
