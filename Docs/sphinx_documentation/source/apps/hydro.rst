.. _sec:apps:hydro:

``exec/hydro`` -- incompressible fluctuating hydrodynamics
==========================================================

Purpose
   Solves the incompressible fluctuating Navier--Stokes (or, in the
   overdamped limit, Stokes) equations on a staggered grid. It is the
   simplest FHDeX application and the best one to start with.

Libraries
   ``src_hydro``, ``src_gmres``, ``src_rng``, ``src_analysis``,
   ``src_common``; ``src_multispec`` is also compiled in.

Algorithm
   Each time step is a predictor--corrector scheme. Viscous terms are
   treated semi-implicitly (Crank--Nicolson) and advection is treated
   explicitly with a trapezoidal correction. Both stages solve a Stokes
   saddle-point system for velocity and pressure with preconditioned GMRES
   from ``src_gmres``. ``algorithm_type = 0`` selects the inertial scheme;
   other values select the overdamped (Stokes) limit, in which the
   deterministic case needs only a single solve per step. Thermal
   fluctuations enter through the stochastic momentum flux
   (``StochMomFlux``), scaled by ``variance_coef_mom``.

Key parameters
   .. list-table::
      :header-rows: 1
      :widths: 30 70

      * - Parameter
        - Meaning
      * - ``prob_type``
        - Initial condition: 0 = fluid at rest, 1 = vortex, 2 and 3 =
          Kelvin--Helmholtz (sine / smooth); see ``exec/hydro/Init.cpp``.
      * - ``visc_coef``, ``visc_type``
        - Viscosity and form of the viscous operator (1 = :math:`\nabla\cdot\beta\nabla`,
          2 = symmetric, 3 = symmetric with bulk viscosity; negative values
          allow spatially varying coefficients).
      * - ``variance_coef_mom``, ``initial_variance_mom``
        - Amplitude of the stochastic momentum flux and of the initial
          velocity fluctuations (0 = deterministic).
      * - ``k_B``, ``T_init``
        - Boltzmann constant and temperature, which set the noise strength.
      * - ``bc_vel_lo/hi``
        - Velocity boundary conditions per direction: -1 periodic, 1 slip,
          2 no-slip.
      * - ``struct_fact_int``, ``n_steps_skip``
        - Accumulate structure factors every ``struct_fact_int`` steps
          after skipping ``n_steps_skip`` steps of equilibration.
      * - ``gmres_*``, ``stag_mg_*``, ``mg_*``
        - Linear solver controls (:ref:`params:gmres`).

Example inputs
   ``inputs_regression_vortex_2d/3d`` (deterministic vortex),
   ``inputs_regression_equil_2d/3d`` and ``inputs_incomp_equil_fluct_*``
   (equilibrium fluctuations), ``inputs_shearwall_*`` (shear flow),
   ``inputs_turbtest_*`` (forced turbulence).

Convergence tests
   ``run_convergence_tests.sh`` runs coarse, medium and fine versions of
   ``inputs_dtrmnistic_convergence_*d``, and ``compute_convergence_rate.sh``
   compares them with the AMReX tool ``DiffSameDomainRefinedStag``. The
   procedure is described in ``exec/hydro/README.CONVERGENCE``.
