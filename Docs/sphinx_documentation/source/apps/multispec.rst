.. _sec:apps:multispec:

``exec/multispec`` -- low-Mach multispecies mixtures
====================================================

Purpose
   Simulates multispecies liquid mixtures in the low-Mach limit,
   with stochastic mass and momentum fluxes. Optional features include
   charged species (electrodiffusion with a Poisson solve), multiphase
   mixtures with surface tension, and gravity. Typical uses are giant
   fluctuations, diffusion-driven instabilities, spinodal decomposition,
   droplets and bubbles, and ionic liquids.

Libraries
   ``src_multispec``, ``src_hydro``, ``src_gmres``, ``src_rng``,
   ``src_analysis``, ``src_common``.

Algorithm
   ``algorithm_type`` selects the time integrator implemented in the
   application directory:

   * ``0`` -- inertial integrator (``AdvanceTimestepInertial.cpp``);
   * ``6`` -- Boussinesq integrator (``AdvanceTimestepBousq.cpp``).

   Other values stop the run with *"algorithm_type not supported"*.

Key parameters
   .. list-table::
      :header-rows: 1
      :widths: 30 70

      * - Parameter
        - Meaning
      * - ``nspecies``, ``molmass``, ``rhobar``, ``rho0``
        - Number of species, mass per molecule (not molar mass), pure
          component densities and reference density.
      * - ``Dbar``
        - Maxwell--Stefan diffusion coefficients. They are given as the
          lower triangle of a symmetric matrix, row by row:
          ``D_12; D_13, D_23; D_14, D_24, D_34; ...``
          (``nspecies*(nspecies-1)/2`` values).
      * - ``is_ideal_mixture``
        - 1 (default) assumes an ideal mixture (thermodynamic factor
          :math:`\Gamma = I`); 0 uses the non-ideal model.
      * - ``c_init_1``, ``c_init_2``
        - Initial mass fractions, interpreted by ``prob_type`` in
          ``Init.cpp``.
      * - ``variance_coef_mass``, ``variance_coef_mom``
        - Amplitude of the stochastic mass and momentum fluxes.
      * - ``bc_mass_lo/hi``
        - -1 periodic, 1 wall, 2 reservoir.
      * - ``use_charged_fluid``, ``charge_per_mass``, ``dielectric_const``,
          ``bc_es_lo/hi``
        - Electrodiffusion for charged species and boundary conditions for
          the electric potential (1 Dirichlet, 2 Neumann).
      * - ``use_multiphase``, ``alpha_gex``, ``n_gex``, ``kc_tension``
        - Multiphase free energy and surface tension.

Example inputs
   ``inputs_regression_equil_*`` (equilibrium fluctuations),
   ``inputs_regression_detbubble_*`` (deterministic bubble),
   ``inputs_regression_spinodal_*`` (charged spinodal decomposition),
   ``inputs_regression_cap_2d``, and many production cases (``inputs_fh_*``
   giant fluctuations, ``inputs_droplet_*``, ``inputs_rtil_*`` ionic
   liquids).

.. note::

   ``MAX_SPEC`` in the ``GNUmakefile`` must be at least ``nspecies``, and
   ``MAX_ELEM`` must equal ``MAX_SPEC*(MAX_SPEC-1)/2``.
