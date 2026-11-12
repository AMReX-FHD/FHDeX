.. _sec:apps:reactdiff:

``exec/reactDiff`` -- stochastic reaction--diffusion
=====================================================

Purpose
   Simulates reaction--diffusion systems of species number densities with
   fluctuations. It implements the fluctuating-hydrodynamics approach to
   reaction--diffusion described by Kim et al., J. Chem. Phys. 146, 124110
   (2017). The inputs files reproduce figures from that paper.

Libraries
   ``src_reactDiff`` (which provides ``main_driver``), ``src_chemistry``,
   ``src_analysis``, ``src_rng``, ``src_common``.

Algorithm
   Reactions and diffusion are combined either by operator splitting or in
   a single unsplit step, selected by ``temporal_integrator``:

   .. list-table::
      :header-rows: 1
      :widths: 15 85

      * - Value
        - Scheme
      * - ``0``
        - D + R (first-order splitting)
      * - ``1``
        - (1/2)R + D + (1/2)R (Strang splitting, option 1)
      * - ``2``
        - (1/2)D + R + (1/2)D (Strang splitting, option 2)
      * - ``-1``
        - unsplit forward Euler
      * - ``-2``
        - unsplit explicit midpoint
      * - ``-3``
        - unsplit multinomial diffusion
      * - ``-4``
        - unsplit implicit midpoint

   For split schemes, ``reactDiff_diffusion_type`` chooses the diffusion
   integrator and ``reactDiff_reaction_type`` the order of the reaction step;
   ``reaction_type`` selects how reactions are sampled (0 deterministic,
   1 CLE, 2 SSA, 3 tau leaping). See :ref:`chap:theory`.

Specifying reactions
   Reactions are given by their reactant and product stoichiometries, one
   pair of lines per reaction, with one entry per species. For example, the
   reaction :math:`U + W \to V + W` in a three-species system is

   .. code-block:: text

      nspecies  = 3
      nreaction = 7
      stoich_1R = 1 0 1     # reactants of reaction 1: U + W
      stoich_1P = 0 1 1     # products of reaction 1:  V + W
      ...
      rate_const = 0.0002 0.0002 1. 0.03666663 4.44444555555 0.00333333 16.66665
      rate_multiplier = 1.
      include_discrete_LMA_correction = 1

   Rates follow the law of mass action. ``rate_multiplier`` scales all of
   them. With ``include_discrete_LMA_correction = 1``, :math:`n^2` and
   :math:`n^3` in the rate expressions become :math:`n(n-1/\Delta V)` and
   :math:`n(n-1/\Delta V)(n-2/\Delta V)`.

Key parameters
   .. list-table::
      :header-rows: 1
      :widths: 30 70

      * - Parameter
        - Meaning
      * - ``D_Fick``
        - Fickian diffusion coefficient for each species.
      * - ``n_init_in_1``, ``n_init_in_2``
        - Initial number densities, interpreted according to ``prob_type``.
      * - ``integer_populations``
        - Start from integer numbers of molecules per cell.
      * - ``variance_coef_mass``, ``initial_variance_mass``
        - Stochastic diffusion on/off (1/0) and initial fluctuations.
      * - ``avg_type``
        - How densities on faces are averaged to weight the stochastic flux
          (1 arithmetic, 2 geometric, 3 harmonic, 10--12 smoothed
          Heaviside variants).
      * - ``bc_mass_lo/hi``
        - -1 periodic, 1 wall (Neumann), 2 reservoir (Dirichlet).
      * - ``cell_depth``
        - Out-of-plane thickness in 2D; sets the cell volume and hence the
          noise strength.

Example inputs
   ``inputs_paper_BPM_2d`` (Baras--Pearson--Mansour limit cycle, Figures
   5--7 of the paper), ``inputs_paper_Lemarchand_3d`` (Figure 8), and
   ``test_Schlogl_2d/`` (Schlögl model, with job scripts).

Output
   Plotfiles ``plt*`` with the number densities of all species every
   ``plot_int`` steps; structure factors when ``struct_fact_int > 0``.
