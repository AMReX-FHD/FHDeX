.. _sec:apps:immersed-ions:

``exec/immersedIons`` -- Brownian ions in a fluctuating fluid
=============================================================

Purpose
   Simulates charged (or neutral) Brownian particles immersed in a
   fluctuating Stokes fluid, with electrostatics, short-range interactions,
   walls and external fields. Typical uses are ionic conductivity,
   electro-osmosis, diffusion and structure factors of electrolytes.

Libraries and build
   ``src_particles``, ``src_hydro``, ``src_gmres``, ``src_electrostatic``,
   ``src_geometry``, ``src_immersed-boundary`` (base kernels),
   ``src_analysis``, ``src_rng``, ``src_common``; built with
   ``USE_PARTICLES = TRUE``, ``USE_FFT = TRUE``, ``DIM = 3``.

Algorithm
   Each time step:

   #. computes short-range (Lennard-Jones) forces and the P3M short-range
      correction (when ``sr_tog != 0`` or ``es_tog = 3``);
   #. spreads charges to the grid and solves the Poisson equation
      (``es_tog`` = 1 or 3), or computes pairwise Coulomb forces
      (``es_tog = 2``);
   #. spreads particle forces to the fluid with the chosen kernel;
   #. solves the fluctuating Stokes equations (``fluid_tog = 1``);
   #. moves the ions with the interpolated fluid velocity plus their dry
      Brownian motion (``move_tog``);
   #. accumulates statistics, distributions, structure factors and output.

Wet and dry diffusion
   The kernel width sets each ion's hydrodynamic ("wet") radius: about
   0.91, 1.255 and 1.481 grid spacings for ``pkernel_fluid`` = 3, 4 and 6.
   For exponential-of-semicircle kernels (``eskernel_fluid``), the radius is
   tabulated in ``wetRad_wN.dat`` files. The total diffusion coefficient
   comes from Stokes--Einstein if ``diameter > 0``, and from ``diff``
   otherwise. The remainder after subtracting the wet part is applied as dry
   diffusion. ``all_dry = 1`` makes all diffusion dry and requires
   ``fluid_tog = 0``.

Key parameters
   .. list-table::
      :header-rows: 1
      :widths: 30 70

      * - Parameter
        - Meaning
      * - ``fluid_tog``
        - 0 no fluid, 1 Stokes (2, low Mach, is not supported).
      * - ``es_tog``
        - 0 off, 1 Poisson, 2 pairwise Coulomb (serial only), 3 P3M.
      * - ``move_tog``
        - Total particle move: 0 off, 1 single step, 2 midpoint.
      * - ``dry_move_tog``
        - Dry particle move: 0 off, 1 on, 2 on without the dry stochastic
          drift correction.
      * - ``sr_tog``, ``eepsilon``, ``sigma``, ``rmin``, ``rmax``
        - Short-range Lennard-Jones interactions between particles
          (``*_wall`` variants for walls).
      * - ``wall_mob``
        - Dry mobility correction near walls: 0 none, 1 infinite plane,
          2 alternative model.
      * - ``pkernel_fluid``, ``pkernel_es``, ``eskernel_fluid``,
          ``eskernel_beta``
        - Kernels for fluid coupling and charge spreading.
      * - ``qval``, ``permittivity``, ``images``
        - Charges, permittivity, and the number of periodic images for
          pairwise Coulomb interactions.
      * - ``eamp``, ``efreq``, ``ephase``
        - Amplitude, frequency and phase of an external electric field.
      * - ``bc_es_lo/hi``, ``potential_lo/hi``
        - Potential boundary conditions: -1 periodic, 1 Dirichlet,
          2 Neumann.
      * - ``particle_count``, ``particle_n0``, ``particle_neff``
        - Number or number density of ions per species.
      * - ``particle_placement``
        - 1 reads ``x y z pinned`` for each particle from ``particles.dat``;
          any other value places particles randomly.
      * - ``msd_int``, ``msd_len``, ``radialdist_int``, ``cartdist_int``
        - Mean-square displacement and pair-distribution diagnostics.

Example inputs
   ``inputs_template`` documents every option. ``regression_inputs/``
   contains conductivity, electro-osmosis and boundary-condition tests.
   ``paper_inputs/`` and ``paper_inputs_2/`` reproduce published results;
   the README maps them to their cases. ``inputs_gpu`` is set up for GPU
   runs.

Output
   Fluid plotfiles ``plt*`` (and ``stagx/y/z*``), particle and charge
   plotfiles ``cplt*``, ``eplt*`` and ``parplt*``, and structure factors
   ``plt_SF_charge`` and ``plt_SF_vel``. Particle diagnostics are written as
   text files, e.g. ``msdEst_<i>``, ``radialDistribution*``,
   ``cartesianDistribution*``, ``currentEst`` and ``bulkFlowEst``.
   ``plot_ascii = 1`` adds ASCII copies of mean fields.
