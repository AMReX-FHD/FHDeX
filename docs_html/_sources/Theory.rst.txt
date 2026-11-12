.. _chap:theory:

Theory and Numerical Methods
============================

This chapter summarises the equations FHDeX solves and the main numerical
ingredients shared by its solvers. Derivations and analysis are in the method
papers listed in :ref:`sec:theory:references`.

Fluctuating hydrodynamics
-------------------------

Fluctuating hydrodynamics augments each dissipative flux of the continuum
equations with a stochastic flux. The noise amplitude is set by the
fluctuation--dissipation balance, so that at equilibrium the fields fluctuate
with the variances predicted by statistical mechanics. A stochastic flux is
written as :math:`\sqrt{2 L k_B T}\,\mathcal{W}`, where :math:`L` is the
corresponding transport coefficient (viscosity, diffusivity, conductivity) and
:math:`\mathcal{W}` is spatio-temporal white noise.

On a grid with cell volume :math:`\Delta V` and time step :math:`\Delta t`,
white noise is approximated by independent standard Gaussian random numbers
:math:`W` scaled by :math:`1/\sqrt{\Delta V\,\Delta t}`. In 2D, the cell
volume includes the out-of-plane thickness ``cell_depth``. The code scales
every stochastic flux by ``variance_coef_mom`` (momentum) or
``variance_coef_mass`` (mass), so setting these to 0 recovers the
deterministic equations.

Incompressible flow
-------------------

The incompressible fluctuating Navier--Stokes equations for velocity
:math:`\mathbf{v}` and pressure :math:`\pi` are

.. math::

   \rho\,\partial_t \mathbf{v} + \nabla\pi
     = -\rho\,\nabla\cdot(\mathbf{v}\mathbf{v}^T) + \nabla\cdot\left(\eta\,\overline{\nabla}\mathbf{v}\right)
       + \nabla\cdot\boldsymbol{\Sigma} + \mathbf{f},
   \qquad \nabla\cdot\mathbf{v} = 0 ,

with stochastic stress :math:`\boldsymbol{\Sigma} = \sqrt{2\eta k_B T}\,\boldsymbol{\mathcal{W}}`.
The form of the viscous operator is selected with ``visc_type``. With
``stoch_stress_form = 1`` (symmetric) the noise tensor is
:math:`(\boldsymbol{\mathcal{W}}+\boldsymbol{\mathcal{W}}^T)/\sqrt{2}`. The
code generates it directly: diagonal entries get variance 2 and
off-diagonal entries variance 1 (``src_hydro/StochMomFlux.cpp``).
``stoch_stress_form = 0`` uses a non-symmetric tensor.

Spatial discretisation
   FHDeX uses a staggered (MAC) grid. Scalars and pressure live at cell
   centres, each velocity component on the faces normal to it, and the
   off-diagonal stochastic stresses on cell edges (nodes in 2D). Discrete
   divergence and gradient are adjoint (:math:`D = -G^T`), which gives the
   discrete fluctuation--dissipation balance that makes equilibrium
   fluctuations correct at all wavenumbers. In the incompressible solvers,
   advection is either centred or uses the BDS scheme (``advection_type``).

Temporal integration
   ``exec/hydro`` uses a predictor--corrector scheme. Viscous terms are
   treated semi-implicitly (Crank--Nicolson) and advection explicitly, with
   trapezoidal averaging in the corrector. Each stage requires the solution
   of a Stokes saddle-point problem (below). ``algorithm_type`` selects
   between the inertial equations and the overdamped (Stokes) limit, where
   the inertia term is dropped.

The Stokes solver
-----------------

Implicit treatment of viscosity together with the incompressibility
constraint leads, at every stage, to the linear system

.. math::

   \begin{pmatrix} A & G \\ -D & 0 \end{pmatrix}
   \begin{pmatrix} \mathbf{v} \\ \pi \end{pmatrix}
   = \begin{pmatrix} \mathbf{b}_u \\ b_p \end{pmatrix},
   \qquad A = \theta_\alpha\,\rho/\Delta t - L_\eta ,

where :math:`L_\eta` is the discrete viscous operator and
:math:`\theta_\alpha` switches the time-derivative (mass) term on or off.
``src_gmres`` solves it with restarted, preconditioned GMRES
(``gmres_max_inner`` sets the restart length). The preconditioner
(``precon_type``) is a block preconditioner:

* ``1`` -- projection preconditioner (default); ``-1`` adds a more expensive
  pressure update;
* ``2``/``-2`` -- lower block-triangular;
* ``3``/``-3`` -- upper block-triangular;
* ``4``/``-4`` -- block diagonal.

The velocity block is approximately inverted by a staggered geometric
multigrid V-cycle (``StagMGSolver``; ``stag_mg_*`` parameters, Jacobi or
coloured Gauss--Seidel smoothing via ``stag_mg_smoother``). The pressure
block uses a cell-centred Poisson solve (MAC projection, ``mg_*``
parameters) built on the AMReX MLMG solvers.

Multispecies mixtures
---------------------

``src_multispec`` evolves the species densities :math:`\rho_k` of an
:math:`N`-species mixture together with the velocity, in the low-Mach
(isothermal, incompressible-like) limit:

.. math::

   \partial_t \rho_k = -\nabla\cdot(\rho_k\mathbf{v})
      - \nabla\cdot\left(\mathbf{F}_k + \widetilde{\mathbf{F}}_k\right),
   \qquad \sum_k \frac{\rho_k}{\bar\rho_k} = 1 ,

where :math:`\mathbf{F}_k` are the deterministic diffusive fluxes from the
Maxwell--Stefan formulation (coefficients ``Dbar``, thermodynamic factor
:math:`\Gamma`, with :math:`\Gamma = I` for ``is_ideal_mixture = 1``) and
:math:`\widetilde{\mathbf{F}}_k` are the stochastic mass fluxes. The
constraint on the pure-component densities ``rhobar`` replaces
:math:`\nabla\cdot\mathbf{v} = 0` and makes the velocity divergence depend on
the diffusive fluxes. Correction fluxes keep :math:`\sum_k \rho_k` consistent
to round-off (``correct_flux``).

Optional physics:

* **Charged species** (``use_charged_fluid``): electrodiffusive fluxes and a
  Lorentz force, with the electric potential from a Poisson equation.
  ``electroneutral = 1`` enforces electroneutrality instead.
* **Multiphase / capillarity** (``use_multiphase``, ``kc_tension``): a
  reversible stress for surface tension (``ComputeDivReversibleStress``).
* **Gravity** (``grav``), optionally in the Boussinesq approximation
  (``algorithm_type = 6``).

The integrators in ``exec/multispec`` are predictor--corrector schemes:
``algorithm_type = 0`` (inertial) and ``algorithm_type = 6`` (Boussinesq).
The stochastic mass flux at the midpoint can be interpreted in the Ito or
Stratonovich sense (``midpoint_stoch_mass_flux_type``), and the densities on
faces used to weight the noise are averaged as set by ``avg_type``.

Compressible flow
-----------------

The compressible solvers evolve the conserved variables

.. math::

   \mathbf{U} = \left(\rho,\; \mathbf{j} = \rho\mathbf{v},\; \rho E,\; \rho_1,\dots,\rho_N\right)

of an ideal-gas mixture,

.. math::

   \partial_t \mathbf{U} = -\nabla\cdot\left(\mathbf{F}_\text{adv} + \mathbf{F}_\text{diss} + \widetilde{\mathbf{F}}\right) + \mathbf{S},

with dissipative fluxes for viscous stress, heat conduction and species
diffusion, and the corresponding stochastic fluxes. The primitive variables
(velocity, temperature, pressure, mass and mole fractions) are recomputed
from :math:`\mathbf{U}` after every stage. Transport coefficients for the
mixture come from kinetic theory (``transport_type``, Giovangigli's
formulation by default).

``src_compressible_stag`` stores momentum on cell faces. It advances
:math:`\mathbf{U}` with the explicit three-stage, third-order Runge--Kutta
scheme of Delong et al. (``RK3stepStag``). Two independent white-noise fields
:math:`W_A, W_B` are generated per step and combined in each stage as
:math:`W_A + w_2 W_B`, with

.. math::

   w_2 = \frac{2\sqrt2 + \sqrt3}{5},\qquad
   \frac{-4\sqrt2 + 3\sqrt3}{5},\qquad
   \frac{\sqrt2 - 2\sqrt3}{10}

in stages 1, 2 and 3. This weighting makes the scheme weakly second-order
accurate for the stochastic terms.

Chemistry and reaction--diffusion
---------------------------------

Chemical reactions are specified by stoichiometric coefficients and rate
constants (``src_chemistry``). The number of reaction events in a cell
during :math:`\Delta t` is sampled according to ``reaction_type``:

* ``0`` -- deterministic rates;
* ``1`` -- chemical Langevin equation (Gaussian noise);
* ``2`` -- SSA (exact stochastic simulation within the cell);
* ``3`` -- tau leaping (Poisson-distributed event counts).

``src_reactDiff`` couples these reaction samplers to diffusion of number
densities. ``temporal_integrator`` selects split or unsplit schemes, and
``reactDiff_diffusion_type`` the diffusion integrator for split schemes:
explicit trapezoidal predictor/corrector (0), Crank--Nicolson (1), explicit
midpoint (2), multinomial diffusion (3) or forward Euler (4). Multinomial
diffusion moves integer numbers of molecules between neighbouring cells and
keeps populations non-negative.

Particle methods
----------------

DSMC
   Direct simulation Monte Carlo represents a dilute gas by simulated
   particles, each standing for ``particle_neff`` real molecules. Every
   step alternates stochastic binary collisions within collision cells with
   free streaming and surface interactions; reservoir surfaces inject
   particles (see :ref:`sec:apps:dsmc`).

Immersed ions
   Ions are Brownian particles coupled to the fluctuating Stokes equations.
   Forces on each ion (electrostatic, short-range, external) are spread to
   the fluid with a regularised kernel (Peskin kernels of width 3, 4 or 6, or
   exponential-of-semicircle kernels), and the fluid velocity is interpolated
   back with the adjoint kernel. Each ion's diffusion is split into a "wet"
   part, carried by the resolved fluctuating flow, and a "dry" part, added as
   independent Brownian motion. Electrostatics use a Poisson solve with a
   P3M short-range correction, or direct pairwise Coulomb sums (see
   :ref:`sec:apps:immersed-ions`).

Random numbers
--------------

All stochastic terms draw from the AMReX random-number generators
(:ref:`amrex:sec:basics:random`) through ``src_rng``. The ``seed``
parameter controls reproducibility: a positive value gives a fixed,
reproducible sequence (offset by MPI rank), and ``0`` seeds from the clock.

.. _sec:theory:references:

References
----------

The methods implemented in FHDeX are described in, among others:

* A. Donev, E. Vanden-Eijnden, A. L. Garcia and J. B. Bell, *On the accuracy
  of finite-volume schemes for fluctuating hydrodynamics*, Commun. Appl.
  Math. Comput. Sci. 5, 149 (2010).
* F. Balboa Usabiaga, J. B. Bell, R. Delgado-Buscalioni, A. Donev,
  T. G. Fai, B. E. Griffith and C. S. Peskin, *Staggered schemes for
  fluctuating hydrodynamics*, Multiscale Model. Simul. 10, 1369 (2012).
* S. Delong, B. E. Griffith, E. Vanden-Eijnden and A. Donev, *Temporal
  integrators for fluctuating hydrodynamics*, Phys. Rev. E 87, 033302 (2013).
* M. Cai, A. J. Nonaka, J. B. Bell, B. E. Griffith and A. Donev, *Efficient
  variable-coefficient finite-volume Stokes solvers*, Commun. Comput. Phys.
  16, 1263 (2014).
* K. Balakrishnan, A. L. Garcia, A. Donev and J. B. Bell, *Fluctuating
  hydrodynamics of multispecies nonreactive mixtures*, Phys. Rev. E 89,
  013017 (2014).
* C. Kim, A. Nonaka, J. B. Bell, A. L. Garcia and A. Donev, *Stochastic
  simulation of reaction-diffusion systems: A fluctuating-hydrodynamics
  approach*, J. Chem. Phys. 146, 124110 (2017).
