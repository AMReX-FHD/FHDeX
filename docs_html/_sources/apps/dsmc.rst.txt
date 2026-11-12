.. _sec:apps:dsmc:

``exec/DSMC`` -- direct simulation Monte Carlo
==============================================

Purpose
   Particle simulation of dilute gas mixtures with direct simulation Monte
   Carlo, including reservoirs and arbitrary planar surfaces. It serves as a
   molecular-level reference for the continuum FHD solvers.

Libraries and build
   ``src_particles`` (through ``MakeDSMC.package``), ``src_geometry``,
   ``src_analysis``, ``src_rng``, ``src_common``; built with
   ``USE_PARTICLES = TRUE``, ``DIM = 3`` and ``MAX_SPEC = 8``.
   Application sources (collisions, particle generation, statistics,
   plotfiles) are in ``exec/DSMC``. Particle data use AMReX particle
   containers (:ref:`amrex:sec:Particles:ParticleContainer`).

Algorithm
   After building the surfaces and initialising particles and collision
   cells, each time step

   #. selects and performs binary collisions in every collision cell
      (``CalcSelections``, ``CollideParticles``);
   #. injects particles at reservoir surfaces (``Source``);
   #. moves particles and applies surface interactions
      (``MoveParticlesCPP``);
   #. accumulates statistics (``EvaluateStats``), which are reset after
      ``n_steps_skip`` steps.

Species
   Each species has a ``mass``, ``diameter`` and ``T_init``. The number of
   simulated particles per species comes from ``particle_count`` if it is
   non-negative, otherwise from the volume fraction ``phi_domain``,
   otherwise from ``rho0`` and the mass fractions ``Yk0``.
   ``particle_neff`` is the number of real molecules per simulated particle.
   With ``particle_input > 0``, initial positions are read from
   ``particles.dat``.

Surfaces
   The six domain faces are always surfaces. Their behaviour follows from the
   boundary-condition parameters (``bc_vel_*``, ``bc_therm_*``,
   ``bc_mass_*``, ``t_lo/hi``, ``n_lo/hi``, ...); for example, reservoir
   faces inject particles at a prescribed state. Additional surfaces are read
   from a file ``paramplanes.dat`` in the run directory, if present. Its first
   entry is the number of extra surfaces, followed by one record per surface:

   .. code-block:: text

      x0 y0 z0          origin
      ux uy uz          first in-plane direction
      vx vy vz          second in-plane direction
      uTop vTop         extent along u and v
      rnx rny rnz       right-side normal
      lnx lny lnz       left-side normal
      porosity specularity temperature density[nspecies] source sink momentumCons   (right side)
      porosity specularity temperature density[nspecies] source sink momentumCons   (left side)
      periodicity area

   The notebook ``exec/DSMC/py/paramPlane.ipynb`` generates this file, and
   ``exec/DSMC/py/paramplanes.dat`` is an example.

Example inputs
   ``inputs_conc``: two species under a concentration gradient between
   mass-conserving reservoirs (``bc_vel = 7`` in x). ``inputs_conc_gpu`` is
   the GPU variant, and ``test_inputs/`` holds periodic-equilibrium, thermal
   and concentration tests.

Output
   ``pltcu*`` (conserved variables, instantaneous and mean), ``pltprim*``
   (primitive variables per species: number density, density, velocity,
   temperature, pressure, energy, concentration), ``pltvar*`` (variances),
   ``pltpart*`` (particles) and ``plt_SF_prim*`` (structure factors).

Related applications
   ``exec/phononDSMC`` and ``exec/phononDetector`` apply the same particle and
   surface machinery to phonon transport. There, all surfaces come from
   ``paramplanes.dat``, with no automatic domain faces.
