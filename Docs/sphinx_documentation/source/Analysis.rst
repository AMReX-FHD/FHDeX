.. _chap:analysis:

Analysis and Output
===================

Fluctuating-hydrodynamics results are statistical, so most analysis happens
while the simulation runs: FHDeX accumulates means, variances and spectra
over many time steps and writes them alongside the instantaneous fields.

Plotfiles
---------

All applications write AMReX plotfiles (:ref:`amrex:sec:IO`) every
``plot_int`` steps, named ``<plot_base_name><step>`` with ``plot_base_name``
defaulting to ``plt``. Staggered (face-centred) fields are written to
separate plotfiles (``stagx``, ``stagy``, ``stagz``) when ``plot_stag = 1``.
Particle applications add particle plotfiles; see the application pages for
the exact names.

Plotfiles can be read by VisIt, ParaView, yt and AMReX's ``amrvis``
(:ref:`amrex:sec:visit`), and compared with the AMReX post-processing tools
(:ref:`amrex:sec:postprocessing`).

Structure factors
-----------------

The static structure factor :math:`S(\mathbf{k}) = \langle \hat{\delta a}(\mathbf{k})\,\hat{\delta b}^*(\mathbf{k}) \rangle`
measures the spectrum of fluctuations and is the main diagnostic for
verifying fluctuating-hydrodynamics solvers against theory. The
``StructFact`` class in ``src_analysis`` computes it with FFTs:

.. code-block:: c++

   StructFact structFact;
   // var_names: one name per field; var_scaling: normalisation per pair
   structFact.define(ba, dmap, var_names, var_scaling);
   // or only selected pairs (a_i, b_i):
   structFact.define(ba, dmap, var_names, var_scaling, s_pairA, s_pairB);

   // every struct_fact_int steps after n_steps_skip:
   structFact.FortStructure(fields);        // add one sample
   structFact.WritePlotFile(step, time, "plt_SF");

Each call to ``FortStructure`` adds one snapshot to a running average. By
default, the class computes all auto- and cross-correlations of the
variables (:math:`N(N+1)/2` pairs for :math:`N` variables), and
``var_scaling`` must have one entry per pair. ``WritePlotFile`` writes two
plotfiles,
``<name>_mag<step>`` (magnitude) and ``<name>_real_imag<step>`` (real and
imaginary parts), with the zero wavenumber shifted to the centre of the
domain.

The relevant inputs are ``struct_fact_int`` (sampling interval; ``<= 0``
disables structure factors) and ``n_steps_skip`` (equilibration steps to
skip before sampling). Some applications also compute structure factors of
quantities averaged along one direction ("flattened", e.g.
``plt_SF_Flattened`` in ``exec/hydro``), or in two slabs separated by
``cross_cell`` (``do_slab_sf``).

Turbulence spectra
   ``IntegratekShells`` and the functions in ``TurbSpectra`` integrate
   spectra over spherical shells in wavenumber space and write the
   shell-averaged kinetic-energy and scalar spectra :math:`E(k)` for
   turbulence studies.

Means, variances and correlations
---------------------------------

``src_common`` provides running statistics (``ComputeBasicStats``),
spatial means and variances (``ComputeSpatialMean``,
``ComputeSpatialVariance``) and plane averages (``WriteHorizontalAverage``).
The compressible staggered solver writes means, variances, covariances and
spatial correlations with the point ``cross_cell`` to its plotfiles, as
selected by ``plot_means``, ``plot_vars``, ``plot_covars`` and
``plot_cross``. Statistics are reset after ``n_steps_skip`` steps so that
transients are excluded.

Checkpoint and restart
----------------------

Checkpoints are written every ``chk_int`` steps to
``<chk_base_name><step>`` (default prefix ``chk``). To restart, set
``restart`` to the step number of the checkpoint. Particle applications
also use ``particle_restart``. Some applications (e.g. ``exec/DSMC``)
discard the accumulated statistics on restart when ``reset_stats = 1``. The hydro checkpoint also stores
the random-number state of each MPI rank, so a restarted run continues the
same random sequence.

Python tools
------------

``python/amrex_to_vtp.py``
   Reads AMReX particle output (``plt*/particle0/Header`` and data files)
   and writes ParaView ``.vtp`` files.

``python/write_pview_file.py``
   Helper functions for writing ParaView files (structured grids and
   particles), used by ``amrex_to_vtp.py``.

``python/convertSeries``
   Shell script that runs ``amrex_to_vtp.py`` on every ``parplt*``
   particle plotfile of an ``exec/immersedIons`` run and collects the
   results as a numbered ``.vtp`` series for ParaView animation.

Several applications also include their own analysis scripts, e.g.
``exec/hydro/PDF`` and ``exec/compressible_stag/TURB_PDFS``.
