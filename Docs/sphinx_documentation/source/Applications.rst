.. _chap:applications:

Applications
============

Each directory under ``exec/`` builds one executable. The pages below
document the main applications; :ref:`sec:apps:overview` lists the
remaining ones.

Every application directory follows the same conventions:

* ``GNUmakefile`` -- build options and the list of ``src_*`` libraries
  (:ref:`chap:getting-started`);
* ``inputs_*`` -- example inputs files. Files named ``inputs_regression_*``
  are short runs (a few time steps) that serve as regression tests and are a
  good starting point for new users;
* optionally, application-specific sources such as ``main_driver.cpp``,
  ``Init.cpp`` (initial conditions selected by ``prob_type``) and
  ``WritePlotFile.cpp``. These override the library versions as described
  in :ref:`chap:code-structure`.

.. toctree::
   :maxdepth: 1

   apps/hydro
   apps/multispec
   apps/compressible_stag
   apps/reactdiff
   apps/dsmc
   apps/immersed_ions

.. _sec:apps:overview:

Other applications
------------------

The following applications are also part of ``exec/``. Each follows the
conventions above; consult its ``inputs_*`` files and ``main_driver.cpp``.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Directory
     - Topic
   * - ``compressible``
     - Collocated compressible FHD (``src_compressible``)
   * - ``compressible_mui``, ``compressible_stag_mui``
     - Compressible FHD coupled to an external particle code through MUI
   * - ``DSMC_granular``, ``phononDSMC``, ``phononDetector``
     - Direct simulation Monte Carlo (``src_particles``, ``src_geometry``)
   * - ``immersed_boundary``, ``cellbdytest_new``
     - Particles and immersed boundaries coupled to the fluid
   * - ``dean_kow``
     - Dean--Kawasaki solvers (``src_deankow``)
   * - ``FPU_FHD_cc``, ``FPU_FHD_stag``, ``FPU_chain``, ``StochasticHeat``,
       ``_1d_dk_hybrid``, ``_1d_staggered_comp``, ``thinFilm``
     - Specialised and one-dimensional models
   * - ``chemistry_testing``, ``structFactTest``, ``Ek_calculator``,
       ``tests``, ``tutorials``, ``lib_hydro``
     - Tests, analysis utilities and tutorials
