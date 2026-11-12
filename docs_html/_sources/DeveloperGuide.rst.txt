.. _chap:developer:

Developer Guide
===============

This chapter collects the conventions to follow when extending FHDeX.
General AMReX programming (``MFIter`` loops, ``ParallelFor``, GPU memory) is
covered in the AMReX documentation; see :ref:`amrex:sec:basics:mfiter`,
:ref:`amrex:sec:basics:parallelfor` and :ref:`amrex:sec:gpu:overview`.

Adding a new application
------------------------

#. Create ``exec/<name>/`` and copy the ``GNUmakefile`` of the closest
   existing application.
#. List the ``src_*`` libraries you need, each with three lines:

   .. code-block:: make

      include ../../src_hydro/Make.package
      VPATH_LOCATIONS   += ../../src_hydro/
      INCLUDE_LOCATIONS += ../../src_hydro/

   Keep ``.`` first in ``VPATH_LOCATIONS``. Order the libraries from most to
   least specific, because Make uses the first file with a given name
   (:ref:`chap:code-structure`).
#. Write ``main_driver.cpp`` in the application directory. It must define
   ``void main_driver(const char* argv)``; ``src_common/main.cpp`` provides
   ``main``. Start by calling the ``Initialize*Namespace()`` functions for
   the parameter namespaces you use.
#. To change the behaviour of a single library file (initial conditions,
   plotfile contents, ...), copy it into the application directory under the
   same name and edit the copy. The library is left untouched.
#. Add an ``inputs_regression_*`` file that runs in a few seconds, so that
   others can check the application still works.

Adding a runtime parameter
--------------------------

#. Declare the variable in ``src_<lib>/<lib>_namespace.H`` with a short
   comment on the same line. The comment becomes its description in
   :ref:`chap:params`:

   .. code-block:: c++

      extern AMREX_GPU_MANAGED amrex::Real my_param;   // what it controls, units

#. Define it in ``<lib>_functions.cpp``.
#. In ``Initialize<Lib>Namespace()``, set its default **before** reading it,
   then read it with ``pp.query`` (optional) or ``pp.get`` (required):

   .. code-block:: c++

      my_param = 1.0;
      pp.query("my_param", my_param);

Following this pattern is all that is needed for the parameter to appear in
the generated documentation. Use ``AMREX_GPU_MANAGED`` for anything that GPU
kernels read, so it is accessible on the device without copies.

GPU support
-----------

FHDeX kernels use ``amrex::ParallelFor`` and AMReX ``Array4`` accessors, so
the same source runs on CPUs and GPUs. GPU builds set ``USE_CUDA = TRUE`` or
``USE_HIP = TRUE`` (:ref:`amrex:sec:gpu:build`). Example build scripts for
NERSC Perlmutter (CUDA) and OLCF Frontier (HIP) are in
``exec/compressible_stag/build_perlmutter.sh`` and ``build_frontier.sh``. In
GPU kernels:

* read parameters only through ``AMREX_GPU_MANAGED`` namespace variables or
  captured local copies;
* do not call host-only functions (``amrex::Print``, ``std::vector``);
* generate random numbers with the AMReX device RNG.

Single precision
----------------

The compressible solvers (``src_compressible``, ``src_compressible_stag``,
``src_chemistry``, ``src_MFsurfchem``) support single-precision builds with
the AMReX option ``PRECISION = FLOAT``. To keep this working:

* write floating-point literals as ``Real(...)`` (e.g. ``Real(0.5)``)
  wherever they meet ``Real`` arithmetic, so that they do not promote
  expressions to double;
* call math functions as ``std::sqrt``, ``std::exp``, ... so that the
  ``float`` overloads are selected, and cast ``std::pow`` exponents to
  ``Real``;
* avoid constants that underflow in single precision (use
  ``std::numeric_limits<Real>::min()`` instead of e.g. ``1.e-60``).

Double-precision builds are unaffected by these rules. Compile with
``-Wdouble-promotion -Wfloat-conversion`` to find violations.

Style
-----

* One ``.cpp`` file per routine or closely related group of routines, with
  declarations in the library's ``<lib>_functions.H``.
* Start every major routine with ``BL_PROFILE_VAR("Name()", Name);`` so that
  it shows up in AMReX profiling (``TINY_PROFILE = TRUE``).
* Use ``amrex::Abort`` with a clear message for invalid parameter
  combinations, as the existing drivers do (e.g. ``nvars`` checks in
  ``src_compressible_stag/main_driver.cpp``).
* Keep new applications out of ``unmaintained/``; move old ones there
  instead of deleting them.

Documentation
-------------

Every new application should get a page under
``Docs/sphinx_documentation/source/apps/``, following the structure of the
existing pages (purpose, libraries, algorithm, key parameters, example
inputs, output). See :ref:`chap:contributing` for how to build and publish
the documentation.
