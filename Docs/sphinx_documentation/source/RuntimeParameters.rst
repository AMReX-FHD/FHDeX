.. _chap:params:

Runtime Parameters
==================

FHDeX reads its runtime parameters with AMReX ``ParmParse``
(:ref:`amrex:sec:basics:parmparse`). They can be given in the inputs file or
on the command line (``./main3d.gnu.MPI.ex inputs max_step=100``). The
command line takes precedence. Parameters are not prefixed: write
``n_cells = 64 64``, not ``common.n_cells``.

The tables below are **generated from the source code** by
``Docs/sphinx_documentation/gen_params.py`` every time the guide is built:

* the parameter list and the parameter names come from the
  ``pp.query``/``pp.get`` calls in each ``Initialize<Name>Namespace()``;
* the type comes from the ``extern`` declaration in
  ``src_*/<name>_namespace.H``. Arrays are written as *element × size*:
  "Real × species" is one value per species (sized by ``MAX_SPEC``),
  "Real × dim" one value per spatial direction, "species pairs" the
  ``nspecies*(nspecies-1)/2`` entries of a symmetric species matrix, and
  "list" a variable-length list;
* the default is the value assigned in the initializer before the parameter
  is read; "(all)" means every array element gets that value. *required* means the parameter is read with ``pp.get`` and the
  run aborts without it. "--" means no default is set in the initializer.
  The value then depends on the variable's initialisation elsewhere, or the
  parameter is only meaningful when given;
* the description is the comment next to the declaration or default.

An empty description means the source has no comment for that parameter.
The way to fill these gaps is to add a comment to the declaration in the
``*_namespace.H`` header; the next build of the guide picks it up.

Which namespaces an application uses depends on the
``Initialize*Namespace()`` calls in its ``main_driver``. ``common`` is always
used. Immersed-boundary parameters (``src_immersed-boundary``) are read
through a Fortran namelist and are not listed here.

.. toctree::
   :maxdepth: 1

   params/common
   params/gmres
   params/multispec
   params/compressible
   params/chemistry
   params/reactDiff
   params/MFsurfchem
   params/dsmc
   params/surfchem_mui
