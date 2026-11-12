.. _chap:contributing:

Contributing to This Documentation
==================================

The documentation sources live in the FHDeX repository under
``Docs/sphinx_documentation/``. The web site and the PDF are both built from
these sources.

Editing a page
--------------

Each page is a reStructuredText file in
``Docs/sphinx_documentation/source/``. Application pages are in
``source/apps/``. The "Edit on GitHub" link at the top of every page opens
the corresponding file. Cross-references into the AMReX documentation use
the ``amrex:`` prefix, e.g. ``:ref:`amrex:sec:basics:multifab```.

Runtime parameters
------------------

The tables in :ref:`chap:params` are generated from the source code by
``Docs/sphinx_documentation/gen_params.py`` on every build and must not be
edited by hand. To improve a description, edit the comment next to the
parameter's declaration in the ``src_*/<name>_namespace.H`` header.

Building locally
----------------

.. code-block:: bash

   pip install sphinx sphinx_rtd_theme
   cd Docs/sphinx_documentation
   make html              # web site in build/html/
   make latexpdf          # PDF in build/latex/fhdex.pdf (needs LaTeX)

Preview the web site with ``python3 -m http.server -d build/html`` and open
http://localhost:8000.

Publishing
----------

Publishing is automatic. On every push to ``main`` of
``AMReX-FHD/FHDeX``, the GitHub Action ``.github/workflows/docs.yml`` runs
``make html`` and deploys ``build/html`` to the ``docs_html/`` folder of the
repository's ``gh-pages`` branch, from where GitHub Pages serves it. Pull
requests also build the documentation, so a broken page shows up as a failed
check before merging.
