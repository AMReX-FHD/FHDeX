FHDeX Documentation
===================

FHDeX is a collection of solvers for **fluctuating hydrodynamics** (FHD) and
related particle methods, built on the `AMReX <https://amrex-codes.github.io/amrex/>`_
software framework. It simulates fluids at mesoscopic scales, where thermal
fluctuations matter: incompressible and compressible flows, multispecies
mixtures, chemical reactions, charged species, and fluids coupled to particles
and immersed boundaries. All solvers run in parallel with MPI and on GPUs.

.. list-table::
   :widths: 50 50
   :class: borderless

   * - **Source code**
     - https://github.com/AMReX-FHD/FHDeX
   * - **AMReX documentation**
     - https://amrex-codes.github.io/amrex/docs_html/

This guide explains how to build and run the FHDeX applications, how the
source code is organised, and which runtime parameters they accept.
FHDeX relies heavily on AMReX, so this guide does not repeat material that the
`AMReX documentation <https://amrex-codes.github.io/amrex/docs_html/>`_
already covers. Instead it links to the relevant AMReX section wherever an
AMReX concept (``MultiFab``, ``ParmParse``, plotfiles, the GNU Make build
system and so on) comes up.

New users should read :ref:`chap:getting-started` first and then the page for
the application they want to run under :ref:`chap:applications`.

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: User Guide

   Introduction
   GettingStarted
   Theory
   CodeStructure
   Applications
   Analysis
   RuntimeParameters
   DeveloperGuide

.. toctree::
   :maxdepth: 1
   :caption: Project

   Contributing
