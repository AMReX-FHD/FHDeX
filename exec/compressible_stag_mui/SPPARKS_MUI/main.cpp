/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "spparks.h"
#include "input.h"

#if defined(MUI)
#include "mui.h"
#include "lib_mpi_split.h"
#elif defined(USE_AMREX_MPMD)
#include <AMReX.H>
#include <AMReX_MPMD.H>
#include <AMReX_ParmParse.H>
#endif

using namespace SPPARKS_NS;

/* ----------------------------------------------------------------------
   main program to drive SPPARKS
------------------------------------------------------------------------- */

int main(int argc, char **argv)
{
#if defined(MUI)

  MPI_Comm comm = mui::mpi_split_by_app(argc,argv);
  mui::uniface2d uniface( "mpi://KMC-side/FHD-KMC-coupling" );
#ifdef SLURM
  // expected args for MUI-SLURM: exec %t %o inputs_file ...
  // since %t and %o were already used by mpi_split_by_app
  // we remove them from args list
  argc -= 2;
  for (int i=1;i<argc;i++) argv[i] = argv[i+2];
#endif
  SPPARKS *spk = new SPPARKS(argc,argv,comm);
  spk->uniface = &uniface;

#elif defined(USE_AMREX_MPMD)

  MPI_Comm comm = amrex::MPMD::Initialize(argc, argv);
  {
      amrex::ParmParse pp("amrex");
      pp.add("the_arena_init_size", 0);
      pp.add("verbose", 0);
  }
  amrex::Initialize(comm);
  SPPARKS *spk = new SPPARKS(argc,argv,comm);

#else

  // non-mui, non-amerx
  MPI_Init(&argc,&argv);
  SPPARKS *spk = new SPPARKS(argc,argv,MPI_COMM_WORLD);

#endif

  spk->input->file();
  delete spk;

#if defined(USE_AMREX_MPMD)
  amrex::Finalize();
  amrex::MPMD::Finalize();
#else
  MPI_Finalize();
#endif
}
