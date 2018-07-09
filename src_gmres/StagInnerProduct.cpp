#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"

#include "AMReX_ArrayLim.H"
#include "AMReX_ParallelDescriptor.H"

using namespace amrex;

void StagInnerProd(const std::array<MultiFab, AMREX_SPACEDIM>& m1,
                   const int& comp1, 
                   const std::array<MultiFab, AMREX_SPACEDIM>& m2,
                   const int& comp2,
                   Vector<Real>& prod_val)
{
  std::fill(prod_val.begin(), prod_val.end(), 0.);

  // Loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
  for (MFIter mfi(m1[0]); mfi.isValid(); ++mfi) {

    // const Box& validBox_cc = enclosedCells(mfi.validbox());

    stag_inner_prod(BL_TO_FORTRAN_N_ANYD(m1[0][mfi],comp1),
		    BL_TO_FORTRAN_N_ANYD(m1[1][mfi],comp1),
#if (AMREX_SPACEDIM == 3)
		    BL_TO_FORTRAN_N_ANYD(m1[2][mfi],comp1),
#endif
		    BL_TO_FORTRAN_N_ANYD(m2[0][mfi],comp2),
		    BL_TO_FORTRAN_N_ANYD(m2[1][mfi],comp2),
#if (AMREX_SPACEDIM == 3)
		    BL_TO_FORTRAN_N_ANYD(m2[2][mfi],comp2),
#endif
		    &prod_val);

  }

  //////////////// Parellel Test /////////////////
  // int myproc = ParallelDescriptor::MyProc();  // Return the rank  
  // int nprocs = ParallelDescriptor::NProcs();  // Return the number of processes
  // ParallelDescriptor::Barrier();
  ////////////////////////////////////////////////

  for (int d = 0; d < AMREX_SPACEDIM; d++) {
      ParallelDescriptor::ReduceRealSum(prod_val.dataPtr(),AMREX_SPACEDIM);
  }

  //////////////// Debug ////////////////////////
  // Print() << prod_val[0];

  // for (MFIter mfi(m1[0]); mfi.isValid(); ++mfi) {
  //   Print() << prod_val << "\n";
  // }

  // double *prod_val_temp = &prod_val[0];
  // Print() << *prod_val_temp << "\n";

  // double prod_val_temp;
  // for (int d = 0; d < AMREX_SPACEDIM; d++) {
  //   prod_val_temp = prod_val[d];
  //   Print() << prod_val_temp;
  // }
  ///////////////////////////////////////////////
  

  //////////// TEST EXAMPLE /////////////////
  // Vector<Real> vect_test(AMREX_SPACEDIM);
  // double *vect_test_pntr;
  // vect_test[0] = 0.5;
  // vect_test[1] = 0.77;
  // vect_test_pntr = &vect_test[0];
  // Print() << *vect_test_pntr << "\n";
  // for (int d = 0; d < AMREX_SPACEDIM; d++) {
  //   Print() << vect_test[d] << "\n";
  // }
  //////////////////////////////////////////
}
