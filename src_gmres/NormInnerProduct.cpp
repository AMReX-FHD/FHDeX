#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"

#include "AMReX_ArrayLim.H"
#include "AMReX_Box.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParallelDescriptor.H"

using namespace amrex;

void SumStag(const std::array<MultiFab, AMREX_SPACEDIM>& m1,
	     const int& comp,
	     amrex::Vector<amrex::Real>& sum,
	     const bool& divide_by_ncells)
{
  // Initialize to zero
  std::fill(sum.begin(), sum.end(), 0.);

  // Loop over boxes
  for (MFIter mfi(m1[0]); mfi.isValid(); ++mfi) {   
    
    // Create cell-centered box from semi-nodal box
    const Box& validBox_cc = enclosedCells(mfi.validbox());

    sum_stag(ARLIM_3D(validBox_cc.loVect()), ARLIM_3D(validBox_cc.hiVect()),
	     BL_TO_FORTRAN_N_ANYD(m1[0][mfi],comp),
	     BL_TO_FORTRAN_N_ANYD(m1[1][mfi],comp),
#if (AMREX_SPACEDIM == 3)
	     BL_TO_FORTRAN_N_ANYD(m1[2][mfi],comp),
#endif
	     sum.dataPtr());
  }

  ParallelDescriptor::ReduceRealSum(sum.dataPtr(),AMREX_SPACEDIM);

  if (divide_by_ncells == 1) {
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      BoxArray ba = m1[0].boxArray();
      ba.enclosedCells();
      long numpts = ba.numPts();
      sum[d] = sum[d]/(double)(numpts);
    }
  }

}

void SumCC(const amrex::MultiFab& m1,
	   const int& comp,
	   amrex::Real& sum,
	   const bool& divide_by_ncells)
{
  sum = 0.;
  sum = m1.MultiFab::sum(comp, false);

  if (divide_by_ncells == 1) {
    BoxArray ba = m1.boxArray();
    long numpts = ba.numPts();
    sum = sum/(double)(numpts);
  }

}

void StagInnerProd(const std::array<MultiFab, AMREX_SPACEDIM>& m1,
                   const int& comp1,
                   const std::array<MultiFab, AMREX_SPACEDIM>& m2,
                   const int& comp2,
                   amrex::Vector<amrex::Real>& prod_val)
{
  std::array<MultiFab, AMREX_SPACEDIM> m1_temp;
  std::array<MultiFab, AMREX_SPACEDIM> m2_temp;

  BoxArray ba;
  DistributionMapping dmap;
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    ba = m1[d].boxArray();
    dmap = m1[d].DistributionMap();
    m1_temp[d].define(ba, dmap, 1, 0);
    m2_temp[d].define(ba, dmap, 1, 0);
  }

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(m1_temp[d],m1[d],comp1,0,1,0);
    MultiFab::Copy(m2_temp[d],m2[d],comp2,0,1,0);
  }

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Multiply(m1_temp[d],m2_temp[d],0,0,1,0);
  }
  
  std::fill(prod_val.begin(), prod_val.end(), 0.);
  SumStag(m1_temp,0,prod_val,false);

}

void CCInnerProd(const amrex::MultiFab& m1,
		 const int& comp1,
		 const amrex::MultiFab& m2,
		 const int& comp2,
		 amrex::Real& prod_val)
{

  amrex::MultiFab m1_temp;
  amrex::MultiFab m2_temp;

  BoxArray ba = m1.boxArray();
  DistributionMapping dmap = m1.DistributionMap();
  m1_temp.define(ba, dmap, 1, 0);
  m2_temp.define(ba, dmap, 1, 0);

  MultiFab::Copy(m1_temp,m1,comp1,0,1,0);
  MultiFab::Copy(m2_temp,m2,comp2,0,1,0);
  MultiFab::Multiply(m1_temp,m2_temp,0,0,1,0);

  SumCC(m1_temp, 0, prod_val,false);

}

void StagL2Norm(const std::array<MultiFab, AMREX_SPACEDIM>& m1,
		const int& comp,
		amrex::Vector<amrex::Real>& norm_l2)
{
  StagInnerProd(m1, comp, m1, comp, norm_l2);
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    norm_l2[d] = sqrt(norm_l2[d]);
  }
}

void CCL2Norm(const amrex::MultiFab& m1,
	      const int& comp,
	      amrex::Real& norm_l2)
{
  CCInnerProd(m1,comp,m1,comp,norm_l2);
  norm_l2 = sqrt(norm_l2);
}
