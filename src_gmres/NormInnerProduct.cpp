#include "gmres_functions.H"

#include "AMReX_ArrayLim.H"
#include "AMReX_Box.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParallelDescriptor.H"

using namespace amrex;

void SumStag(const Geometry& geom,
             const std::array<MultiFab, AMREX_SPACEDIM>& m1,
	     const int& comp,
	     amrex::Vector<amrex::Real>& sum,
	     const bool& divide_by_ncells)
{
  BL_PROFILE_VAR("SumStag()",SumStag);

  // Initialize to zero
  std::fill(sum.begin(), sum.end(), 0.);

  ReduceOps<ReduceOpSum> reduce_op;

  //////// x-faces

  ReduceData<Real> reduce_datax(reduce_op);
  using ReduceTuple = typename decltype(reduce_datax)::Type;

  for (MFIter mfi(m1[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
      const Box& bx = mfi.tilebox();
      auto const& fab = m1[0].array(mfi);

      int xlo = bx.smallEnd(0);
      int xhi = bx.bigEnd(0);

      reduce_op.eval(bx, reduce_datax,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
          Real weight = (i>xlo && i<xhi) ? 1.0 : 0.5;
          return {fab(i,j,k)*weight};
      });
  }

  sum[0] = amrex::get<0>(reduce_datax.value());
  ParallelDescriptor::ReduceRealSum(sum[0]);

  //////// y-faces

  ReduceData<Real> reduce_datay(reduce_op);

  for (MFIter mfi(m1[1],TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
      const Box& bx = mfi.tilebox();
      auto const& fab = m1[1].array(mfi);

      int ylo = bx.smallEnd(1);
      int yhi = bx.bigEnd(1);

      reduce_op.eval(bx, reduce_datay,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
          Real weight = (j>ylo && j<yhi) ? 1.0 : 0.5;
          return {fab(i,j,k)*weight};
      });
  }

  sum[1] = amrex::get<0>(reduce_datay.value());
  ParallelDescriptor::ReduceRealSum(sum[1]);

#if (AMREX_SPACEDIM == 3)

  //////// z-faces

  ReduceData<Real> reduce_dataz(reduce_op);

  for (MFIter mfi(m1[2],TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
      const Box& bx = mfi.tilebox();
      auto const& fab = m1[2].array(mfi);

      int zlo = bx.smallEnd(2);
      int zhi = bx.bigEnd(2);

      reduce_op.eval(bx, reduce_dataz,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
          Real weight = (k>zlo && k<zhi) ? 1.0 : 0.5;
          return {fab(i,j,k)*weight};
      });
  }

  sum[2] = amrex::get<0>(reduce_dataz.value());
  ParallelDescriptor::ReduceRealSum(sum[2]);

#endif

  if (divide_by_ncells == true) {
    BoxArray ba_temp = m1[0].boxArray();
    ba_temp.enclosedCells();
    long numpts = ba_temp.numPts();
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      sum[d] = sum[d]/(double)(numpts);
    }
  }
}

void SumCC(const amrex::MultiFab& m1,
	   const int& comp,
	   amrex::Real& sum,
	   const bool& divide_by_ncells)
{
  BL_PROFILE_VAR("SumCC()",SumCC);

  sum = 0.;
  sum = m1.MultiFab::sum(comp, false);

  if (divide_by_ncells == 1) {
    BoxArray ba_temp = m1.boxArray();
    long numpts = ba_temp.numPts();
    sum = sum/(double)(numpts);
  }
}

void StagInnerProd(const Geometry& geom,
                   const std::array<MultiFab, AMREX_SPACEDIM>& m1,
                   const int& comp1,
                   const std::array<MultiFab, AMREX_SPACEDIM>& m2,
                   const int& comp2,
                   std::array<MultiFab, AMREX_SPACEDIM>& mscr,
                   amrex::Vector<amrex::Real>& prod_val)
{
  BL_PROFILE_VAR("StagInnerProd()",StagInnerProd);

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(mscr[d],m1[d],comp1,0,1,0);
    MultiFab::Multiply(mscr[d],m2[d],comp2,0,1,0);
  }

  std::fill(prod_val.begin(), prod_val.end(), 0.);
  SumStag(geom,mscr,0,prod_val,false);
}

void CCInnerProd(const amrex::MultiFab& m1,
		 const int& comp1,
		 const amrex::MultiFab& m2,
		 const int& comp2,
                 amrex::MultiFab& mscr,
		 amrex::Real& prod_val)
{

  BL_PROFILE_VAR("CCInnerProd()",CCInnerProd);

  MultiFab::Copy(mscr,m1,comp1,0,1,0);
  MultiFab::Multiply(mscr,m2,comp2,0,1,0);

  prod_val = 0.;
  SumCC(mscr,0,prod_val,false);
}

void StagL2Norm(const Geometry& geom,
                const std::array<MultiFab, AMREX_SPACEDIM>& m1,
		const int& comp,
                std::array<MultiFab, AMREX_SPACEDIM>& mscr,
		Real& norm_l2)
{

    BL_PROFILE_VAR("StagL2Norm()",StagL2Norm);

    Vector<Real> inner_prod(AMREX_SPACEDIM);

    StagInnerProd(geom, m1, comp, m1, comp, mscr, inner_prod);
    norm_l2 = sqrt(std::accumulate(inner_prod.begin(), inner_prod.end(), 0.));
}

void CCL2Norm(const amrex::MultiFab& m1,
	      const int& comp,
              amrex::MultiFab& mscr,
	      amrex::Real& norm_l2)
{

  BL_PROFILE_VAR("CCL2Norm()",CCL2Norm);

  norm_l2 = 0.;
  CCInnerProd(m1,comp,m1,comp,mscr,norm_l2);
  norm_l2 = sqrt(norm_l2);
}
