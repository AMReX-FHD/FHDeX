#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"

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

  Box domain_bx(geom.Domain());
  
  MultiFab AMREX_D_DECL(xmf,ymf,zmf);

  // these face-centered masks are 1 on all faces, except for faces at
  // grid boundaries where they equal 2
  // we use this since we want half-weighting toward the sum at all grid boundaries
  AMREX_D_TERM(auto xmask = m1[0].OverlapMask(Periodicity(domain_bx.size()));,
               auto ymask = m1[1].OverlapMask(Periodicity(domain_bx.size()));,
               auto zmask = m1[2].OverlapMask(Periodicity(domain_bx.size()));)
  
  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<Real> reduce_data(reduce_op);
  using ReduceTuple = typename decltype(reduce_data)::Type;

  for (MFIter mfi(m1[0]); mfi.isValid(); ++mfi)
  {
      const Box& bx = mfi.validbox();
      auto const& fab = m1[0].array(mfi);
      auto const& msk = xmask->array(mfi);
      reduce_op.eval(bx, reduce_data,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
          return {fab(i,j,k) / msk(i,j,k)};
      });
  }

  sum[0] = amrex::get<0>(reduce_data.value());
  ParallelDescriptor::ReduceRealSum(sum[0]);

  for (MFIter mfi(m1[1]); mfi.isValid(); ++mfi)
  {
      const Box& bx = mfi.validbox();
      auto const& fab = m1[1].array(mfi);
      auto const& msk = ymask->array(mfi);
      reduce_op.eval(bx, reduce_data,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
          return {fab(i,j,k) / msk(i,j,k)};
      });
  }

  sum[1] = amrex::get<0>(reduce_data.value());
  ParallelDescriptor::ReduceRealSum(sum[1]);

#if (AMREX_SPACEDIM == 3)


  for (MFIter mfi(m1[2]); mfi.isValid(); ++mfi)
  {
      const Box& bx = mfi.validbox();
      auto const& fab = m1[2].array(mfi);
      auto const& msk = zmask->array(mfi);
      reduce_op.eval(bx, reduce_data,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
      {
          return {fab(i,j,k) / msk(i,j,k)};
      });
  }

  sum[2] = amrex::get<0>(reduce_data.value());
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
                   amrex::Vector<amrex::Real>& prod_val)
{
  BL_PROFILE_VAR("StagInnerProd()",StagInnerProd);

  std::array<MultiFab, AMREX_SPACEDIM> prod_temp;

  DistributionMapping dmap = m1[0].DistributionMap();
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    prod_temp[d].define(m1[d].boxArray(), dmap, 1, 0);
    MultiFab::Copy(prod_temp[d],m1[d],comp1,0,1,0);
    MultiFab::Multiply(prod_temp[d],m2[d],comp2,0,1,0);
  }
  
  std::fill(prod_val.begin(), prod_val.end(), 0.);
  SumStag(geom,prod_temp,0,prod_val,false);
}

void CCInnerProd(const amrex::MultiFab& m1,
		 const int& comp1,
		 const amrex::MultiFab& m2,
		 const int& comp2,
		 amrex::Real& prod_val)
{
  
  BL_PROFILE_VAR("CCInnerProd()",CCInnerProd);

  amrex::MultiFab prod_temp;
  prod_temp.define(m1.boxArray(), m1.DistributionMap(), 1, 0);

  MultiFab::Copy(prod_temp,m1,comp1,0,1,0);
  MultiFab::Multiply(prod_temp,m2,comp2,0,1,0);
  
  prod_val = 0.;
  SumCC(prod_temp,0,prod_val,false);
}

void StagL2Norm(const Geometry& geom,
                const std::array<MultiFab, AMREX_SPACEDIM>& m1,
		const int& comp,
		Real& norm_l2)
{

    BL_PROFILE_VAR("StagL2Norm()",StagL2Norm);

    Vector<Real> inner_prod(AMREX_SPACEDIM);
    StagInnerProd(geom, m1, comp, m1, comp, inner_prod);
    norm_l2 = sqrt(std::accumulate(inner_prod.begin(), inner_prod.end(), 0.));
}

void CCL2Norm(const amrex::MultiFab& m1,
	      const int& comp,
	      amrex::Real& norm_l2)
{
  
  BL_PROFILE_VAR("CCL2Norm()",CCL2Norm);

  norm_l2 = 0.;
  CCInnerProd(m1,comp,m1,comp,norm_l2);
  norm_l2 = sqrt(norm_l2);
}
