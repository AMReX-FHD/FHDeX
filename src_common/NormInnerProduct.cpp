#include "common_functions.H"

void SumStag(const std::array<MultiFab, AMREX_SPACEDIM>& m1,
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
    const Box& bx_grid = mfi.validbox();

    auto const& fab = m1[0].array(mfi);

    int xlo = bx_grid.smallEnd(0);
    int xhi = bx_grid.bigEnd(0);

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
    const Box& bx_grid = mfi.validbox();

    auto const& fab = m1[1].array(mfi);

    int ylo = bx_grid.smallEnd(1);
    int yhi = bx_grid.bigEnd(1);

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
    const Box& bx_grid = mfi.validbox();

    auto const& fab = m1[2].array(mfi);

    int zlo = bx_grid.smallEnd(2);
    int zhi = bx_grid.bigEnd(2);

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

  // divide by the number of cells
  if (divide_by_ncells == true) {
    BoxArray ba_temp = m1[0].boxArray();
    ba_temp.enclosedCells();
    long numpts = ba_temp.numPts();
    for (int d=0; d<AMREX_SPACEDIM; d++) {
      sum[d] = sum[d]/(double)(numpts);
    }
  }
}

void SumEdge(const std::array<MultiFab, NUM_EDGE>& m1,
  amrex::Vector<amrex::Real>& sum,
  const bool& divide_by_ncells)
{
  BL_PROFILE_VAR("SumEdge()",SumEdge);

  // Initialize to zero
  std::fill(sum.begin(), sum.end(), 0.);

  ReduceOps<ReduceOpSum> reduce_op;

  //////// xy-edges

  ReduceData<Real> reduce_dataxy(reduce_op);
  using ReduceTuple = typename decltype(reduce_dataxy)::Type;

  for (MFIter mfi(m1[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    const Box& bx_grid = mfi.validbox();

    auto const& fab = m1[0].array(mfi);

    int xlo = bx_grid.smallEnd(0);
    int xhi = bx_grid.bigEnd(0);

    int ylo = bx_grid.smallEnd(1);
    int yhi = bx_grid.bigEnd(1);

    reduce_op.eval(bx, reduce_dataxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
    {
      Real weight;
      if (i>xlo && i<xhi && j>ylo && j<yhi) {
        weight = 1.;
      } else if ( (i==xlo || i==xhi) && (j==ylo || j==yhi) ) {
        weight = 0.25;
      } else {
        weight = 0.5;
      }
      return {fab(i,j,k)*weight};
    });
  }

  sum[0] = amrex::get<0>(reduce_dataxy.value());
  ParallelDescriptor::ReduceRealSum(sum[0]);

  if (AMREX_SPACEDIM == 2) {
    return;
  }

  //////// xz-edges

  ReduceData<Real> reduce_dataxz(reduce_op);
  using ReduceTuple = typename decltype(reduce_dataxz)::Type;

  for (MFIter mfi(m1[1],TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    const Box& bx_grid = mfi.validbox();

    auto const& fab = m1[1].array(mfi);

    int xlo = bx_grid.smallEnd(0);
    int xhi = bx_grid.bigEnd(0);

    int zlo = bx_grid.smallEnd(2);
    int zhi = bx_grid.bigEnd(2);

    reduce_op.eval(bx, reduce_dataxz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
    {
      Real weight;
      if (i>xlo && i<xhi && k>zlo && k<zhi) {
        weight = 1.;
      } else if ( (i==xlo || i==xhi) && (k==zlo || k==zhi) ) {
        weight = 0.25;
      } else {
        weight = 0.5;
      }
      return {fab(i,j,k)*weight};
    });
  }

  sum[1] = amrex::get<0>(reduce_dataxz.value());
  ParallelDescriptor::ReduceRealSum(sum[1]);

  //////// yz-edges

  ReduceData<Real> reduce_datayz(reduce_op);
  using ReduceTuple = typename decltype(reduce_datayz)::Type;

  for (MFIter mfi(m1[2],TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    const Box& bx_grid = mfi.validbox();

    auto const& fab = m1[2].array(mfi);

    int ylo = bx_grid.smallEnd(1);
    int yhi = bx_grid.bigEnd(1);

    int zlo = bx_grid.smallEnd(2);
    int zhi = bx_grid.bigEnd(2);

    reduce_op.eval(bx, reduce_datayz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
    {
      Real weight;
      if (j>ylo && j<yhi && k>zlo && k<zhi) {
        weight = 1.;
      } else if ( (j==ylo || j==yhi) && (k==zlo || k==zhi) ) {
        weight = 0.25;
      } else {
        weight = 0.5;
      }
      return {fab(i,j,k)*weight};
    });
  }

  sum[2] = amrex::get<0>(reduce_datayz.value());
  ParallelDescriptor::ReduceRealSum(sum[2]);

  // divide by the number of cells
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

void StagInnerProd(const std::array<MultiFab, AMREX_SPACEDIM>& m1,
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
  SumStag(mscr,prod_val,false);
}

void EdgeInnerProd(const std::array<MultiFab, NUM_EDGE>& m1,
  const int& comp1,
  const std::array<MultiFab, NUM_EDGE>& m2,
  const int& comp2,
  std::array<MultiFab, NUM_EDGE>& mscr,
  amrex::Vector<amrex::Real>& prod_val)
{
  BL_PROFILE_VAR("EdgeInnerProd()",EdgeInnerProd);

  for (int d=0; d<NUM_EDGE; d++) {
    MultiFab::Copy(mscr[d],m1[d],comp1,0,1,0);
    MultiFab::Multiply(mscr[d],m2[d],comp2,0,1,0);
  }

  std::fill(prod_val.begin(), prod_val.end(), 0.);
  SumEdge(mscr,prod_val,false);
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

void CCMoments(const amrex::MultiFab& m1,
  const int& comp1,
  amrex::MultiFab& mscr,
  const int& power,
  amrex::Real& prod_val)
{

  BL_PROFILE_VAR("CCMoments()",CCMoments);

  MultiFab::Copy(mscr,m1,comp1,0,1,0);
  for(int i=1; i<power; i++){
    MultiFab::Multiply(mscr,m1,comp1,0,1,0);
  }

  prod_val = 0.;
  SumCC(mscr,0,prod_val,false);
}

void FCMoments(const std::array<MultiFab, AMREX_SPACEDIM>& m1,
  const amrex::Vector<int>& comps,
  std::array<MultiFab, AMREX_SPACEDIM>&  mscr,
  const int& power,
  amrex::Vector<amrex::Real>& prod_val)
{

  BL_PROFILE_VAR("FCMoments()",FCMoments);

  if (comps.size() != AMREX_SPACEDIM) amrex::Abort("FCMoments:: Vector of comps needs to same size as AMREX_SPACEDIM");
  if (prod_val.size() != AMREX_SPACEDIM) amrex::Abort("FCMoments:: Vector of prod_val needs to same size as AMREX_SPACEDIM");

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    MultiFab::Copy(mscr[d],m1[d],comps[d],0,1,0);
    for(int i=1; i<power; i++){
      MultiFab::Multiply(mscr[d],m1[d],comps[d],0,1,0);
    }
  }
  SumStag(mscr,prod_val,false);
}


void StagL2Norm(const std::array<MultiFab, AMREX_SPACEDIM>& m1,
  const int& comp,
  std::array<MultiFab, AMREX_SPACEDIM>& mscr,
  Real& norm_l2)
{

  BL_PROFILE_VAR("StagL2Norm()",StagL2Norm);

  Vector<Real> inner_prod(AMREX_SPACEDIM);

  StagInnerProd(m1, comp, m1, comp, mscr, inner_prod);
  norm_l2 = sqrt(std::accumulate(inner_prod.begin(), inner_prod.end(), 0.));
}

void EdgeL2Norm(const std::array<MultiFab, NUM_EDGE>& m1,
  const int& comp,
  std::array<MultiFab, NUM_EDGE>& mscr,
  Real& norm_l2)
{

  BL_PROFILE_VAR("EdgeL2Norm()",EdgeL2Norm);

  Vector<Real> inner_prod(NUM_EDGE);

  EdgeInnerProd(m1, comp, m1, comp, mscr, inner_prod);
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