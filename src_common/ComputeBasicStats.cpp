#include "common_functions.H"
#include "RunningMean.H"

#include <AMReX_ParReduce.H>

// Sum of component comp over the valid region, accumulated in double regardless of Real.
// MultiFab::sum accumulates in Real, which in single precision loses ~eps*sqrt(N) (or worse)
// relative accuracy over N cells and can overflow for large summands.
double SumDouble (const MultiFab& mf, int comp, bool local)
{
    BL_PROFILE_VAR("SumDouble()",SumDouble);

    auto const& ma = mf.const_arrays();
    double r = ParReduce(TypeList<ReduceOpSum>{}, TypeList<double>{}, mf, IntVect(0),
                         [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<double>
                         {
                             return { static_cast<double>(ma[b](i,j,k,comp)) };
                         });
    if (!local) {
        ParallelDescriptor::ReduceRealSum(r);
    }
    return r;
}

// Sum of (mf(comp) - shift)^power over the valid region; the power is formed in double inside
// the reduction, so e.g. 4th moments of large gradients do not overflow float.
double SumPowDouble (const MultiFab& mf, int comp, int power, double shift, bool local)
{
    BL_PROFILE_VAR("SumPowDouble()",SumPowDouble);

    auto const& ma = mf.const_arrays();
    double r = ParReduce(TypeList<ReduceOpSum>{}, TypeList<double>{}, mf, IntVect(0),
                         [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<double>
                         {
                             double x = static_cast<double>(ma[b](i,j,k,comp)) - shift;
                             double p = 1.0;
                             for (int n=0; n<power; ++n) { p *= x; }
                             return { p };
                         });
    if (!local) {
        ParallelDescriptor::ReduceRealSum(r);
    }
    return r;
}

Real ComputeSpatialMean(MultiFab& mf, const int& incomp)
{
    BL_PROFILE_VAR("ComputeSpatialMean()",ComputeSpatialMean);

    Long npts = mf.boxArray().numPts();

    Real average = Real(SumDouble(mf, incomp) / static_cast<double>(npts));

    return average;

}

Real ComputeSpatialVariance(MultiFab& mf, const int& incomp)
{
    BL_PROFILE_VAR("ComputeSpatialVariance()",ComputeSpatialVariance);

    Long npts = mf.boxArray().numPts();

    // two-pass variance, both passes accumulated in double
    double average = SumDouble(mf, incomp) / static_cast<double>(npts);

    Real variance = Real(SumPowDouble(mf, incomp, 2, average) / static_cast<double>(npts-1));

    return variance;
}

void ComputeBasicStats(MultiFab & instant, MultiFab & means,
                       const int incomp, const int outcomp, const int steps)
{
    BL_PROFILE_VAR("ComputeBasicStats()",ComputeBasicStats);

    // increment-form running mean (see RunningMean.H); (mean*(steps-1) + x)/steps swamps
    // the new sample in single precision
    const Real stepsInv = Real(1.0/static_cast<double>(steps));

    for ( MFIter mfi(instant); mfi.isValid(); ++mfi ) {

        Box tile_box  = mfi.tilebox();

        Array4<Real> means_data = means.array(mfi);
        Array4<Real> instant_data = instant.array(mfi);

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            runningMean(means_data(i,j,k,outcomp), instant_data(i,j,k,incomp), stepsInv);
        });

    }
}

void OutputVolumeMean(const MultiFab & instant, const int comp, const Real domainVol, std::string filename, const Geometry geom)
{
    BL_PROFILE_VAR("OutputVolumeMean()",OutputVolumeMean);

    Real dV = AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));
#if (AMREX_SPACEDIM == 2)
    dV *= cell_depth;
#endif

    Real result = MaskedSum(instant, comp, geom.periodicity())*dV/domainVol;

    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream ofs(filename, std::ofstream::app);
        ofs << result << "\n";
        ofs.close();
    }

}

Real MaskedSum(const MultiFab & inFab,int comp, const Periodicity& period)
{
    BL_PROFILE_VAR("MaskedSum()",MaskedSum);

    MultiFab tmpmf(inFab.boxArray(), inFab.DistributionMap(), 1, 0,
                   MFInfo(), inFab.Factory());

    MultiFab::Copy(tmpmf, inFab, comp, 0, 1, 0);

//#ifdef AMREX_USE_EB
//    if ( this -> hasEBFabFactory() && set_covered )
//        EB_set_covered( tmpmf, 0.0 );
//#endif

    auto mask = tmpmf.OverlapMask(period);
    MultiFab::Divide(tmpmf, *mask, 0, 0, 1, 0);

    return Real(SumDouble(tmpmf, 0));
}
