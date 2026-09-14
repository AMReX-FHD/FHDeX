#include "common_functions.H"

Real ComputeSpatialMean(MultiFab& mf, const int& incomp)
{
    BL_PROFILE_VAR("ComputeSpatialMean()",ComputeSpatialMean);

    int npts = (AMREX_SPACEDIM == 2) ? n_cells[0]*n_cells[1] : n_cells[0]*n_cells[1]*n_cells[2];

    Real average = mf.sum(incomp) / npts;

    return average;

}

Real ComputeSpatialVariance(MultiFab& mf, const int& incomp)
{
    BL_PROFILE_VAR("ComputeSpatialVariance()",ComputeSpatialVariance);

    int npts = (AMREX_SPACEDIM == 2) ? n_cells[0]*n_cells[1] : n_cells[0]*n_cells[1]*n_cells[2];

    Real average = mf.sum(incomp) / npts;

    BoxArray ba = mf.boxArray();
    DistributionMapping dmap = mf.DistributionMap();

    // MultiFab with one component and no ghost cells
    MultiFab temp(ba, dmap, 1, 0);

    // set temp to the average
    temp.setVal(average);

    // subtract mf from temp; "temp = temp - mf"
    MultiFab::Subtract(temp,mf,incomp,0,1,0);

    // square the contents of temp
    MultiFab::Multiply(temp,temp,0,0,1,0);

    // compute the variance
    Real variance = temp.sum(0) / (npts-1);

    return variance;
}

void ComputeBasicStats(MultiFab & instant, MultiFab & means,
                       const int incomp, const int outcomp, const int steps)
{
    BL_PROFILE_VAR("ComputeBasicStats()",ComputeBasicStats);

    const Real stepsInv = 1.0/steps;
    const int stepsMinusOne = steps-1;

    for ( MFIter mfi(instant); mfi.isValid(); ++mfi ) {

        Box tile_box  = mfi.tilebox();

        Array4<Real> means_data = means.array(mfi);
        Array4<Real> instant_data = instant.array(mfi);

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            means_data(i,j,k,outcomp) = (means_data(i,j,k,outcomp)*stepsMinusOne + instant_data(i,j,k,incomp))*stepsInv;
        });

    }
}

void OutputVolumeMean(const MultiFab & instant, const int comp, const Real domainVol, std::string filename, const Geometry geom)
{
    BL_PROFILE_VAR("ComputeBasicStats()",ComputeBasicStats);

    Real result = (MaskedSum(instant, comp, geom.periodicity())*geom.CellSize()[0]*geom.CellSize()[1]*geom.CellSize()[2])/domainVol;

    if(ParallelDescriptor::MyProc() == 0) {
        std::ofstream ofs(filename, std::ofstream::app);
        ofs << result/domainVol << "\n";
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

    return tmpmf.sum(0, 0);
}
