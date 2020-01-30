#include "common_functions.H"
#include "common_functions_F.H"



void ComputeBasicStats(const MultiFab & instant, MultiFab & means, MultiFab & vars,
                       const int incomp, const int outcomp, const int steps)
{
    BL_PROFILE_VAR("ComputeBasicStats()",ComputeBasicStats);

    for ( MFIter mfi(instant); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        compute_means(BL_TO_FORTRAN_FAB(instant[mfi]),
                      BL_TO_FORTRAN_FAB(means[mfi]),
                      &incomp, &outcomp, &steps);

        compute_vars(BL_TO_FORTRAN_FAB(instant[mfi]),
                     BL_TO_FORTRAN_FAB(means[mfi]),
                     BL_TO_FORTRAN_FAB(vars[mfi]),
                     &incomp, &outcomp, &outcomp, &steps);

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

//Note comp indexes from 1.
Real SumFab(const MultiFab & in, const int ng, const int comp)
{
    BL_PROFILE_VAR("SumFab()",SumFab);

    Real total;

    for ( MFIter mfi(in); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        sum_fab(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                BL_TO_FORTRAN_FAB(in[mfi]),
                &ng, &total, &comp);
    }

    return total;    
}

//Note comp indexes from 1.
void XMeanFab(const MultiFab & in, MultiFab & out, const int ng)
{
    BL_PROFILE_VAR("SumFab()",SumFab);

    for ( MFIter mfi(in); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        x_mean_fab(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                BL_TO_FORTRAN_FAB(in[mfi]), BL_TO_FORTRAN_FAB(out[mfi]),
                &ng);
    }
 
}

void MaxSpeed(std::array< MultiFab, AMREX_SPACEDIM >& umac)
{
    BL_PROFILE_VAR("MaxSpeed()",MaxSpeed);

    Real maxspeed;

    for ( MFIter mfi(umac[1]); mfi.isValid(); ++mfi ) {
        const Box& bx = enclosedCells(mfi.validbox());
#if(AMREX_SPACEDIM == 3)
        max_speed(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                BL_TO_FORTRAN_3D(umac[1][mfi]),
                BL_TO_FORTRAN_3D(umac[2][mfi]),
                BL_TO_FORTRAN_3D(umac[3][mfi]),
                &maxspeed);
#endif
#if(AMREX_SPACEDIM == 2)
        max_speed(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                BL_TO_FORTRAN_3D(umac[1][mfi]),
                BL_TO_FORTRAN_3D(umac[2][mfi]),
                &maxspeed);
#endif
    }

    ParallelDescriptor::ReduceRealMax(maxspeed);

    Print() << "Max speed: " << maxspeed << std::endl;
 
}
