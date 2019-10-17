#include "common_functions.H"
#include "common_functions_F.H"


//Computes divergence at cell centres from velcocities at cell faces
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
