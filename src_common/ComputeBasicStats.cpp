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


