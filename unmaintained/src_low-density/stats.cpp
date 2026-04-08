#include "compressible_functions.H"





void evaluateStats(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar, const int steps, const amrex::Real* dx)
{


    double totalMass;

    // Loop over boxes
    for ( MFIter mfi(cons); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        evaluate_means(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                       cons[mfi].dataPtr(),
                       consMean[mfi].dataPtr(),
                       &steps);

    }



    //Print() << "Total mass: " << totalMass << "\n";

    for ( MFIter mfi(cons); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        evaluate_corrs(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                       cons[mfi].dataPtr(),
                       consMean[mfi].dataPtr(),
                       consVar[mfi].dataPtr(),
                       &steps);
    }

}

