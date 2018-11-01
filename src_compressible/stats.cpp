#include "compressible_functions.H"
#include "compressible_functions_F.H"


void evaluateStats(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar, const MultiFab& prim, MultiFab& primMean, MultiFab& primVar, MultiFab& spatialCross, const int steps, const amrex::Real* dx)
{
    double del1;
    double del2;
    double del3;
    // Loop over boxes
    for ( MFIter mfi(prim); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        evaluate_means(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       cons[mfi].dataPtr(),  
                       consMean[mfi].dataPtr(),
                       prim[mfi].dataPtr(),
                       primMean[mfi].dataPtr(), &steps, &del1, &del2, &del3);

    }

    ParallelDescriptor::ReduceRealSum(del1);
    ParallelDescriptor::ReduceRealSum(del2);
    ParallelDescriptor::ReduceRealSum(del3);

    for ( MFIter mfi(prim); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        evaluate_corrs(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       cons[mfi].dataPtr(),  
                       consMean[mfi].dataPtr(),
                       consVar[mfi].dataPtr(),
                       prim[mfi].dataPtr(),
                       primMean[mfi].dataPtr(),
                       primVar[mfi].dataPtr(),
                       spatialCross[mfi].dataPtr(),
                       &steps, &del1, &del2, &del3);

    }
}
