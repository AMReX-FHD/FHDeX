#include "compressible_functions.H"
#include "compressible_functions_F.H"

#include "common_namespace.H"


using namespace common;


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

void yzAverage(const MultiFab& consMean, const MultiFab& consVar, const MultiFab& primMean, const MultiFab& primVar, const MultiFab& spatialCross, MultiFab& consMeanAv, MultiFab& consVarAv, MultiFab& primMeanAv, MultiFab& primVarAv, MultiFab& spatialCrossAv)
{

    // Loop over boxes

    int three = 3; //Look up how to avoid this later?
    int primVarsPlusFive = nprimvars + 5;

    for ( MFIter mfi(consMean); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        multifab_yzav(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       consMean[mfi].dataPtr(),
                       consMeanAv[mfi].dataPtr(), &nvars);

        multifab_yzav(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       consVar[mfi].dataPtr(),
                       consVarAv[mfi].dataPtr(), &nvars);

        multifab_yzav(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       primMean[mfi].dataPtr(),
                       primMeanAv[mfi].dataPtr(), &nprimvars);

        multifab_yzav(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       primVar[mfi].dataPtr(),
                       primVarAv[mfi].dataPtr(), &nprimvars);

        multifab_yzav(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       primVar[mfi].dataPtr(),
                       primVarAv[mfi].dataPtr(), &primVarsPlusFive);

        multifab_yzav(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       spatialCross[mfi].dataPtr(),
                       spatialCrossAv[mfi].dataPtr(), &three);

    }

}
