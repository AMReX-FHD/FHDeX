#include "compressible_functions.H"
#include "compressible_functions_F.H"

#include "common_namespace.H"


using namespace common;


void evaluateStats(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar, const MultiFab& prim, MultiFab& primMean, MultiFab& primVar, MultiFab& eta, MultiFab& etaMean, MultiFab& kappa, MultiFab& kappaMean, const int steps, const amrex::Real* dx)
{
    double totalMass;

    // Loop over boxes
    for ( MFIter mfi(prim); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        evaluate_means(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       cons[mfi].dataPtr(),  
                       consMean[mfi].dataPtr(),
                       prim[mfi].dataPtr(),
                       primMean[mfi].dataPtr(), &steps, &totalMass);

    }

    ParallelDescriptor::ReduceRealSum(totalMass);

    //Print() << "Total mass: " << totalMass << "\n";

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
                       &steps);
    }

}

void yzAverage(const MultiFab& consMean, const MultiFab& consVar, const MultiFab& primMean, const MultiFab& primVar, const MultiFab& spatialCross, const MultiFab& etaMean, const MultiFab& kappaMean, MultiFab& consMeanAv, MultiFab& consVarAv, MultiFab& primMeanAv, MultiFab& primVarAv, MultiFab& spatialCrossAv, MultiFab& etaMeanAv, MultiFab& kappaMeanAv)
{  

    // Loop over boxes

    int six = 6; //Look up how to avoid this later?
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
                       spatialCross[mfi].dataPtr(),
                       spatialCrossAv[mfi].dataPtr(), &six);

    }

}
