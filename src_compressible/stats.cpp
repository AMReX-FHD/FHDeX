#include "compressible_functions.H"
#include "compressible_functions_F.H"

#include "common_namespace.H"


using namespace common;


void evaluateStats(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar, const MultiFab& prim, MultiFab& primMean, MultiFab& primVar, MultiFab& spatialCross, Real* delHolder1, Real* delHolder2, Real* delHolder3, Real* delHolder4, Real* delHolder5, Real* delHolder6, const int steps, const amrex::Real* dx)
{
    double del1;
    double del2;
    double del3;
    double del4;
    double del5;
    double del6;

    // Loop over boxes
    for ( MFIter mfi(prim); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        evaluate_means(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       cons[mfi].dataPtr(),  
                       consMean[mfi].dataPtr(),
                       prim[mfi].dataPtr(),
                       primMean[mfi].dataPtr(), &steps, delHolder1, delHolder2, delHolder3, delHolder4, delHolder5, delHolder6);

    }

    //this is a bit of a hack? Reduce real sum should work with vectors
    for(int i=0;i<(n_cells[1]*n_cells[2]);i++)
    {
        //Fix to directly address array elements

        del1 = delHolder1[i];
        del2 = delHolder2[i];
        del3 = delHolder3[i];
        del4 = delHolder4[i];
        del5 = delHolder5[i];
        del6 = delHolder6[i];

        ParallelDescriptor::ReduceRealSum(del1);
        ParallelDescriptor::ReduceRealSum(del2);
        ParallelDescriptor::ReduceRealSum(del3);
        ParallelDescriptor::ReduceRealSum(del4);
        ParallelDescriptor::ReduceRealSum(del5);
        ParallelDescriptor::ReduceRealSum(del6);

        delHolder1[i] = del1;
        delHolder2[i] = del2;
        delHolder3[i] = del3;
        delHolder4[i] = del4;
        delHolder5[i] = del5;
        delHolder6[i] = del6;
    }

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
                       &steps, delHolder1, delHolder2, delHolder3, delHolder4, delHolder5, delHolder6);
    }

}

void yzAverage(const MultiFab& consMean, const MultiFab& consVar, const MultiFab& primMean, const MultiFab& primVar, const MultiFab& spatialCross, MultiFab& consMeanAv, MultiFab& consVarAv, MultiFab& primMeanAv, MultiFab& primVarAv, MultiFab& spatialCrossAv)
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
                       primVar[mfi].dataPtr(),
                       primVarAv[mfi].dataPtr(), &primVarsPlusFive);

        multifab_yzav(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       spatialCross[mfi].dataPtr(),
                       spatialCrossAv[mfi].dataPtr(), &six);

    }

}
