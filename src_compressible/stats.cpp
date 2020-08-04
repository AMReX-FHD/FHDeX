#include "compressible_functions.H"

#include "common_functions.H"

void evaluateStats(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                   const MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
                   MultiFab& spatialCross, MultiFab& miscStats, Real* miscVals,
                   const int steps, const amrex::Real* dx)
{
    BL_PROFILE_VAR("evaluateStats()",evaluateStats);
    
    double totalMass = 0.;
    double stepsminusone = steps - 1.;
    double stepsinv = 1./steps;

    Real fracvec[nspecies];
    Real massvec[nspecies];

    /*
      0  = mean xmom
      1  = instant xmom
      2  = mean xvel
      3  = mean rho
      4  = instant rho
      5  = instant xvel
      6  = instant energy
      7  = mean energy
      8  = instant ymom
      9  = mean ymom
      10 = instant zmom
      11 = mean zmom
      12 = mean cv
      13 = mean temperature
      14 = mean yvel
      15 = mean zvel
      16 = instant temperature
    */
    for (int i=0; i<20; ++i) {
        miscVals[i] = 0.;
    }

    int counter = 0;
    
    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<const Real> cu        = cons.array(mfi);
        const Array4<      Real> cumeans   = consMean.array(mfi);
        const Array4<const Real> prim      = prim_in.array(mfi);
        const Array4<      Real> primmeans = primMean.array(mfi);
        const Array4<      Real> miscstats = miscStats.array(mfi);

        // on host, not gpu
        for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {
                for (auto i = lo.x; i <= hi.x; ++i) {

                    for (int l=0; l<nvars; ++l) {
                        cumeans(i,j,k,l) = (cumeans(i,j,k,l)*stepsminusone + cu(i,j,k,l))*stepsinv;
                    }
                    
/*
            fracvec = cumeans(i,j,k,6:nvars)/cumeans(i,j,k,1)
            massvec = cumeans(i,j,k,6:nvars)

            densitymeaninv = 1.0/cumeans(i,j,k,1)

            primmeans(i,j,k,1) = cumeans(i,j,k,1)
            primmeans(i,j,k,2) = cumeans(i,j,k,2)*densitymeaninv
            primmeans(i,j,k,3) = cumeans(i,j,k,3)*densitymeaninv
            primmeans(i,j,k,4) = cumeans(i,j,k,4)*densitymeaninv

            call get_temperature(cumeans(i,j,k,5), massvec, primmeans(i,j,k,5))
            call get_pressure_gas(primmeans(i,j,k,6), fracvec, cumeans(i,j,k,1),cumeans(i,j,k,5))

            totalmass = totalmass + cu(i,j,k,1)
*/
                    
                }
            }
        }
    }

    // parallel reduce sum totalmass
    
    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);
    

        // cross_cell check FIXME
        for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {
                for (auto i = lo.x; i <= hi.x; ++i) {
                    

                    
                }
            }
        }
    }
    
    // parallel reduce sum miscvals and counter


    for (int i=0; i<17; ++i) {
        miscVals[i] /= counter;
    }

    /*
    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.validbox();
        
        evaluate_means(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       cons[mfi].dataPtr(),  
                       consMean[mfi].dataPtr(),
                       prim_in[mfi].dataPtr(),
                       primMean[mfi].dataPtr(), &steps, 
                       miscStats[mfi].dataPtr(), miscVals, &totalMass);
    }
    */


    for(int i=0;i<20;i++)
    {
        //Fix to directly address array elements

        Real temp = miscVals[i];
        ParallelDescriptor::ReduceRealSum(temp);
        miscVals[i] = temp;
    }

    ParallelDescriptor::ReduceRealSum(totalMass);

    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.validbox();

        evaluate_corrs(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),  
                       cons[mfi].dataPtr(),  
                       consMean[mfi].dataPtr(),
                       consVar[mfi].dataPtr(),
                       prim_in[mfi].dataPtr(),
                       primMean[mfi].dataPtr(),
                       primVar[mfi].dataPtr(),
                       spatialCross[mfi].dataPtr(),
                       &steps, miscStats[mfi].dataPtr(), miscVals);
    }
}

void yzAverage(const MultiFab& consMean,
               const MultiFab& consVar,
               const MultiFab& primMean,
               const MultiFab& primVar,
               const MultiFab& spatialCross,
               MultiFab& consMeanAv,
               MultiFab& consVarAv,
               MultiFab& primMeanAv,
               MultiFab& primVarAv,
               MultiFab& spatialCrossAv)
{
    BL_PROFILE_VAR("yzAverage()",yzAverage);

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
