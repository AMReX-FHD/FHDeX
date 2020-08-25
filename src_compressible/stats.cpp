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

    GpuArray<Real,MAX_SPECIES> fracvec;
    GpuArray<Real,MAX_SPECIES> massvec;

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

#if 1
    
    // from namelist
    GpuArray<Real,MAX_SPECIES> molmass_gpu;
    for (int n=0; n<nspecies; ++n) {
        molmass_gpu[n] = molmass[n];
    }
    GpuArray<Real,MAX_SPECIES> hcv_gpu;
    for (int n=0; n<nspecies; ++n) {
        hcv_gpu[n] = hcv[n];
    }
    
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

            for (int l=5; l<nvars; ++l) {
                fracvec[l-5] = cumeans(i,j,k,l)/cumeans(i,j,k,0);
                massvec[l-5] = cumeans(i,j,k,l);;
            }

            Real densitymeaninv = 1.0/cumeans(i,j,k,0);

            primmeans(i,j,k,0) = cumeans(i,j,k,0);
            primmeans(i,j,k,1) = cumeans(i,j,k,1)*densitymeaninv;
            primmeans(i,j,k,2) = cumeans(i,j,k,2)*densitymeaninv;
            primmeans(i,j,k,3) = cumeans(i,j,k,3)*densitymeaninv;

            GetTemperature(cumeans(i,j,k,4), massvec, primmeans(i,j,k,4), nspecies, hcv_gpu);
            GetPressureGas(primmeans(i,j,k,5), fracvec, cumeans(i,j,k,0), cumeans(i,j,k,4),
                           nspecies, Runiv, molmass_gpu);

            totalMass = totalMass + cu(i,j,k,0);
                    
        }
        }
        }
    }

    // parallel reduce sum totalMass
    ParallelDescriptor::ReduceRealSum(totalMass);
    
    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<const Real> cu        = cons.array(mfi);
        const Array4<      Real> cumeans   = consMean.array(mfi);
        const Array4<const Real> prim      = prim_in.array(mfi);
        const Array4<      Real> primmeans = primMean.array(mfi);

        if (cross_cell >= lo.x && cross_cell <= hi.x) {
            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {

                miscVals[0] = miscVals[0] + cumeans(cross_cell,j,k,1);   //slice average of mean x momentum
                miscVals[1] = miscVals[1] + cu(cross_cell,j,k,1);        //slice average of instant x momentum
                miscVals[2] = miscVals[2] + primmeans(cross_cell,j,k,1); //slice average of mean x velocity
                miscVals[3] = miscVals[3] + cumeans(cross_cell,j,k,0);   //slice average of mean rho
                miscVals[4] = miscVals[4] + cu(cross_cell,j,k,0);        //slice average of instant rho
                miscVals[5] = miscVals[5] + prim(cross_cell,j,k,1);      //slice average of instant x velocity
                miscVals[6] = miscVals[6] + cu(cross_cell,j,k,4);        //slice average of instant energy
                miscVals[7] = miscVals[7] + cumeans(cross_cell,j,k,4);   //slice average of mean energy
                
                miscVals[8] = miscVals[8] + cu(cross_cell,j,k,2);        //slice average of instant y momentum
                miscVals[9] = miscVals[9] + cumeans(cross_cell,j,k,2);   //slice average of mean y momentum
                
                miscVals[10] = miscVals[10] + cu(cross_cell,j,k,3);      //slice average of instant z momentum
                miscVals[11] = miscVals[11] + cumeans(cross_cell,j,k,3); //slice average of mean z momentum

                Real cv = 0;
                for (int l=0; l<nspecies; ++l) {
                    cv = cv + hcv[l]*cumeans(cross_cell,j,k,5+l)/cumeans(cross_cell,j,k,0);
                }

                miscVals[12] = miscVals[12] + cv; //slice average mean cv
                
                miscVals[13] = miscVals[13] + primmeans(cross_cell,j,k,4); //slice average of mean temperature
                
                miscVals[14] = miscVals[14] + primmeans(cross_cell,j,k,2); //slice average of mean y velocity
                miscVals[15] = miscVals[15] + primmeans(cross_cell,j,k,3); //slice average of mean z velocity
                
                miscVals[16] = miscVals[16] + prim(cross_cell,j,k,4);      //slice average of instant temperature

                counter = counter + 1;
            }
            }
        }
    }

    // parallel reduce sum miscVals and counter
    ParallelDescriptor::ReduceRealSum(miscVals,20);
    ParallelDescriptor::ReduceIntSum(counter);

    for (int i=0; i<17; ++i) {
        miscVals[i] /= counter;
    }

#else
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
#endif

    // parallel reduce sum miscVals and counter
    ParallelDescriptor::ReduceRealSum(miscVals,20);
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
