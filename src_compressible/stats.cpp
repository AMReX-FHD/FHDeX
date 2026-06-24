#include "compressible_functions.H"

#include "common_functions.H"

void evaluateStats(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                   const MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
                   MultiFab& spatialCross, MultiFab& miscStats, Real* miscVals,
                   const int steps, const amrex::Real* /*dx*/)
{
    BL_PROFILE_VAR("evaluateStats()",evaluateStats);

    double totalMass = 0.;
    double stepsminusone = steps - 1.;
    double stepsinv = 1./steps;

    GpuArray<Real,MAX_SPECIES> fracvec;

    int n_cells_yz = n_cells[1]*n_cells[2];

    /* miscVals
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

    //////////////////
    // evaluate_means
    //////////////////

    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<const Real> cu        = cons.array(mfi);
        const Array4<      Real> cumeans   = consMean.array(mfi);
        const Array4<      Real> primmeans = primMean.array(mfi);

        // on host, not gpu
        for (auto k = lo.z; k <= hi.z; ++k) {
        for (auto j = lo.y; j <= hi.y; ++j) {
        for (auto i = lo.x; i <= hi.x; ++i) {

            for (int l=0; l<nvars; ++l) {
                cumeans(i,j,k,l) = (cumeans(i,j,k,l)*stepsminusone + cu(i,j,k,l))*stepsinv;
            }

            Real densitymeaninv = 1.0/cumeans(i,j,k,0);

            for (int l=5; l<nvars; ++l) {
                fracvec[l-5] = cumeans(i,j,k,l) * densitymeaninv;
            }

            primmeans(i,j,k,0) = cumeans(i,j,k,0);
            primmeans(i,j,k,1) = cumeans(i,j,k,1)*densitymeaninv;
            primmeans(i,j,k,2) = cumeans(i,j,k,2)*densitymeaninv;
            primmeans(i,j,k,3) = cumeans(i,j,k,3)*densitymeaninv;

            Real vsqr = primmeans(i,j,k,1)*primmeans(i,j,k,1) +
                        primmeans(i,j,k,2)*primmeans(i,j,k,2) +
                        primmeans(i,j,k,3)*primmeans(i,j,k,3);

            Real intenergy = cumeans(i,j,k,4)/cumeans(i,j,k,0) - 0.5*vsqr;

            GetTemperature(intenergy, fracvec, primmeans(i,j,k,4));
            GetPressureGas(primmeans(i,j,k,5), fracvec, cumeans(i,j,k,0), primmeans(i,j,k,4));

            totalMass = totalMass + cu(i,j,k,0);

        }
        }
        }
    }

    // parallel reduce sum totalMass
    ParallelDescriptor::ReduceRealSum(totalMass);

    for (int i=0; i<20; ++i) {
        miscVals[i] = 0.;
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
            }
            }
        }
    }

    // parallel reduce sum miscVals
    ParallelDescriptor::ReduceRealSum(miscVals,20);

    // compute the mean value of miscVals at the cross_cell slice
    for (int i=0; i<20; ++i) {
        miscVals[i] /= n_cells_yz;
    }

    //////////////////
    // evaluate_corrs
    //////////////////

    int nstats = 19;
    Vector<Real>  yzAvMeans(n_cells[0]*nstats, 0.0); // yz-average at all x

    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<const Real> cu        = cons.array(mfi);
        const Array4<const Real> cumeans   = consMean.array(mfi);
        const Array4<const Real> prim      = prim_in.array(mfi);
        const Array4<const Real> primmeans = primMean.array(mfi);

        for (auto k = lo.z; k <= hi.z; ++k) {
        for (auto j = lo.y; j <= hi.y; ++j) {
        for (auto i = lo.x; i <= hi.x; ++i) {

            yzAvMeans[i*nstats+0] += cu(i,j,k,0); // rho instant slices
            yzAvMeans[i*nstats+1] += cumeans(i,j,k,0); // rho mean slices
            yzAvMeans[i*nstats+2] += cu(i,j,k,4); // energy instant slices
            yzAvMeans[i*nstats+3] += cumeans(i,j,k,4); // energy mean slices

            yzAvMeans[i*nstats+4] += cu(i,j,k,1); // x momentum instant slices
            yzAvMeans[i*nstats+5] += cumeans(i,j,k,1); // x momentum mean slices

            yzAvMeans[i*nstats+6] += cu(i,j,k,2); // y momentum instant slices
            yzAvMeans[i*nstats+7] += cumeans(i,j,k,2); // y momentum mean slices

            yzAvMeans[i*nstats+8] += cu(i,j,k,3); // z momentum instant slices
            yzAvMeans[i*nstats+9] += cumeans(i,j,k,3); // z momentum mean slices

            yzAvMeans[i*nstats+10] += prim(i,j,k,1); // x vel instant slices
            yzAvMeans[i*nstats+11] += primmeans(i,j,k,1); // x vel mean slices

            yzAvMeans[i*nstats+12] += prim(i,j,k,2); // y vel instant slices
            yzAvMeans[i*nstats+13] += primmeans(i,j,k,2); // y vel mean slices

            yzAvMeans[i*nstats+14] +=  prim(i,j,k,3); // z vel instant slices
            yzAvMeans[i*nstats+15] += primmeans(i,j,k,3); // z vel mean slices

            Real cv = 0;
            for (int l=0; l<nspecies; ++l) {
                cv = cv + hcv[l]*cumeans(i,j,k,5+l)/cumeans(i,j,k,0);
            }

            yzAvMeans[i*nstats+16] += cv; // cv mean slices
            yzAvMeans[i*nstats+17] += prim(i,j,k,4); // temperature instant slices
            yzAvMeans[i*nstats+18] += primmeans(i,j,k,4); // temperature mean slices
        }
        }
        }
    }

    // parallel reduce yzAvMeans
    ParallelDescriptor::ReduceRealSum(yzAvMeans.dataPtr(),n_cells[0]*nstats);

    // compute mean over each slice in i for each variable
    for (auto n = 0; n<n_cells[0]*nstats; ++n) {
        yzAvMeans[n] /= n_cells_yz;
    }

    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<const Real> cu        = cons.array(mfi);
        const Array4<const Real> cumeans   = consMean.array(mfi);
        const Array4<      Real> cuvars    = consVar.array(mfi);
        const Array4<const Real> prim      = prim_in.array(mfi);
        const Array4<const Real> primmeans = primMean.array(mfi);
        const Array4<      Real> primvars  = primVar.array(mfi);
        const Array4<      Real> spatialcross = spatialCross.array(mfi);
        const Array4<      Real> miscstats = miscStats.array(mfi);

        for (auto k = lo.z; k <= hi.z; ++k) {
        for (auto j = lo.y; j <= hi.y; ++j) {
        for (auto i = lo.x; i <= hi.x; ++i) {

            Real cv = 0.;
            for (int l=0; l<nspecies; ++l) {
                cv = cv + hcv[l]*cumeans(i,j,k,5+l)/cumeans(i,j,k,0);
            }

            // Real cvinv = 1.0/cv;
            // Real cvinvS = 1.0/yzAvMeans[i*nstats+16];
            // Real cvinvSstar = 1.0/miscVals[12];

            // Vars
#if 0
            Real qmean = cv*primmeans(i,j,k,4)-0.5*(  primmeans(i,j,k,1)*primmeans(i,j,k,1)
                                                    + primmeans(i,j,k,2)*primmeans(i,j,k,2)
                                                    + primmeans(i,j,k,3)*primmeans(i,j,k,3));

            Real qmeanS = yzAvMeans[i*nstats+16]*yzAvMeans[i*nstats+18]-0.5*(  yzAvMeans[i*nstats+11]*yzAvMeans[i*nstats+11]
                                                                 + yzAvMeans[i*nstats+13]*yzAvMeans[i*nstats+13]
                                                                 + yzAvMeans[i*nstats+15]*yzAvMeans[i*nstats+15]);

            Real qmeanSstar =  miscVals[12]*yzAvMeans[i*nstats+18]-0.5*(  miscVals[2]*miscVals[2]
                                                                          + miscVals[14]*miscVals[14]
                                                                          + miscVals[15]*miscVals[15]);
#endif

            Real densitymeaninv = 1.0/cumeans(i,j,k,0);
//            Real densitymeaninvS = 1.0/yzAvMeans[i*nstats+1];
//            Real densitymeaninvSstar = 1.0/miscVals[3];

            Real delrho = cu(i,j,k,0) - cumeans(i,j,k,0);
            Real delpx = cu(i,j,k,1) - cumeans(i,j,k,1);
            Real delpy = cu(i,j,k,2) - cumeans(i,j,k,2);
            Real delpz = cu(i,j,k,3) - cumeans(i,j,k,3);
            Real delenergy = cu(i,j,k,4) - cumeans(i,j,k,4);

            Real delrhoS = yzAvMeans[i*nstats+0] - yzAvMeans[i*nstats+1]; // rho(x) - <rho(x)>, sliced
            Real delES = yzAvMeans[i*nstats+3] - yzAvMeans[i*nstats+2]; // E(x) - <E(x)>, sliced
//            Real delpxS = yzAvMeans[i*nstats+5] - yzAvMeans[i*nstats+4];
//            Real delpyS = yzAvMeans[i*nstats+7] - yzAvMeans[i*nstats+6];
//            Real delpzS = yzAvMeans[i*nstats+9] - yzAvMeans[i*nstats+8];

            Real delrhoSstar = miscVals[4] - miscVals[3];
            // Real delESstar = miscVals[6] - miscVals[7];
            Real delpxSstar = miscVals[1] - miscVals[0];
//            Real delpySstar = miscVals[8] - miscVals[9];
//            Real delpzSstar = miscVals[10] - miscVals[11];

            cuvars(i,j,k,0) = (cuvars(i,j,k,0)*stepsminusone + delrho*delrho)*stepsinv;
            cuvars(i,j,k,1) = (cuvars(i,j,k,1)*stepsminusone + delpx*delpx)*stepsinv;
            cuvars(i,j,k,2) = (cuvars(i,j,k,2)*stepsminusone + delpy*delpy)*stepsinv;
            cuvars(i,j,k,3) = (cuvars(i,j,k,3)*stepsminusone + delpz*delpz)*stepsinv;
            cuvars(i,j,k,4) = (cuvars(i,j,k,4)*stepsminusone + delenergy*delenergy)*stepsinv;

            Real delvelx = (delpx - primmeans(i,j,k,1)*delrho)*densitymeaninv;
            Real delvely = (delpy - primmeans(i,j,k,2)*delrho)*densitymeaninv;
            Real delvelz = (delpz - primmeans(i,j,k,3)*delrho)*densitymeaninv;

            primvars(i,j,k,0) = cuvars(i,j,k,0);
            primvars(i,j,k,1) = (primvars(i,j,k,1)*stepsminusone + delvelx*delvelx)*stepsinv;
            primvars(i,j,k,2) = (primvars(i,j,k,2)*stepsminusone + delvely*delvely)*stepsinv;
            primvars(i,j,k,3) = (primvars(i,j,k,3)*stepsminusone + delvelz*delvelz)*stepsinv;

            Real delg = primmeans(i,j,k,1)*delpx + primmeans(i,j,k,2)*delpy + primmeans(i,j,k,3)*delpz;

            // Real delgS = yzAvMeans[i*nstats+11]*delpxS + yzAvMeans[i*nstats+13]*delpyS + yzAvMeans[i*nstats+15]*delpzS;

            // Real delgSstar = miscVals[2]*delpxSstar + miscVals[14]*delpySstar + miscVals[15]*delpzSstar;

            primvars(i,j,k,nprimvars) = (primvars(i,j,k,nprimvars)*stepsminusone + delg*delg)*stepsinv; // gvar

            primvars(i,j,k,nprimvars+1) = (primvars(i,j,k,nprimvars+1)*stepsminusone + delg*delenergy)*stepsinv; // kgcross
            primvars(i,j,k,nprimvars+2) = (primvars(i,j,k,nprimvars+2)*stepsminusone + delrho*delenergy)*stepsinv; // krcross
            primvars(i,j,k,nprimvars+3) = (primvars(i,j,k,nprimvars+3)*stepsminusone + delrho*delg)*stepsinv; // rgcross

            Real delT = prim(i,j,k,4) - primmeans(i,j,k,4);
            primvars(i,j,k,4) =  (primvars(i,j,k,4)*stepsminusone + delT*delT)*stepsinv;

            /*
            primvars(i,j,k,4) = (primvars(i,j,k,4)*stepsminusone + cvinv*cvinv*densitymeaninv*densitymeaninv*
                                 (cuvars(i,j,k,4) + primvars(i,j,k,nprimvars) - 2*primvars(i,j,k,nprimvars+1)
                                  + qmean*(qmean*cuvars(i,j,k,0) - 2*primvars(i,j,k,nprimvars+2) + 2*primvars(i,j,k,nprimvars+3))))*stepsinv;
            */

            // Real deltemp = (delenergy - delg - qmean*delrho)*cvinv*densitymeaninv;

            // Real deltempS = (delES - delgS - qmeanS*delrhoS)*cvinvS*densitymeaninvS;

            // Real deltempSstar = (delESstar - delgSstar - qmeanSstar*delrhoSstar)*cvinvSstar*densitymeaninvSstar;

            miscstats(i,j,k,0) = (miscstats(i,j,k,0)*stepsminusone + miscVals[1]*yzAvMeans[i*nstats+0])*stepsinv; // <p(x*)rho(x)>, sliced

            Real delpdelrho = miscstats(i,j,k,0) - miscVals[0]*cumeans(i,j,k,0); // <p(x*)rho(x)> - <p(x*)><rho(x)>, sliced

            miscstats(i,j,k,1) = (miscstats(i,j,k,1)*stepsminusone + delrhoS*delrhoSstar)*stepsinv; // <(rho(x*)-<rho(x*)>)(rho(x)-<rho(x)>)>, sliced

            miscstats(i,j,k,2) = (miscstats(i,j,k,2)*stepsminusone + miscVals[16]*yzAvMeans[i*nstats+17])*stepsinv; // <(T(x*)T(x))>
            miscstats(i,j,k,3) = (miscstats(i,j,k,3)*stepsminusone + miscVals[16]*yzAvMeans[i*nstats+0])*stepsinv; // <(T(x*)rho(x))>

            miscstats(i,j,k,4) = (miscstats(i,j,k,4)*stepsminusone + delrhoS*delpxSstar)*stepsinv; // <(jx(x*)-<jx(x*)>)(rho(x)-<rho(x)>)>, sliced
            miscstats(i,j,k,5) = (miscstats(i,j,k,5)*stepsminusone + delES*delrhoSstar)*stepsinv; // <(rho(x*)-<rho(x*)>)(rhoE(x)-<rhoE(x)>)>, sliced

            spatialcross(i,j,k,0) = miscVals[13];
            spatialcross(i,j,k,1) = yzAvMeans[i*nstats+18];
            spatialcross(i,j,k,2) = miscstats(i,j,k,2);

            spatialcross(i,j,k,3) = miscstats(i,j,k,2) - yzAvMeans[i*nstats+18]*miscVals[13];
            spatialcross(i,j,k,4) = miscstats(i,j,k,3) - yzAvMeans[i*nstats+1]*miscVals[13];

            if (miscVals[3] == 0.) {
                spatialcross(i,j,k,5) = 0.;
            } else {
                spatialcross(i,j,k,5) = (delpdelrho - miscVals[2]*miscstats(i,j,k,1))/miscVals[3];
            }
            spatialcross(i,j,k,6) = miscstats(i,j,k,4);
            spatialcross(i,j,k,7) = miscstats(i,j,k,5);
        }
        }
        }

    } // end MFIter
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

    WriteHorizontalAverageToMF(consMean, consMeanAv,
                               0, 0, consMean.nComp());
    WriteHorizontalAverageToMF(consVar, consVarAv,
                               0, 0, consVar.nComp());
    WriteHorizontalAverageToMF(primMean, primMeanAv,
                               0, 0, primMean.nComp());
    WriteHorizontalAverageToMF(primVar, primVarAv,
                               0, 0, primVar.nComp());
    WriteHorizontalAverageToMF(spatialCross, spatialCrossAv,
                               0, 0, spatialCrossAv.nComp());

}
