#include "compressible_functions.H"

#include "common_functions.H"

// runningMean / runningMeanKahan: increment-form running averages (Kahan-compensated
// for accumulators with nonzero expectation); see RunningMean.H
#include "RunningMean.H"

// Component offsets into statsKahan (Kahan compensation terms for the running statistics).
// statsKahan must have statsKahanNComp() components (see compressible_functions.H).
//   [0, nvars)                         : consMean
//   [nvars, nvars+5)                   : consVar 0-4
//   [nvars+5, nvars+5+nprimvars+5)     : primVar
//   [nvars+nprimvars+10, ... +10)      : miscStats
int statsKahanNComp () { return nvars + 5 + (nprimvars+5) + 10; }

void evaluateStats(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                   const MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
                   MultiFab& spatialCross, MultiFab& miscStats, Real* miscVals,
                   MultiFab& statsKahan,
                   const int steps, const amrex::Real* /*dx*/)
{
    BL_PROFILE_VAR("evaluateStats()",evaluateStats);

    const Real stepsinv = Real(1.0/static_cast<double>(steps));

    const int kcu   = 0;
    const int kcuv  = nvars;
    const int kprv  = nvars + 5;
    const int kmisc = nvars + 5 + nprimvars + 5;

    GpuArray<Real,MAX_SPECIES> fracvec;

#ifdef AMREX_USE_GPU
    // The loops below run on the host over Array4 views of the MultiFab data, which is only
    // valid if that data lives in managed memory (e.g., amrex.the_arena_is_managed = 1).
    if (!cons.arena()->isManaged()) {
        amrex::Abort("evaluateStats (src_compressible/stats.cpp) runs on the host and is not ported to GPU; "
                     "set amrex.the_arena_is_managed = 1 or stats_int <= 0");
    }
#endif

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
        const Array4<      Real> sK        = statsKahan.array(mfi);

        // on host, not gpu
        for (auto k = lo.z; k <= hi.z; ++k) {
        for (auto j = lo.y; j <= hi.y; ++j) {
        for (auto i = lo.x; i <= hi.x; ++i) {

            runningMeanKahan(cumeans(i,j,k,0), sK(i,j,k,kcu+0), cu(i,j,k,0), stepsinv); // rho
            for (int l=1; l<4; ++l) {
                runningMean(cumeans(i,j,k,l), cu(i,j,k,l), stepsinv); // momentum (zero expectation)
            }
            for (int l=4; l<nvars; ++l) {
                runningMeanKahan(cumeans(i,j,k,l), sK(i,j,k,kcu+l), cu(i,j,k,l), stepsinv); // rhoE, rhoYk
            }

            Real densitymeaninv = Real(1.0)/cumeans(i,j,k,0);

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

            Real intenergy = cumeans(i,j,k,4)/cumeans(i,j,k,0) - Real(0.5)*vsqr;

            GetTemperature(intenergy, fracvec, primmeans(i,j,k,4));
            GetPressureGas(primmeans(i,j,k,5), fracvec, cumeans(i,j,k,0), primmeans(i,j,k,4));

        }
        }
        }
    }

    // Slice sums over Ny*Nz cells are accumulated in double: the slice-averaged fluctuations are
    // O(sigma/sqrt(Ny*Nz)), below single-precision summation error for O(1) values.
    // Slots 17-19 hold slice averages of per-cell fluctuations, so they need no differencing.
    double miscValsD[20];
    for (int i=0; i<20; ++i) {
        miscValsD[i] = 0.;
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

                miscValsD[0] += cumeans(cross_cell,j,k,1);   //slice average of mean x momentum
                miscValsD[1] += cu(cross_cell,j,k,1);        //slice average of instant x momentum
                miscValsD[2] += primmeans(cross_cell,j,k,1); //slice average of mean x velocity
                miscValsD[3] += cumeans(cross_cell,j,k,0);   //slice average of mean rho
                miscValsD[4] += cu(cross_cell,j,k,0);        //slice average of instant rho
                miscValsD[5] += prim(cross_cell,j,k,1);      //slice average of instant x velocity
                miscValsD[6] += cu(cross_cell,j,k,4);        //slice average of instant energy
                miscValsD[7] += cumeans(cross_cell,j,k,4);   //slice average of mean energy

                miscValsD[8] += cu(cross_cell,j,k,2);        //slice average of instant y momentum
                miscValsD[9] += cumeans(cross_cell,j,k,2);   //slice average of mean y momentum

                miscValsD[10] += cu(cross_cell,j,k,3);      //slice average of instant z momentum
                miscValsD[11] += cumeans(cross_cell,j,k,3); //slice average of mean z momentum

                Real cv = 0;
                for (int l=0; l<nspecies; ++l) {
                    cv = cv + hcv[l]*cumeans(cross_cell,j,k,5+l)/cumeans(cross_cell,j,k,0);
                }

                miscValsD[12] += cv; //slice average mean cv

                miscValsD[13] += primmeans(cross_cell,j,k,4); //slice average of mean temperature

                miscValsD[14] += primmeans(cross_cell,j,k,2); //slice average of mean y velocity
                miscValsD[15] += primmeans(cross_cell,j,k,3); //slice average of mean z velocity

                miscValsD[16] += prim(cross_cell,j,k,4);      //slice average of instant temperature

                miscValsD[17] += cu(cross_cell,j,k,0) - cumeans(cross_cell,j,k,0);     //slice average of rho fluctuation
                miscValsD[18] += cu(cross_cell,j,k,1) - cumeans(cross_cell,j,k,1);     //slice average of x momentum fluctuation
                miscValsD[19] += prim(cross_cell,j,k,4) - primmeans(cross_cell,j,k,4); //slice average of temperature fluctuation
            }
            }
        }
    }

    // parallel reduce sum miscVals
    ParallelDescriptor::ReduceRealSum(miscValsD,20);

    // compute the mean value of miscVals at the cross_cell slice
    const double n_cells_yz_inv = 1.0/static_cast<double>(n_cells_yz);
    for (int i=0; i<20; ++i) {
        miscVals[i] = Real(miscValsD[i]*n_cells_yz_inv);
    }

    //////////////////
    // evaluate_corrs
    //////////////////

    // slots 19-21 hold slice averages of per-cell fluctuations (see miscVals 17-19)
    int nstats = 22;
    Vector<double> yzAvMeansD(n_cells[0]*nstats, 0.0); // yz-sums at all x, accumulated in double

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

            yzAvMeansD[i*nstats+0] += cu(i,j,k,0); // rho instant slices
            yzAvMeansD[i*nstats+1] += cumeans(i,j,k,0); // rho mean slices
            yzAvMeansD[i*nstats+2] += cu(i,j,k,4); // energy instant slices
            yzAvMeansD[i*nstats+3] += cumeans(i,j,k,4); // energy mean slices

            yzAvMeansD[i*nstats+4] += cu(i,j,k,1); // x momentum instant slices
            yzAvMeansD[i*nstats+5] += cumeans(i,j,k,1); // x momentum mean slices

            yzAvMeansD[i*nstats+6] += cu(i,j,k,2); // y momentum instant slices
            yzAvMeansD[i*nstats+7] += cumeans(i,j,k,2); // y momentum mean slices

            yzAvMeansD[i*nstats+8] += cu(i,j,k,3); // z momentum instant slices
            yzAvMeansD[i*nstats+9] += cumeans(i,j,k,3); // z momentum mean slices

            yzAvMeansD[i*nstats+10] += prim(i,j,k,1); // x vel instant slices
            yzAvMeansD[i*nstats+11] += primmeans(i,j,k,1); // x vel mean slices

            yzAvMeansD[i*nstats+12] += prim(i,j,k,2); // y vel instant slices
            yzAvMeansD[i*nstats+13] += primmeans(i,j,k,2); // y vel mean slices

            yzAvMeansD[i*nstats+14] +=  prim(i,j,k,3); // z vel instant slices
            yzAvMeansD[i*nstats+15] += primmeans(i,j,k,3); // z vel mean slices

            Real cv = 0;
            for (int l=0; l<nspecies; ++l) {
                cv = cv + hcv[l]*cumeans(i,j,k,5+l)/cumeans(i,j,k,0);
            }

            yzAvMeansD[i*nstats+16] += cv; // cv mean slices
            yzAvMeansD[i*nstats+17] += prim(i,j,k,4); // temperature instant slices
            yzAvMeansD[i*nstats+18] += primmeans(i,j,k,4); // temperature mean slices

            yzAvMeansD[i*nstats+19] += cu(i,j,k,0) - cumeans(i,j,k,0);     // rho fluctuation slices
            yzAvMeansD[i*nstats+20] += cu(i,j,k,4) - cumeans(i,j,k,4);     // energy fluctuation slices
            yzAvMeansD[i*nstats+21] += prim(i,j,k,4) - primmeans(i,j,k,4); // temperature fluctuation slices
        }
        }
        }
    }

    // parallel reduce yzAvMeans
    ParallelDescriptor::ReduceRealSum(yzAvMeansD.dataPtr(),n_cells[0]*nstats);

    // compute mean over each slice in i for each variable
    Vector<Real> yzAvMeans(n_cells[0]*nstats);
    for (auto n = 0; n<n_cells[0]*nstats; ++n) {
        yzAvMeans[n] = Real(yzAvMeansD[n]*n_cells_yz_inv);
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
        const Array4<      Real> sK        = statsKahan.array(mfi);

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

            Real densitymeaninv = Real(1.0)/cumeans(i,j,k,0);
//            Real densitymeaninvS = 1.0/yzAvMeans[i*nstats+1];
//            Real densitymeaninvSstar = 1.0/miscVals[3];

            Real delrho = cu(i,j,k,0) - cumeans(i,j,k,0);
            Real delpx = cu(i,j,k,1) - cumeans(i,j,k,1);
            Real delpy = cu(i,j,k,2) - cumeans(i,j,k,2);
            Real delpz = cu(i,j,k,3) - cumeans(i,j,k,3);
            Real delenergy = cu(i,j,k,4) - cumeans(i,j,k,4);

            // slice-averaged fluctuations, accumulated directly (not as instant minus mean slice averages)
            Real delrhoS = yzAvMeans[i*nstats+19]; // rho(x) - <rho(x)>, sliced
            Real delES = -yzAvMeans[i*nstats+20]; // <E(x)> - E(x), sliced (sign as in the original yzAv[3]-yzAv[2])
            Real delTS = yzAvMeans[i*nstats+21];   // T(x) - <T(x)>, sliced
//            Real delpxS = yzAvMeans[i*nstats+5] - yzAvMeans[i*nstats+4];
//            Real delpyS = yzAvMeans[i*nstats+7] - yzAvMeans[i*nstats+6];
//            Real delpzS = yzAvMeans[i*nstats+9] - yzAvMeans[i*nstats+8];

            Real delrhoSstar = miscVals[17]; // rho(x*) - <rho(x*)>, sliced
            // Real delESstar = miscVals[6] - miscVals[7];
            Real delpxSstar = miscVals[18];  // jx(x*) - <jx(x*)>, sliced
            Real delTSstar = miscVals[19];   // T(x*) - <T(x*)>, sliced
//            Real delpySstar = miscVals[8] - miscVals[9];
//            Real delpzSstar = miscVals[10] - miscVals[11];

            runningMeanKahan(cuvars(i,j,k,0), sK(i,j,k,kcuv+0), delrho*delrho, stepsinv);
            runningMeanKahan(cuvars(i,j,k,1), sK(i,j,k,kcuv+1), delpx*delpx, stepsinv);
            runningMeanKahan(cuvars(i,j,k,2), sK(i,j,k,kcuv+2), delpy*delpy, stepsinv);
            runningMeanKahan(cuvars(i,j,k,3), sK(i,j,k,kcuv+3), delpz*delpz, stepsinv);
            runningMeanKahan(cuvars(i,j,k,4), sK(i,j,k,kcuv+4), delenergy*delenergy, stepsinv);

            Real delvelx = (delpx - primmeans(i,j,k,1)*delrho)*densitymeaninv;
            Real delvely = (delpy - primmeans(i,j,k,2)*delrho)*densitymeaninv;
            Real delvelz = (delpz - primmeans(i,j,k,3)*delrho)*densitymeaninv;

            primvars(i,j,k,0) = cuvars(i,j,k,0);
            runningMeanKahan(primvars(i,j,k,1), sK(i,j,k,kprv+1), delvelx*delvelx, stepsinv);
            runningMeanKahan(primvars(i,j,k,2), sK(i,j,k,kprv+2), delvely*delvely, stepsinv);
            runningMeanKahan(primvars(i,j,k,3), sK(i,j,k,kprv+3), delvelz*delvelz, stepsinv);

            Real delg = primmeans(i,j,k,1)*delpx + primmeans(i,j,k,2)*delpy + primmeans(i,j,k,3)*delpz;

            // Real delgS = yzAvMeans[i*nstats+11]*delpxS + yzAvMeans[i*nstats+13]*delpyS + yzAvMeans[i*nstats+15]*delpzS;

            // Real delgSstar = miscVals[2]*delpxSstar + miscVals[14]*delpySstar + miscVals[15]*delpzSstar;

            runningMeanKahan(primvars(i,j,k,nprimvars), sK(i,j,k,kprv+nprimvars), delg*delg, stepsinv); // gvar

            runningMean(primvars(i,j,k,nprimvars+1), delg*delenergy, stepsinv); // kgcross
            runningMeanKahan(primvars(i,j,k,nprimvars+2), sK(i,j,k,kprv+nprimvars+2), delrho*delenergy, stepsinv); // krcross
            runningMean(primvars(i,j,k,nprimvars+3), delrho*delg, stepsinv); // rgcross

            Real delT = prim(i,j,k,4) - primmeans(i,j,k,4);
            runningMeanKahan(primvars(i,j,k,4), sK(i,j,k,kprv+4), delT*delT, stepsinv);

            /*
            primvars(i,j,k,4) = (primvars(i,j,k,4)*stepsminusone + cvinv*cvinv*densitymeaninv*densitymeaninv*
                                 (cuvars(i,j,k,4) + primvars(i,j,k,nprimvars) - 2*primvars(i,j,k,nprimvars+1)
                                  + qmean*(qmean*cuvars(i,j,k,0) - 2*primvars(i,j,k,nprimvars+2) + 2*primvars(i,j,k,nprimvars+3))))*stepsinv;
            */

            // Real deltemp = (delenergy - delg - qmean*delrho)*cvinv*densitymeaninv;

            // Real deltempS = (delES - delgS - qmeanS*delrhoS)*cvinvS*densitymeaninvS;

            // Real deltempSstar = (delESstar - delgSstar - qmeanSstar*delrhoSstar)*cvinvSstar*densitymeaninvSstar;

            runningMean(miscstats(i,j,k,0), miscVals[1]*yzAvMeans[i*nstats+0], stepsinv); // <p(x*)rho(x)>, sliced

            runningMeanKahan(miscstats(i,j,k,1), sK(i,j,k,kmisc+1), delrhoS*delrhoSstar, stepsinv); // <(rho(x*)-<rho(x*)>)(rho(x)-<rho(x)>)>, sliced

            runningMeanKahan(miscstats(i,j,k,2), sK(i,j,k,kmisc+2), miscVals[16]*yzAvMeans[i*nstats+17], stepsinv); // <(T(x*)T(x))>
            runningMeanKahan(miscstats(i,j,k,3), sK(i,j,k,kmisc+3), miscVals[16]*yzAvMeans[i*nstats+0], stepsinv); // <(T(x*)rho(x))>

            runningMean(miscstats(i,j,k,4), delrhoS*delpxSstar, stepsinv); // <(jx(x*)-<jx(x*)>)(rho(x)-<rho(x)>)>, sliced

            // <p(x*)rho(x)> - <p(x*)><rho(x)>, sliced: use the centred accumulator (miscstats 4) rather than
            // differencing miscstats 0 and <p(x*)><rho(x)>
            Real delpdelrho = miscstats(i,j,k,4);
            runningMeanKahan(miscstats(i,j,k,5), sK(i,j,k,kmisc+5), delES*delrhoSstar, stepsinv); // <(rho(x*)-<rho(x*)>)(rhoE(x)-<rhoE(x)>)>, sliced

            // centred temperature correlations; <T*T> - <T><T*> is the difference of two ~T^2 numbers
            // for an O(dT^2/(Ny*Nz)) result, which is pure roundoff in single precision
            runningMeanKahan(miscstats(i,j,k,6), sK(i,j,k,kmisc+6), delTS*delTSstar, stepsinv);   // <(T(x*)-<T(x*)>)(T(x)-<T(x)>)>, sliced
            runningMeanKahan(miscstats(i,j,k,7), sK(i,j,k,kmisc+7), delrhoS*delTSstar, stepsinv); // <(T(x*)-<T(x*)>)(rho(x)-<rho(x)>)>, sliced

            spatialcross(i,j,k,0) = miscVals[13];
            spatialcross(i,j,k,1) = yzAvMeans[i*nstats+18];
            spatialcross(i,j,k,2) = miscstats(i,j,k,2);

            spatialcross(i,j,k,3) = miscstats(i,j,k,6);
            spatialcross(i,j,k,4) = miscstats(i,j,k,7);

            if (miscVals[3] == Real(0.)) {
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
