#include "compressible_functions.H"
#include "compressible_functions_stag.H"

#include "common_functions.H"

// runningMean / runningMeanKahan: increment-form running averages (Kahan-compensated
// for accumulators with nonzero expectation); see RunningMean.H
#include "RunningMean.H"

namespace {
// Sum a device vector across MPI ranks through a host copy (MPI cannot be assumed to accept
// device pointers).
void ReduceDeviceVectorSum (amrex::Gpu::DeviceVector<Real>& v, int n)
{
    amrex::Gpu::HostVector<Real> h(n);
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, v.begin(), v.begin()+n, h.begin());
    amrex::Gpu::streamSynchronize();
    ParallelDescriptor::ReduceRealSum(h.data(), n);
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h.begin(), h.end(), v.begin());
    amrex::Gpu::streamSynchronize();
}
}

///////////////////////////////////////////
// Evaluate Stats for the 3D case /////////
/// ///////////////////////////////////////
void evaluateStatsStag3D(MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                         MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
                         const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                         std::array<MultiFab, AMREX_SPACEDIM>& velMean,
                         std::array<MultiFab, AMREX_SPACEDIM>& velVar,
                         const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                         std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                         std::array<MultiFab, AMREX_SPACEDIM>& cumomVar,
                         MultiFab& coVar,
                         MultiFab& mom3,
                         MultiFab& mom4,
                         MultiFab& theta, MultiFab& thetaMean, MultiFab& thetaVar, MultiFab& thetacoVar,
                         Vector<Real>& dataSliceMeans_xcross,
                         Vector<Real>& spatialCross3D, const int ncross,
                         const amrex::Box& domain,
                         StatsKahan& kahan,
                         const int steps,
                         const Geometry& geom)
{
    BL_PROFILE_VAR("evaluateStatsStag3D()",evaluateStatsStag3D);

    //// Evaluate Means
    if (plot_means) {
        EvaluateStatsMeans(cons,consMean,prim_in,primMean,vel,velMean,cumom,cumomMean,theta,thetaMean,kahan,steps);
        consMean.FillBoundary(geom.periodicity());
        primMean.FillBoundary(geom.periodicity());
    }

    //// Evaluate Variances and Covariances
    if ((plot_vars) or (plot_covars)) {
        if (!plot_means) Abort("plot_vars and plot_covars require plot_means");
        EvaluateVarsCoVarsMom3(cons,consMean,consVar,prim_in,primMean,primVar,velMean,velVar,
                           cumom,cumomMean,cumomVar,coVar,mom3,mom4,theta,thetaMean,thetaVar,thetacoVar,kahan,steps);
        consVar.FillBoundary(geom.periodicity());
        primVar.FillBoundary(geom.periodicity());
    }

    //// Evaluate Spatial Correlations

    // contains yz-averaged running & instantaneous averages of conserved variables (2*nvars) + primitive variables [vx, vy, vz, T, Yk]: 2*4 + 2*nspecies
    int nstats = 2*nvars+8+2*nspecies;

    amrex::Gpu::HostVector<Real> cu_avg;
    amrex::Gpu::HostVector<Real> cumeans_avg;
    amrex::Gpu::HostVector<Real> prim_avg;
    amrex::Gpu::HostVector<Real> primmeans_avg;

    cu_avg = sumToLine(cons,0,nvars,domain,0,false);
    cumeans_avg = sumToLine(consMean,0,nvars,domain,0,false);
    prim_avg = sumToLine(prim_in,1,nspecies+4,domain,0,false);
    primmeans_avg = sumToLine(primMean,1,nspecies+4,domain,0,false);

    const Real nyzinv = Real(1.0/(static_cast<double>(n_cells[1])*static_cast<double>(n_cells[2])));
    for (int i=0; i<nvars*domain.length(0); ++i) {
        cu_avg[i] *= nyzinv;
        cumeans_avg[i] *= nyzinv;
    }
    for (int i=0; i<(nspecies+4)*domain.length(0); ++i) {
        prim_avg[i] *= nyzinv;
        primmeans_avg[i] *= nyzinv;
    }

    // Update Spatial Correlations
    if (plot_cross) {
        // Plane-averaged fluctuations are O(sigma/sqrt(Ny*Nz)) while the plane averages themselves
        // are O(1); differencing two separately rounded plane sums (cu_avg - cumeans_avg) leaves
        // mostly roundoff in single precision. Average the per-cell fluctuations instead.
        // Same component layout as cu_avg and prim_avg.
        MultiFab delcons(cons.boxArray(), cons.DistributionMap(), nvars, 0);
        MultiFab::Copy    (delcons, cons,     0, 0, nvars, 0);
        MultiFab::Subtract(delcons, consMean, 0, 0, nvars, 0);
        MultiFab delprim(prim_in.boxArray(), prim_in.DistributionMap(), nspecies+4, 0);
        MultiFab::Copy    (delprim, prim_in,  1, 0, nspecies+4, 0);
        MultiFab::Subtract(delprim, primMean, 1, 0, nspecies+4, 0);

        amrex::Gpu::HostVector<Real> delcu_avg   = sumToLine(delcons,0,nvars,domain,0,false);
        amrex::Gpu::HostVector<Real> delprim_avg = sumToLine(delprim,0,nspecies+4,domain,0,false);
        for (int i=0; i<nvars*domain.length(0); ++i) {
            delcu_avg[i] *= nyzinv;
        }
        for (int i=0; i<(nspecies+4)*domain.length(0); ++i) {
            delprim_avg[i] *= nyzinv;
        }

        EvaluateSpatialCorrelations3D(spatialCross3D,kahan.spatialCrossVec,dataSliceMeans_xcross,cu_avg,cumeans_avg,prim_avg,primmeans_avg,
                                      delcu_avg,delprim_avg,steps,nstats,ncross);
    }
}


//////////////////////////////////////////////
//// Evaluate Stats for the 2D case (do later)
///// ////////////////////////////////////////
void evaluateStatsStag2D(MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                         MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
                         const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                         std::array<MultiFab, AMREX_SPACEDIM>& velMean,
                         std::array<MultiFab, AMREX_SPACEDIM>& velVar,
                         const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                         std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                         std::array<MultiFab, AMREX_SPACEDIM>& cumomVar,
                         MultiFab& coVar,
                         MultiFab& mom3,
                         MultiFab& mom4,
                         MultiFab& theta, MultiFab& thetaMean, MultiFab& thetaVar, MultiFab& thetacoVar,
                         MultiFab& /*spatialCross2D*/, const int /*ncross*/,
                         StatsKahan& kahan,
                         const int steps,
                         const Geometry& geom)
{
    BL_PROFILE_VAR("evaluateStatsStag2D()",evaluateStatsStag2D);

    //// Evaluate Means
    if (plot_means) {
        EvaluateStatsMeans(cons,consMean,prim_in,primMean,vel,velMean,cumom,cumomMean,theta,thetaMean,kahan,steps);
        consMean.FillBoundary(geom.periodicity());
        primMean.FillBoundary(geom.periodicity());
    }

    //// Evaluate Variances and Covariances
    if ((plot_vars) or (plot_covars)) {
        if (!plot_means) Abort("plot_vars and plot_covars require plot_means");
        EvaluateVarsCoVarsMom3(cons,consMean,consVar,prim_in,primMean,primVar,velMean,velVar,
                           cumom,cumomMean,cumomVar,coVar,mom3,mom4,theta,thetaMean,thetaVar,thetacoVar,kahan,steps);
        consVar.FillBoundary(geom.periodicity());
        primVar.FillBoundary(geom.periodicity());
    }

    //// Evaluate Spatial Correlations (do later)

    //// contains yz-averaged running & instantaneous averages of conserved variables (2*nvars) + primitive variables [vx, vy, vz, T, Yk]: 2*4 + 2*nspecies
    //int nstats = 2*nvars+8+2*nspecies;

    //// Get all nstats at xcross for all j and k, and store in data_xcross
    //amrex::Gpu::DeviceVector<Real> data_xcross(nstats*n_cells[1]*n_cells[2], 0.0); // values at x* for a given y and z
    //GetPencilCross(data_xcross,consMean,primMean,prim_in,cons,nstats);
    //ParallelDescriptor::ReduceRealSum(data_xcross.data(),nstats*n_cells[1]*n_cells[2]);

    //// Update Spatial Correlations
    //EvaluateSpatialCorrelations2D(spatialCross2D,data_xcross,consMean,primMean,prim_in,cons,vel,velMean,cumom,cumomMean,steps,nstats);
}


///////////////////////////////////////////
// Evaluate Stats for the 1D case /////////
/// ///////////////////////////////////////
void evaluateStatsStag1D(MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                         MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
                         const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                         std::array<MultiFab, AMREX_SPACEDIM>& velMean,
                         std::array<MultiFab, AMREX_SPACEDIM>& velVar,
                         const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                         std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                         std::array<MultiFab, AMREX_SPACEDIM>& cumomVar,
                         MultiFab& coVar,
                         MultiFab& mom3,
                         MultiFab& mom4,
                         MultiFab& theta, MultiFab& thetaMean, MultiFab& thetaVar, MultiFab& thetacoVar,
                         MultiFab& spatialCross1D, const int ncross,
                         StatsKahan& kahan,
                         const int steps,
                         const Geometry& geom)
{
    BL_PROFILE_VAR("evaluateStatsStag1D()",evaluateStatsStag1D);

    //// Evaluate Means
    if (plot_means) {
        EvaluateStatsMeans(cons,consMean,prim_in,primMean,vel,velMean,cumom,cumomMean,theta,thetaMean,kahan,steps);
        consMean.FillBoundary(geom.periodicity());
        primMean.FillBoundary(geom.periodicity());
    }

    //// Evaluate Variances and Covariances
    if ((plot_vars) or (plot_covars)) {
        if (!plot_means) Abort("plot_vars and plot_covars require plot_means");
        EvaluateVarsCoVarsMom3(cons,consMean,consVar,prim_in,primMean,primVar,velMean,velVar,
                           cumom,cumomMean,cumomVar,coVar,mom3,mom4,theta,thetaMean,thetaVar,thetacoVar,kahan,steps);
        consVar.FillBoundary(geom.periodicity());
        primVar.FillBoundary(geom.periodicity());
    }

    //// Evaluate Spatial Correlations

    // contains yz-averaged running & instantaneous averages of conserved variables (2*nvars) + primitive variables [vx, vy, vz, T, Yk]: 2*4 + 2*nspecies
    int nstats = 2*nvars+8+2*nspecies;

    // Get all nstats at xcross for all j and k, and store in data_xcross
    if (plot_cross) {
        if (all_correl == 0) { // for a specified cross_cell
            amrex::Gpu::DeviceVector<Real> data_xcross(nstats*n_cells[1]*n_cells[2], 0.0); // values at x* for a given y and z
            GetPencilCross(data_xcross,consMean,primMean,prim_in,cons,nstats,cross_cell);
            ReduceDeviceVectorSum(data_xcross,nstats*n_cells[1]*n_cells[2]);

            // Update Spatial Correlations
            EvaluateSpatialCorrelations1D(spatialCross1D,kahan.spatialCrossMF,data_xcross,consMean,primMean,prim_in,cons,vel,velMean,cumom,cumomMean,steps,nstats,ncross,0);
        }
        else { // at five equi-distant x*
            GpuArray<int, 5 > x_star;
            x_star[0] = 0;
            x_star[1] = (int)amrex::Math::floor(1.0*n_cells[0]/4.0);
            x_star[2] = (int)amrex::Math::floor(2.0*n_cells[0]/4.0);
            x_star[3] = (int)amrex::Math::floor(3.0*n_cells[0]/4.0);
            x_star[4] = n_cells[0] - 1;

            for (int i=0; i<5; ++i) {
                amrex::Gpu::DeviceVector<Real> data_xcross(nstats*n_cells[1]*n_cells[2], 0.0); // values at x* for a given y and z
                GetPencilCross(data_xcross,consMean,primMean,prim_in,cons,nstats,x_star[i]);
                ReduceDeviceVectorSum(data_xcross,nstats*n_cells[1]*n_cells[2]);
                EvaluateSpatialCorrelations1D(spatialCross1D,kahan.spatialCrossMF,data_xcross,consMean,primMean,prim_in,cons,vel,velMean,cumom,cumomMean,steps,nstats,ncross,i);
            }
        }
    }
}


///////////////////////
// evaluate_means
///////////////////////
void EvaluateStatsMeans(MultiFab& cons, MultiFab& consMean,
                        MultiFab& prim_in, MultiFab& primMean,
                        const std::array<MultiFab, AMREX_SPACEDIM>& /*vel*/,
                        std::array<MultiFab, AMREX_SPACEDIM>& velMean,
                        const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                        std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                        MultiFab& theta, MultiFab& thetaMean,
                        StatsKahan& kahan,
                        const int steps)
{
    BL_PROFILE_VAR("EvaluateStatsMeans()",EvaluateStatsMeans);

    Real stepsinv = Real(1.0/static_cast<double>(steps));

    // Loop over boxes
    for ( MFIter mfi(prim_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<const Real> momx      = cumom[0].array(mfi);
        const Array4<const Real> momy      = cumom[1].array(mfi);
        const Array4<const Real> momz      = cumom[2].array(mfi);
        const Array4<      Real> momxmeans = cumomMean[0].array(mfi);
        const Array4<      Real> momymeans = cumomMean[1].array(mfi);
        const Array4<      Real> momzmeans = cumomMean[2].array(mfi);

        // update mean momentum
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
             runningMean(momxmeans(i,j,k), momx(i,j,k), stepsinv);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            runningMean(momymeans(i,j,k), momy(i,j,k), stepsinv);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            runningMean(momzmeans(i,j,k), momz(i,j,k), stepsinv);
        });

    }

    for ( MFIter mfi(prim_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bxg = mfi.growntilebox(ngc[0]);
        const Array4<      Real> cu        = cons.array(mfi);
        const Array4<      Real> cumeans   = consMean.array(mfi);
        const Array4<      Real> cumeansK  = kahan.cuMeans.array(mfi);

        // update mean density (required in the ghost region as well)
        amrex::ParallelFor(bxg, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            runningMeanKahan(cumeans(i,j,k,0), cumeansK(i,j,k,0), cu(i,j,k,0), stepsinv);
        });
    }

    for ( MFIter mfi(prim_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx  = mfi.tilebox();

        const Array4<      Real> cu        = cons.array(mfi);
        const Array4<      Real> cumeans   = consMean.array(mfi);
        const Array4<      Real> prim      = prim_in.array(mfi);
        const Array4<      Real> primmeans = primMean.array(mfi);

        const Array4<const Real> momxmeans = cumomMean[0].array(mfi);
        const Array4<const Real> momymeans = cumomMean[1].array(mfi);
        const Array4<const Real> momzmeans = cumomMean[2].array(mfi);

        const Array4<      Real> surfcov      = (nspec_surfcov>0) ? theta.array(mfi) : cons.array(mfi);
        const Array4<      Real> surfcovmeans = (nspec_surfcov>0) ? thetaMean.array(mfi) : consMean.array(mfi);

        const Array4<      Real> cumeansK      = kahan.cuMeans.array(mfi);
        const Array4<      Real> primmeansK    = kahan.primMeans.array(mfi);
        const Array4<      Real> surfcovmeansK = (nspec_surfcov>0) ? kahan.surfcovMeans.array(mfi) : Array4<Real>{};

        // update mean other values (primitive and conserved)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            GpuArray<Real,MAX_SPECIES> fracvec;
            cumeans(i,j,k,1) = Real(0.5)*(momxmeans(i,j,k) + momxmeans(i+1,j,k)); // jxmeans on CC
            cumeans(i,j,k,2) = Real(0.5)*(momymeans(i,j,k) + momymeans(i,j+1,k)); // jymeans on CC
            cumeans(i,j,k,3) = Real(0.5)*(momzmeans(i,j,k) + momzmeans(i,j,k+1)); // jzmeans on CC

            runningMeanKahan(cumeans(i,j,k,4), cumeansK(i,j,k,4), cu(i,j,k,4), stepsinv); //rhoEmeans

            for (int l=5; l<nvars; ++l) {
                runningMeanKahan(cumeans(i,j,k,l), cumeansK(i,j,k,l), cu(i,j,k,l), stepsinv); //rhoYkmeans
                fracvec[l-5] = cumeans(i,j,k,l)/cumeans(i,j,k,0); // Ykmeans
                primmeans(i,j,k,l+1) = fracvec[l-5]; // Ykmeans
            }
            for (int l=0; l<nspecies; ++l) {
                runningMeanKahan(primmeans(i,j,k,6+nspecies+l), primmeansK(i,j,k,6+nspecies+l), prim(i,j,k,6+nspecies+l), stepsinv); // Xkmeans
            }

            Real densitymeaninv = Real(1.0)/cumeans(i,j,k,0);
            primmeans(i,j,k,1) = densitymeaninv*Real(0.5)*(momxmeans(i,j,k) + momxmeans(i+1,j,k)); // velxmeans on CC
            primmeans(i,j,k,2) = densitymeaninv*Real(0.5)*(momymeans(i,j,k) + momymeans(i,j+1,k)); // velymeans on CC
            primmeans(i,j,k,3) = densitymeaninv*Real(0.5)*(momzmeans(i,j,k) + momzmeans(i,j,k+1)); // velzmeans on CC

            // mean COM velocity
            runningMean(primmeans(i,j,k,nprimvars+0), prim(i,j,k,1), stepsinv);
            runningMean(primmeans(i,j,k,nprimvars+1), prim(i,j,k,2), stepsinv);
            runningMean(primmeans(i,j,k,nprimvars+2), prim(i,j,k,3), stepsinv);

            primmeans(i,j,k,0) = cumeans(i,j,k,0); //rhomeans

            Real kinenergy = 0.;
            kinenergy += (momxmeans(i+1,j,k) + momxmeans(i,j,k))*(momxmeans(i+1,j,k) + momxmeans(i,j,k));
            kinenergy += (momymeans(i,j+1,k) + momymeans(i,j,k))*(momymeans(i,j+1,k) + momymeans(i,j,k));
            kinenergy += (momzmeans(i,j,k+1) + momzmeans(i,j,k))*(momzmeans(i,j,k+1) + momzmeans(i,j,k));
            kinenergy *= (Real(0.125)/cumeans(i,j,k,0));

            Real intenergy = (cumeans(i,j,k,4)-kinenergy)/cumeans(i,j,k,0);

            runningMeanKahan(primmeans(i,j,k,4), primmeansK(i,j,k,4), prim(i,j,k,4), stepsinv); // Tmean
            //GetTemperature(intenergy, fracvec, primmeans(i,j,k,4)); // Tmean
            //primmeans(i,j,k,5) = (primmeans(i,j,k,5)*stepsminusone + prim(i,j,k,5))*stepsinv; // Pmean
            GetPressureGas(primmeans(i,j,k,5), fracvec, cumeans(i,j,k,0), primmeans(i,j,k,4)); // Pmean

            if (nspec_surfcov>0) {
                for (int m=0; m<nspec_surfcov; ++m) {
                    runningMeanKahan(surfcovmeans(i,j,k,m), surfcovmeansK(i,j,k,m), surfcov(i,j,k,m), stepsinv);
                }
            }
        });
    }

    for ( MFIter mfi(prim_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<      Real> velxmeans = velMean[0].array(mfi);
        const Array4<      Real> velymeans = velMean[1].array(mfi);
        const Array4<      Real> velzmeans = velMean[2].array(mfi);

        const Array4<const Real> momxmeans = cumomMean[0].array(mfi);
        const Array4<const Real> momymeans = cumomMean[1].array(mfi);
        const Array4<const Real> momzmeans = cumomMean[2].array(mfi);

        const Array4<      Real> cumeans   = consMean.array(mfi);

        // update mean velocities
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real densitymeaninv = Real(2.0)/(cumeans(i-1,j,k,0)+cumeans(i,j,k,0));
            velxmeans(i,j,k) = momxmeans(i,j,k)*densitymeaninv;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real densitymeaninv = Real(2.0)/(cumeans(i,j-1,k,0)+cumeans(i,j,k,0));
            velymeans(i,j,k) = momymeans(i,j,k)*densitymeaninv;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real densitymeaninv = Real(2.0)/(cumeans(i,j,k-1,0)+cumeans(i,j,k,0));
            velzmeans(i,j,k) = momzmeans(i,j,k)*densitymeaninv;
        });

    } // end MFIter
}


/////////////////////////////////////
// evaluate variances and covariances
/////////////////////////////////////
void EvaluateVarsCoVarsMom3(const MultiFab& cons, const MultiFab& consMean, MultiFab& consVar,
                            const MultiFab& prim_in, const MultiFab& primMean, MultiFab& primVar,
                            const std::array<MultiFab, AMREX_SPACEDIM>& velMean,
                            std::array<MultiFab, AMREX_SPACEDIM>& velVar,
                            const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                            const std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                            std::array<MultiFab, AMREX_SPACEDIM>& cumomVar,
                            MultiFab& coVar,
                            MultiFab& mom3,
                            MultiFab& mom4,
                            const MultiFab& theta, const MultiFab& thetaMean,
                            MultiFab& thetaVar, MultiFab& thetacoVar,
                            StatsKahan& kahan,
                            const int steps)
{
    BL_PROFILE_VAR("EvaluateVarsCoVarsMom3()",EvaluateVarsCoVarsMom3);

    Real stepsinv = Real(1.0/static_cast<double>(steps));

    // Loop over boxes
    for ( MFIter mfi(cons,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<const Real> cu        = cons.array(mfi);
        const Array4<const Real> cumeans   = consMean.array(mfi);

        const Array4<const Real> velxmeans = velMean[0].array(mfi);
        const Array4<const Real> velymeans = velMean[1].array(mfi);
        const Array4<const Real> velzmeans = velMean[2].array(mfi);
        const Array4<      Real> velxvars  = velVar[0].array(mfi);
        const Array4<      Real> velyvars  = velVar[1].array(mfi);
        const Array4<      Real> velzvars  = velVar[2].array(mfi);

        const Array4<const Real> momx      = cumom[0].array(mfi);
        const Array4<const Real> momy      = cumom[1].array(mfi);
        const Array4<const Real> momz      = cumom[2].array(mfi);
        const Array4<const Real> momxmeans = cumomMean[0].array(mfi);
        const Array4<const Real> momymeans = cumomMean[1].array(mfi);
        const Array4<const Real> momzmeans = cumomMean[2].array(mfi);
        const Array4<      Real> momxvars  = cumomVar[0].array(mfi);
        const Array4<      Real> momyvars  = cumomVar[1].array(mfi);
        const Array4<      Real> momzvars  = cumomVar[2].array(mfi);

        const Array4<const Real> surfcov      = (nspec_surfcov>0) ? theta.array(mfi) : cons.array(mfi);
        const Array4<const Real> surfcovmeans = (nspec_surfcov>0) ? thetaMean.array(mfi) : consMean.array(mfi);
        const Array4<      Real> surfcovvars  = (nspec_surfcov>0) ? thetaVar.array(mfi) : consVar.array(mfi);
        const Array4<      Real> surfcovcoVars= (nspec_surfcov>0) ? thetacoVar.array(mfi) : consVar.array(mfi);

        const Array4<      Real> velxvarsK = kahan.velVars[0].array(mfi);
        const Array4<      Real> velyvarsK = kahan.velVars[1].array(mfi);
        const Array4<      Real> velzvarsK = kahan.velVars[2].array(mfi);
        const Array4<      Real> momxvarsK = kahan.cumomVars[0].array(mfi);
        const Array4<      Real> momyvarsK = kahan.cumomVars[1].array(mfi);
        const Array4<      Real> momzvarsK = kahan.cumomVars[2].array(mfi);

        // update momentum and velocity variances
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real deljx = momx(i,j,k) - momxmeans(i,j,k);
            runningMeanKahan(momxvars(i,j,k), momxvarsK(i,j,k), deljx*deljx, stepsinv); // <jx jx>

            Real densitymeaninv = Real(2.0)/(cumeans(i-1,j,k,0)+cumeans(i,j,k,0));
            Real delrho = Real(0.5)*(cu(i-1,j,k,0) + cu(i,j,k,0)) - Real(0.5)*(cumeans(i-1,j,k,0) + cumeans(i,j,k,0));
            Real delvelx = (deljx - velxmeans(i,j,k)*delrho)*densitymeaninv;
            runningMeanKahan(velxvars(i,j,k), velxvarsK(i,j,k), delvelx*delvelx, stepsinv); // <vx vx>
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real deljy = momy(i,j,k) - momymeans(i,j,k);
            runningMeanKahan(momyvars(i,j,k), momyvarsK(i,j,k), deljy*deljy, stepsinv); // <jy jy>

            Real densitymeaninv = Real(2.0)/(cumeans(i,j-1,k,0)+cumeans(i,j,k,0));
            Real delrho = Real(0.5)*(cu(i,j-1,k,0) + cu(i,j,k,0)) - Real(0.5)*(cumeans(i,j-1,k,0) + cumeans(i,j,k,0));
            Real delvely = (deljy - velymeans(i,j,k)*delrho)*densitymeaninv;
            runningMeanKahan(velyvars(i,j,k), velyvarsK(i,j,k), delvely*delvely, stepsinv); // <vx vx>
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            Real deljz = momz(i,j,k) - momzmeans(i,j,k);
            runningMeanKahan(momzvars(i,j,k), momzvarsK(i,j,k), deljz*deljz, stepsinv); // <jz jz>

            Real densitymeaninv = Real(2.0)/(cumeans(i,j,k-1,0)+cumeans(i,j,k,0));
            Real delrho = Real(0.5)*(cu(i,j,k-1,0) + cu(i,j,k,0)) - Real(0.5)*(cumeans(i,j,k-1,0) + cumeans(i,j,k,0));
            Real delvelz = (deljz - velzmeans(i,j,k)*delrho)*densitymeaninv;
            runningMeanKahan(velzvars(i,j,k), velzvarsK(i,j,k), delvelz*delvelz, stepsinv); // <vz vz>
        });

    }

    // Loop over boxes
    for ( MFIter mfi(prim_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<const Real> cu        = cons.array(mfi);
        const Array4<const Real> cumeans   = consMean.array(mfi);
        const Array4<      Real> cuvars    = consVar.array(mfi);
        const Array4<const Real> prim      = prim_in.array(mfi);
        const Array4<const Real> primmeans = primMean.array(mfi);
        const Array4<      Real> primvars  = primVar.array(mfi);

        const Array4<      Real> covars    = coVar.array(mfi);

        const Array4<const Real> velxmeans = velMean[0].array(mfi);
        const Array4<const Real> velymeans = velMean[1].array(mfi);
        const Array4<const Real> velzmeans = velMean[2].array(mfi);

        const Array4<const Real> momx      = cumom[0].array(mfi);
        const Array4<const Real> momy      = cumom[1].array(mfi);
        const Array4<const Real> momz      = cumom[2].array(mfi);
        const Array4<const Real> momxmeans = cumomMean[0].array(mfi);
        const Array4<const Real> momymeans = cumomMean[1].array(mfi);
        const Array4<const Real> momzmeans = cumomMean[2].array(mfi);

        const Array4<const Real> surfcov      = (nspec_surfcov>0) ? theta.array(mfi) : cons.array(mfi);
        const Array4<const Real> surfcovmeans = (nspec_surfcov>0) ? thetaMean.array(mfi) : consMean.array(mfi);
        const Array4<      Real> surfcovvars  = (nspec_surfcov>0) ? thetaVar.array(mfi) : consVar.array(mfi);
        const Array4<      Real> surfcovcoVars= (nspec_surfcov>0) ? thetacoVar.array(mfi) : consVar.array(mfi);

        const Array4<      Real> cuvarsK      = kahan.cuVars.array(mfi);
        const Array4<      Real> primvarsK    = kahan.primVars.array(mfi);
        const Array4<      Real> covarsK      = kahan.coVars.array(mfi);
        const Array4<      Real> surfcovvarsK = (nspec_surfcov>0) ? kahan.surfcovVars.array(mfi) : Array4<Real>{};

        // other cell-centered variances and covariances (temperature fluctuation and variance)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            // conserved variable variances (rho, rhoE, rhoYk)
            Real delrho = cu(i,j,k,0) - cumeans(i,j,k,0);
            runningMeanKahan(cuvars(i,j,k,0), cuvarsK(i,j,k,0), delrho*delrho, stepsinv); // <rho rho>

            Real delenergy = cu(i,j,k,4) - cumeans(i,j,k,4);
            runningMeanKahan(cuvars(i,j,k,4), cuvarsK(i,j,k,4), delenergy*delenergy, stepsinv); // <rhoE rhoE>

            for (int l=5; l<nvars; ++l) {
                Real delmassvec = cu(i,j,k,l) - cumeans(i,j,k,l);
                runningMeanKahan(cuvars(i,j,k,l), cuvarsK(i,j,k,l), delmassvec*delmassvec, stepsinv); // <rhoYk rhoYk>
            }

            // del rhoYk --> delYk
            GpuArray<Real,MAX_SPECIES> delrhoYk;
            GpuArray<Real,MAX_SPECIES> Ykmean;
            GpuArray<Real,MAX_SPECIES> delYk;
            for (int ns=0; ns<nspecies; ++ns) {
                delrhoYk[ns] = cu(i,j,k,5+ns) - cumeans(i,j,k,5+ns); // delrhoYk
                Ykmean[ns] = primmeans(i,j,k,6+ns); // Ykmean
                delYk[ns] = (delrhoYk[ns] - Ykmean[ns]*delrho)/cumeans(i,j,k,0); // delYk = (delrhoYk - Ykmean*delrho)/rhomean
                runningMeanKahan(primvars(i,j,k,6+ns), primvarsK(i,j,k,6+ns), delYk[ns]*delYk[ns], stepsinv);
            }

            // primitive variable variances (rho)
            primvars(i,j,k,0) = cuvars(i,j,k,0);

            Real vx = Real(0.5)*(velxmeans(i,j,k) + velxmeans(i+1,j,k));
            Real vy = Real(0.5)*(velymeans(i,j,k) + velymeans(i,j+1,k));
            Real vz = Real(0.5)*(velzmeans(i,j,k) + velzmeans(i,j,k+1));
            Real T = primmeans(i,j,k,4);
            Real densitymeaninv = Real(1.0)/cumeans(i,j,k,0);

            Real cv = 0.;
            for (int l=0; l<nspecies; ++l) {
                cv = cv + hcv[l]*cumeans(i,j,k,5+l)/cumeans(i,j,k,0);
            }
            Real cvinv = Real(1.0)/cv;

            Real qmean = cv*T-Real(0.5)*(vx*vx + vy*vy + vz*vz);

            Real deljx = Real(0.5)*(momx(i,j,k) + momx(i+1,j,k)) - Real(0.5)*(momxmeans(i,j,k) + momxmeans(i+1,j,k));
            Real deljy = Real(0.5)*(momy(i,j,k) + momy(i,j+1,k)) - Real(0.5)*(momymeans(i,j,k) + momymeans(i,j+1,k));
            Real deljz = Real(0.5)*(momz(i,j,k) + momz(i,j,k+1)) - Real(0.5)*(momzmeans(i,j,k) + momzmeans(i,j,k+1));

            runningMeanKahan(cuvars(i,j,k,1), cuvarsK(i,j,k,1), deljx*deljx, stepsinv); // <jx jx> on CC
            runningMeanKahan(cuvars(i,j,k,2), cuvarsK(i,j,k,2), deljy*deljy, stepsinv); // <jy jy>  on CC
            runningMeanKahan(cuvars(i,j,k,3), cuvarsK(i,j,k,3), deljz*deljz, stepsinv); // <jz jz> on CC

            Real delvelx = (deljx - vx*delrho)*densitymeaninv;
            Real delvely = (deljy - vy*delrho)*densitymeaninv;
            Real delvelz = (deljz - vz*delrho)*densitymeaninv;

            runningMeanKahan(primvars(i,j,k,1), primvarsK(i,j,k,1), delvelx*delvelx, stepsinv); // <vx vx> on CC
            runningMeanKahan(primvars(i,j,k,2), primvarsK(i,j,k,2), delvely*delvely, stepsinv); // <vy vy>  on CC
            runningMeanKahan(primvars(i,j,k,3), primvarsK(i,j,k,3), delvelz*delvelz, stepsinv); // <vz vz> on CC

            Real delg = vx*deljx + vy*deljy + vz*deljz;

            runningMeanKahan(primvars(i,j,k,nprimvars+0), primvarsK(i,j,k,nprimvars+0), delg*delg, stepsinv); // gvar
            runningMean(primvars(i,j,k,nprimvars+1), delg*delenergy, stepsinv); // kgcross
            runningMeanKahan(primvars(i,j,k,nprimvars+2), primvarsK(i,j,k,nprimvars+2), delrho*delenergy, stepsinv); // krcross
            runningMean(primvars(i,j,k,nprimvars+3), delrho*delg, stepsinv); // rgcross

            runningMeanKahan(primvars(i,j,k,nprimvars+4), primvarsK(i,j,k,nprimvars+4), cvinv*cvinv*densitymeaninv*densitymeaninv*
                             (cuvars(i,j,k,4) + primvars(i,j,k,nprimvars+0) - 2*primvars(i,j,k,nprimvars+1)
                              + qmean*(qmean*cuvars(i,j,k,0) - 2*primvars(i,j,k,nprimvars+2) + 2*primvars(i,j,k,nprimvars+3))), stepsinv); // <T T>

            // use instantaneous value (the above presumably does not work for multispecies)
            Real delT = prim(i,j,k,4) - primmeans(i,j,k,4);
            runningMeanKahan(primvars(i,j,k,4), primvarsK(i,j,k,4), delT*delT, stepsinv); // <T T> -- direct

            //Real deltemp = (delenergy - delg - qmean*delrho)*cvinv*densitymeaninv;
            Real deltemp = delT;

            runningMean(covars(i,j,k,0), delrho*deljx, stepsinv); // <rho jx>
            runningMean(covars(i,j,k,1), delrho*deljy, stepsinv); // <rho jy>
            runningMean(covars(i,j,k,2), delrho*deljz, stepsinv); // <rho jz>
            runningMean(covars(i,j,k,3), deljx*deljy, stepsinv); // <jx jy>
            runningMean(covars(i,j,k,4), deljy*deljz, stepsinv); // <jy jz>
            runningMean(covars(i,j,k,5), deljx*deljz, stepsinv); // <jx jz>
            runningMeanKahan(covars(i,j,k,6), covarsK(i,j,k,6), delrho*delenergy, stepsinv); // <rho rhoE>
            runningMean(covars(i,j,k,7), deljx*delenergy, stepsinv); // <jx rhoE>
            runningMean(covars(i,j,k,8), deljy*delenergy, stepsinv); // <jy rhoE>
            runningMean(covars(i,j,k,9), deljz*delenergy, stepsinv); // <jz rhoE>
            runningMeanKahan(covars(i,j,k,10), covarsK(i,j,k,10), delrhoYk[0]*delrhoYk[nspecies-1], stepsinv); // <rhoYklightest rhoYkheaviest>
            runningMean(covars(i,j,k,11), delrho*delvelx, stepsinv); // <rho vx>
            runningMean(covars(i,j,k,12), delrho*delvely, stepsinv); // <rho vy>
            runningMean(covars(i,j,k,13), delrho*delvelz, stepsinv); // <rho vz>
            runningMean(covars(i,j,k,14), delvelx*delvely, stepsinv); // <vx vy>
            runningMean(covars(i,j,k,15), delvely*delvelz, stepsinv); // <vy vz>
            runningMean(covars(i,j,k,16), delvelx*delvelz, stepsinv); // <vx vz>
            runningMean(covars(i,j,k,17), delrho*deltemp, stepsinv); // <rho T>
            runningMean(covars(i,j,k,18), delvelx*deltemp, stepsinv); // <vx T>
            runningMean(covars(i,j,k,19), delvely*deltemp, stepsinv); // <vy T>
            runningMean(covars(i,j,k,20), delvelz*deltemp, stepsinv); // <vz T>
            runningMeanKahan(covars(i,j,k,21), covarsK(i,j,k,21), delYk[0]*delYk[nspecies-1], stepsinv); // <Yklightest Ykheaviest>
            runningMean(covars(i,j,k,22), delYk[0]*delvelx, stepsinv); // <Yklightest velx>
            runningMean(covars(i,j,k,23), delYk[nspecies-1]*delvelx, stepsinv); // <Ykheaviest velx>
            runningMean(covars(i,j,k,24), delrhoYk[0]*delvelx, stepsinv); // <rhoYklightest velx>
            runningMean(covars(i,j,k,25), delrhoYk[nspecies-1]*delvelx, stepsinv); // <rhoYkheaviest velx>

            if (nspec_surfcov>0) {
                for (int m=0; m<nspec_surfcov; ++m) {
                    Real delsurfcov = surfcov(i,j,k,m) - surfcovmeans(i,j,k,m);
                    runningMeanKahan(surfcovvars(i,j,k,m), surfcovvarsK(i,j,k,m), delsurfcov*delsurfcov, stepsinv);
                    runningMean(surfcovcoVars(i,j,k,6*m), delrhoYk[m]*deltemp, stepsinv);        // <rhoYk T>
                    runningMean(surfcovcoVars(i,j,k,6*m+1), delsurfcov*delrhoYk[m], stepsinv); // <theta rhoYk>
                    runningMean(surfcovcoVars(i,j,k,6*m+2), delsurfcov*delvelx, stepsinv);     // <theta vx>
                    runningMean(surfcovcoVars(i,j,k,6*m+3), delsurfcov*delvely, stepsinv);     // <theta vy>
                    runningMean(surfcovcoVars(i,j,k,6*m+4), delsurfcov*delvelz, stepsinv);     // <theta vz>
                    runningMean(surfcovcoVars(i,j,k,6*m+5), delsurfcov*deltemp, stepsinv);     // <theta T>
                }
            }
        });

    } // end MFIter

    // Loop again for third moment calculations
    if (plot_mom3 || plot_mom4) {

        // Loop over boxes
        for ( MFIter mfi(prim_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            const Array4<const Real> cu        = cons.array(mfi);
            const Array4<const Real> cumeans   = consMean.array(mfi);
            const Array4<const Real> prim      = prim_in.array(mfi);
            const Array4<const Real> primmeans = primMean.array(mfi);

            const Array4<      Real> m3        = mom3.array(mfi);
            const Array4<      Real> m4        = mom4.array(mfi);
            const Array4<      Real> m4K       = (plot_mom4) ? kahan.mom4.array(mfi) : Array4<Real>{};

            const Array4<const Real> velxmeans = velMean[0].array(mfi);
            const Array4<const Real> velymeans = velMean[1].array(mfi);
            const Array4<const Real> velzmeans = velMean[2].array(mfi);

            const Array4<const Real> momx      = cumom[0].array(mfi);
            const Array4<const Real> momy      = cumom[1].array(mfi);
            const Array4<const Real> momz      = cumom[2].array(mfi);
            const Array4<const Real> momxmeans = cumomMean[0].array(mfi);
            const Array4<const Real> momymeans = cumomMean[1].array(mfi);
            const Array4<const Real> momzmeans = cumomMean[2].array(mfi);

            // other cell-centered variances and covariances (temperature fluctuation and variance)
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // conserved variable variances (rho, rhoE, rhoYk)
                Real delrho = cu(i,j,k,0) - cumeans(i,j,k,0);
                if (plot_mom3) runningMean(m3(i,j,k,0), delrho*delrho*delrho, stepsinv); // <rho rho rho>
                if (plot_mom4) runningMeanKahan(m4(i,j,k,0), m4K(i,j,k,0), delrho*delrho*delrho*delrho, stepsinv); // <rho rho rho>

                Real delenergy = cu(i,j,k,4) - cumeans(i,j,k,4);
                if (plot_mom3) runningMean(m3(i,j,k,4), delenergy*delenergy*delenergy, stepsinv); // <rhoE rhoE rhoE>
                if (plot_mom4) runningMeanKahan(m4(i,j,k,4), m4K(i,j,k,4), delenergy*delenergy*delenergy*delenergy, stepsinv); // <rhoE rhoE rhoE>

                // del rhoYk --> delYk
                GpuArray<Real,MAX_SPECIES> delrhoYk;
                GpuArray<Real,MAX_SPECIES> Ykmean;
                GpuArray<Real,MAX_SPECIES> delYk;
                for (int ns=0; ns<nspecies; ++ns) {
                    delrhoYk[ns] = cu(i,j,k,5+ns) - cumeans(i,j,k,5+ns); // delrhoYk
                    Ykmean[ns] = primmeans(i,j,k,6+ns); // Ykmean
                    delYk[ns] = (delrhoYk[ns] - Ykmean[ns]*delrho)/cumeans(i,j,k,0); // delYk = (delrhoYk - Ykmean*delrho)/rhomean
                    if (plot_mom3) runningMean(m3(i,j,k,5+ns), delYk[ns]*delYk[ns]*delYk[ns], stepsinv);
                    if (plot_mom4) runningMeanKahan(m4(i,j,k,5+ns), m4K(i,j,k,5+ns), delYk[ns]*delYk[ns]*delYk[ns]*delYk[ns], stepsinv);
                }

                Real vx = Real(0.5)*(velxmeans(i,j,k) + velxmeans(i+1,j,k));
                Real vy = Real(0.5)*(velymeans(i,j,k) + velymeans(i,j+1,k));
                Real vz = Real(0.5)*(velzmeans(i,j,k) + velzmeans(i,j,k+1));
                Real T = primmeans(i,j,k,4);
                Real densitymeaninv = Real(1.0)/cumeans(i,j,k,0);

                Real deljx = Real(0.5)*(momx(i,j,k) + momx(i+1,j,k)) - Real(0.5)*(momxmeans(i,j,k) + momxmeans(i+1,j,k));
                Real deljy = Real(0.5)*(momy(i,j,k) + momy(i,j+1,k)) - Real(0.5)*(momymeans(i,j,k) + momymeans(i,j+1,k));
                Real deljz = Real(0.5)*(momz(i,j,k) + momz(i,j,k+1)) - Real(0.5)*(momzmeans(i,j,k) + momzmeans(i,j,k+1));

                Real delvelx = (deljx - vx*delrho)*densitymeaninv;
                Real delvely = (deljy - vy*delrho)*densitymeaninv;
                Real delvelz = (deljz - vz*delrho)*densitymeaninv;

                if (plot_mom3) {
                    runningMean(m3(i,j,k,1), delvelx*delvelx*delvelx, stepsinv);
                    runningMean(m3(i,j,k,2), delvely*delvely*delvely, stepsinv);
                    runningMean(m3(i,j,k,3), delvelz*delvelz*delvelz, stepsinv);
                }

                if (plot_mom4) {
                    runningMeanKahan(m4(i,j,k,1), m4K(i,j,k,1), delvelx*delvelx*delvelx*delvelx, stepsinv);
                    runningMeanKahan(m4(i,j,k,2), m4K(i,j,k,2), delvely*delvely*delvely*delvely, stepsinv);
                    runningMeanKahan(m4(i,j,k,3), m4K(i,j,k,3), delvelz*delvelz*delvelz*delvelz, stepsinv);
                }

                // use instantaneous value (the above presumably does not work for multispecies)
                Real delT = prim(i,j,k,4) - primmeans(i,j,k,4);
                if (plot_mom3) runningMean(m3(i,j,k,5+nspecies), delT*delT*delT, stepsinv);
                if (plot_mom4) runningMeanKahan(m4(i,j,k,5+nspecies), m4K(i,j,k,5+nspecies), delT*delT*delT*delT, stepsinv);
            });
        } // end MFIter
    } // end plot_mom3 || plot_mom4
}


//////////////////////////////////////////
// Get Pencil values at x* ///////////////
// for all j and k ///////////////////////
// ///////////////////////////////////////
void GetPencilCross(amrex::Gpu::DeviceVector<Real>& data_xcross_in,
                    const MultiFab& consMean,
                    const MultiFab& primMean,
                    const MultiFab& prim_in,
                    const MultiFab& cons,
                    const int nstats,
                    const int x_star)
{
    BL_PROFILE_VAR("GetPencilCross()",GetPencilCross);

    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const Array4<const Real> cumeans   = consMean.array(mfi);
        const Array4<const Real> primmeans = primMean.array(mfi);
        const Array4<const Real> prim      = prim_in.array(mfi);
        const Array4<const Real> cu        = cons.array(mfi);

        Real* data_xcross = data_xcross_in.data();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i==x_star) {
                int index = k*n_cells[1]*nstats + j*nstats;
                data_xcross[index + 0]  = cu(i,j,k,0);          // rho-instant
                data_xcross[index + 1]  = cumeans(i,j,k,0);     // rho-mean
                data_xcross[index + 2]  = cu(i,j,k,4);          // energy-instant
                data_xcross[index + 3]  = cumeans(i,j,k,4);     // energy-mean
                data_xcross[index + 4]  = cu(i,j,k,1);          // jx-instant
                data_xcross[index + 5]  = cumeans(i,j,k,1);     // jx-mean
                data_xcross[index + 6]  = cu(i,j,k,2);          // jy-instant
                data_xcross[index + 7]  = cumeans(i,j,k,2);     // jy-mean
                data_xcross[index + 8]  = cu(i,j,k,3);          // jz-instant
                data_xcross[index + 9]  = cumeans(i,j,k,3);     // jz-mean
                data_xcross[index + 10] = prim(i,j,k,1);        // velx-instant
                data_xcross[index + 11] = primmeans(i,j,k,1);   // velx-mean
                data_xcross[index + 12] = prim(i,j,k,2);        // vely-instant
                data_xcross[index + 13] = primmeans(i,j,k,2);   // vely-mean
                data_xcross[index + 14] = prim(i,j,k,3);        // velz-instant
                data_xcross[index + 15] = primmeans(i,j,k,3);   // velz-mean
                data_xcross[index + 16] = prim(i,j,k,4);        // T-instant
                data_xcross[index + 17] = primmeans(i,j,k,4);   // T-mean
                for (int ns=0; ns<nspecies; ++ns) {
                    data_xcross[index + 18+4*ns+0]  = cu(i,j,k,5+ns);        // rhoYk-instant
                    data_xcross[index + 18+4*ns+1]  = cumeans(i,j,k,5+ns);   // rhoYk-mean
                    data_xcross[index + 18+4*ns+2]  = prim(i,j,k,6+ns);      // Yk-instant
                    data_xcross[index + 18+4*ns+3]  = primmeans(i,j,k,6+ns); // Yk-mean
                }
           }

        }); // end MFITer
    }
}


//////////////////////////////////////////
// Update Spatial Correlations ///////////
// ///////////////////////////////////////
void EvaluateSpatialCorrelations3D(Vector<Real>& spatialCross,
                                   Vector<Real>& spatialCrossK,
                                   Vector<Real>& data_xcross,
                                   amrex::Gpu::HostVector<Real>& cu_avg,
                                   amrex::Gpu::HostVector<Real>& cumeans_avg,
                                   amrex::Gpu::HostVector<Real>& prim_avg,
                                   amrex::Gpu::HostVector<Real>& primmeans_avg,
                                   amrex::Gpu::HostVector<Real>& delcu_avg,
                                   amrex::Gpu::HostVector<Real>& delprim_avg,
                                   const int steps,
                                   const int /*nstats*/,
                                   const int ncross)
{

    BL_PROFILE_VAR("EvaluateSpatialCorrelations3D()",EvaluateSpatialCorrelations3D);

    // delcu_avg / delprim_avg are plane averages of the per-cell fluctuations (cons - consMean,
    // prim - primMean), with the layout of cu_avg / prim_avg. Fluctuations are taken from them
    // rather than as differences of instantaneous and mean plane averages (cancellation in FP32).

    Real stepsinv = Real(1.0/static_cast<double>(steps));

    int nprims = nspecies + 4;

    // Fill in data_xcross
    data_xcross[0]  = cu_avg[0+nvars*cross_cell];                                 // rho-instant
    data_xcross[1]  = cumeans_avg[0+nvars*cross_cell];                            // rho-mean
    data_xcross[2]  = cu_avg[4+nvars*cross_cell];                                 // energy-instant
    data_xcross[3]  = cumeans_avg[4+nvars*cross_cell];                            // energy-mean
    data_xcross[4]  = cu_avg[1+nvars*cross_cell];                                 // jx-instant
    data_xcross[5]  = cumeans_avg[1+nvars*cross_cell];                            // jx-mean
    data_xcross[6]  = cu_avg[2+nvars*cross_cell];                                 // jy-instant
    data_xcross[7]  = cumeans_avg[2+nvars*cross_cell];                            // jy-mean
    data_xcross[8]  = cu_avg[3+nvars*cross_cell];                                 // jz-instant
    data_xcross[9]  = cumeans_avg[3+nvars*cross_cell];                            // jz-mean
    data_xcross[10] = prim_avg[0+nprims*cross_cell];                              // velx-instant
    data_xcross[11] = primmeans_avg[0+nprims*cross_cell];                         // velx-mean
    data_xcross[12] = prim_avg[1+nprims*cross_cell];                              // vely-instant
    data_xcross[13] = primmeans_avg[1+nprims*cross_cell];                         // vely-mean
    data_xcross[14] = prim_avg[2+nprims*cross_cell];                              // velz-instant
    data_xcross[15] = primmeans_avg[2+nprims*cross_cell];                         // velz-mean
    data_xcross[16] = prim_avg[3+nprims*cross_cell];                              // T-instant
    data_xcross[17] = primmeans_avg[3+nprims*cross_cell];                         // T-mean
    for (int ns=0; ns<nspecies; ++ns) {
        data_xcross[18+4*ns+0] = cu_avg[5+ns+nvars*cross_cell];                   // rhoYk-instant
        data_xcross[18+4*ns+1] = cumeans_avg[5+ns+nvars*cross_cell];              // rhoYk-mean
        data_xcross[18+4*ns+2] = prim_avg[4+ns+nprims*cross_cell];                // Yk-instant
        data_xcross[18+4*ns+3] = primmeans_avg[4+ns+nprims*cross_cell];           // Yk-mean
    }

    // Get mean values
    Real meanrhocross = data_xcross[1];
    Vector<Real>  meanYkcross(nspecies, 0.0);
    for (int ns=0; ns<nspecies; ++ns) {
        meanYkcross[ns] =  data_xcross[18+4*ns+3];
    }
    Real meanuxcross = data_xcross[11];

    // Get fluctuations of the conserved variables at the cross cell
    Real delrhocross = delcu_avg[0+nvars*cross_cell];
    Real delKcross   = delcu_avg[4+nvars*cross_cell];
    Real deljxcross  = delcu_avg[1+nvars*cross_cell];
    Real deljycross  = delcu_avg[2+nvars*cross_cell];
    Real deljzcross  = delcu_avg[3+nvars*cross_cell];
    Vector<Real>  delrhoYkcross(nspecies, 0.0);
    for (int ns=0; ns<nspecies; ++ns) {
        delrhoYkcross[ns] =  delcu_avg[5+ns+nvars*cross_cell];
    }

    // Get fluctuations of some primitive variables (for direct fluctuation calculations)
    Real delTcross = delprim_avg[3+nprims*cross_cell];
    Real delvxcross = delprim_avg[0+nprims*cross_cell];
    Vector<Real>  delYkcross(nspecies, 0.0);
    for (int ns=0; ns<nspecies; ++ns) {
        delYkcross[ns] =  delprim_avg[4+ns+nprims*cross_cell];
    }

    // evaluate heat stuff at the cross cell
    Real cvcross = 0.;
    for (int l=0; l<nspecies; ++l) {
        cvcross = cvcross + hcv[l]*data_xcross[18+4*l+1]/data_xcross[1];
    }
    Real cvinvcross = Real(1.0)/cvcross;
    Real qmeancross = cvcross*data_xcross[17] -
                     Real(0.5)*(data_xcross[11]*data_xcross[11] + data_xcross[13]*data_xcross[13] + data_xcross[15]*data_xcross[15]);

    // Get fluctuations of derived hydrodynamic quantities at the cross cell
    // delG = \vec{v}\cdot\vec{\deltaj}
    Real delGcross = data_xcross[11]*deljxcross + data_xcross[13]*deljycross +
                    data_xcross[15]*deljzcross;

    /////////////////////////////////////////////////////////////
    // evaluate x-spatial correlations
    /////////////////////////////////////////////////////////////
    // int ncross = 37+nspecies+2; check main_drive.cpp for latest
    for (int i=0; i<n_cells[0]; ++i) {

        // Get mean values
        Real meanrho = cumeans_avg[0+nvars*i];
        Vector<Real>  meanYk(nspecies, 0.0);
        for (int ns=0; ns<nspecies; ++ns) {
            meanYk[ns] = primmeans_avg[4+ns+nprims*i];
        }

        // Get fluctuations of the conserved variables
        Real delrho = delcu_avg[0+nvars*i];
        Real delK   = delcu_avg[4+nvars*i];
        Real deljx  = delcu_avg[1+nvars*i];
        Real deljy  = delcu_avg[2+nvars*i];
        Real deljz  = delcu_avg[3+nvars*i];
        Vector<Real>  delrhoYk(nspecies, 0.0);
        for (int ns=0; ns<nspecies; ++ns) {
            delrhoYk[ns] = delcu_avg[5+ns+nvars*i];
        }

        // Get fluctuations of some primitive variables (for direct fluctuation calculations)
        Real delT = delprim_avg[3+nprims*i];
        Real delvx = delprim_avg[0+nprims*i];
        Vector<Real>  delYk(nspecies, 0.0);
        for (int ns=0; ns<nspecies; ++ns) {
            delYk[ns] = delprim_avg[4+ns+nprims*i];
        }

        // evaluate heat stuff at the cross cell
        Real cv = 0.;
        for (int l=0; l<nspecies; ++l) {
            cv = cv + hcv[l]*cumeans_avg[5+l+nvars*i]/cumeans_avg[0+nvars*i];
        }
        Real cvinv = Real(1.0)/cv;
        Real qmean = cv*primmeans_avg[3+nprims*i] -
                         Real(0.5)*(primmeans_avg[0+nprims*i]*primmeans_avg[0+nprims*i] +
                              primmeans_avg[1+nprims*i]*primmeans_avg[1+nprims*i] +
                              primmeans_avg[2+nprims*i]*primmeans_avg[2+nprims*i]);

        // Get fluctuations of derived hydrodynamic quantities
        // delG = \vec{v}\cdot\vec{\deltaj}
        Real delG = primmeans_avg[0+nprims*i]*deljx + primmeans_avg[1+nprims*i]*deljy +
                    primmeans_avg[2+nprims*i]*deljz;

        // First update correlations of conserved quantities (we will do rhoYk later)
        runningMeanKahan(spatialCross[i*ncross+0], spatialCrossK[i*ncross+0], delrhocross*delrho, stepsinv); // <delrho(x*)delrho(x)>
        runningMeanKahan(spatialCross[i*ncross+1], spatialCrossK[i*ncross+1], delKcross*delK, stepsinv);     // <delK(x*)delK(x)>
        runningMeanKahan(spatialCross[i*ncross+2], spatialCrossK[i*ncross+2], deljxcross*deljx, stepsinv);   // <deljx(x*)deljx(x)>
        runningMeanKahan(spatialCross[i*ncross+3], spatialCrossK[i*ncross+3], deljycross*deljy, stepsinv);   // <deljy(x*)deljy(x)>
        runningMeanKahan(spatialCross[i*ncross+4], spatialCrossK[i*ncross+4], deljzcross*deljz, stepsinv);   // <deljz(x*)deljz(x)>
        runningMean(spatialCross[i*ncross+5], deljxcross*delrho, stepsinv);  // <deljx(x*)delrho(x)>
        runningMean(spatialCross[i*ncross+6], deljxcross*delrhoYk[0], stepsinv);  // <deljx(x*)delrhoYkL(x)>
        runningMean(spatialCross[i*ncross+7], deljxcross*delrhoYk[nspecies-1], stepsinv);  // <deljx(x*)delrhoYkH(x)>
        runningMeanKahan(spatialCross[i*ncross+8], spatialCrossK[i*ncross+8], delrhocross*delrhoYk[0], stepsinv); // <delrho(x*)delrhoYkL(x)>
        runningMeanKahan(spatialCross[i*ncross+9], spatialCrossK[i*ncross+9], delrhocross*delrhoYk[nspecies-1], stepsinv); // <delrho(x*)delrhoYkH(x)>
        runningMeanKahan(spatialCross[i*ncross+10], spatialCrossK[i*ncross+10], delrhoYkcross[0]*delrho, stepsinv); // <delrhoYkL(x*)delrho(x)>
        runningMeanKahan(spatialCross[i*ncross+11], spatialCrossK[i*ncross+11], delrhoYkcross[nspecies-1]*delrho, stepsinv); // <delrhoYkH(x*)delrho(x)>

        // Some more cross-correlations for hydrodynamical variables later
        runningMeanKahan(spatialCross[i*ncross+12], spatialCrossK[i*ncross+12], delGcross*delG, stepsinv); // <delG(x*)delG(x)>
        runningMean(spatialCross[i*ncross+13], delGcross*delK, stepsinv); // <delG(x*)delK(x)>
        runningMean(spatialCross[i*ncross+14], delKcross*delG, stepsinv); // <delK(x*)delG(x)>
        runningMeanKahan(spatialCross[i*ncross+15], spatialCrossK[i*ncross+15], delrhocross*delK, stepsinv); // <delrho(x*)delK(x)>
        runningMeanKahan(spatialCross[i*ncross+16], spatialCrossK[i*ncross+16], delKcross*delrho, stepsinv); // <delK(x*)delrho(x)>
        runningMean(spatialCross[i*ncross+17], delrhocross*delG, stepsinv); // <delrho(x*)delG(x)>
        runningMean(spatialCross[i*ncross+18], delGcross*delrho, stepsinv); // <delG(x*)delrho(x)>

        // Next we do cross-correlations with and between hydrodynamical variables
        // <delT(x*)delT(x)> = (1/cv*/cv/<rho(x)>/<rho(x*)>)(<delK*delK> + <delG*delG> - <delG*delK> - <delK*delG>
        //                      + <Q><Q*><delrho*delrho> - <Q*><delrho*delK> - <Q><delK*delrho> + <Q*><delrho*delG> + <Q><delG*delrho>)
        spatialCross[i*ncross+19] = (cvinvcross*cvinv/(meanrhocross*meanrho))*
                                    (spatialCross[i*ncross+1] + spatialCross[i*ncross+12] - spatialCross[i*ncross+13] - spatialCross[i*ncross+14]
                                     + qmean*qmeancross*spatialCross[i*ncross+0] - qmeancross*spatialCross[i*ncross+15] - qmean*spatialCross[i*ncross+16]
                                     + qmeancross*spatialCross[i*ncross+17] + qmean*spatialCross[i*ncross+18]);

        // <delT(x*)delrho(x)> = (1/cv/<rho(x*)>)*(<delK*delrho> - <delG*delrho> - <Q*><delrhodelrho*>)
        spatialCross[i*ncross+20] = (cvinvcross/meanrhocross)*(spatialCross[i*ncross+16] - spatialCross[i*ncross+18] - qmeancross*spatialCross[i*ncross+0]);

        // <delu(x*)delrho> = (1/<rho(x*)>)*(<deljx(x*)delrho(x)> - <u(x*)><<delrho(x*)delrho(x)>)
        spatialCross[i*ncross+21] = (Real(1.0)/meanrhocross)*(spatialCross[i*ncross+5] - meanuxcross*spatialCross[i*ncross+0]);

        // <delu(x*)del(rhoYkL)> = (1/<rho(x*)>)*(<deljx(x*)del(rhoYkL)> - <u(x*)><delrho(x*)del(rhoYkL)>)
        spatialCross[i*ncross+22] = (Real(1.0)/meanrhocross)*(spatialCross[i*ncross+6] - meanuxcross*spatialCross[i*ncross+8]);

        // <delu(x*)del(rhoYkH)> = (1/<rho(x*)>)*(<deljx(x*)del(rhoYkH)> - <u(x*)><delrho(x*)del(rhoYkH)>)
        spatialCross[i*ncross+23] = (Real(1.0)/meanrhocross)*(spatialCross[i*ncross+7] - meanuxcross*spatialCross[i*ncross+9]);

        // <delu(x*)del(YkL)> = (1/<rho(x*)>/<rho(x)>)*(<deljx(x*)del(rhoYkL) - <u(x*)><delrho(x*)del(rhoYkL)>
        //                      - <YkL(x)><deljx(x*)delrho(x)> + <u(x*)><YkL(x)><delrho(x*)delrho(x)>)
        spatialCross[i*ncross+24] = (Real(1.0)/(meanrho*meanrhocross))*(spatialCross[i*ncross+6] - meanuxcross*spatialCross[i*ncross+8]
                                                                     - meanYk[0]*spatialCross[i*ncross+5] + meanuxcross*meanYk[0]*spatialCross[i*ncross+0]);

        // <delu(x*)del(YkH)> = (1/<rho(x*)>/<rho(x)>)*(<deljx(x*)del(rhoYkH) - <u(x*)><delrho(x*)del(rhoYkH)>
        //                      - <YkH(x)><deljx(x*)delrho(x)> + <u(x*)><YkH(x)><delrho(x*)delrho(x)>)
        spatialCross[i*ncross+25] = (Real(1.0)/(meanrho*meanrhocross))*(spatialCross[i*ncross+7] - meanuxcross*spatialCross[i*ncross+9]
                                                                     - meanYk[nspecies-1]*spatialCross[i*ncross+5] + meanuxcross*meanYk[nspecies-1]*spatialCross[i*ncross+0]);

        // Direct -- <delT(x*)delT(x)>
        runningMeanKahan(spatialCross[i*ncross+26], spatialCrossK[i*ncross+26], delTcross*delT, stepsinv);

        // Direct -- <delT(x*)delrho(x)>
        runningMean(spatialCross[i*ncross+27], delTcross*delrho, stepsinv);

        // Direct -- <delT(x*)delu(x)>
        runningMean(spatialCross[i*ncross+28], delTcross*delvx, stepsinv);

        // Direct -- <delu(x*)delrho>
        runningMean(spatialCross[i*ncross+29], delvxcross*delrho, stepsinv);

        // Direct -- <delu(x*) delYkL(x)>
        runningMean(spatialCross[i*ncross+30], delvxcross*delYk[0], stepsinv);

        // Direct -- <delu(x*) delYkH(x)>
        runningMean(spatialCross[i*ncross+31], delvxcross*delYk[nspecies-1], stepsinv);

        // Direct -- <delu(x*) delrhoYkL(x)>
        runningMean(spatialCross[i*ncross+32], delvxcross*delrhoYk[0], stepsinv);

        // Direct -- <delu(x*) delrhoYkH(x)>
        runningMean(spatialCross[i*ncross+33], delvxcross*delrhoYk[nspecies-1], stepsinv);

        // Direct <delYkL(x*)delYkL(x)>
        runningMeanKahan(spatialCross[i*ncross+34], spatialCrossK[i*ncross+34], delYkcross[0]*delYk[0], stepsinv);

        // Direct <delYkH(x*)delYkH(x)>
        runningMeanKahan(spatialCross[i*ncross+35], spatialCrossK[i*ncross+35], delYkcross[nspecies-1]*delYk[nspecies-1], stepsinv);

        // Direct <delYkL(x*)delYkH(x)>
        runningMeanKahan(spatialCross[i*ncross+36], spatialCrossK[i*ncross+36], delYkcross[0]*delYk[nspecies-1], stepsinv);

        // Last we rhoYk for species
        for (int ns=0; ns<nspecies; ++ns) {
            runningMeanKahan(spatialCross[i*ncross+37+ns], spatialCrossK[i*ncross+37+ns], delrhoYkcross[ns]*delrhoYk[ns], stepsinv); // <delrhoYk(x*)delrhoYk(x)>
        }

        // <delrhoYkL(x)delrhoYkH(x*)
        runningMeanKahan(spatialCross[i*ncross+37+nspecies+0], spatialCrossK[i*ncross+37+nspecies+0], delrhoYkcross[nspecies-1]*delrhoYk[0], stepsinv);

        // <delYkL(x*)delYkL(x)> = (1/<rho(x*)>/<rho(x)>)*(<delrhoYkL(x*)delrhoYkL> - <YkL(x*)><delrho(x*)delrhoYkL(x)>
        //                                                 - <YkL(x)><delrhoYkL(x*)delrho(x) + <YkL(x*)><YkL(x)><delrho(x*)delrho(x)>)
        Real delrhoYkdelrhoYk = spatialCross[i*ncross+37] + (delrhoYkcross[0]*delrhoYk[0] - spatialCross[i*ncross+37])*stepsinv;
        spatialCross[i*ncross+37+nspecies+1] = (Real(1.0)/(meanrho*meanrhocross))*(delrhoYkdelrhoYk - meanYkcross[0]*spatialCross[i*ncross+8]
                                                                    - meanYk[0]*spatialCross[i*ncross+10] + meanYkcross[0]*meanYk[0]*spatialCross[i*ncross+0]);

        // <delYkH(x*)delYkH(x)> = (1/<rho(x*)>/<rho(x)>)*(<delrhoYkH(x*)delrhoYkH> - <YkH(x*)><delrho(x*)delrhoYkH(x)>
        //                                                 - <YkH(x)><delrhoYkH(x*)delrho(x) + <YkH(x*)><YkH(x)><delrho(x*)delrho(x)>)
        delrhoYkdelrhoYk = spatialCross[i*ncross+37+nspecies-1] + (delrhoYkcross[nspecies-1]*delrhoYk[nspecies-1] - spatialCross[i*ncross+37+nspecies-1])*stepsinv;
        spatialCross[i*ncross+37+nspecies+2] = (Real(1.0)/(meanrho*meanrhocross))*(delrhoYkdelrhoYk - meanYkcross[nspecies-1]*spatialCross[i*ncross+9]
                                                    - meanYk[nspecies-1]*spatialCross[i*ncross+11] + meanYkcross[nspecies-1]*meanYk[nspecies-1]*spatialCross[i*ncross+0]);
    }
}

//    OLDER VERSION
//    BL_PROFILE_VAR("EvaluateSpatialCorrelations3D()",EvaluateSpatialCorrelations3D);
//
//    Real stepsminusone = steps - 1.;
//    Real stepsinv = 1./steps;
//
//    // Get mean values
//    Real meanrhocross = data_xcross[1];
//    Vector<Real>  meanYkcross(nspecies, 0.0);
//    for (int ns=0; ns<nspecies; ++ns) {
//        meanYkcross[ns] =  data_xcross[18+4*ns+3];
//    }
//    Real meanuxcross = data_xcross[11];
//
//    // Get fluctuations of the conserved variables at the cross cell
//    Real delrhocross = data_xcross[0] - data_xcross[1];
//    Real delKcross   = data_xcross[2] - data_xcross[3];
//    Real deljxcross  = data_xcross[4] - data_xcross[5];
//    Real deljycross  = data_xcross[6] - data_xcross[7];
//    Real deljzcross  = data_xcross[8] - data_xcross[9];
//    Vector<Real>  delrhoYkcross(nspecies, 0.0);
//    for (int ns=0; ns<nspecies; ++ns) {
//        delrhoYkcross[ns] =  data_xcross[18+4*ns+0] - data_xcross[18+4*ns+1];
//    }
//
//    // Get fluctuations of some primitive variables (for direct fluctuation calculations)
//    Real delTcross = data_xcross[16] - data_xcross[17];
//    Real delvxcross = data_xcross[10] - data_xcross[11];
//    Vector<Real>  delYkcross(nspecies, 0.0);
//    for (int ns=0; ns<nspecies; ++ns) {
//        delYkcross[ns] =  data_xcross[18+4*ns+2] - data_xcross[18+4*ns+3];
//    }
//
//    // evaluate heat stuff at the cross cell
//    Real cvcross = 0.;
//    for (int l=0; l<nspecies; ++l) {
//        cvcross = cvcross + hcv[l]*data_xcross[18+4*l+1]/data_xcross[1];
//    }
//    Real cvinvcross = 1.0/cvcross;
//    Real qmeancross = cvcross*data_xcross[17] -
//                     0.5*(data_xcross[11]*data_xcross[11] + data_xcross[13]*data_xcross[13] + data_xcross[15]*data_xcross[15]);
//
//    // Get fluctuations of derived hydrodynamic quantities at the cross cell
//    // delG = \vec{v}\cdot\vec{\deltaj}
//    Real delGcross = data_xcross[11]*(data_xcross[4]-data_xcross[5]) + data_xcross[13]*(data_xcross[6]-data_xcross[7]) +
//                    data_xcross[15]*(data_xcross[8]-data_xcross[9]);
//
//    /////////////////////////////////////////////////////////////
//    // evaluate x-spatial correlations
//    /////////////////////////////////////////////////////////////
//    // int ncross = 37+nspecies+2; check main_drive.cpp for latest
//    for (int i=0; i<n_cells[0]; ++i) {
//
//        // Get mean values
//        Real meanrho = data_x[i*nstats+1];
//        Vector<Real>  meanYk(nspecies, 0.0);
//        for (int ns=0; ns<nspecies; ++ns) {
//            meanYk[ns] = data_x[i*nstats+18+4*ns+3];
//        }
//
//        // Get fluctuations of the conserved variables
//        Real delrho = data_x[i*nstats+0] - data_x[i*nstats+1];
//        Real delK   = data_x[i*nstats+2] - data_x[i*nstats+3];
//        Real deljx  = data_x[i*nstats+4] - data_x[i*nstats+5];
//        Real deljy  = data_x[i*nstats+6] - data_x[i*nstats+7];
//        Real deljz  = data_x[i*nstats+8] - data_x[i*nstats+9];
//        Vector<Real>  delrhoYk(nspecies, 0.0);
//        for (int ns=0; ns<nspecies; ++ns) {
//            delrhoYk[ns] = data_x[i*nstats+18+4*ns+0] - data_x[i*nstats+18+4*ns+1];
//        }
//
//        // Get fluctuations of some primitive variables (for direct fluctuation calculations)
//        Real delT = data_x[i*nstats+16] - data_x[i*nstats+17];
//        Real delvx = data_x[i*nstats+10] - data_x[i*nstats+11];
//        Vector<Real>  delYk(nspecies, 0.0);
//        for (int ns=0; ns<nspecies; ++ns) {
//            delYk[ns] = data_x[i*nstats+18+4*ns+2] - data_x[i*nstats+18+4*ns+3];
//        }
//
//        // evaluate heat stuff at the cross cell
//        Real cv = 0.;
//        for (int l=0; l<nspecies; ++l) {
//            cv = cv + hcv[l]*data_xcross[18+4*l+1]/data_xcross[1];
//        }
//        Real cvinv = 1.0/cv;
//        Real qmean = cv*data_x[i*nstats+17] -
//                         0.5*(data_x[i*nstats+11]*data_x[i*nstats+11] + data_x[i*nstats+13]*data_x[i*nstats+13] + data_x[i*nstats+15]*data_x[i*nstats+15]);
//
//        // Get fluctuations of derived hydrodynamic quantities
//        // delG = \vec{v}\cdot\vec{\deltaj}
//        Real delG = data_x[i*nstats+11]*(data_x[i*nstats+4] - data_x[i*nstats+5]) + data_x[i*nstats+13]*(data_x[i*nstats+6] - data_x[i*nstats+7]) +
//                    data_x[i*nstats+15]*(data_x[i*nstats+8] - data_x[i*nstats+9]);
//
//        // First update correlations of conserved quantities (we will do rhoYk later)
//        spatialCross[i*ncross+0]  = (spatialCross[i*ncross+0]*stepsminusone + delrhocross*delrho)*stepsinv; // <delrho(x*)delrho(x)>
//        spatialCross[i*ncross+1]  = (spatialCross[i*ncross+1]*stepsminusone + delKcross*delK)*stepsinv;     // <delK(x*)delK(x)>
//        spatialCross[i*ncross+2]  = (spatialCross[i*ncross+2]*stepsminusone + deljxcross*deljx)*stepsinv;   // <deljx(x*)deljx(x)>
//        spatialCross[i*ncross+3]  = (spatialCross[i*ncross+3]*stepsminusone + deljycross*deljy)*stepsinv;   // <deljy(x*)deljy(x)>
//        spatialCross[i*ncross+4]  = (spatialCross[i*ncross+4]*stepsminusone + deljzcross*deljz)*stepsinv;   // <deljz(x*)deljz(x)>
//        spatialCross[i*ncross+5]  = (spatialCross[i*ncross+5]*stepsminusone + deljxcross*delrho)*stepsinv;  // <deljx(x*)delrho(x)>
//        spatialCross[i*ncross+6]  = (spatialCross[i*ncross+6]*stepsminusone + deljxcross*delrhoYk[0])*stepsinv;  // <deljx(x*)delrhoYkL(x)>
//        spatialCross[i*ncross+7]  = (spatialCross[i*ncross+7]*stepsminusone + deljxcross*delrhoYk[nspecies-1])*stepsinv;  // <deljx(x*)delrhoYkH(x)>
//        spatialCross[i*ncross+8]  = (spatialCross[i*ncross+8]*stepsminusone + delrhocross*delrhoYk[0])*stepsinv; // <delrho(x*)delrhoYkL(x)>
//        spatialCross[i*ncross+9]  = (spatialCross[i*ncross+9]*stepsminusone + delrhocross*delrhoYk[nspecies-1])*stepsinv; // <delrho(x*)delrhoYkH(x)>
//        spatialCross[i*ncross+10] = (spatialCross[i*ncross+10]*stepsminusone + delrhoYkcross[0]*delrho)*stepsinv; // <delrhoYkL(x*)delrho(x)>
//        spatialCross[i*ncross+11] = (spatialCross[i*ncross+11]*stepsminusone + delrhoYkcross[nspecies-1]*delrho)*stepsinv; // <delrhoYkH(x*)delrho(x)>
//
//        // Some more cross-correlations for hydrodynamical variables later
//        spatialCross[i*ncross+12] = (spatialCross[i*ncross+12]*stepsminusone + delGcross*delG)*stepsinv; // <delG(x*)delG(x)>
//        spatialCross[i*ncross+13] = (spatialCross[i*ncross+13]*stepsminusone + delGcross*delK)*stepsinv; // <delG(x*)delK(x)>
//        spatialCross[i*ncross+14] = (spatialCross[i*ncross+14]*stepsminusone + delKcross*delG)*stepsinv; // <delK(x*)delG(x)>
//        spatialCross[i*ncross+15] = (spatialCross[i*ncross+15]*stepsminusone + delrhocross*delK)*stepsinv; // <delrho(x*)delK(x)>
//        spatialCross[i*ncross+16] = (spatialCross[i*ncross+16]*stepsminusone + delKcross*delrho)*stepsinv; // <delK(x*)delrho(x)>
//        spatialCross[i*ncross+17] = (spatialCross[i*ncross+17]*stepsminusone + delrhocross*delG)*stepsinv; // <delrho(x*)delG(x)>
//        spatialCross[i*ncross+18] = (spatialCross[i*ncross+18]*stepsminusone + delGcross*delrho)*stepsinv; // <delG(x*)delrho(x)>
//
//        // Next we do cross-correlations with and between hydrodynamical variables
//        // <delT(x*)delT(x)> = (1/cv*/cv/<rho(x)>/<rho(x*)>)(<delK*delK> + <delG*delG> - <delG*delK> - <delK*delG>
//        //                      + <Q><Q*><delrho*delrho> - <Q*><delrho*delK> - <Q><delK*delrho> + <Q*><delrho*delG> + <Q><delG*delrho>)
//        spatialCross[i*ncross+19] = (cvinvcross*cvinv/(meanrhocross*meanrho))*
//                                    (spatialCross[i*ncross+1] + spatialCross[i*ncross+12] - spatialCross[i*ncross+13] - spatialCross[i*ncross+14]
//                                     + qmean*qmeancross*spatialCross[i*ncross+0] - qmeancross*spatialCross[i*ncross+15] - qmean*spatialCross[i*ncross+16]
//                                     + qmeancross*spatialCross[i*ncross+17] + qmean*spatialCross[i*ncross+18]);
//
//        // <delT(x*)delrho(x)> = (1/cv/<rho(x*)>)*(<delK*delrho> - <delG*delrho> - <Q*><delrhodelrho*>)
//        spatialCross[i*ncross+20] = (cvinvcross*meanrhocross)*(spatialCross[i*ncross+16] - spatialCross[i*ncross+18] - qmeancross*spatialCross[i*ncross+0]);
//
//        // <delu(x*)delrho> = (1/<rho(x*)>)*(<deljx(x*)delrho(x)> - <u(x*)><<delrho(x*)delrho(x)>)
//        spatialCross[i*ncross+21] = (1.0/meanrhocross)*(spatialCross[i*ncross+5] - meanuxcross*spatialCross[i*ncross+0]);
//
//        // <delu(x*)del(rhoYkL)> = (1/<rho(x*)>)*(<deljx(x*)del(rhoYkL)> - <u(x*)><delrho(x*)del(rhoYkL)>)
//        spatialCross[i*ncross+22] = (1.0/meanrhocross)*(spatialCross[i*ncross+6] - meanuxcross*spatialCross[i*ncross+8]);
//
//        // <delu(x*)del(rhoYkH)> = (1/<rho(x*)>)*(<deljx(x*)del(rhoYkH)> - <u(x*)><delrho(x*)del(rhoYkH)>)
//        spatialCross[i*ncross+23] = (1.0/meanrhocross)*(spatialCross[i*ncross+7] - meanuxcross*spatialCross[i*ncross+9]);
//
//        // <delu(x*)del(YkL)> = (1/<rho(x*)>/<rho(x)>)*(<deljx(x*)del(rhoYkL) - <u(x*)><delrho(x*)del(rhoYkL)>
//        //                      - <YkL(x)><deljx(x*)delrho(x)> + <u(x*)><YkL(x)><delrho(x*)delrho(x)>)
//        spatialCross[i*ncross+24] = (1.0/(meanrho*meanrhocross))*(spatialCross[i*ncross+6] - meanuxcross*spatialCross[i*ncross+8]
//                                                                 - meanYk[0]*spatialCross[i*ncross+5] + meanuxcross*meanYk[0]*spatialCross[i*ncross+0]);
//
//        // <delu(x*)del(YkH)> = (1/<rho(x*)>/<rho(x)>)*(<deljx(x*)del(rhoYkH) - <u(x*)><delrho(x*)del(rhoYkH)>
//        //                      - <YkH(x)><deljx(x*)delrho(x)> + <u(x*)><YkH(x)><delrho(x*)delrho(x)>)
//        spatialCross[i*ncross+25] = (1.0/(meanrho*meanrhocross))*(spatialCross[i*ncross+7] - meanuxcross*spatialCross[i*ncross+9]
//                                                                 - meanYk[nspecies-1]*spatialCross[i*ncross+5] + meanuxcross*meanYk[nspecies-1]*spatialCross[i*ncross+0]);
//
//        // Direct -- <delT(x*)delT(x)>
//        spatialCross[i*ncross+26] = (spatialCross[i*ncross+26]*stepsminusone + delTcross*delT)*stepsinv;
//
//        // Direct -- <delT(x*)delrho(x)>
//        spatialCross[i*ncross+27] = (spatialCross[i*ncross+27]*stepsminusone + delTcross*delrho)*stepsinv;
//
//        // Direct -- <delT(x*)delu(x)>
//        spatialCross[i*ncross+28] = (spatialCross[i*ncross+28]*stepsminusone + delTcross*delvx)*stepsinv;
//
//        // Direct -- <delu(x*)delrho>
//        spatialCross[i*ncross+29] = (spatialCross[i*ncross+29]*stepsminusone + delvxcross*delrho)*stepsinv;
//
//        // Direct -- <delu(x*)del(rhoYkL)
//        spatialCross[i*ncross+30] = (spatialCross[i*ncross+30]*stepsminusone + delvxcross*delYk[0])*stepsinv;
//
//        // Direct -- <delu(x*)del(rhoYkH)
//        spatialCross[i*ncross+31] = (spatialCross[i*ncross+31]*stepsminusone + delvxcross*delYk[nspecies-1])*stepsinv;
//
//        // Direct -- <delu(x*)del(YkL)
//        spatialCross[i*ncross+32] = (spatialCross[i*ncross+32]*stepsminusone + delvxcross*delrhoYk[0])*stepsinv;
//
//        // Direct -- <delu(x*)del(YkH)
//        spatialCross[i*ncross+33] = (spatialCross[i*ncross+33]*stepsminusone + delvxcross*delrhoYk[nspecies-1])*stepsinv;
//
//        // Direct <delYkL(x*)delYkL(x)>
//        spatialCross[i*ncross+34] = (spatialCross[i*ncross+34]*stepsminusone + delYkcross[0]*delYk[0])*stepsinv;
//
//        // Direct <delYkH(x*)delYkH(x)>
//        spatialCross[i*ncross+35] = (spatialCross[i*ncross+35]*stepsminusone + delYkcross[nspecies-1]*delYk[nspecies-1])*stepsinv;
//
//        // Direct <delYkL(x*)delYkH(x)>
//        spatialCross[i*ncross+36] = (spatialCross[i*ncross+36]*stepsminusone + delYkcross[0]*delYk[nspecies-1])*stepsinv;
//
//        // Last we rhoYk for species
//        for (int ns=0; ns<nspecies; ++ns) {
//            spatialCross[i*ncross+37+ns] = (spatialCross[i*ncross+37+ns]*stepsminusone + delrhoYkcross[ns]*delrhoYk[ns])*stepsinv; // <delrhoYk(x*)delrhoYk(x)>
//        }
//
//        // <delYkL(x*)delYkL(x)> = (1/<rho(x*)>/<rho(x)>)*(<delrhoYkL(x*)delrhoYkL> - <YkL(x*)><delrho(x*)delrhoYkL(x)>
//        //                                                 - <YkL(x)><delrhoYkL(x*)delrho(x) + <YkL(x*)><YkL(x)><delrho(x*)delrho(x)>)
//        Real delrhoYkdelrhoYk = (spatialCross[i*ncross+37]*stepsminusone + delrhoYkcross[0]*delrhoYk[0])*stepsinv;
//        spatialCross[i*ncross+37+nspecies] = (1.0/(meanrho*meanrhocross))*(delrhoYkdelrhoYk - meanYkcross[0]*spatialCross[i*ncross+8]
//                                                                - meanYk[0]*spatialCross[i*ncross+10] + meanYkcross[0]*meanYk[0]*spatialCross[i*ncross+0]);
//
//        // <delYkL(x*)delYkH(x)> = (1/<rho(x*)>/<rho(x)>)*(<delrhoYkH(x*)delrhoYkH> - <YkH(x*)><delrho(x*)delrhoYkH(x)>
//        //                                                 - <YkH(x)><delrhoYkH(x*)delrho(x) + <YkH(x*)><YkH(x)><delrho(x*)delrho(x)>)
//        delrhoYkdelrhoYk = (spatialCross[i*ncross+37+nspecies-1]*stepsminusone + delrhoYkcross[nspecies-1]*delrhoYk[nspecies-1])*stepsinv;
//        spatialCross[i*ncross+37+nspecies+1] = (1.0/(meanrho*meanrhocross))*(delrhoYkdelrhoYk - meanYkcross[nspecies-1]*spatialCross[i*ncross+9]
//                                                - meanYk[nspecies-1]*spatialCross[i*ncross+11] + meanYkcross[nspecies-1]*meanYk[nspecies-1]*spatialCross[i*ncross+0]);
//    }
//}


//////////////////////////////////////////
// Update Spatial Correlations ///////////
// ///////////////////////////////////////
void EvaluateSpatialCorrelations1D(MultiFab& spatialCross1D,
                                   MultiFab& spatialCross1DK,
                                   amrex::Gpu::DeviceVector<Real>& data_xcross_in,
                                   const MultiFab& consMean,
                                   const MultiFab& primMean,
                                   const MultiFab& prim_in,
                                   const MultiFab& cons,
                                   const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                                   const std::array<MultiFab, AMREX_SPACEDIM>& velMean,
                                   const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                                   const std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                                   const int steps,
                                   const int nstats,
                                   const int ncross,
                                   const int star_index)
{

    BL_PROFILE_VAR("EvaluateSpatialCorrelations1D()",EvaluateSpatialCorrelations1D);

    Real stepsinv = Real(1.0/static_cast<double>(steps));

    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const Array4<const Real> cumeans   = consMean.array(mfi);
        const Array4<const Real> primmeans = primMean.array(mfi);
        const Array4<const Real> prim      = prim_in.array(mfi);
        const Array4<const Real> cu        = cons.array(mfi);

        const Array4<const Real> velx      = vel[0].array(mfi);
        const Array4<const Real> velxmeans = velMean[0].array(mfi);
        const Array4<const Real> velymeans = velMean[1].array(mfi);
        const Array4<const Real> velzmeans = velMean[2].array(mfi);

        const Array4<const Real> momx      = cumom[0].array(mfi);
        const Array4<const Real> momy      = cumom[1].array(mfi);
        const Array4<const Real> momz      = cumom[2].array(mfi);
        const Array4<const Real> momxmeans = cumomMean[0].array(mfi);
        const Array4<const Real> momymeans = cumomMean[1].array(mfi);
        const Array4<const Real> momzmeans = cumomMean[2].array(mfi);

        const Array4<      Real> spatialCross = spatialCross1D.array(mfi);
        const Array4<      Real> spatialCrossK = spatialCross1DK.array(mfi);

        Real* data_xcross = data_xcross_in.data();

        //for (auto k = lo.z; k <= hi.z; ++k) {
        //for (auto j = lo.y; j <= hi.y; ++j) {
        //for (auto i = lo.x; i <= hi.x; ++i) {

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            //////////////////////////////////////////////
            // Get fluctuations at xcross for this j and k
            //////////////////////////////////////////////
            int index = k*n_cells[1]*nstats + j*nstats;

            Real meanrhocross = data_xcross[index + 1];
            Real meanuxcross  = data_xcross[index + 11];

            GpuArray<Real,MAX_SPECIES> meanYkcross;
            GpuArray<Real,MAX_SPECIES> delYkcross;
            GpuArray<Real,MAX_SPECIES> delrhoYkcross;
            for (int ns=0; ns<nspecies; ++ns) {
                meanYkcross[ns] =  data_xcross[index + 18+4*ns+3];
                delYkcross[ns]  =  data_xcross[index + 18+4*ns+2] - data_xcross[index + 18+4*ns+3];
                delrhoYkcross[ns] =  data_xcross[index + 18+4*ns+0] - data_xcross[index + 18+4*ns+1];
            }

            Real delrhocross = data_xcross[index + 0] - data_xcross[index + 1];
            Real delKcross   = data_xcross[index + 2] - data_xcross[index + 3];
            Real deljxcross  = data_xcross[index + 4] - data_xcross[index + 5];
            Real deljycross  = data_xcross[index + 6] - data_xcross[index + 7];
            Real deljzcross  = data_xcross[index + 8] - data_xcross[index + 9];

            Real delTcross = data_xcross[index + 16] - data_xcross[index + 17];
            Real delvxcross = data_xcross[index + 10] - data_xcross[index + 11];

            Real cvcross = 0.;
            for (int l=0; l<nspecies; ++l) {
                cvcross = cvcross + hcv[l]*data_xcross[index + 18+4*l+1]/data_xcross[index + 1];
            }
            Real cvinvcross = Real(1.0)/cvcross;
            Real vxmeancross = data_xcross[index + 11];
            Real vymeancross = data_xcross[index + 13];
            Real vzmeancross = data_xcross[index + 15];
            Real qmeancross = cvcross*data_xcross[index + 17] -
                             Real(0.5)*(vxmeancross*vxmeancross + vymeancross*vymeancross + vzmeancross*vzmeancross);

            // delG = \vec{v}\cdot\vec{\deltaj}
            Real delGcross = vxmeancross*deljxcross + vymeancross*deljycross + vzmeancross*deljzcross;


            ////////////////////////////////////////
            // Get fluctuations at this cell (i,j,k)
            ////////////////////////////////////////

            Real meanrho = cumeans(i,j,k,0);

            GpuArray<Real,MAX_SPECIES>  meanYk;
            GpuArray<Real,MAX_SPECIES> delrhoYk;
            GpuArray<Real,MAX_SPECIES> delYk;
            for (int ns=0; ns<nspecies; ++ns) {
                meanYk[ns] = primmeans(i,j,k,6+ns);;
                delrhoYk[ns] = cu(i,j,k,5+ns) - cumeans(i,j,k,5+ns);
                delYk[ns]    = prim(i,j,k,6+ns) - primmeans(i,j,k,6+ns);
            }

            Real delrho = cu(i,j,k,0) - cumeans(i,j,k,0);
            Real delK   = cu(i,j,k,4) - cumeans(i,j,k,4);
            Real deljx  = Real(0.5)*(momx(i,j,k) + momx(i+1,j,k)) - Real(0.5)*(momxmeans(i,j,k) + momxmeans(i+1,j,k));
            Real deljy  = Real(0.5)*(momy(i,j,k) + momy(i,j+1,k)) - Real(0.5)*(momymeans(i,j,k) + momymeans(i,j+1,k));
            Real deljz  = Real(0.5)*(momz(i,j,k) + momz(i,j,k+1)) - Real(0.5)*(momzmeans(i,j,k) + momzmeans(i,j,k+1));

            Real delT = prim(i,j,k,4) - primmeans(i,j,k,4);
            Real delvx = Real(0.5)*(velx(i,j,k) + velx(i+1,j,k)) - Real(0.5)*(velxmeans(i,j,k) + velxmeans(i+1,j,k));

            Real cv = 0.;
            for (int l=0; l<nspecies; ++l) {
                cv = cv + hcv[l]*cumeans(i,j,k,5+l)/cumeans(i,j,k,0);
            }
            Real cvinv = Real(1.0)/cv;
            Real vxmean = Real(0.5)*(velxmeans(i,j,k) + velxmeans(i+1,j,k));
            Real vymean = Real(0.5)*(velymeans(i,j,k) + velymeans(i,j+1,k));
            Real vzmean = Real(0.5)*(velzmeans(i,j,k) + velzmeans(i,j,k+1));

            Real qmean = cv*primmeans(i,j,k,4) - Real(0.5)*(vxmean*vxmean + vymean*vymean + vzmean*vzmean);

            // delG = \vec{v}\cdot\vec{\deltaj}
            Real delG = vxmean*deljx + vymean*deljy + vzmean*deljz;

            // Spatial Correlation Calculations

            int MFindex = star_index*ncross;

            runningMeanKahan(spatialCross(i,j,k,MFindex+0), spatialCrossK(i,j,k,MFindex+0), delrhocross*delrho, stepsinv); // <delrho(x*)delrho(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+1), spatialCrossK(i,j,k,MFindex+1), delKcross*delK, stepsinv);     // <delK(x*)delK(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+2), spatialCrossK(i,j,k,MFindex+2), deljxcross*deljx, stepsinv);   // <deljx(x*)deljx(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+3), spatialCrossK(i,j,k,MFindex+3), deljycross*deljy, stepsinv);   // <deljy(x*)deljy(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+4), spatialCrossK(i,j,k,MFindex+4), deljzcross*deljz, stepsinv);   // <deljz(x*)deljz(x)>
            runningMean(spatialCross(i,j,k,MFindex+5), deljxcross*delrho, stepsinv);  // <deljx(x*)delrho(x)>
            runningMean(spatialCross(i,j,k,MFindex+6), deljxcross*delrhoYk[0], stepsinv);  // <deljx(x*)delrhoYkL(x)>
            runningMean(spatialCross(i,j,k,MFindex+7), deljxcross*delrhoYk[nspecies-1], stepsinv);  // <deljx(x*)delrhoYkH(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+8), spatialCrossK(i,j,k,MFindex+8), delrhocross*delrhoYk[0], stepsinv); // <delrho(x*)delrhoYkL(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+9), spatialCrossK(i,j,k,MFindex+9), delrhocross*delrhoYk[nspecies-1], stepsinv); // <delrho(x*)delrhoYkH(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+10), spatialCrossK(i,j,k,MFindex+10), delrhoYkcross[0]*delrho, stepsinv); // <delrhoYkL(x*)delrho(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+11), spatialCrossK(i,j,k,MFindex+11), delrhoYkcross[nspecies-1]*delrho, stepsinv); // <delrhoYkH(x*)delrho(x)>

            runningMeanKahan(spatialCross(i,j,k,MFindex+12), spatialCrossK(i,j,k,MFindex+12), delGcross*delG, stepsinv); // <delG(x*)delG(x)>
            runningMean(spatialCross(i,j,k,MFindex+13), delGcross*delK, stepsinv); // <delG(x*)delK(x)>
            runningMean(spatialCross(i,j,k,MFindex+14), delKcross*delG, stepsinv); // <delK(x*)delG(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+15), spatialCrossK(i,j,k,MFindex+15), delrhocross*delK, stepsinv); // <delrho(x*)delK(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+16), spatialCrossK(i,j,k,MFindex+16), delKcross*delrho, stepsinv); // <delK(x*)delrho(x)>
            runningMean(spatialCross(i,j,k,MFindex+17), delrhocross*delG, stepsinv); // <delrho(x*)delG(x)>
            runningMean(spatialCross(i,j,k,MFindex+18), delGcross*delrho, stepsinv); // <delG(x*)delrho(x)>

            // <delT(x*)delT(x)> = (1/cv*/cv/<rho(x)>/<rho(x*)>)(<delK*delK> + <delG*delG> - <delG*delK> - <delK*delG>
            //                      + <Q><Q*><delrho*delrho> - <Q*><delrho*delK> - <Q><delK*delrho> + <Q*><delrho*delG> + <Q><delG*delrho>)
            spatialCross(i,j,k,MFindex+19) = (cvinvcross*cvinv/(meanrhocross*meanrho))*
                                        (spatialCross(i,j,k,MFindex+1) + spatialCross(i,j,k,MFindex+12) - spatialCross(i,j,k,MFindex+13) - spatialCross(i,j,k,MFindex+14)
                                         + qmean*qmeancross*spatialCross(i,j,k,MFindex+0) - qmeancross*spatialCross(i,j,k,MFindex+15) - qmean*spatialCross(i,j,k,MFindex+16)
                                         + qmeancross*spatialCross(i,j,k,MFindex+17) + qmean*spatialCross(i,j,k,MFindex+18));

            // <delT(x*)delrho(x)> = (1/cv/<rho(x*)>)*(<delK*delrho> - <delG*delrho> - <Q*><delrhodelrho*>)
            spatialCross(i,j,k,MFindex+20) = (cvinvcross/meanrhocross)*(spatialCross(i,j,k,MFindex+16) - spatialCross(i,j,k,MFindex+18) - qmeancross*spatialCross(i,j,k,MFindex+0));

            // <delu(x*)delrho> = (1/<rho(x*)>)*(<deljx(x*)delrho(x)> - <u(x*)><<delrho(x*)delrho(x)>)
            spatialCross(i,j,k,MFindex+21) = (Real(1.0)/meanrhocross)*(spatialCross(i,j,k,MFindex+5) - meanuxcross*spatialCross(i,j,k,MFindex+0));

            // <delu(x*)del(rhoYkL)> = (1/<rho(x*)>)*(<deljx(x*)del(rhoYkL)> - <u(x*)><delrho(x*)del(rhoYkL)>)
            spatialCross(i,j,k,MFindex+22) = (Real(1.0)/meanrhocross)*(spatialCross(i,j,k,MFindex+6) - meanuxcross*spatialCross(i,j,k,MFindex+8));

            // <delu(x*)del(rhoYkH)> = (1/<rho(x*)>)*(<deljx(x*)del(rhoYkH)> - <u(x*)><delrho(x*)del(rhoYkH)>)
            spatialCross(i,j,k,MFindex+23) = (Real(1.0)/meanrhocross)*(spatialCross(i,j,k,MFindex+7) - meanuxcross*spatialCross(i,j,k,MFindex+9));

            // <delu(x*)del(YkL)> = (1/<rho(x*)>/<rho(x)>)*(<deljx(x*)del(rhoYkL) - <u(x*)><delrho(x*)del(rhoYkL)>
            //                      - <YkL(x)><deljx(x*)delrho(x)> + <u(x*)><YkL(x)><delrho(x*)delrho(x)>)
            spatialCross(i,j,k,MFindex+24) = (Real(1.0)/(meanrho*meanrhocross))*(spatialCross(i,j,k,MFindex+6) - meanuxcross*spatialCross(i,j,k,MFindex+8)
                                                                     - meanYk[0]*spatialCross(i,j,k,MFindex+5) + meanuxcross*meanYk[0]*spatialCross(i,j,k,MFindex+0));

            // <delu(x*)del(YkH)> = (1/<rho(x*)>/<rho(x)>)*(<deljx(x*)del(rhoYkH) - <u(x*)><delrho(x*)del(rhoYkH)>
            //                      - <YkH(x)><deljx(x*)delrho(x)> + <u(x*)><YkH(x)><delrho(x*)delrho(x)>)
            spatialCross(i,j,k,MFindex+25) = (Real(1.0)/(meanrho*meanrhocross))*(spatialCross(i,j,k,MFindex+7) - meanuxcross*spatialCross(i,j,k,MFindex+9)
                                                                     - meanYk[nspecies-1]*spatialCross(i,j,k,MFindex+5) +
                                                                     meanuxcross*meanYk[nspecies-1]*spatialCross(i,j,k,MFindex+0));

            // Direct -- <delT(x*)delT(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+26), spatialCrossK(i,j,k,MFindex+26), delTcross*delT, stepsinv);

            // Direct -- <delT(x*)delrho(x)>
            runningMean(spatialCross(i,j,k,MFindex+27), delTcross*delrho, stepsinv);

            // Direct -- <delT(x*)delu(x)>
            runningMean(spatialCross(i,j,k,MFindex+28), delTcross*delvx, stepsinv);

            // Direct -- <delu(x*)delrho>
            runningMean(spatialCross(i,j,k,MFindex+29), delvxcross*delrho, stepsinv);

            // Direct -- <delu(x*) delYkL(x)>
            runningMean(spatialCross(i,j,k,MFindex+30), delvxcross*delYk[0], stepsinv);

            // Direct -- <delu(x*) delYkH(x)>
            runningMean(spatialCross(i,j,k,MFindex+31), delvxcross*delYk[nspecies-1], stepsinv);

            // Direct -- <delu(x*) delrhoYkL(x)>
            runningMean(spatialCross(i,j,k,MFindex+32), delvxcross*delrhoYk[0], stepsinv);

            // Direct -- <delu(x*) delrhoYkH(x)>
            runningMean(spatialCross(i,j,k,MFindex+33), delvxcross*delrhoYk[nspecies-1], stepsinv);

            // Direct <delYkL(x*)delYkL(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+34), spatialCrossK(i,j,k,MFindex+34), delYkcross[0]*delYk[0], stepsinv);

            // Direct <delYkH(x*)delYkH(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+35), spatialCrossK(i,j,k,MFindex+35), delYkcross[nspecies-1]*delYk[nspecies-1], stepsinv);

            // Direct <delYkL(x*)delYkH(x)>
            runningMeanKahan(spatialCross(i,j,k,MFindex+36), spatialCrossK(i,j,k,MFindex+36), delYkcross[0]*delYk[nspecies-1], stepsinv);

            // Last we rhoYk for species
            for (int ns=0; ns<nspecies; ++ns) {
                runningMeanKahan(spatialCross(i,j,k,MFindex+37+ns), spatialCrossK(i,j,k,MFindex+37+ns), delrhoYkcross[ns]*delrhoYk[ns], stepsinv); // <delrhoYk(x*)delrhoYk(x)>
            }

            // <delrhoYkL(x)delrhoYkH(x*)
            runningMeanKahan(spatialCross(i,j,k,MFindex+37+nspecies+0), spatialCrossK(i,j,k,MFindex+37+nspecies+0), delrhoYkcross[nspecies-1]*delrhoYk[0], stepsinv);

            // <delYkL(x*)delYkL(x)> = (1/<rho(x*)>/<rho(x)>)*(<delrhoYkL(x*)delrhoYkL> - <YkL(x*)><delrho(x*)delrhoYkL(x)>
            //                                                 - <YkL(x)><delrhoYkL(x*)delrho(x) + <YkL(x*)><YkL(x)><delrho(x*)delrho(x)>)
            Real delrhoYkdelrhoYk = spatialCross(i,j,k,MFindex+37) + (delrhoYkcross[0]*delrhoYk[0] - spatialCross(i,j,k,MFindex+37))*stepsinv;
            spatialCross(i,j,k,MFindex+37+nspecies+1) = (Real(1.0)/(meanrho*meanrhocross))*(delrhoYkdelrhoYk - meanYkcross[0]*spatialCross(i,j,k,MFindex+8)
                                                                    - meanYk[0]*spatialCross(i,j,k,MFindex+10) + meanYkcross[0]*meanYk[0]*spatialCross(i,j,k,MFindex+0));

            // <delYkH(x*)delYkH(x)> = (1/<rho(x*)>/<rho(x)>)*(<delrhoYkH(x*)delrhoYkH> - <YkH(x*)><delrho(x*)delrhoYkH(x)>
            //                                                 - <YkH(x)><delrhoYkH(x*)delrho(x) + <YkH(x*)><YkH(x)><delrho(x*)delrho(x)>)
            delrhoYkdelrhoYk = spatialCross(i,j,k,MFindex+37+nspecies-1) + (delrhoYkcross[nspecies-1]*delrhoYk[nspecies-1] - spatialCross(i,j,k,MFindex+37+nspecies-1))*stepsinv;
            spatialCross(i,j,k,MFindex+37+nspecies+2) = (Real(1.0)/(meanrho*meanrhocross))*(delrhoYkdelrhoYk - meanYkcross[nspecies-1]*spatialCross(i,j,k,MFindex+9)
                                                    - meanYk[nspecies-1]*spatialCross(i,j,k,MFindex+11) + meanYkcross[nspecies-1]*meanYk[nspecies-1]*spatialCross(i,j,k,MFindex+0));

        });
       // }
       // }
       // }
    }
}

