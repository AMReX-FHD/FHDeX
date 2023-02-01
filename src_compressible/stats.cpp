#include "compressible_functions.H"

#include "common_functions.H"


///////////////////////////////////////////
// Evaluate Stats for the 3D case /////////
/// ///////////////////////////////////////
void evaluateStats3D(MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                     MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
                     MultiFab& coVar, 
                     Vector<Real>& spatialCross3D, const int ncross,
                     const amrex::Box& domain,
                     const int steps,
                     const Geometry& geom)
{
    BL_PROFILE_VAR("evaluateStats3D()",evaluateStats3D);
    
    //// Evaluate Means
    if (plot_means) {
        EvaluateStatsMeans(cons,consMean,prim_in,primMean,steps);
        consMean.FillBoundary(geom.periodicity());
        primMean.FillBoundary(geom.periodicity());
    }

    //// Evaluate Variances and Covariances
    if ((plot_vars) or (plot_covars)) {
        EvaluateVarsCoVars(cons,consMean,consVar,prim_in,primMean,primVar,coVar,steps);
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

    for (int i=0; i<nvars*domain.length(0); ++i) {
        cu_avg[i] /= n_cells[1]*n_cells[2];
        cumeans_avg[i] /= n_cells[1]*n_cells[2];
    }
    for (int i=0; i<(nspecies+4)*domain.length(0); ++i) {
        prim_avg[i] /= n_cells[1]*n_cells[2];
        primmeans_avg[i] /= n_cells[1]*n_cells[2];
    }

    // Update Spatial Correlations
    if (plot_cross) EvaluateSpatialCorrelations3D(spatialCross3D,cu_avg,cumeans_avg,prim_avg,primmeans_avg,steps,nstats,ncross);
}


////////////////////////////////////////////////
////// Evaluate Stats for the 2D case (do later)
/////// ////////////////////////////////////////
//void evaluateStatsStag2D(MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
//                         MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
//                         const std::array<MultiFab, AMREX_SPACEDIM>& vel, 
//                         std::array<MultiFab, AMREX_SPACEDIM>& velMean, 
//                         std::array<MultiFab, AMREX_SPACEDIM>& velVar, 
//                         const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
//                         std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
//                         std::array<MultiFab, AMREX_SPACEDIM>& cumomVar,
//                         MultiFab& coVar, 
//                         MultiFab& theta, MultiFab& thetaMean, MultiFab& thetaVar,
//                         MultiFab& /*spatialCross2D*/, const int /*ncross*/,
//                         const int steps,
//                         const Geometry& geom)
//{
//    BL_PROFILE_VAR("evaluateStatsStag2D()",evaluateStatsStag2D);
//    
//    //// Evaluate Means
//    if (plot_means) {
//        EvaluateStatsMeans(cons,consMean,prim_in,primMean,vel,velMean,cumom,cumomMean,theta,thetaMean,steps);
//        consMean.FillBoundary(geom.periodicity());
//        primMean.FillBoundary(geom.periodicity());
//    }
//
//    //// Evaluate Variances and Covariances
//    if ((plot_vars) or (plot_covars)) {
//        EvaluateVarsCoVars(cons,consMean,consVar,prim_in,primMean,primVar,velMean,velVar,
//                           cumom,cumomMean,cumomVar,coVar,theta,thetaMean,thetaVar,steps);
//        consVar.FillBoundary(geom.periodicity());
//        primVar.FillBoundary(geom.periodicity());
//    }
//
//    //// Evaluate Spatial Correlations (do later)
//    
//    //// contains yz-averaged running & instantaneous averages of conserved variables (2*nvars) + primitive variables [vx, vy, vz, T, Yk]: 2*4 + 2*nspecies 
//    //int nstats = 2*nvars+8+2*nspecies;
//
//    //// Get all nstats at xcross for all j and k, and store in data_xcross
//    //amrex::Gpu::DeviceVector<Real> data_xcross(nstats*n_cells[1]*n_cells[2], 0.0); // values at x* for a given y and z
//    //GetPencilCross(data_xcross,consMean,primMean,prim_in,cons,nstats);
//    //ParallelDescriptor::ReduceRealSum(data_xcross.data(),nstats*n_cells[1]*n_cells[2]);
//
//    //// Update Spatial Correlations
//    //EvaluateSpatialCorrelations2D(spatialCross2D,data_xcross,consMean,primMean,prim_in,cons,vel,velMean,cumom,cumomMean,steps,nstats);
//}


///////////////////////////////////////////
// Evaluate Stats for the 1D case /////////
/// ///////////////////////////////////////
//void evaluateStatsStag1D(MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
//                         MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
//                         const std::array<MultiFab, AMREX_SPACEDIM>& vel, 
//                         std::array<MultiFab, AMREX_SPACEDIM>& velMean, 
//                         std::array<MultiFab, AMREX_SPACEDIM>& velVar, 
//                         const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
//                         std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
//                         std::array<MultiFab, AMREX_SPACEDIM>& cumomVar,
//                         MultiFab& coVar, 
//                         MultiFab& theta, MultiFab& thetaMean, MultiFab& thetaVar,
//                         MultiFab& spatialCross1D, const int ncross,
//                         const int steps,
//                         const Geometry& geom)
//{
//    BL_PROFILE_VAR("evaluateStatsStag1D()",evaluateStatsStag1D);
//    
//    //// Evaluate Means
//    if (plot_means) {
//        EvaluateStatsMeans(cons,consMean,prim_in,primMean,vel,velMean,cumom,cumomMean,theta,thetaMean,steps);
//        consMean.FillBoundary(geom.periodicity());
//        primMean.FillBoundary(geom.periodicity());
//    }
//
//    //// Evaluate Variances and Covariances
//    if ((plot_vars) or (plot_covars)) {
//        EvaluateVarsCoVars(cons,consMean,consVar,prim_in,primMean,primVar,velMean,velVar,
//                           cumom,cumomMean,cumomVar,coVar,theta,thetaMean,thetaVar,steps);
//        consVar.FillBoundary(geom.periodicity());
//        primVar.FillBoundary(geom.periodicity());
//    }
//
//    //// Evaluate Spatial Correlations
//    
//    // contains yz-averaged running & instantaneous averages of conserved variables (2*nvars) + primitive variables [vx, vy, vz, T, Yk]: 2*4 + 2*nspecies 
//    int nstats = 2*nvars+8+2*nspecies;
//
//    // Get all nstats at xcross for all j and k, and store in data_xcross
//    if (plot_cross) {
//        if (all_correl == 0) { // for a specified cross_cell
//            amrex::Gpu::DeviceVector<Real> data_xcross(nstats*n_cells[1]*n_cells[2], 0.0); // values at x* for a given y and z
//            GetPencilCross(data_xcross,consMean,primMean,prim_in,cons,nstats,cross_cell);
//            ParallelDescriptor::ReduceRealSum(data_xcross.data(),nstats*n_cells[1]*n_cells[2]);
//
//            // Update Spatial Correlations
//            EvaluateSpatialCorrelations1D(spatialCross1D,data_xcross,consMean,primMean,prim_in,cons,vel,velMean,cumom,cumomMean,steps,nstats,ncross,0);
//        }
//        else { // at five equi-distant x*
//            GpuArray<int, 5 > x_star;
//            x_star[0] = 0;
//            x_star[1] = (int)amrex::Math::floor(1.0*n_cells[0]/4.0);
//            x_star[2] = (int)amrex::Math::floor(2.0*n_cells[0]/4.0);
//            x_star[3] = (int)amrex::Math::floor(3.0*n_cells[0]/4.0);
//            x_star[4] = n_cells[0] - 1;
//
//            amrex::Gpu::DeviceVector<Real> data_xcross(nstats*n_cells[1]*n_cells[2], 0.0); // values at x* for a given y and z
//            for (int i=0; i<5; ++i) {
//                GetPencilCross(data_xcross,consMean,primMean,prim_in,cons,nstats,x_star[i]);
//                EvaluateSpatialCorrelations1D(spatialCross1D,data_xcross,consMean,primMean,prim_in,cons,vel,velMean,cumom,cumomMean,steps,nstats,ncross,i);
//            }
//        }
//    }
//}


///////////////////////
// evaluate_means
///////////////////////
void EvaluateStatsMeans(MultiFab& cons, MultiFab& consMean,
                        MultiFab& prim_in, MultiFab& primMean,
                        const int steps)
{
    BL_PROFILE_VAR("EvaluateStatsMeans()",EvaluateStatsMeans);
    
    double stepsminusone = steps - 1.;
    double stepsinv = 1./steps;

    // Loop over boxes
    for ( MFIter mfi(prim_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx  = mfi.tilebox();

        const Array4<      Real> cu        = cons.array(mfi);
        const Array4<      Real> cumeans   = consMean.array(mfi);
        const Array4<      Real> prim      = prim_in.array(mfi);
        const Array4<      Real> primmeans = primMean.array(mfi);

        // update mean other values (primitive and conserved) 
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {    
            GpuArray<Real,MAX_SPECIES> fracvec;
            cumeans(i,j,k,0) = (cumeans(i,j,k,0)*stepsminusone + cu(i,j,k,0))*stepsinv;
            cumeans(i,j,k,1) = (cumeans(i,j,k,1)*stepsminusone + cu(i,j,k,1))*stepsinv;
            cumeans(i,j,k,2) = (cumeans(i,j,k,2)*stepsminusone + cu(i,j,k,2))*stepsinv;
            cumeans(i,j,k,3) = (cumeans(i,j,k,3)*stepsminusone + cu(i,j,k,3))*stepsinv;
            cumeans(i,j,k,4) = (cumeans(i,j,k,4)*stepsminusone + cu(i,j,k,4))*stepsinv; //rhoEmeans

            for (int l=5; l<nvars; ++l) {
                cumeans(i,j,k,l) = (cumeans(i,j,k,l)*stepsminusone + cu(i,j,k,l))*stepsinv; //rhoYkmeans
                fracvec[l-5] = cumeans(i,j,k,l)/cumeans(i,j,k,0); // Ykmeans
                primmeans(i,j,k,l+1) = fracvec[l-5]; // Ykmeans
            }
            for (int l=0; l<nspecies; ++l) {
                primmeans(i,j,k,6+nspecies+l) = (primmeans(i,j,k,6+nspecies+l)*stepsminusone + prim(i,j,k,6+nspecies+l))*stepsinv; // Xkmeans
            }

            Real densitymeaninv = 1.0/cumeans(i,j,k,0);
            primmeans(i,j,k,0) = cumeans(i,j,k,0); //rhomeans
            primmeans(i,j,k,1) = densitymeaninv*cumeans(i,j,k,1);
            primmeans(i,j,k,2) = densitymeaninv*cumeans(i,j,k,2);
            primmeans(i,j,k,3) = densitymeaninv*cumeans(i,j,k,3);

            // mean COM velocity
            primmeans(i,j,k,nprimvars+0) = (primmeans(i,j,k,nprimvars+0)*stepsminusone + prim(i,j,k,1))*stepsinv;
            primmeans(i,j,k,nprimvars+1) = (primmeans(i,j,k,nprimvars+1)*stepsminusone + prim(i,j,k,2))*stepsinv;
            primmeans(i,j,k,nprimvars+2) = (primmeans(i,j,k,nprimvars+2)*stepsminusone + prim(i,j,k,3))*stepsinv;

            Real kinenergy = 0.;
            kinenergy += cumeans(i,j,k,1)*cumeans(i,j,k,1);
            kinenergy += cumeans(i,j,k,2)*cumeans(i,j,k,2);
            kinenergy += cumeans(i,j,k,3)*cumeans(i,j,k,3);
            kinenergy *= (0.5/cumeans(i,j,k,0));

            Real intenergy = (cumeans(i,j,k,4)-kinenergy)/cumeans(i,j,k,0);

            //primmeans(i,j,k,4) = (primmeans(i,j,k,4)*stepsminusone + prim(i,j,k,4))*stepsinv; // Tmean
            GetTemperature(intenergy, fracvec, primmeans(i,j,k,4)); // Tmean
            //primmeans(i,j,k,5) = (primmeans(i,j,k,5)*stepsminusone + prim(i,j,k,5))*stepsinv; // Pmean
            GetPressureGas(primmeans(i,j,k,5), fracvec, cumeans(i,j,k,0), primmeans(i,j,k,4)); // Pmean

        });

    } // end MFIter
}   

/////////////////////////////////////
// evaluate variances and covariances
/////////////////////////////////////
void EvaluateVarsCoVars(const MultiFab& cons, const MultiFab& consMean, MultiFab& consVar,
                        const MultiFab& prim_in, const MultiFab& primMean, MultiFab& primVar,
                        MultiFab& coVar,const int steps)
{
    BL_PROFILE_VAR("EvaluateVarsCoVars()",EvaluateVarsCoVars);
    
    double stepsminusone = steps - 1.;
    double stepsinv = 1./steps;

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

        // other cell-centered variances and covariances (temperature fluctuation and variance)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {    

            // conserved variable variances (rho, rhoE, rhoYk)
            Real delrho = cu(i,j,k,0) - cumeans(i,j,k,0);
            cuvars(i,j,k,0) = (cuvars(i,j,k,0)*stepsminusone + delrho*delrho)*stepsinv; // <rho rho>
            
            Real delenergy = cu(i,j,k,4) - cumeans(i,j,k,4);
            cuvars(i,j,k,4) = (cuvars(i,j,k,4)*stepsminusone + delenergy*delenergy)*stepsinv; // <rhoE rhoE>
            
            for (int l=5; l<nvars; ++l) {
                Real delmassvec = cu(i,j,k,l) - cumeans(i,j,k,l);
                cuvars(i,j,k,l) = (cuvars(i,j,k,l)*stepsminusone + delmassvec*delmassvec)*stepsinv; // <rhoYk rhoYk>
            }

            // del rhoYk --> delYk
            GpuArray<Real,MAX_SPECIES> delrhoYk;
            GpuArray<Real,MAX_SPECIES> Ykmean;
            GpuArray<Real,MAX_SPECIES> delYk;
            for (int ns=0; ns<nspecies; ++ns) {
                delrhoYk[ns] = cu(i,j,k,5+ns) - cumeans(i,j,k,5+ns); // delrhoYk
                Ykmean[ns] = primmeans(i,j,k,6+ns); // Ykmean
                delYk[ns] = (delrhoYk[ns] - Ykmean[ns]*delrho)/cumeans(i,j,k,0); // delYk = (delrhoYk - Ykmean*delrho)/rhomean
                primvars(i,j,k,6+ns) = (primvars(i,j,k,6+ns)*stepsminusone + delYk[ns]*delYk[ns])*stepsinv;
            }
            
            // primitive variable variances (rho)
            primvars(i,j,k,0) = cuvars(i,j,k,0);

            Real vx = primmeans(i,j,k,1);
            Real vy = primmeans(i,j,k,2);
            Real vz = primmeans(i,j,k,3);
            Real T  = primmeans(i,j,k,4);
            Real densitymeaninv = 1.0/cumeans(i,j,k,0);

            Real cv = 0.;
            for (int l=0; l<nspecies; ++l) {
                cv = cv + hcv[l]*cumeans(i,j,k,5+l)/cumeans(i,j,k,0);
            }
            Real cvinv = 1.0/cv;

            Real qmean = cv*T-0.5*(vx*vx + vy*vy + vz*vz);

            Real deljx = cu(i,j,k,1) - cumeans(i,j,k,1);
            Real deljy = cu(i,j,k,2) - cumeans(i,j,k,2);
            Real deljz = cu(i,j,k,3) - cumeans(i,j,k,3);

            cuvars(i,j,k,1) = (cuvars(i,j,k,1)*stepsminusone + deljx*deljx)*stepsinv; 
            cuvars(i,j,k,2) = (cuvars(i,j,k,2)*stepsminusone + deljy*deljy)*stepsinv; 
            cuvars(i,j,k,3) = (cuvars(i,j,k,3)*stepsminusone + deljz*deljz)*stepsinv; 

            Real delvelx = (deljx - vx*delrho)*densitymeaninv;
            Real delvely = (deljy - vy*delrho)*densitymeaninv;
            Real delvelz = (deljz - vz*delrho)*densitymeaninv;

            primvars(i,j,k,1) = (primvars(i,j,k,1)*stepsminusone + delvelx*delvelx)*stepsinv; 
            primvars(i,j,k,2) = (primvars(i,j,k,2)*stepsminusone + delvely*delvely)*stepsinv; 
            primvars(i,j,k,3) = (primvars(i,j,k,3)*stepsminusone + delvelz*delvelz)*stepsinv; 

            Real delg = vx*deljx + vy*deljy + vz*deljz;

            primvars(i,j,k,nprimvars+0) = (primvars(i,j,k,nprimvars+0)*stepsminusone + delg*delg)*stepsinv; // gvar
            primvars(i,j,k,nprimvars+1) = (primvars(i,j,k,nprimvars+1)*stepsminusone + delg*delenergy)*stepsinv; // kgcross
            primvars(i,j,k,nprimvars+2) = (primvars(i,j,k,nprimvars+2)*stepsminusone + delrho*delenergy)*stepsinv; // krcross
            primvars(i,j,k,nprimvars+3) = (primvars(i,j,k,nprimvars+3)*stepsminusone + delrho*delg)*stepsinv; // rgcross

            primvars(i,j,k,nprimvars+4) = (primvars(i,j,k,nprimvars+4)*stepsminusone + cvinv*cvinv*densitymeaninv*densitymeaninv*
                                 (cuvars(i,j,k,4) + primvars(i,j,k,nprimvars+0) - 2*primvars(i,j,k,nprimvars+1)
                                  + qmean*(qmean*cuvars(i,j,k,0) - 2*primvars(i,j,k,nprimvars+2) + 2*primvars(i,j,k,nprimvars+3))))*stepsinv; // <T T>
            
            // use instantaneous value (the above presumably does not work for multispecies)
            Real delT = prim(i,j,k,4) - primmeans(i,j,k,4);
            primvars(i,j,k,4)   = (primvars(i,j,k,4)*stepsminusone + delT*delT)*stepsinv; // <T T> -- direct

            //Real deltemp = (delenergy - delg - qmean*delrho)*cvinv*densitymeaninv;
            Real deltemp = delT; 

            covars(i,j,k,0)  = (covars(i,j,k,0)*stepsminusone + delrho*deljx)*stepsinv; // <rho jx>
            covars(i,j,k,1)  = (covars(i,j,k,1)*stepsminusone + delrho*deljy)*stepsinv; // <rho jy>
            covars(i,j,k,2)  = (covars(i,j,k,2)*stepsminusone + delrho*deljz)*stepsinv; // <rho jz>
            covars(i,j,k,3)  = (covars(i,j,k,3)*stepsminusone + deljx*deljy)*stepsinv; // <jx jy>
            covars(i,j,k,4)  = (covars(i,j,k,4)*stepsminusone + deljy*deljz)*stepsinv; // <jy jz>
            covars(i,j,k,5)  = (covars(i,j,k,5)*stepsminusone + deljx*deljz)*stepsinv; // <jx jz>
            covars(i,j,k,6)  = (covars(i,j,k,6)*stepsminusone + delrho*delenergy)*stepsinv; // <rho rhoE>
            covars(i,j,k,7)  = (covars(i,j,k,7)*stepsminusone + deljx*delenergy)*stepsinv; // <jx rhoE>
            covars(i,j,k,8)  = (covars(i,j,k,8)*stepsminusone + deljy*delenergy)*stepsinv; // <jy rhoE>
            covars(i,j,k,9)  = (covars(i,j,k,9)*stepsminusone + deljz*delenergy)*stepsinv; // <jz rhoE>
            covars(i,j,k,10) = (covars(i,j,k,10)*stepsminusone + delrhoYk[0]*delrhoYk[nspecies-1])*stepsinv; // <rhoYklightest rhoYkheaviest>
            covars(i,j,k,11) = (covars(i,j,k,11)*stepsminusone + delrho*delvelx)*stepsinv; // <rho vx>
            covars(i,j,k,12) = (covars(i,j,k,12)*stepsminusone + delrho*delvely)*stepsinv; // <rho vy>
            covars(i,j,k,13) = (covars(i,j,k,13)*stepsminusone + delrho*delvelz)*stepsinv; // <rho vz>
            covars(i,j,k,14) = (covars(i,j,k,14)*stepsminusone + delvelx*delvely)*stepsinv; // <vx vy>
            covars(i,j,k,15) = (covars(i,j,k,15)*stepsminusone + delvely*delvelz)*stepsinv; // <vy vz>
            covars(i,j,k,16) = (covars(i,j,k,16)*stepsminusone + delvelx*delvelz)*stepsinv; // <vx vz>
            covars(i,j,k,17) = (covars(i,j,k,17)*stepsminusone + delrho*deltemp)*stepsinv; // <rho T>
            covars(i,j,k,18) = (covars(i,j,k,18)*stepsminusone + delvelx*deltemp)*stepsinv; // <vx T>
            covars(i,j,k,19) = (covars(i,j,k,19)*stepsminusone + delvely*deltemp)*stepsinv; // <vy T>
            covars(i,j,k,20) = (covars(i,j,k,20)*stepsminusone + delvelz*deltemp)*stepsinv; // <vz T>
            covars(i,j,k,21) = (covars(i,j,k,21)*stepsminusone + delYk[0]*delYk[nspecies-1])*stepsinv; // <Yklightest Ykheaviest>
            covars(i,j,k,22) = (covars(i,j,k,22)*stepsminusone + delYk[0]*delvelx)*stepsinv; // <Yklightest velx>
            covars(i,j,k,23) = (covars(i,j,k,23)*stepsminusone + delYk[nspecies-1]*delvelx)*stepsinv; // <Ykheaviest velx>
            covars(i,j,k,24) = (covars(i,j,k,24)*stepsminusone + delrhoYk[0]*delvelx)*stepsinv; // <rhoYklightest velx>
            covars(i,j,k,25) = (covars(i,j,k,25)*stepsminusone + delrhoYk[nspecies-1]*delvelx)*stepsinv; // <rhoYkheaviest velx>

        });

    } // end MFIter
}

///////////////////////////////////////////
// Get Slice Average at x and x* //////////
/// ///////////////////////////////////////
//void GetSliceAverageCross(Vector<Real>& dataAvMeans_x,
//                         Vector<Real>& dataAvMeans_xcross,
//                         const MultiFab& consMean,
//                         const MultiFab& primMean,
//                         const MultiFab& prim_in,
//                         const MultiFab& cons,
//                         const std::array<MultiFab, AMREX_SPACEDIM>& vel,
//                         const std::array<MultiFab, AMREX_SPACEDIM>& velMean,
//                         const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
//                         const std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
//                         const int nstats)
//{
//    BL_PROFILE_VAR("GetSliceAverageCross()",GetSliceAverageCross);
//
//    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {
//
//        const Box& bx = mfi.validbox();
//
//        const auto lo = amrex::lbound(bx);
//        const auto hi = amrex::ubound(bx);
//
//        const Array4<const Real> cumeans   = consMean.array(mfi);
//        const Array4<const Real> primmeans = primMean.array(mfi);
//        const Array4<const Real> prim      = prim_in.array(mfi);
//        const Array4<const Real> cu        = cons.array(mfi);
//
//        const Array4<const Real> velx      = vel[0].array(mfi);
//        const Array4<const Real> vely      = vel[1].array(mfi);
//        const Array4<const Real> velz      = vel[2].array(mfi);
//        const Array4<const Real> velxmeans = velMean[0].array(mfi);
//        const Array4<const Real> velymeans = velMean[1].array(mfi);
//        const Array4<const Real> velzmeans = velMean[2].array(mfi);
//
//        const Array4<const Real> momx      = cumom[0].array(mfi);
//        const Array4<const Real> momy      = cumom[1].array(mfi);
//        const Array4<const Real> momz      = cumom[2].array(mfi);
//        const Array4<const Real> momxmeans = cumomMean[0].array(mfi);
//        const Array4<const Real> momymeans = cumomMean[1].array(mfi);
//        const Array4<const Real> momzmeans = cumomMean[2].array(mfi);
//
//        for (int i=0; i<nstats; ++i) {
//            dataAvMeans_xcross[i] = 0.;
//        }
//        
//        for (auto k = lo.z; k <= hi.z; ++k) {
//        for (auto j = lo.y; j <= hi.y; ++j) {
//        for (auto i = lo.x; i <= hi.x; ++i) {
//            if (i==cross_cell) {
//                dataAvMeans_xcross[0]  += cu(i,j,k,0);                                 // rho-instant
//                dataAvMeans_xcross[1]  += cumeans(i,j,k,0);                            // rho-mean
//                dataAvMeans_xcross[2]  += cu(i,j,k,4);                                 // energy-instant
//                dataAvMeans_xcross[3]  += cumeans(i,j,k,4);                            // energy-mean
//                dataAvMeans_xcross[4]  += 0.5*(momx(i,j,k) + momx(i+1,j,k));           // jx-instant
//                dataAvMeans_xcross[5]  += 0.5*(momxmeans(i,j,k) + momxmeans(i+1,j,k)); // jx-mean
//                dataAvMeans_xcross[6]  += 0.5*(momy(i,j,k) + momy(i,j+1,k));           // jy-instant
//                dataAvMeans_xcross[7]  += 0.5*(momymeans(i,j,k) + momymeans(i,j+1,k)); // jy-mean
//                dataAvMeans_xcross[8]  += 0.5*(momz(i,j,k) + momz(i,j,k+1));           // jz-instant
//                dataAvMeans_xcross[9]  += 0.5*(momzmeans(i,j,k) + momzmeans(i,j,k+1)); // jz-mean
//                dataAvMeans_xcross[10] += 0.5*(velx(i,j,k) + velx(i+1,j,k));           // velx-instant
//                dataAvMeans_xcross[11] += 0.5*(velxmeans(i,j,k) + velxmeans(i+1,j,k)); // velx-mean
//                dataAvMeans_xcross[12] += 0.5*(vely(i,j,k) + vely(i,j+1,k));           // vely-instant
//                dataAvMeans_xcross[13] += 0.5*(velymeans(i,j,k) + velymeans(i,j+1,k)); // vely-mean
//                dataAvMeans_xcross[14] += 0.5*(velz(i,j,k) + velz(i,j,k+1));           // velz-instant
//                dataAvMeans_xcross[15] += 0.5*(velzmeans(i,j,k) + velzmeans(i,j,k+1)); // velz-mean
//                dataAvMeans_xcross[16] += prim(i,j,k,4);                               // T-instant
//                dataAvMeans_xcross[17] += primmeans(i,j,k,4);                          // T-mean
//                for (int ns=0; ns<nspecies; ++ns) {
//                    dataAvMeans_xcross[18+4*ns+0]   += cu(i,j,k,5+ns);                 // rhoYk-instant
//                    dataAvMeans_xcross[18+4*ns+1]   += cumeans(i,j,k,5+ns);            // rhoYk-mean
//                    dataAvMeans_xcross[18+4*ns+2]   += prim(i,j,k,6+ns);               // Yk-instant
//                    dataAvMeans_xcross[18+4*ns+3]   += primmeans(i,j,k,6+ns);          // Yk-mean
//                }
//            }
//            dataAvMeans_x[i*nstats+0]  += cu(i,j,k,0);                                 // rho-instant
//            dataAvMeans_x[i*nstats+1]  += cumeans(i,j,k,0);                            // rho-mean
//            dataAvMeans_x[i*nstats+2]  += cu(i,j,k,4);                                 // energy-instant
//            dataAvMeans_x[i*nstats+3]  += cumeans(i,j,k,4);                            // energy-mean
//            dataAvMeans_x[i*nstats+4]  += 0.5*(momx(i,j,k) + momx(i+1,j,k));           // jx-instant
//            dataAvMeans_x[i*nstats+5]  += 0.5*(momxmeans(i,j,k) + momxmeans(i+1,j,k)); // jx-mean
//            dataAvMeans_x[i*nstats+6]  += 0.5*(momy(i,j,k) + momy(i,j+1,k));           // jy-instant
//            dataAvMeans_x[i*nstats+7]  += 0.5*(momymeans(i,j,k) + momymeans(i,j+1,k)); // jy-mean
//            dataAvMeans_x[i*nstats+8]  += 0.5*(momz(i,j,k) + momz(i,j,k+1));           // jz-instant
//            dataAvMeans_x[i*nstats+9]  += 0.5*(momzmeans(i,j,k) + momzmeans(i,j,k+1)); // jz-mean
//            dataAvMeans_x[i*nstats+10] += 0.5*(velx(i,j,k) + velx(i+1,j,k));           // velx-instant
//            dataAvMeans_x[i*nstats+11] += 0.5*(velxmeans(i,j,k) + velxmeans(i+1,j,k)); // velx-mean
//            dataAvMeans_x[i*nstats+12] += 0.5*(vely(i,j,k) + vely(i,j+1,k));           // vely-instant
//            dataAvMeans_x[i*nstats+13] += 0.5*(velymeans(i,j,k) + velymeans(i,j+1,k)); // vely-mean
//            dataAvMeans_x[i*nstats+14] += 0.5*(velz(i,j,k) + velz(i,j,k+1));           // velz-instant
//            dataAvMeans_x[i*nstats+15] += 0.5*(velzmeans(i,j,k) + velzmeans(i,j,k+1)); // velz-mean
//            dataAvMeans_x[i*nstats+16] += prim(i,j,k,4);                               // T-instant
//            dataAvMeans_x[i*nstats+17] += primmeans(i,j,k,4);                          // T-mean
//            for (int ns=0; ns<nspecies; ++ns) {
//                dataAvMeans_x[i*nstats+18+4*ns+0]   += cu(i,j,k,5+ns);                 // rhoYk-instant
//                dataAvMeans_x[i*nstats+18+4*ns+1]   += cumeans(i,j,k,5+ns);            // rhoYk-mean
//                dataAvMeans_x[i*nstats+18+4*ns+2]   += prim(i,j,k,6+ns);               // Yk-instant
//                dataAvMeans_x[i*nstats+18+4*ns+3]   += primmeans(i,j,k,6+ns);          // Yk-mean
//            }
//        }
//        }
//        }
//
//        for (int i=0; i<n_cells[0]*nstats; ++i) {
//            dataAvMeans_x[i] /= (n_cells[1]*n_cells[2]);
//        }
//        for (int i=0; i<nstats; ++i) {
//            dataAvMeans_xcross[i] /= (n_cells[1]*n_cells[2]);
//        }
//
//    } // end MFITer
//}

//////////////////////////////////////////
// Get Pencil values at x* ///////////////
// for all j and k ///////////////////////
// ///////////////////////////////////////
//void GetPencilCross(amrex::Gpu::DeviceVector<Real>& data_xcross_in,
//                    const MultiFab& consMean,
//                    const MultiFab& primMean,
//                    const MultiFab& prim_in,
//                    const MultiFab& cons,
//                    const int nstats,
//                    const int x_star)
//{
//    BL_PROFILE_VAR("GetPencilCross()",GetPencilCross);
//
//    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {
//
//        const Box& bx = mfi.validbox();
//
//        const Array4<const Real> cumeans   = consMean.array(mfi);
//        const Array4<const Real> primmeans = primMean.array(mfi);
//        const Array4<const Real> prim      = prim_in.array(mfi);
//        const Array4<const Real> cu        = cons.array(mfi);
//        
//        Real* data_xcross = data_xcross_in.data();
//        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//        {
//            if (i==x_star) {
//                int index = k*n_cells[1]*nstats + j*nstats;
//                data_xcross[index + 0]  = cu(i,j,k,0);          // rho-instant
//                data_xcross[index + 1]  = cumeans(i,j,k,0);     // rho-mean
//                data_xcross[index + 2]  = cu(i,j,k,4);          // energy-instant
//                data_xcross[index + 3]  = cumeans(i,j,k,4);     // energy-mean
//                data_xcross[index + 4]  = cu(i,j,k,1);          // jx-instant
//                data_xcross[index + 5]  = cumeans(i,j,k,1);     // jx-mean
//                data_xcross[index + 6]  = cu(i,j,k,2);          // jy-instant
//                data_xcross[index + 7]  = cumeans(i,j,k,2);     // jy-mean
//                data_xcross[index + 8]  = cu(i,j,k,3);          // jz-instant
//                data_xcross[index + 9]  = cumeans(i,j,k,3);     // jz-mean
//                data_xcross[index + 10] = prim(i,j,k,1);        // velx-instant
//                data_xcross[index + 11] = primmeans(i,j,k,1);   // velx-mean
//                data_xcross[index + 12] = prim(i,j,k,2);        // vely-instant
//                data_xcross[index + 13] = primmeans(i,j,k,2);   // vely-mean
//                data_xcross[index + 14] = prim(i,j,k,3);        // velz-instant
//                data_xcross[index + 15] = primmeans(i,j,k,3);   // velz-mean
//                data_xcross[index + 16] = prim(i,j,k,4);        // T-instant
//                data_xcross[index + 17] = primmeans(i,j,k,4);   // T-mean
//                for (int ns=0; ns<nspecies; ++ns) {
//                    data_xcross[index + 18+4*ns+0]  = cu(i,j,k,5+ns);        // rhoYk-instant
//                    data_xcross[index + 18+4*ns+1]  = cumeans(i,j,k,5+ns);   // rhoYk-mean
//                    data_xcross[index + 18+4*ns+2]  = prim(i,j,k,6+ns);      // Yk-instant
//                    data_xcross[index + 18+4*ns+3]  = primmeans(i,j,k,6+ns); // Yk-mean
//                }
//           }
//
//        }); // end MFITer
//    }
//}

//////////////////////////////////////////
// Update Spatial Correlations ///////////
// ///////////////////////////////////////
void EvaluateSpatialCorrelations3D(Vector<Real>& spatialCross,
                                   amrex::Gpu::HostVector<Real>& cu_avg, 
                                   amrex::Gpu::HostVector<Real>& cumeans_avg, 
                                   amrex::Gpu::HostVector<Real>& prim_avg, 
                                   amrex::Gpu::HostVector<Real>& primmeans_avg, 
                                   const int steps,
                                   const int /*nstats*/,
                                   const int ncross)
{
    
    BL_PROFILE_VAR("EvaluateSpatialCorrelations3D()",EvaluateSpatialCorrelations3D);
    
    double stepsminusone = steps - 1.;
    double stepsinv = 1./steps;

    int nprims = nspecies + 4;

    Vector<Real> data_xcross(2*nvars+8+2*nspecies, 0.0);

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
    Real meanuxcross = data_xcross[11];

    GpuArray<Real,MAX_SPECIES> meanYkcross;
    GpuArray<Real,MAX_SPECIES> delYkcross;
    GpuArray<Real,MAX_SPECIES> delrhoYkcross;
    for (int ns=0; ns<nspecies; ++ns) {
        meanYkcross[ns] =  data_xcross[18+4*ns+3];
        delrhoYkcross[ns] =  data_xcross[18+4*ns+0] - data_xcross[18+4*ns+1];
        delYkcross[ns] =  data_xcross[18+4*ns+2] - data_xcross[18+4*ns+3];
    }

    // Get fluctuations of the conserved variables at the cross cell
    Real delrhocross = data_xcross[0] - data_xcross[1];
    Real delKcross   = data_xcross[2] - data_xcross[3];
    Real deljxcross  = data_xcross[4] - data_xcross[5];
    Real deljycross  = data_xcross[6] - data_xcross[7];
    Real deljzcross  = data_xcross[8] - data_xcross[9];

    // Get fluctuations of some primitive variables (for direct fluctuation calculations)
    Real delTcross = data_xcross[16] - data_xcross[17];
    Real delvxcross = data_xcross[10] - data_xcross[11];
    
    // evaluate heat stuff at the cross cell
    Real cvcross = 0.;
    for (int l=0; l<nspecies; ++l) {
        cvcross = cvcross + hcv[l]*data_xcross[18+4*l+1]/data_xcross[1];
    }
    Real cvinvcross = 1.0/cvcross;
    Real vxmeancross = data_xcross[11];
    Real vymeancross = data_xcross[13];
    Real vzmeancross = data_xcross[15];
    Real qmeancross = cvcross*data_xcross[17] - 
                     0.5*(vxmeancross*vxmeancross + vymeancross*vymeancross + vzmeancross*vzmeancross); 

    // Get fluctuations of derived hydrodynamic quantities at the cross cell
    // delG = \vec{v}\cdot\vec{\deltaj}
    Real delGcross = vxmeancross*deljxcross + vymeancross*deljycross + vzmeancross*deljzcross;

    /////////////////////////////////////////////////////////////
    // evaluate x-spatial correlations
    /////////////////////////////////////////////////////////////
    // int ncross = 37+nspecies+2; check main_drive.cpp for latest
    for (int i=0; i<n_cells[0]; ++i) {

        // Get mean values
        Real meanrho = cumeans_avg[0+nvars*i];
        GpuArray<Real,MAX_SPECIES>  meanYk;
        GpuArray<Real,MAX_SPECIES> delrhoYk;
        GpuArray<Real,MAX_SPECIES> delYk;
        for (int ns=0; ns<nspecies; ++ns) {
            meanYk[ns] = primmeans_avg[4+ns+nprims*i];
            delrhoYk[ns] = cu_avg[5+ns+nvars*i] - cumeans_avg[5+ns+nvars*i];
            delYk[ns] = prim_avg[4+ns+nprims*i] - primmeans_avg[4+ns+nprims*i];
        }

        // Get fluctuations of the conserved variables
        Real delrho = cu_avg[0+nvars*i] - cumeans_avg[0+nvars*i];
        Real delK   = cu_avg[4+nvars*i] - cumeans_avg[4+nvars*i];
        Real deljx  = cu_avg[1+nvars*i] - cumeans_avg[1+nvars*i];
        Real deljy  = cu_avg[2+nvars*i] - cumeans_avg[2+nvars*i];
        Real deljz  = cu_avg[3+nvars*i] - cumeans_avg[3+nvars*i];

        // Get fluctuations of some primitive variables (for direct fluctuation calculations)
        Real delT = prim_avg[3+nprims*i] - primmeans_avg[3+nprims*i];
        Real delvx = prim_avg[0+nprims*i] - primmeans_avg[0+nprims*i];
    
        // evaluate heat stuff at the cross cell
        Real cv = 0.;
        for (int l=0; l<nspecies; ++l) {
            cv = cv + hcv[l]*cumeans_avg[5+l+nvars*i]/cumeans_avg[0+nvars*i];
        }
        Real cvinv = 1.0/cv;
        Real qmean = cv*primmeans_avg[3+nprims*i] - 
                         0.5*(primmeans_avg[0+nprims*i]*primmeans_avg[0+nprims*i] +
                              primmeans_avg[1+nprims*i]*primmeans_avg[1+nprims*i] + 
                              primmeans_avg[2+nprims*i]*primmeans_avg[2+nprims*i]);

        // Get fluctuations of derived hydrodynamic quantities
        // delG = \vec{v}\cdot\vec{\deltaj}
        Real delG = primmeans_avg[0+nprims*i]*deljx + primmeans_avg[1+nprims*i]*deljy +
                    primmeans_avg[2+nprims*i]*deljz;
    
        // First update correlations of conserved quantities (we will do rhoYk later)
        spatialCross[i*ncross+0]  = (spatialCross[i*ncross+0]*stepsminusone + delrhocross*delrho)*stepsinv; // <delrho(x*)delrho(x)>
        spatialCross[i*ncross+1]  = (spatialCross[i*ncross+1]*stepsminusone + delKcross*delK)*stepsinv;     // <delK(x*)delK(x)>
        spatialCross[i*ncross+2]  = (spatialCross[i*ncross+2]*stepsminusone + deljxcross*deljx)*stepsinv;   // <deljx(x*)deljx(x)>
        spatialCross[i*ncross+3]  = (spatialCross[i*ncross+3]*stepsminusone + deljycross*deljy)*stepsinv;   // <deljy(x*)deljy(x)>
        spatialCross[i*ncross+4]  = (spatialCross[i*ncross+4]*stepsminusone + deljzcross*deljz)*stepsinv;   // <deljz(x*)deljz(x)>
        spatialCross[i*ncross+5]  = (spatialCross[i*ncross+5]*stepsminusone + deljxcross*delrho)*stepsinv;  // <deljx(x*)delrho(x)>
        spatialCross[i*ncross+6]  = (spatialCross[i*ncross+6]*stepsminusone + deljxcross*delrhoYk[0])*stepsinv;  // <deljx(x*)delrhoYkL(x)>
        spatialCross[i*ncross+7]  = (spatialCross[i*ncross+7]*stepsminusone + deljxcross*delrhoYk[nspecies-1])*stepsinv;  // <deljx(x*)delrhoYkH(x)>
        spatialCross[i*ncross+8]  = (spatialCross[i*ncross+8]*stepsminusone + delrhocross*delrhoYk[0])*stepsinv; // <delrho(x*)delrhoYkL(x)>
        spatialCross[i*ncross+9]  = (spatialCross[i*ncross+9]*stepsminusone + delrhocross*delrhoYk[nspecies-1])*stepsinv; // <delrho(x*)delrhoYkH(x)>
        spatialCross[i*ncross+10] = (spatialCross[i*ncross+10]*stepsminusone + delrhoYkcross[0]*delrho)*stepsinv; // <delrhoYkL(x*)delrho(x)>
        spatialCross[i*ncross+11] = (spatialCross[i*ncross+11]*stepsminusone + delrhoYkcross[nspecies-1]*delrho)*stepsinv; // <delrhoYkH(x*)delrho(x)>
        
        // Some more cross-correlations for hydrodynamical variables later
        spatialCross[i*ncross+12] = (spatialCross[i*ncross+12]*stepsminusone + delGcross*delG)*stepsinv; // <delG(x*)delG(x)>
        spatialCross[i*ncross+13] = (spatialCross[i*ncross+13]*stepsminusone + delGcross*delK)*stepsinv; // <delG(x*)delK(x)>
        spatialCross[i*ncross+14] = (spatialCross[i*ncross+14]*stepsminusone + delKcross*delG)*stepsinv; // <delK(x*)delG(x)>
        spatialCross[i*ncross+15] = (spatialCross[i*ncross+15]*stepsminusone + delrhocross*delK)*stepsinv; // <delrho(x*)delK(x)>
        spatialCross[i*ncross+16] = (spatialCross[i*ncross+16]*stepsminusone + delKcross*delrho)*stepsinv; // <delK(x*)delrho(x)>
        spatialCross[i*ncross+17] = (spatialCross[i*ncross+17]*stepsminusone + delrhocross*delG)*stepsinv; // <delrho(x*)delG(x)>
        spatialCross[i*ncross+18] = (spatialCross[i*ncross+18]*stepsminusone + delGcross*delrho)*stepsinv; // <delG(x*)delrho(x)>

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
        spatialCross[i*ncross+21] = (1.0/meanrhocross)*(spatialCross[i*ncross+5] - meanuxcross*spatialCross[i*ncross+0]);  

        // <delu(x*)del(rhoYkL)> = (1/<rho(x*)>)*(<deljx(x*)del(rhoYkL)> - <u(x*)><delrho(x*)del(rhoYkL)>)
        spatialCross[i*ncross+22] = (1.0/meanrhocross)*(spatialCross[i*ncross+6] - meanuxcross*spatialCross[i*ncross+8]);  

        // <delu(x*)del(rhoYkH)> = (1/<rho(x*)>)*(<deljx(x*)del(rhoYkH)> - <u(x*)><delrho(x*)del(rhoYkH)>)
        spatialCross[i*ncross+23] = (1.0/meanrhocross)*(spatialCross[i*ncross+7] - meanuxcross*spatialCross[i*ncross+9]);  

        // <delu(x*)del(YkL)> = (1/<rho(x*)>/<rho(x)>)*(<deljx(x*)del(rhoYkL) - <u(x*)><delrho(x*)del(rhoYkL)> 
        //                      - <YkL(x)><deljx(x*)delrho(x)> + <u(x*)><YkL(x)><delrho(x*)delrho(x)>)
        spatialCross[i*ncross+24] = (1.0/(meanrho*meanrhocross))*(spatialCross[i*ncross+6] - meanuxcross*spatialCross[i*ncross+8] 
                                                                 - meanYk[0]*spatialCross[i*ncross+5] + meanuxcross*meanYk[0]*spatialCross[i*ncross+0]);

        // <delu(x*)del(YkH)> = (1/<rho(x*)>/<rho(x)>)*(<deljx(x*)del(rhoYkH) - <u(x*)><delrho(x*)del(rhoYkH)> 
        //                      - <YkH(x)><deljx(x*)delrho(x)> + <u(x*)><YkH(x)><delrho(x*)delrho(x)>)
        spatialCross[i*ncross+25] = (1.0/(meanrho*meanrhocross))*(spatialCross[i*ncross+7] - meanuxcross*spatialCross[i*ncross+9] 
                                                                 - meanYk[nspecies-1]*spatialCross[i*ncross+5] + meanuxcross*meanYk[nspecies-1]*spatialCross[i*ncross+0]);

        // Direct -- <delT(x*)delT(x)>
        spatialCross[i*ncross+26] = (spatialCross[i*ncross+26]*stepsminusone + delTcross*delT)*stepsinv;

        // Direct -- <delT(x*)delrho(x)>
        spatialCross[i*ncross+27] = (spatialCross[i*ncross+27]*stepsminusone + delTcross*delrho)*stepsinv;

        // Direct -- <delT(x*)delu(x)>
        spatialCross[i*ncross+28] = (spatialCross[i*ncross+28]*stepsminusone + delTcross*delvx)*stepsinv;

        // Direct -- <delu(x*)delrho>
        spatialCross[i*ncross+29] = (spatialCross[i*ncross+29]*stepsminusone + delvxcross*delrho)*stepsinv;

        // Direct -- <delu(x*)del(rhoYkL)
        spatialCross[i*ncross+30] = (spatialCross[i*ncross+30]*stepsminusone + delvxcross*delYk[0])*stepsinv;

        // Direct -- <delu(x*)del(rhoYkH)
        spatialCross[i*ncross+31] = (spatialCross[i*ncross+31]*stepsminusone + delvxcross*delYk[nspecies-1])*stepsinv;

        // Direct -- <delu(x*)del(rhoYkL)
        spatialCross[i*ncross+32] = (spatialCross[i*ncross+32]*stepsminusone + delvxcross*delrhoYk[0])*stepsinv;

        // Direct -- <delu(x*)del(rhoYkH)
        spatialCross[i*ncross+33] = (spatialCross[i*ncross+33]*stepsminusone + delvxcross*delrhoYk[nspecies-1])*stepsinv;

        // Direct <delYkL(x*)delYkL(x)>
        spatialCross[i*ncross+34] = (spatialCross[i*ncross+34]*stepsminusone + delYkcross[0]*delYk[0])*stepsinv;

        // Direct <delYkH(x*)delYkH(x)>
        spatialCross[i*ncross+35] = (spatialCross[i*ncross+35]*stepsminusone + delYkcross[nspecies-1]*delYk[nspecies-1])*stepsinv;

        // Direct <delYkL(x*)delYkH(x)>
        spatialCross[i*ncross+36] = (spatialCross[i*ncross+36]*stepsminusone + delYkcross[0]*delYk[nspecies-1])*stepsinv;

        // Last we rhoYk for species
        for (int ns=0; ns<nspecies; ++ns) {
            spatialCross[i*ncross+37+ns] = (spatialCross[i*ncross+37+ns]*stepsminusone + delrhoYkcross[ns]*delrhoYk[ns])*stepsinv; // <delrhoYk(x*)delrhoYk(x)>
        }

        // <delrhoYkL(x)delrhoYkH(x*)
        spatialCross[i*ncross+37+nspecies+0] = (spatialCross[i*ncross+37+nspecies+0]*stepsminusone + delrhoYkcross[nspecies-1]*delrhoYk[0])*stepsinv;

        // <delYkL(x*)delYkL(x)> = (1/<rho(x*)>/<rho(x)>)*(<delrhoYkL(x*)delrhoYkL> - <YkL(x*)><delrho(x*)delrhoYkL(x)>
        //                                                 - <YkL(x)><delrhoYkL(x*)delrho(x) + <YkL(x*)><YkL(x)><delrho(x*)delrho(x)>)
        Real delrhoYkdelrhoYk = (spatialCross[i*ncross+37]*stepsminusone + delrhoYkcross[0]*delrhoYk[0])*stepsinv;
        spatialCross[i*ncross+37+nspecies+1] = (1.0/(meanrho*meanrhocross))*(delrhoYkdelrhoYk - meanYkcross[0]*spatialCross[i*ncross+8]
                                                                - meanYk[0]*spatialCross[i*ncross+10] + meanYkcross[0]*meanYk[0]*spatialCross[i*ncross+0]);

        // <delYkL(x*)delYkH(x)> = (1/<rho(x*)>/<rho(x)>)*(<delrhoYkH(x*)delrhoYkH> - <YkH(x*)><delrho(x*)delrhoYkH(x)>
        //                                                 - <YkH(x)><delrhoYkH(x*)delrho(x) + <YkH(x*)><YkH(x)><delrho(x*)delrho(x)>)
        delrhoYkdelrhoYk = (spatialCross[i*ncross+37+nspecies-1]*stepsminusone + delrhoYkcross[nspecies-1]*delrhoYk[nspecies-1])*stepsinv;
        spatialCross[i*ncross+37+nspecies+2] = (1.0/(meanrho*meanrhocross))*(delrhoYkdelrhoYk - meanYkcross[nspecies-1]*spatialCross[i*ncross+9]
                                                - meanYk[nspecies-1]*spatialCross[i*ncross+11] + meanYkcross[nspecies-1]*meanYk[nspecies-1]*spatialCross[i*ncross+0]);
    }
}



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
