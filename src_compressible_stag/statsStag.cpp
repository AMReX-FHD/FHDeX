#include "compressible_functions.H"
#include "compressible_functions_stag.H"

#include "common_functions.H"

void evaluateStatsStag(const MultiFab& cons, MultiFab& consMean, MultiFab& consVar,
                       const MultiFab& prim_in, MultiFab& primMean, MultiFab& primVar,
                       const std::array<MultiFab, AMREX_SPACEDIM>& vel, 
                       std::array<MultiFab, AMREX_SPACEDIM>& velMean, 
                       std::array<MultiFab, AMREX_SPACEDIM>& velVar, 
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                       std::array<MultiFab, AMREX_SPACEDIM>& cumomMean,
                       std::array<MultiFab, AMREX_SPACEDIM>& cumomVar,
                       MultiFab& coVar, 
                       Vector<Real>& yzAvMeans_cross,
                       Vector<Real>& spatialCross,
                       const int steps, const amrex::Real* dx)
{
    BL_PROFILE_VAR("evaluateStatsStag()",evaluateStatsStag);
    
    double stepsminusone = steps - 1.;
    double stepsinv = 1./steps;

    //////////////////
    // evaluate_means
    //////////////////
    
    // Loop over boxes
    for ( MFIter mfi(prim_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx  = mfi.tilebox();
        const Box& bxg = mfi.growntilebox(ngc[0]);
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<const Real> cu        = cons.array(mfi);
        const Array4<      Real> cumeans   = consMean.array(mfi);
        const Array4<      Real> primmeans = primMean.array(mfi);

        const Array4<      Real> velxmeans = velMean[0].array(mfi);
        const Array4<      Real> velymeans = velMean[1].array(mfi);
        const Array4<      Real> velzmeans = velMean[2].array(mfi);

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
             momxmeans(i,j,k) = (momxmeans(i,j,k)*stepsminusone + momx(i,j,k))*stepsinv;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) 
        {
            momymeans(i,j,k) = (momymeans(i,j,k)*stepsminusone + momy(i,j,k))*stepsinv;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) 
        {
            momzmeans(i,j,k) = (momzmeans(i,j,k)*stepsminusone + momz(i,j,k))*stepsinv;
        });

        // update mean density (required in the ghost region as well)
        amrex::ParallelFor(bxg, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cumeans(i,j,k,0) = (cumeans(i,j,k,0)*stepsminusone + cu(i,j,k,0))*stepsinv;
        });

        // update mean other values (primitive and conserved) 
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {    
            GpuArray<Real,MAX_SPECIES> fracvec;
            cumeans(i,j,k,1) = 0.5*(momxmeans(i,j,k) + momxmeans(i+1,j,k)); // jxmeans on CC
            cumeans(i,j,k,2) = 0.5*(momymeans(i,j,k) + momymeans(i,j+1,k)); // jymeans on CC
            cumeans(i,j,k,3) = 0.5*(momzmeans(i,j,k) + momzmeans(i,j,k+1)); // jzmeans on CC
            cumeans(i,j,k,4) = (cumeans(i,j,k,4)*stepsminusone + cu(i,j,k,4))*stepsinv; //rhoEmeans

            for (int l=5; l<nvars; ++l) {
                cumeans(i,j,k,l) = (cumeans(i,j,k,l)*stepsminusone + cu(i,j,k,l))*stepsinv; //rhoYkmeans
                fracvec[l-5] = cumeans(i,j,k,l)/cumeans(i,j,k,0); // Ykmeans
                primmeans(i,j,k,l+1) = fracvec[l-5]; // Ykmeans
            }

            Real densitymeaninv = 1.0/cumeans(i,j,k,0);
            primmeans(i,j,k,1) = densitymeaninv*0.5*(momxmeans(i,j,k) + momxmeans(i+1,j,k)); // velxmeans on CC
            primmeans(i,j,k,2) = densitymeaninv*0.5*(momymeans(i,j,k) + momymeans(i,j+1,k)); // velymeans on CC
            primmeans(i,j,k,3) = densitymeaninv*0.5*(momzmeans(i,j,k) + momzmeans(i,j,k+1)); // velzmeans on CC

            primmeans(i,j,k,0) = cumeans(i,j,k,0); //rhomeans

            Real kinenergy = 0.;
            kinenergy += (momxmeans(i+1,j,k) + momxmeans(i,j,k))*(momxmeans(i+1,j,k) + momxmeans(i,j,k));
            kinenergy += (momymeans(i,j+1,k) + momymeans(i,j,k))*(momymeans(i,j+1,k) + momymeans(i,j,k));
            kinenergy += (momzmeans(i,j,k+1) + momzmeans(i,j,k))*(momzmeans(i,j,k+1) + momzmeans(i,j,k));
            kinenergy *= (0.125/cumeans(i,j,k,0));

            Real intenergy = (cumeans(i,j,k,4)-kinenergy)/cumeans(i,j,k,0);

            GetTemperature(intenergy, fracvec, primmeans(i,j,k,4)); // Tmean
            GetPressureGas(primmeans(i,j,k,5), fracvec, cumeans(i,j,k,0), primmeans(i,j,k,4)); // Pmean

        });

        // update mean velocities
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) 
        {
            Real densitymeaninv = 2.0/(cumeans(i-1,j,k,0)+cumeans(i,j,k,0));
            velxmeans(i,j,k) = momxmeans(i,j,k)*densitymeaninv;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) 
        {
            Real densitymeaninv = 2.0/(cumeans(i,j-1,k,0)+cumeans(i,j,k,0));
            velymeans(i,j,k) = momymeans(i,j,k)*densitymeaninv;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) 
        {
            Real densitymeaninv = 2.0/(cumeans(i,j,k-1,0)+cumeans(i,j,k,0));
            velzmeans(i,j,k) = momzmeans(i,j,k)*densitymeaninv;
        });

    } // end MFIter
    
    /////////////////////////////////////
    // evaluate variances and covariances
    /////////////////////////////////////
    
    // Loop over boxes

    for ( MFIter mfi(prim_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.tilebox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Array4<const Real> cu        = cons.array(mfi);
        const Array4<      Real> cumeans   = consMean.array(mfi);
        const Array4<      Real> cuvars    = consVar.array(mfi);
        const Array4<      Real> primmeans = primMean.array(mfi);
        const Array4<      Real> primvars  = primVar.array(mfi);

        const Array4<      Real> covars    = coVar.array(mfi);

        const Array4<      Real> velxmeans = velMean[0].array(mfi);
        const Array4<      Real> velymeans = velMean[1].array(mfi);
        const Array4<      Real> velzmeans = velMean[2].array(mfi);
        const Array4<      Real> velxvars  = velVar[0].array(mfi);
        const Array4<      Real> velyvars  = velVar[1].array(mfi);
        const Array4<      Real> velzvars  = velVar[2].array(mfi);

        const Array4<const Real> momx      = cumom[0].array(mfi);
        const Array4<const Real> momy      = cumom[1].array(mfi);
        const Array4<const Real> momz      = cumom[2].array(mfi);
        const Array4<      Real> momxmeans = cumomMean[0].array(mfi);
        const Array4<      Real> momymeans = cumomMean[1].array(mfi);
        const Array4<      Real> momzmeans = cumomMean[2].array(mfi);
        const Array4<      Real> momxvars  = cumomVar[0].array(mfi);
        const Array4<      Real> momyvars  = cumomVar[1].array(mfi);
        const Array4<      Real> momzvars  = cumomVar[2].array(mfi);

        // update momentum and velocity variances
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) 
        {
            Real deljx = momx(i,j,k) - momxmeans(i,j,k);
            momxvars(i,j,k) = (momxvars(i,j,k)*stepsminusone + deljx*deljx)*stepsinv; // <jx jx>

            Real densitymeaninv = 2.0/(cumeans(i-1,j,k,0)+cumeans(i,j,k,0));
            Real delrho = 0.5*(cu(i-1,j,k,0) + cu(i,j,k,0)) - 0.5*(cumeans(i-1,j,k,0) + cumeans(i,j,k,0));
            Real delvelx = (deljx - velxmeans(i,j,k)*delrho)*densitymeaninv;
            velxvars(i,j,k) = (velxvars(i,j,k)*stepsminusone + delvelx*delvelx)*stepsinv; // <vx vx>
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) 
        {
            Real deljy = momy(i,j,k) - momymeans(i,j,k);
            momyvars(i,j,k) = (momyvars(i,j,k)*stepsminusone + deljy*deljy)*stepsinv; // <jy jy>

            Real densitymeaninv = 2.0/(cumeans(i,j-1,k,0)+cumeans(i,j,k,0));
            Real delrho = 0.5*(cu(i,j-1,k,0) + cu(i,j,k,0)) - 0.5*(cumeans(i,j-1,k,0) + cumeans(i,j,k,0));
            Real delvely = (deljy - velymeans(i,j,k)*delrho)*densitymeaninv;
            velyvars(i,j,k) = (velyvars(i,j,k)*stepsminusone + delvely*delvely)*stepsinv; // <vx vx>
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) 
        {
            Real deljz = momz(i,j,k) - momzmeans(i,j,k);
            momzvars(i,j,k) = (momzvars(i,j,k)*stepsminusone + deljz*deljz)*stepsinv; // <jz jz>

            Real densitymeaninv = 2.0/(cumeans(i,j,k-1,0)+cumeans(i,j,k,0));
            Real delrho = 0.5*(cu(i,j,k-1,0) + cu(i,j,k,0)) - 0.5*(cumeans(i,j,k-1,0) + cumeans(i,j,k,0));
            Real delvelz = (deljz - velzmeans(i,j,k)*delrho)*densitymeaninv;
            velzvars(i,j,k) = (velzvars(i,j,k)*stepsminusone + delvelz*delvelz)*stepsinv; // <vz vz>
        });

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

            Real vx = 0.5*(velxmeans(i,j,k) + velxmeans(i+1,j,k));
            Real vy = 0.5*(velymeans(i,j,k) + velymeans(i,j+1,k));
            Real vz = 0.5*(velzmeans(i,j,k) + velzmeans(i,j,k+1));
            Real T = primmeans(i,j,k,4);
            Real densitymeaninv = 1.0/cumeans(i,j,k,0);

            Real cv = 0.;
            for (int l=0; l<nspecies; ++l) {
                cv = cv + hcv[l]*cumeans(i,j,k,5+l)/cumeans(i,j,k,0);
            }
            Real cvinv = 1.0/cv;

            Real qmean = cv*T-0.5*(vx*vx + vy*vy + vz*vz);

            Real deljx = 0.5*(momx(i,j,k) + momx(i+1,j,k)) - 0.5*(momxmeans(i,j,k) + momxmeans(i+1,j,k));
            Real deljy = 0.5*(momy(i,j,k) + momy(i,j+1,k)) - 0.5*(momymeans(i,j,k) + momymeans(i,j+1,k));
            Real deljz = 0.5*(momz(i,j,k) + momz(i,j,k+1)) - 0.5*(momzmeans(i,j,k) + momzmeans(i,j,k+1));

            cuvars(i,j,k,1) = (cuvars(i,j,k,1)*stepsminusone + deljx*deljx)*stepsinv; // <jx jx> on CC
            cuvars(i,j,k,2) = (cuvars(i,j,k,2)*stepsminusone + deljy*deljy)*stepsinv; // <jy jy>  on CC
            cuvars(i,j,k,3) = (cuvars(i,j,k,3)*stepsminusone + deljz*deljz)*stepsinv; // <jz jz> on CC

            Real delvelx = (deljx - vx*delrho)*densitymeaninv;
            Real delvely = (deljy - vy*delrho)*densitymeaninv;
            Real delvelz = (deljz - vz*delrho)*densitymeaninv;

            primvars(i,j,k,1) = (primvars(i,j,k,1)*stepsminusone + delvelx*delvelx)*stepsinv; // <vx vx> on CC
            primvars(i,j,k,2) = (primvars(i,j,k,2)*stepsminusone + delvely*delvely)*stepsinv; // <vy vy>  on CC
            primvars(i,j,k,3) = (primvars(i,j,k,3)*stepsminusone + delvelz*delvelz)*stepsinv; // <vz vz> on CC

            Real delg = vx*deljx + vy*deljy + vz*deljz;

            primvars(i,j,k,nprimvars)   = (primvars(i,j,k,nprimvars)*stepsminusone + delg*delg)*stepsinv; // gvar
            primvars(i,j,k,nprimvars+1) = (primvars(i,j,k,nprimvars+1)*stepsminusone + delg*delenergy)*stepsinv; // kgcross
            primvars(i,j,k,nprimvars+2) = (primvars(i,j,k,nprimvars+2)*stepsminusone + delrho*delenergy)*stepsinv; // krcross
            primvars(i,j,k,nprimvars+3) = (primvars(i,j,k,nprimvars+3)*stepsminusone + delrho*delg)*stepsinv; // rgcross

            primvars(i,j,k,4) = (primvars(i,j,k,4)*stepsminusone + cvinv*cvinv*densitymeaninv*densitymeaninv*
                                 (cuvars(i,j,k,4) + primvars(i,j,k,nprimvars) - 2*primvars(i,j,k,nprimvars+1)
                                  + qmean*(qmean*cuvars(i,j,k,0) - 2*primvars(i,j,k,nprimvars+2) + 2*primvars(i,j,k,nprimvars+3))))*stepsinv; // <T T> 

            Real deltemp = (delenergy - delg - qmean*delrho)*cvinv*densitymeaninv;

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

    /////////////////////////////////////////////////////////////
    // evaluate slice average at the cross cell and at every cell
    /////////////////////////////////////////////////////////////

    // contains yz-averaged running & instantaneous averages of conserved variables (2*nvars) + primitive variables [vx, vy, vz, T, Yk]: 2*4 + 2*nspecies 
    int nstats = 2*nvars+8+2*nspecies;
    Vector<Real>  yzAvMeans(n_cells[0]*nstats, 0.0); // yz-average at all x

    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<const Real> cumeans   = consMean.array(mfi);
        const Array4<const Real> primmeans = primMean.array(mfi);
        const Array4<const Real> prim      = prim_in.array(mfi);
        const Array4<const Real> cu        = cons.array(mfi);

        const Array4<const Real> velx      = vel[0].array(mfi);
        const Array4<const Real> vely      = vel[1].array(mfi);
        const Array4<const Real> velz      = vel[2].array(mfi);
        const Array4<      Real> velxmeans = velMean[0].array(mfi);
        const Array4<      Real> velymeans = velMean[1].array(mfi);
        const Array4<      Real> velzmeans = velMean[2].array(mfi);

        const Array4<const Real> momx      = cumom[0].array(mfi);
        const Array4<const Real> momy      = cumom[1].array(mfi);
        const Array4<const Real> momz      = cumom[2].array(mfi);
        const Array4<      Real> momxmeans = cumomMean[0].array(mfi);
        const Array4<      Real> momymeans = cumomMean[1].array(mfi);
        const Array4<      Real> momzmeans = cumomMean[2].array(mfi);

        for (int i=0; i<nstats; ++i) {
            yzAvMeans_cross[i] = 0.;
        }
        
        for (auto k = lo.z; k <= hi.z; ++k) {
        for (auto j = lo.y; j <= hi.y; ++j) {
        for (auto i = lo.x; i <= hi.x; ++i) {
            if (i==cross_cell) {
                yzAvMeans_cross[0]  += cu(i,j,k,0);                                 // rho-instant
                yzAvMeans_cross[1]  += cumeans(i,j,k,0);                            // rho-mean
                yzAvMeans_cross[2]  += cu(i,j,k,4);                                 // energy-instant
                yzAvMeans_cross[3]  += cumeans(i,j,k,4);                            // energy-mean
                yzAvMeans_cross[4]  += 0.5*(momx(i,j,k) + momx(i+1,j,k));           // jx-instant
                yzAvMeans_cross[5]  += 0.5*(momxmeans(i,j,k) + momxmeans(i+1,j,k)); // jx-mean
                yzAvMeans_cross[6]  += 0.5*(momy(i,j,k) + momy(i,j+1,k));           // jy-instant
                yzAvMeans_cross[7]  += 0.5*(momymeans(i,j,k) + momymeans(i,j+1,k)); // jy-mean
                yzAvMeans_cross[8]  += 0.5*(momz(i,j,k) + momz(i,j,k+1));           // jz-instant
                yzAvMeans_cross[9]  += 0.5*(momzmeans(i,j,k) + momzmeans(i,j,k+1)); // jz-mean
                yzAvMeans_cross[10] += 0.5*(velx(i,j,k) + velx(i+1,j,k));           // velx-instant
                yzAvMeans_cross[11] += 0.5*(velxmeans(i,j,k) + velxmeans(i+1,j,k)); // velx-mean
                yzAvMeans_cross[12] += 0.5*(vely(i,j,k) + vely(i,j+1,k));           // vely-instant
                yzAvMeans_cross[13] += 0.5*(velymeans(i,j,k) + velymeans(i,j+1,k)); // vely-mean
                yzAvMeans_cross[14] += 0.5*(velz(i,j,k) + velz(i,j,k+1));           // velz-instant
                yzAvMeans_cross[15] += 0.5*(velzmeans(i,j,k) + velzmeans(i,j,k+1)); // velz-mean
                yzAvMeans_cross[16] += prim(i,j,k,4);                               // T-instant
                yzAvMeans_cross[17] += primmeans(i,j,k,4);                          // T-mean
                for (int ns=0; ns<nspecies; ++ns) {
                    yzAvMeans_cross[18+4*ns+0]   += cu(i,j,k,5+ns);                 // rhoYk-instant
                    yzAvMeans_cross[18+4*ns+1]   += cumeans(i,j,k,5+ns);            // rhoYk-mean
                    yzAvMeans_cross[18+4*ns+2]   += prim(i,j,k,6+ns);               // Yk-instant
                    yzAvMeans_cross[18+4*ns+3]   += primmeans(i,j,k,6+ns);          // Yk-mean
                }
            }
            yzAvMeans[i*nstats+0]  += cu(i,j,k,0);                                 // rho-instant
            yzAvMeans[i*nstats+1]  += cumeans(i,j,k,0);                            // rho-mean
            yzAvMeans[i*nstats+2]  += cu(i,j,k,4);                                 // energy-instant
            yzAvMeans[i*nstats+3]  += cumeans(i,j,k,4);                            // energy-mean
            yzAvMeans[i*nstats+4]  += 0.5*(momx(i,j,k) + momx(i+1,j,k));           // jx-instant
            yzAvMeans[i*nstats+5]  += 0.5*(momxmeans(i,j,k) + momxmeans(i+1,j,k)); // jx-mean
            yzAvMeans[i*nstats+6]  += 0.5*(momy(i,j,k) + momy(i,j+1,k));           // jy-instant
            yzAvMeans[i*nstats+7]  += 0.5*(momymeans(i,j,k) + momymeans(i,j+1,k)); // jy-mean
            yzAvMeans[i*nstats+8]  += 0.5*(momz(i,j,k) + momz(i,j,k+1));           // jz-instant
            yzAvMeans[i*nstats+9]  += 0.5*(momzmeans(i,j,k) + momzmeans(i,j,k+1)); // jz-mean
            yzAvMeans[i*nstats+10] += 0.5*(velx(i,j,k) + velx(i+1,j,k));           // velx-instant
            yzAvMeans[i*nstats+11] += 0.5*(velxmeans(i,j,k) + velxmeans(i+1,j,k)); // velx-mean
            yzAvMeans[i*nstats+12] += 0.5*(vely(i,j,k) + vely(i,j+1,k));           // vely-instant
            yzAvMeans[i*nstats+13] += 0.5*(velymeans(i,j,k) + velymeans(i,j+1,k)); // vely-mean
            yzAvMeans[i*nstats+14] += 0.5*(velz(i,j,k) + velz(i,j,k+1));           // velz-instant
            yzAvMeans[i*nstats+15] += 0.5*(velzmeans(i,j,k) + velzmeans(i,j,k+1)); // velz-mean
            yzAvMeans[i*nstats+16] += prim(i,j,k,4);                               // T-instant
            yzAvMeans[i*nstats+17] += primmeans(i,j,k,4);                          // T-mean
            for (int ns=0; ns<nspecies; ++ns) {
                yzAvMeans[i*nstats+18+4*ns+0]   += cu(i,j,k,5+ns);                 // rhoYk-instant
                yzAvMeans[i*nstats+18+4*ns+1]   += cumeans(i,j,k,5+ns);            // rhoYk-mean
                yzAvMeans[i*nstats+18+4*ns+2]   += prim(i,j,k,6+ns);               // Yk-instant
                yzAvMeans[i*nstats+18+4*ns+3]   += primmeans(i,j,k,6+ns);          // Yk-mean
            }
        }
        }
        }

        for (int i=0; i<n_cells[0]*nstats; ++i) {
            yzAvMeans[i] /= (n_cells[1]*n_cells[2]);
        }
        for (int i=0; i<nstats; ++i) {
            yzAvMeans_cross[i] /= (n_cells[1]*n_cells[2]);
        }

    } // end MFITer

    // parallel reduce sum yzAvMeans and yzAvMeans_cross
    ParallelDescriptor::ReduceRealSum(yzAvMeans.dataPtr(),n_cells[0]*nstats);
    ParallelDescriptor::ReduceRealSum(yzAvMeans_cross.dataPtr(),nstats);

    // Get mean values
    Real meanrhostar = yzAvMeans_cross[1];
    Vector<Real>  meanYkstar(nspecies, 0.0);
    for (int ns=0; ns<nspecies; ++ns) {
        meanYkstar[ns] =  yzAvMeans_cross[18+4*ns+3];
    }
    Real meanuxstar = yzAvMeans_cross[11];

    // Get fluctuations of the conserved variables at the cross cell
    Real delrhostar = yzAvMeans_cross[0] - yzAvMeans_cross[1];
    Real delKstar   = yzAvMeans_cross[2] - yzAvMeans_cross[3];
    Real deljxstar  = yzAvMeans_cross[4] - yzAvMeans_cross[5];
    Real deljystar  = yzAvMeans_cross[6] - yzAvMeans_cross[7];
    Real deljzstar  = yzAvMeans_cross[8] - yzAvMeans_cross[9];
    Vector<Real>  delrhoYkstar(nspecies, 0.0);
    for (int ns=0; ns<nspecies; ++ns) {
        delrhoYkstar[ns] =  yzAvMeans_cross[18+4*ns+0] - yzAvMeans_cross[18+4*ns+1];
    }
    
    // evaluate heat stuff at the cross cell
    Real cv = 0.;
    for (int l=0; l<nspecies; ++l) {
        cv = cv + hcv[l]*yzAvMeans_cross[18+4*l+1]/yzAvMeans_cross[1];
    }
    Real cvinv = 1.0/cv;
    Real qmeanstar = cv*yzAvMeans_cross[17] - 
                     0.5*(yzAvMeans_cross[11]*yzAvMeans_cross[11] + yzAvMeans_cross[13]*yzAvMeans_cross[13] + yzAvMeans_cross[15]*yzAvMeans_cross[15]);

    // Get fluctuations of derived hydrodynamic quantities at the cross cell
    // delG = \vec{v}\cdot\vec{\deltaj}
    Real delGstar = yzAvMeans_cross[11]*(yzAvMeans_cross[4]-yzAvMeans_cross[5]) + yzAvMeans_cross[13]*(yzAvMeans_cross[6]-yzAvMeans_cross[7]) + 
                    yzAvMeans_cross[15]*(yzAvMeans_cross[8]-yzAvMeans_cross[9]);

    /////////////////////////////////////////////////////////////
    // evaluate x-spatial correlations
    /////////////////////////////////////////////////////////////
    int ncross = 28+nspecies;
    for (int i=0; i<n_cells[0]; ++i) {

        // Get mean values
        Real meanrho = yzAvMeans[i*nstats+1];
        Vector<Real>  meanYk(nspecies, 0.0);
        for (int ns=0; ns<nspecies; ++ns) {
            meanYk[ns] = yzAvMeans[i*nstats+18+4*ns+3];
        }

        // Get fluctuations of the conserved variables
        Real delrho = yzAvMeans[i*nstats+0] - yzAvMeans[i*nstats+1];
        Real delK   = yzAvMeans[i*nstats+2] - yzAvMeans[i*nstats+3];
        Real deljx  = yzAvMeans[i*nstats+4] - yzAvMeans[i*nstats+5];
        Real deljy  = yzAvMeans[i*nstats+6] - yzAvMeans[i*nstats+7];
        Real deljz  = yzAvMeans[i*nstats+8] - yzAvMeans[i*nstats+9];
        Vector<Real>  delrhoYk(nspecies, 0.0);
        for (int ns=0; ns<nspecies; ++ns) {
            delrhoYk[ns] = yzAvMeans[i*nstats+18+4*ns+0] - yzAvMeans[i*nstats+18+4*ns+1];
        }
    
        Real qmean = cv*yzAvMeans[i*nstats+17] - 
                         0.5*(yzAvMeans[i*nstats+11]*yzAvMeans[i*nstats+11] + yzAvMeans[i*nstats+13]*yzAvMeans[i*nstats+13] + yzAvMeans[i*nstats+15]*yzAvMeans[i*nstats+15]);

        // Get fluctuations of derived hydrodynamic quantities
        // delG = \vec{v}\cdot\vec{\deltaj}
        Real delG = yzAvMeans[i*nstats+11]*(yzAvMeans[i*nstats+4] - yzAvMeans[i*nstats+5]) + yzAvMeans[i*nstats+13]*(yzAvMeans[i*nstats+6] - yzAvMeans[i*nstats+7]) +
                    yzAvMeans[i*nstats+15]*(yzAvMeans[i*nstats+8] - yzAvMeans[i*nstats+9]);
    
        // First update correlations of conserved quantities (we will do rhoYk later)
        spatialCross[i*ncross+0]  = (spatialCross[i*ncross+0]*stepsminusone + delrhostar*delrho)*stepsinv; // <delrho(x*)delrho(x)>
        spatialCross[i*ncross+1]  = (spatialCross[i*ncross+1]*stepsminusone + delKstar*delK)*stepsinv;     // <delK(x*)delK(x)>
        spatialCross[i*ncross+2]  = (spatialCross[i*ncross+2]*stepsminusone + deljxstar*deljx)*stepsinv;   // <deljx(x*)deljx(x)>
        spatialCross[i*ncross+3]  = (spatialCross[i*ncross+3]*stepsminusone + deljystar*deljy)*stepsinv;   // <deljy(x*)deljy(x)>
        spatialCross[i*ncross+4]  = (spatialCross[i*ncross+4]*stepsminusone + deljzstar*deljz)*stepsinv;   // <deljz(x*)deljz(x)>
        spatialCross[i*ncross+5]  = (spatialCross[i*ncross+5]*stepsminusone + deljxstar*delrho)*stepsinv;  // <deljx(x*)delrho(x)>
        spatialCross[i*ncross+6]  = (spatialCross[i*ncross+6]*stepsminusone + deljxstar*delrhoYk[0])*stepsinv;  // <deljx(x*)delrhoYkL(x)>
        spatialCross[i*ncross+7]  = (spatialCross[i*ncross+7]*stepsminusone + deljxstar*delrhoYk[nspecies-1])*stepsinv;  // <deljx(x*)delrhoYkH(x)>
        spatialCross[i*ncross+8]  = (spatialCross[i*ncross+8]*stepsminusone + delrhostar*delrhoYk[0])*stepsinv; // <delrho(x*)delrhoYkL(x)>
        spatialCross[i*ncross+9]  = (spatialCross[i*ncross+9]*stepsminusone + delrhostar*delrhoYk[nspecies-1])*stepsinv; // <delrho(x*)delrhoYkH(x)>
        spatialCross[i*ncross+10] = (spatialCross[i*ncross+10]*stepsminusone + delrhoYkstar[0]*delrho)*stepsinv; // <delrhoYkL(x*)delrho(x)>
        spatialCross[i*ncross+11] = (spatialCross[i*ncross+11]*stepsminusone + delrhoYkstar[nspecies-1]*delrho)*stepsinv; // <delrhoYkH(x*)delrho(x)>
        
        // Some more cross-correlations for hydrodynamical variables later
        spatialCross[i*ncross+12] = (spatialCross[i*ncross+12]*stepsminusone + delGstar*delG)*stepsinv; // <delG(x*)delG(x)>
        spatialCross[i*ncross+13] = (spatialCross[i*ncross+13]*stepsminusone + delGstar*delK)*stepsinv; // <delG(x*)delK(x)>
        spatialCross[i*ncross+14] = (spatialCross[i*ncross+14]*stepsminusone + delKstar*delG)*stepsinv; // <delK(x*)delG(x)>
        spatialCross[i*ncross+15] = (spatialCross[i*ncross+15]*stepsminusone + delrhostar*delK)*stepsinv; // <delrho(x*)delK(x)>
        spatialCross[i*ncross+16] = (spatialCross[i*ncross+16]*stepsminusone + delKstar*delrho)*stepsinv; // <delK(x*)delrho(x)>
        spatialCross[i*ncross+17] = (spatialCross[i*ncross+17]*stepsminusone + delrhostar*delG)*stepsinv; // <delrho(x*)delG(x)>
        spatialCross[i*ncross+18] = (spatialCross[i*ncross+18]*stepsminusone + delGstar*delrho)*stepsinv; // <delG(x*)delrho(x)>

        // Next we do cross-correlations with and between hydrodynamical variables
        // <delT(x*)delT(x)> = (1/cv*/cv/<rho(x)>/<rho(x*)>)(<delK*delK> + <delG*delG> - <delG*delK> - <delK*delG> 
        //                      + <Q><Q*><delrho*delrho> - <Q*><delrho*delK> - <Q><delK*delrho> + <Q*><delrho*delG> + <Q><delG*delrho>)
        spatialCross[i*ncross+19] = (cvinv*cvinv/(meanrhostar*meanrho))*
                                    (spatialCross[i*ncross+1] + spatialCross[i*ncross+12] - spatialCross[i*ncross+13] - spatialCross[i*ncross+14]
                                     + qmean*qmeanstar*spatialCross[i*ncross+0] - qmeanstar*spatialCross[i*ncross+15] - qmean*spatialCross[i*ncross+16]
                                     + qmeanstar*spatialCross[i*ncross+17] + qmean*spatialCross[i*ncross+18]);

        // <delT(x*)delrho(x)> = (1/cv/<rho(x*)>)*(<delK*delrho> - <delG*delrho> - <Q*><delrhodelrho*>)
        spatialCross[i*ncross+20] = (cvinv*meanrhostar)*(spatialCross[i*ncross+16] - spatialCross[i*ncross+18] - qmeanstar*spatialCross[i*ncross+0]);

        // <delu(x*)delrho> = (1/<rho(x*)>)*(<deljx(x*)delrho(x)> - <u(x*)><<delrho(x*)delrho(x)>) 
        spatialCross[i*ncross+21] = (1.0/meanrhostar)*(spatialCross[i*ncross+5] - meanuxstar*spatialCross[i*ncross+0]);  

        // <delu(x*)del(rhoYkL)> = (1/<rho(x*)>)*(<deljx(x*)del(rhoYkL)> - <u(x*)><delrho(x*)del(rhoYkL)>)
        spatialCross[i*ncross+22] = (1.0/meanrhostar)*(spatialCross[i*ncross+6] - meanuxstar*spatialCross[i*ncross+8]);  

        // <delu(x*)del(rhoYkH)> = (1/<rho(x*)>)*(<deljx(x*)del(rhoYkH)> - <u(x*)><delrho(x*)del(rhoYkH)>)
        spatialCross[i*ncross+23] = (1.0/meanrhostar)*(spatialCross[i*ncross+7] - meanuxstar*spatialCross[i*ncross+9]);  

        // <delu(x*)del(YkL)> = (1/<rho(x*)>/<rho(x)>)*(<deljx(x*)del(rhoYkL) - <u(x*)><delrho(x*)del(rhoYkL)> 
        //                      - <YkL(x)><deljx(x*)delrho(x)> + <u(x*)><YkL(x)><delrho(x*)delrho(x)>)
        spatialCross[i*ncross+24] = (1.0/(meanrho*meanrhostar))*(spatialCross[i*ncross+6] - meanuxstar*spatialCross[i*ncross+8] 
                                                                 - meanYk[0]*spatialCross[i*ncross+5] + meanuxstar*meanYk[0]*spatialCross[i*ncross+0]);

        // <delu(x*)del(YkH)> = (1/<rho(x*)>/<rho(x)>)*(<deljx(x*)del(rhoYkH) - <u(x*)><delrho(x*)del(rhoYkH)> 
        //                      - <YkH(x)><deljx(x*)delrho(x)> + <u(x*)><YkH(x)><delrho(x*)delrho(x)>)
        spatialCross[i*ncross+25] = (1.0/(meanrho*meanrhostar))*(spatialCross[i*ncross+7] - meanuxstar*spatialCross[i*ncross+9] 
                                                                 - meanYk[nspecies-1]*spatialCross[i*ncross+5] + meanuxstar*meanYk[nspecies-1]*spatialCross[i*ncross+0]);

        // <delYkL(x*)delYkL(x)> = (1/<rho(x*)>/<rho(x)>)*(<delrhoYkL(x*)delrhoYkL> - <YkL(x*)><delrho(x*)delrhoYkL(x)>
        //                                                 - <YkL(x)><delrhoYkL(x*)delrho(x) + <YkL(x*)><YkL(x)><delrho(x*)delrho(x)>)
        Real delrhoYkdelrhoYk = (spatialCross[i*ncross+28]*stepsminusone + delrhoYkstar[0]*delrhoYk[0])*stepsinv;
        spatialCross[i*ncross+26] = (1.0/(meanrho*meanrhostar))*(delrhoYkdelrhoYk - meanYkstar[0]*spatialCross[i*ncross+8]
                                                                - meanYk[0]*spatialCross[i*ncross+10] + meanYkstar[0]*meanYk[0]*spatialCross[i*ncross+0]);

        // <delYkL(x*)delYkH(x)> = (1/<rho(x*)>/<rho(x)>)*(<delrhoYkH(x*)delrhoYkH> - <YkH(x*)><delrho(x*)delrhoYkH(x)>
        //                                                 - <YkH(x)><delrhoYkH(x*)delrho(x) + <YkH(x*)><YkH(x)><delrho(x*)delrho(x)>)
        delrhoYkdelrhoYk = (spatialCross[i*ncross+28+nspecies-1]*stepsminusone + delrhoYkstar[nspecies-1]*delrhoYk[nspecies-1])*stepsinv;
        spatialCross[i*ncross+27] = (1.0/(meanrho*meanrhostar))*(delrhoYkdelrhoYk - meanYkstar[nspecies-1]*spatialCross[i*ncross+9]
                                                                - meanYk[nspecies-1]*spatialCross[i*ncross+11] + meanYkstar[nspecies-1]*meanYk[nspecies-1]*spatialCross[i*ncross+0]);

        // Last we rhoYk for species
        for (int ns=0; ns<nspecies; ++ns) {
            spatialCross[i*ncross+28+ns] = (spatialCross[i*ncross+28+ns]*stepsminusone + delrhoYkstar[ns]*delrhoYk[ns])*stepsinv; // <delrhoYk(x*)delrhoYk(x)>
        }
    }
    
    //Vector<Real> yzAv(n_cells[0]*(nvars+2),0.); 
    //Vector<Real> yzAvCross(n_cells[0]*(nvars+2)*(nvars+2),0.); 
    //Vector<Real> yzAv_cross(nvars+2,0.); 

    //for (MFIter mfi(cons, false); mfi.isValid(); ++mfi) {

    //    const Box& bx = mfi.validbox();
    //    const auto lo = amrex::lbound(bx);
    //    const auto hi = amrex::ubound(bx);

    //    const Array4<const Real> primmeans = primMean.array(mfi);
    //    const Array4<const Real> cumeans   = consMean.array(mfi);

    //    for (auto i = lo.x; i <= hi.x; ++i) {
    //    for (auto k = lo.z; k <= hi.z; ++k) {
    //    for (auto j = lo.y; j <= hi.y; ++j) {
    //        for (auto l = 0; l < nvars; ++l) {
    //            yzAv[i*(nvars+2) + l] += cumeans(i,j,k,l); 
    //        }
    //        yzAv[i*(nvars+2) + nvars + 0] += primmeans(i,j,k,1); // <vx>
    //        yzAv[i*(nvars+2) + nvars + 1] += primmeans(i,j,k,4); // <T>
    //        for (auto n = 0; 
    //    }
    //    }
    //    }

    //} // end MFiter

    //// sum over all processors
    //ParallelDescriptor::ReduceRealSum(yzAv.dataPtr(),n_cells[0]*(nvars+2));

    //// divide by the number of yz-face cells
    //int n_face_cells = n_cells[1]*n_cells[2];
    //for (auto i=0; i<n_cells[0]; ++i) {
    //    for (auto l=0; l<nvars+2; ++l) {
    //        yzAv[i*(nvars+2) + l] /= n_face_cells;
    //    }
    //}

    //// copy the cross_cell yz-averaged value in the vector
    //for (auto i=0; i<n_cells[0]; ++i) {
    //    if (i == cross_cell) {
    //        for (auto l=0; l<nvars+2; ++l) {
    //            yzAv_cross[l] = yzAv[i*(nvars+2) + l];
    //        }
    //    }
    //}

    /////////////////////////////////////////////////
    //// increment yz average of conserved variables
    //// to the running mean and compute fluctuation
    /////////////////////////////////////////////////
    //
    //Vector<Real> delyzAv(n_cells[0]*(nvars+2),0.); // yz-average of the variance of conserved variables at current snapshot + two primitive variables [vx, T]
    //Vector<Real> delyzAv_cross(nvars+2,0.); // yz-averaged of the variance conserved variable at current snapshot at the cross cell + two primitive variables [vx, T]

    //for (auto i=0; i<n_cells[0]; ++i) {
    //    for (auto l=0; l<nvars; ++l) {
    //        yzAvMeans[i*nvars + l] = (yzAvMeans[i*nvars + l]*stepsminusone + yzAv[i*nvars + l])*stepsinv;
    //        delyzAv[i*nvars + l] = yzAv[i*nvars + l] - yzAvMeans[i*nvars + l];
    //    }
    //}

    //for (auto l=0; l<nvars; ++l) {
    //    yzAvMeans_cross[l] = (yzAvMeans_cross[l]*stepsminusone + yzAv_cross[l])*stepsinv;
    //    delyzAv_cross[l] = yzAv_cross[l] - yzAvMeans_cross[l];
    //}

    /////////////////////////////////////////////////
    //// increment yz running average of spatial
    //// correlation of fluctuations of conserved
    //// variables, i.e., <delA(x)delB(x*)>, where
    //// x* is the cross_cell, and A and B are all
    //// combinations of two conserved variables
    /////////////////////////////////////////////////
    //
    //int flag;
    //for (auto i=0; i<n_cells[0]; ++i) {
    //    flag = 0;
    //    for (auto n=0; n<nvars; ++n) {
    //        for (auto m=0; m<nvars; ++m) {

    //            spatialCross[i*nvars*nvars + flag] = (spatialCross[i*nvars*nvars + flag]*stepsminusone + 
    //                                                         delyzAv[i*nvars + n]*delyzAv_cross[m])*stepsinv;
    //            flag += 1;
    //        }
    //    }
    //}
                
}
