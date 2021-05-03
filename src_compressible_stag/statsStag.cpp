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
        const Array4<const Real> prim      = prim_in.array(mfi);
        const Array4<      Real> primmeans = primMean.array(mfi);

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
        const Array4<const Real> prim      = prim_in.array(mfi);
        const Array4<      Real> primmeans = primMean.array(mfi);
        const Array4<      Real> primvars  = primVar.array(mfi);

        const Array4<      Real> covars    = coVar.array(mfi);

        const Array4<const Real> velx      = vel[0].array(mfi);
        const Array4<const Real> vely      = vel[1].array(mfi);
        const Array4<const Real> velz      = vel[2].array(mfi);
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

            // del rhoYk
            GpuArray<Real,MAX_SPECIES> delrhoYk;
            for (int ns=0; ns<nspecies; ++ns) {
                delrhoYk[ns] = cu(i,j,k,5+ns) - cumeans(i,j,k,5+ns);
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

        });

    } // end MFIter

    /////////////////////////////////////////////////////////////
    // evaluate slice average at the cross cell and at every cell
    /////////////////////////////////////////////////////////////

    // contains yz-averaged running & instantaneous averages of conserved variables at every x + four primitive variables [vx, vy, vz, T]: 2*nvars + 2*4
    Vector<Real>  yzAvMeans(n_cells[0]*(2*nvars+8), 0.0); // yz-average at all x

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

        for (int i=0; i<(2*nvars+8); ++i) {
            yzAvMeans_cross[i] = 0.;
        }
        
        int counter = 0; // number of y * number of z cells
        
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
                yzAvMeans_cross[17] += primmeans(i,j,k,0);                          // T-mean
                for (int ns=0; ns<nspecies; ++ns) {
                    yzAvMeans_cross[18+2*ns+0]   += cu(i,j,k,5+ns);                 // rhoYk-instant
                    yzAvMeans_cross[18+2*ns+1]   += cumeans(i,j,k,5+ns);            // rhoYk-mean
                }
                counter = counter + 1;
            }
            yzAvMeans[i*(2*nvars+8)+0]  += cu(i,j,k,0);                                 // rho-instant
            yzAvMeans[i*(2*nvars+8)+1]  += cumeans(i,j,k,0);                            // rho-mean
            yzAvMeans[i*(2*nvars+8)+2]  += cu(i,j,k,4);                                 // energy-instant
            yzAvMeans[i*(2*nvars+8)+3]  += cumeans(i,j,k,4);                            // energy-mean
            yzAvMeans[i*(2*nvars+8)+4]  += 0.5*(momx(i,j,k) + momx(i+1,j,k));           // jx-instant
            yzAvMeans[i*(2*nvars+8)+5]  += 0.5*(momxmeans(i,j,k) + momxmeans(i+1,j,k)); // jx-mean
            yzAvMeans[i*(2*nvars+8)+6]  += 0.5*(momy(i,j,k) + momy(i,j+1,k));           // jy-instant
            yzAvMeans[i*(2*nvars+8)+7]  += 0.5*(momymeans(i,j,k) + momymeans(i,j+1,k)); // jy-mean
            yzAvMeans[i*(2*nvars+8)+8]  += 0.5*(momz(i,j,k) + momz(i,j,k+1));           // jz-instant
            yzAvMeans[i*(2*nvars+8)+9]  += 0.5*(momzmeans(i,j,k) + momzmeans(i,j,k+1)); // jz-mean
            yzAvMeans[i*(2*nvars+8)+10] += 0.5*(velx(i,j,k) + velx(i+1,j,k));           // velx-instant
            yzAvMeans[i*(2*nvars+8)+11] += 0.5*(velxmeans(i,j,k) + velxmeans(i+1,j,k)); // velx-mean
            yzAvMeans[i*(2*nvars+8)+12] += 0.5*(vely(i,j,k) + vely(i,j+1,k));           // vely-instant
            yzAvMeans[i*(2*nvars+8)+13] += 0.5*(velymeans(i,j,k) + velymeans(i,j+1,k)); // vely-mean
            yzAvMeans[i*(2*nvars+8)+14] += 0.5*(velz(i,j,k) + velz(i,j,k+1));           // velz-instant
            yzAvMeans[i*(2*nvars+8)+15] += 0.5*(velzmeans(i,j,k) + velzmeans(i,j,k+1)); // velz-mean
            yzAvMeans[i*(2*nvars+8)+16] += prim(i,j,k,4);                               // T-instant
            yzAvMeans[i*(2*nvars+8)+17] += primmeans(i,j,k,0);                          // T-mean
            for (int ns=0; ns<nspecies; ++ns) {
                yzAvMeans[i*(2*nvars+8)+18+2*ns+0]   += cu(i,j,k,5+ns);                 // rhoYk-instant
                yzAvMeans[i*(2*nvars+8)+18+2*ns+1]   += cumeans(i,j,k,5+ns);            // rhoYk-mean
            }
        }
        }
        }

        for (int i=0; i<n_cells[0]*(2*nvars+8); ++i) {
            yzAvMeans[i] /= counter;
        }
        for (int i=0; i<(2*nvars+8); ++i) {
            yzAvMeans_cross[i] /= counter;
        }

    } // end MFITer

    // parallel reduce sum yzAvMeans and yzAvMeans_cross
    ParallelDescriptor::ReduceRealSum(yzAvMeans.dataPtr(),n_cells[0]*(2*nvars+8));
    ParallelDescriptor::ReduceRealSum(yzAvMeans_cross.dataPtr(),2*nvars+8);

    /////////////////////////////////////////////////////////////
    // evaluate x-spatial correlations
    /////////////////////////////////////////////////////////////
    int ncross = 11+nspecies;
    for (int i=0; i<n_cells[0]; ++i) {

        spatialCross[i*ncross+0] = (spatialCross[i*ncross+0]*stepsminusone + yzAvMeans_cross[16]*yzAvMeans[i*(2*nvars+8)+16])*stepsinv; // <T(x*)T(x)>
        spatialCross[i*ncross+1] = spatialCross[i*ncross+0] - yzAvMeans_cross[17]*yzAvMeans[i*(2*nvars+8)+17]; // <T(x*)T(x)> - <T(x*)><T(x)>
        spatialCross[i*ncross+2] = (spatialCross[i*ncross+2]*stepsminusone + yzAvMeans_cross[16]*yzAvMeans[i*(2*nvars+8)+0])*stepsinv; // <T(x*)rho(x)>
        spatialCross[i*ncross+3] = spatialCross[i*ncross+2] - yzAvMeans_cross[17]*yzAvMeans[i*(2*nvars+8)+1]; // <T(x*)rho(x)> - <T(x*)><rho(x)>
        
        Real delrhostar = yzAvMeans_cross[0] - yzAvMeans_cross[1]; // <rho(x*) - <rho(x*)>>
        Real delrho     = yzAvMeans[i*(2*nvars+8)+0] - yzAvMeans[i*(2*nvars+8)+1]; // <rho(x) - <rho(x)>>
        spatialCross[i*ncross+4] = (spatialCross[i*ncross+4]*stepsminusone + delrhostar*delrho)*stepsinv; // <delrho(x*) delrho(x)>
        spatialCross[i*ncross+5] = (spatialCross[i*ncross+5]*stepsminusone + yzAvMeans_cross[6]*yzAvMeans[i*(2*nvars+8)+0])*stepsinv; // <jx(x*)rho(x)>
        
        Real deljxdelrho = spatialCross[i*ncross+5] - yzAvMeans_cross[5]*yzAvMeans[i*(2*nvars+8)+1]; // <jx(x*)rho(x)> - <jx(x*)><rho(x)> = <deljx(x*)delrho(x)>
        // <delu(x*)delrho> = (<deljx(x*)delrho(x)> - <u(x*)><<delrho(x*) delrho(x)>)/<rho(x*)> -- see Garcia 2007
        spatialCross[i*ncross+6] = (deljxdelrho - yzAvMeans_cross[11]*spatialCross[i*ncross+4])/yzAvMeans_cross[1];  
        
        Real delrhoEstar = yzAvMeans_cross[2] - yzAvMeans_cross[3]; // <rhoE(x*) - <rhoE(x*)>>
        Real delrhoE     = yzAvMeans[i*(2*nvars+8)+2] - yzAvMeans[i*(2*nvars+8)+3]; // <rhoE(x) - <rhoE(x)>>
        spatialCross[i*ncross+7] = (spatialCross[i*ncross+7]*stepsminusone + delrhoEstar*delrhoE)*stepsinv; // <delrhoE(x*) delrhoE(x)>
        
        Real deljxstar = yzAvMeans_cross[4] - yzAvMeans_cross[5]; // <jx(x*) - <jx(x*)>>
        Real deljx     = yzAvMeans[i*(2*nvars+8)+4] - yzAvMeans[i*(2*nvars+8)+5]; // <jx(x) - <jx(x)>>
        spatialCross[i*ncross+8] = (spatialCross[i*ncross+8]*stepsminusone + deljxstar*deljx)*stepsinv; // <deljx(x*) deljx(x)>
        
        Real deljystar = yzAvMeans_cross[6] - yzAvMeans_cross[7]; // <jy(x*) - <jy(x*)>>
        Real deljy     = yzAvMeans[i*(2*nvars+8)+6] - yzAvMeans[i*(2*nvars+8)+7]; // <jy(x) - <jy(x)>>
        spatialCross[i*ncross+9] = (spatialCross[i*ncross+9]*stepsminusone + deljystar*deljy)*stepsinv; // <deljy(x*) deljy(x)>
        
        Real deljzstar = yzAvMeans_cross[8] - yzAvMeans_cross[9]; // <jz(x*) - <jz(x*)>>
        Real deljz     = yzAvMeans[i*(2*nvars+8)+8] - yzAvMeans[i*(2*nvars+8)+9]; // <jz(x) - <jz(x)>>
        spatialCross[i*ncross+10] = (spatialCross[i*ncross+10]*stepsminusone + deljzstar*deljz)*stepsinv; // <deljz(x*) deljz(x)>

        for (int ns=0; ns<nspecies; ++ns) {
            Real delrhoykstar = yzAvMeans_cross[18+2*ns+0] - yzAvMeans_cross[18+2*ns+1]; // <rhoyk(x*) - <rhoyk(x*)>>
            Real delrhoyk     = yzAvMeans[i*(2*nvars+8)+18+2*ns+0] - yzAvMeans[i*(2*nvars+8)+18+2*ns+1]; // <rhoyk(x) - <rhoyk(x)>>
            spatialCross[i*ncross+11+ns] = (spatialCross[i*ncross+11+ns]*stepsminusone + delrhoykstar*delrhoyk)*stepsinv; // <delrhoyk(x*) delrhoyk(x)>
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

//void yzAverage(const MultiFab& consMean,
//               const MultiFab& consVar,
//               const MultiFab& primMean,
//               const MultiFab& primVar,
//               const MultiFab& spatialCross,
//               MultiFab& consMeanAv,
//               MultiFab& consVarAv,
//               MultiFab& primMeanAv,
//               MultiFab& primVarAv,
//               MultiFab& spatialCrossAv)
//{
//    BL_PROFILE_VAR("yzAverage()",yzAverage);
//
//    // Loop over boxes
//
//    int six = 6; //Look up how to avoid this later?
//    int primVarsPlusFive = nprimvars + 5;
//
//    for ( MFIter mfi(consMean); mfi.isValid(); ++mfi)
//    {
//        const Box& bx = mfi.validbox();
//
//        multifab_yzav(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
//                       consMean[mfi].dataPtr(),
//                       consMeanAv[mfi].dataPtr(), &nvars);
//
//        multifab_yzav(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
//                       consVar[mfi].dataPtr(),
//                       consVarAv[mfi].dataPtr(), &nvars);
//
//        multifab_yzav(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
//                       primMean[mfi].dataPtr(),
//                       primMeanAv[mfi].dataPtr(), &nprimvars);
//
//        multifab_yzav(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
//                       primVar[mfi].dataPtr(),
//                       primVarAv[mfi].dataPtr(), &nprimvars);
//
//        multifab_yzav(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
//                       spatialCross[mfi].dataPtr(),
//                       spatialCrossAv[mfi].dataPtr(), &six);
//
//    }
//
//}
