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
                       Vector<Real>& cuyzAvMeans, Vector<Real>& cuyzAvMeans_cross,
                       Vector<Real>& spatialCross,
                       const int steps, const amrex::Real* dx)
{
    BL_PROFILE_VAR("evaluateStatsStag()",evaluateStatsStag);
    
    double stepsminusone = steps - 1.;
    double stepsinv = 1./steps;

    GpuArray<Real,MAX_SPECIES> fracvec;
    //GpuArray<Real,MAX_SPECIES> massvec;

    //////////////////
    // evaluate_means
    //////////////////
    
    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.validbox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);
        const auto lox = amrex::lbound(tbx);
        const auto hix = amrex::ubound(tbx);
        const auto loy = amrex::lbound(tby);
        const auto hiy = amrex::ubound(tby);
        const auto loz = amrex::lbound(tbz);
        const auto hiz = amrex::ubound(tbz);

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

        // x-face (mean jx)
        // on host, not gpu
        for (auto k = lox.z; k <= hix.z; ++k) {
        for (auto j = lox.y; j <= hix.y; ++j) {
        for (auto i = lox.x; i <= hix.x; ++i) {
            momxmeans(i,j,k) = (momxmeans(i,j,k)*stepsminusone + momx(i,j,k))*stepsinv;
        }
        }
        }

        // y-face (mean jy)
        // on host, not gpu
        for (auto k = loy.z; k <= hiy.z; ++k) {
        for (auto j = loy.y; j <= hiy.y; ++j) {
        for (auto i = loy.x; i <= hiy.x; ++i) {
            momymeans(i,j,k) = (momymeans(i,j,k)*stepsminusone + momy(i,j,k))*stepsinv;
        }
        }
        }

        // z-face (mean jz)
        // on host, not gpu
        for (auto k = loz.z; k <= hiz.z; ++k) {
        for (auto j = loz.y; j <= hiz.y; ++j) {
        for (auto i = loz.x; i <= hiz.x; ++i) {
            momzmeans(i,j,k) = (momzmeans(i,j,k)*stepsminusone + momz(i,j,k))*stepsinv;
        }
        }
        }

        // cell centers (mean rho, rhoE, rhoYk, Yk, P, T) 
        // also jx, velx, jy, vely, jz, velz on CC
        // on host, not gpu
        for (auto k = lo.z-ngc[2]; k <= hi.z+ngc[2]; ++k) {
        for (auto j = lo.y-ngc[1]; j <= hi.y+ngc[1]; ++j) {
        for (auto i = lo.x-ngc[0]; i <= hi.x+ngc[0]; ++i) {
            cumeans(i,j,k,0) = (cumeans(i,j,k,0)*stepsminusone + cu(i,j,k,0))*stepsinv; //rhomeans
        }
        }
        }

        for (auto k = lo.z; k <= hi.z; ++k) {
        for (auto j = lo.y; j <= hi.y; ++j) {
        for (auto i = lo.x; i <= hi.x; ++i) {
            
        //    cumeans(i,j,k,0) = (cumeans(i,j,k,0)*stepsminusone + cu(i,j,k,0))*stepsinv; //rhomeans
            cumeans(i,j,k,1) = 0.5*(momxmeans(i,j,k) + momxmeans(i+1,j,k)); // jxmeans on CC
            cumeans(i,j,k,2) = 0.5*(momymeans(i,j,k) + momymeans(i,j+1,k)); // jymeans on CC
            cumeans(i,j,k,3) = 0.5*(momzmeans(i,j,k) + momzmeans(i,j,k+1)); // jzmeans on CC
            cumeans(i,j,k,4) = (cumeans(i,j,k,4)*stepsminusone + cu(i,j,k,4))*stepsinv; //rhoEmeans

            for (int l=5; l<nvars; ++l) {
                cumeans(i,j,k,l) = (cumeans(i,j,k,l)*stepsminusone + cu(i,j,k,l))*stepsinv; //rhoYkmeans
                fracvec[l-5] = cumeans(i,j,k,l)/cumeans(i,j,k,0); // Ykmeans
                //massvec[l-5] = cumeans(i,j,k,l); //rhoYk
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
        }
        }
        }

        // x-face (mean velx)
        // on host, not gpu
        for (auto k = lox.z; k <= hix.z; ++k) {
        for (auto j = lox.y; j <= hix.y; ++j) {
        for (auto i = lox.x; i <= hix.x; ++i) {
            Real densitymeaninv = 2.0/(cumeans(i-1,j,k,0)+cumeans(i,j,k,0));
            velxmeans(i,j,k) = momxmeans(i,j,k)*densitymeaninv;
        }
        }
        }

        // y-face (mean vely)
        // on host, not gpu
        for (auto k = loy.z; k <= hiy.z; ++k) {
        for (auto j = loy.y; j <= hiy.y; ++j) {
        for (auto i = loy.x; i <= hiy.x; ++i) {
            Real densitymeaninv = 2.0/(cumeans(i,j-1,k,0)+cumeans(i,j,k,0));
            velymeans(i,j,k) = momymeans(i,j,k)*densitymeaninv;
        }
        }
        }

        // z-face (mean velz)
        // on host, not gpu
        for (auto k = loz.z; k <= hiz.z; ++k) {
        for (auto j = loz.y; j <= hiz.y; ++j) {
        for (auto i = loz.x; i <= hiz.x; ++i) {
            Real densitymeaninv = 2.0/(cumeans(i,j,k-1,0)+cumeans(i,j,k,0));
            velzmeans(i,j,k) = momzmeans(i,j,k)*densitymeaninv;
        }
        }
        }
    }

    /////////////////////////////////////////////
    // evaluate yz average of conserved variables
    /////////////////////////////////////////////

    Vector<Real> cuyzAv(n_cells[0]*nvars,0.); // yz-average of conserved variable at current snapshot
    Vector<Real> cuyzAv_cross(nvars,0.); // yz-averaged of conserved variable at current snapshot at the cross cell

    for (MFIter mfi(cons, false); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<const Real> cu = cons.array(mfi);

        for (auto i = lo.x; i <= hi.x; ++i) {
        for (auto l = 0;    l < nvars; ++l) {
        for (auto k = lo.z; k <= hi.z; ++k) {
        for (auto j = lo.y; j <= hi.y; ++j) {
            cuyzAv[i*nvars + l] += cu(i,j,k,l); 
        }
        }
        }
        }

    } // end MFiter

    // sum over all processors
    ParallelDescriptor::ReduceRealSum(cuyzAv.dataPtr(),n_cells[0]*nvars);

    // divide by the number of yz-face cells
    int n_face_cells = n_cells[1]*n_cells[2];
    for (auto i=0; i<n_cells[0]; ++i) {
        for (auto l=0; l<nvars; ++l) {
            cuyzAv[i*nvars + l] /= n_face_cells;
        }
    }

    // copy the cross_cell yz-averaged value in the vector
    for (auto i=0; i<n_cells[0]; ++i) {
        if (i == cross_cell) {
            for (auto l=0; l<nvars; ++l) {
                cuyzAv_cross[l] = cuyzAv[i*nvars + l];
            }
        }
    }

    ///////////////////////////////////////////////
    // increment yz average of conserved variables
    // to the running mean and compute fluctuation
    ///////////////////////////////////////////////
    
    Vector<Real> delcuyzAv(n_cells[0]*nvars,0.); // yz-average of conserved variable at current snapshot 
    Vector<Real> delcuyzAv_cross(nvars,0.); // yz-averaged of conserved variable at current snapshot at the cross cell

    for (auto i=0; i<n_cells[0]; ++i) {
        for (auto l=0; l<nvars; ++l) {
            cuyzAvMeans[i*nvars + l] = (cuyzAvMeans[i*nvars + l]*stepsminusone + cuyzAv[i*nvars + l])*stepsinv;
            delcuyzAv[i*nvars + l] = cuyzAv[i*nvars + l] - cuyzAvMeans[i*nvars + l];
        }
    }

    for (auto l=0; l<nvars; ++l) {
        cuyzAvMeans_cross[l] = (cuyzAvMeans_cross[l]*stepsminusone + cuyzAv_cross[l])*stepsinv;
        delcuyzAv_cross[l] = cuyzAv_cross[l] - cuyzAvMeans_cross[l];
    }

    ///////////////////////////////////////////////
    // increment yz running average of spatial
    // correlation of fluctuations of conserved
    // variables, i.e., <delA(x)delB(x*)>, where
    // x* is the cross_cell, and A and B are all
    // combinations of two conserved variables
    ///////////////////////////////////////////////
    
    int flag;
    for (auto i=0; i<n_cells[0]; ++i) {
        flag = 0;
        for (auto n=0; n<nvars; ++n) {
            for (auto m=0; m<nvars; ++m) {

                spatialCross[i*nvars*nvars + flag] = (spatialCross[i*nvars*nvars + flag]*stepsminusone + 
                                                             delcuyzAv[i*nvars + n]*delcuyzAv_cross[m])*stepsinv;
                flag += 1;
            }
        }
    }
                
    
    /////////////////////////////////////
    // evaluate variances and covariances
    /////////////////////////////////////
    
    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.validbox();
        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);
        const auto lox = amrex::lbound(tbx);
        const auto hix = amrex::ubound(tbx);
        const auto loy = amrex::lbound(tby);
        const auto hiy = amrex::ubound(tby);
        const auto loz = amrex::lbound(tbz);
        const auto hiz = amrex::ubound(tbz);

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

        // x-face: jx and velx
        // on host, not gpu
        for (auto k = lox.z; k <= hix.z; ++k) {
        for (auto j = lox.y; j <= hix.y; ++j) {
        for (auto i = lox.x; i <= hix.x; ++i) {
            
            Real deljx = momx(i,j,k) - momxmeans(i,j,k);
            momxvars(i,j,k) = (momxvars(i,j,k)*stepsminusone + deljx*deljx)*stepsinv; // <jx jx>

            Real densitymeaninv = 2.0/(cumeans(i-1,j,k,0)+cumeans(i,j,k,0));
            Real delrho = 0.5*(cu(i-1,j,k,0) + cu(i,j,k,0)) - 0.5*(cumeans(i-1,j,k,0) + cumeans(i,j,k,0));
            Real delvelx = (deljx - velxmeans(i,j,k)*delrho)*densitymeaninv;
            velxvars(i,j,k) = (velxvars(i,j,k)*stepsminusone + delvelx*delvelx)*stepsinv; // <vx vx>
        }
        }
        }

        // y-face: jy and vely
        // on host, not gpu
        for (auto k = loy.z; k <= hiy.z; ++k) {
        for (auto j = loy.y; j <= hiy.y; ++j) {
        for (auto i = loy.x; i <= hiy.x; ++i) {
            
            Real deljy = momy(i,j,k) - momymeans(i,j,k);
            momyvars(i,j,k) = (momyvars(i,j,k)*stepsminusone + deljy*deljy)*stepsinv; // <jy jy>

            Real densitymeaninv = 2.0/(cumeans(i,j-1,k,0)+cumeans(i,j,k,0));
            Real delrho = 0.5*(cu(i,j-1,k,0) + cu(i,j,k,0)) - 0.5*(cumeans(i,j-1,k,0) + cumeans(i,j,k,0));
            Real delvely = (deljy - velymeans(i,j,k)*delrho)*densitymeaninv;
            velyvars(i,j,k) = (velyvars(i,j,k)*stepsminusone + delvely*delvely)*stepsinv; // <vx vx>
            std::cout << i << " " << j << " " << k << " " <<  stepsinv << " " << velyvars(i,j,k) << "\n";
        }
        }
        }

        // z-face: jz and velz
        // on host, not gpu
        for (auto k = loz.z; k <= hiz.z; ++k) {
        for (auto j = loz.y; j <= hiz.y; ++j) {
        for (auto i = loz.x; i <= hiz.x; ++i) {
            
            Real deljz = momz(i,j,k) - momzmeans(i,j,k);
            momzvars(i,j,k) = (momzvars(i,j,k)*stepsminusone + deljz*deljz)*stepsinv; // <jz jz>

            Real densitymeaninv = 2.0/(cumeans(i,j,k-1,0)+cumeans(i,j,k,0));
            Real delrho = 0.5*(cu(i,j,k-1,0) + cu(i,j,k,0)) - 0.5*(cumeans(i,j,k-1,0) + cumeans(i,j,k,0));
            Real delvelz = (deljz - velzmeans(i,j,k)*delrho)*densitymeaninv;
            velzvars(i,j,k) = (velzvars(i,j,k)*stepsminusone + delvelz*delvelz)*stepsinv; // <vz vz>
        }
        }
        }

        // cell centers (temperature fluctuation and variance)
        // on host, not gpu
        for (auto k = lo.z; k <= hi.z; ++k) {
        for (auto j = lo.y; j <= hi.y; ++j) {
        for (auto i = lo.x; i <= hi.x; ++i) {

            // conserved variable variances (rho, rhoE, rhoYk)
            Real delrho = cu(i,j,k,0) - cumeans(i,j,k,0);
            cuvars(i,j,k,0) = (cuvars(i,j,k,0)*stepsminusone + delrho*delrho)*stepsinv; // <rho rho>
            
            Real delenergy = cu(i,j,k,4) - cumeans(i,j,k,4);
            cuvars(i,j,k,4) = (cuvars(i,j,k,4)*stepsminusone + delenergy*delenergy)*stepsinv; // <rhoE rhoE>
            
            for (int l=5; l<nvars; ++l) {
                Real delmassvec = cu(i,j,k,l) - cumeans(i,j,k,l);
                cuvars(i,j,k,l) = (cuvars(i,j,k,l)*stepsminusone + delmassvec*delmassvec)*stepsinv; // <rhoYk rhoYk>
            }
            Real delrho1 = cu(i,j,k,5) - cumeans(i,j,k,5); 
            Real delrho4 = cu(i,j,k,nvars-1) - cumeans(i,j,k,nvars-1);
            
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
            covars(i,j,k,10) = (covars(i,j,k,10)*stepsminusone + delrho1*delrho4)*stepsinv; // <rhoYk1 rhoYk4>
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
        }
        }
        }
    } // end MFIter
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
