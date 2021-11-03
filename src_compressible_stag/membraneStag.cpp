#include "compressible_functions_stag.H"
#include "compressible_functions.H"
#include "common_functions.H"
#include "rng_functions.H"
#include <math.h>

void doMembraneStag(MultiFab& cons, 
                    std::array< MultiFab, AMREX_SPACEDIM >& cumom,
                    MultiFab& prim, 
                    std::array< MultiFab, AMREX_SPACEDIM >& vel,
                    std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                    std::array< MultiFab, 2 >& edgeflux_x,
                    const amrex::Geometry geom, const amrex::Real dt)
{
    BL_PROFILE_VAR("doMembraneStag()",doMembraneStag);

    faceflux[0].setVal(0.0); // set mass & energy flux to zero
    edgeflux_x[0].setVal(0.0); // set flux of y-momentum to zero
    edgeflux_x[1].setVal(0.0); // set flux of z-momentum to zero
    
    doLangevin(cons,prim,faceflux,edgeflux_x,vel,geom,dt);

    faceflux[0].OverrideSync(geom.periodicity());
    edgeflux_x[0].OverrideSync(geom.periodicity());
    edgeflux_x[1].OverrideSync(geom.periodicity());

    applyEffusion(faceflux,edgeflux_x,cons,cumom);

    cons.FillBoundary(geom.periodicity()); // need correct periodicity for density to calculate velocities below
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        cumom[d].FillBoundary(geom.periodicity());
    }
    conservedToPrimitiveStag(prim, vel, cons, cumom);
    cons.FillBoundary(geom.periodicity()); // cell-centered momentum is filled in the line above 
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        vel[d].FillBoundary(geom.periodicity());
        cumom[d].FillBoundary(geom.periodicity());
    }
    prim.FillBoundary(geom.periodicity());
}

void doLangevin(MultiFab& cons_in, MultiFab& prim_in,
                std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                std::array< MultiFab, 2 >& edgeflux_x,
                std::array<MultiFab, AMREX_SPACEDIM>& vel,
                const amrex::Geometry geom,
                const amrex::Real dt)
{

    BL_PROFILE_VAR("doLangevin()",doLangevin);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    Real vol = dx[0]*dx[1]*dx[2];
    Real area = dx[1]*dx[2];
    
    GpuArray<Real,MAX_SPECIES> mass;
    GpuArray<Real,MAX_SPECIES> alpha;
    for (int l=0;l<nspecies;++l) {
        mass[l] = molmass[l]/(Runiv/k_B);
        alpha[l] = area*dt*transmission[l]*sqrt(k_B/(2.0*3.142*mass[l]));
    }
    
    for ( MFIter mfi(faceflux[0],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const Array4<Real>& prim = prim_in.array(mfi);
        const Array4<Real>& cons = cons_in.array(mfi);
        const Array4<Real>& xflux = faceflux[0].array(mfi);
        const Array4<Real>& edgex_v = edgeflux_x[0].array(mfi);
        const Array4<Real>& edgex_w = edgeflux_x[1].array(mfi);
        const Array4<Real>& vely = vel[1].array(mfi);
        const Array4<Real>& velz = vel[2].array(mfi);

//        amrex::Print() << bx.smallEnd(0) << " " << bx.bigEnd(0) << " " << membrane_cell << "\n";
        if ((bx.smallEnd(0) == membrane_cell) or (bx.bigEnd(0) == membrane_cell)) { // not really required

            amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
            {
                if (i == membrane_cell) {

                    Real TL = prim(i-1,j,k,4);
                    Real TR = prim(i,j,k,4); 
                    Real sqrtTL = sqrt(TL);
                    Real sqrtTR = sqrt(TR);
//                    amrex::Print() << i << " " << j << " " << k << " " << TL << " " << TR << "\n";

                    Real vyL = vely(i-1,j,k);
                    Real vyR = vely(i,j,k);
                    Real vzL = velz(i-1,j,k);
                    Real vzR = velz(i,j,k);

                    GpuArray<Real,MAX_SPECIES> rhoL;
                    GpuArray<Real,MAX_SPECIES> rhoR;

                    for (int l=0;l<nspecies;++l) {
                        rhoL[l] = cons(i-1,j,k,5+l);
                        rhoR[l] = cons(i,j,k,5+l);
                    }

                    Real EL, ER, pyL, pyR, pzL, pzR, meL, meR, pyeL, pyeR, pzeL, pzeR;

                    // Compute effusive flux per species
                    GpuArray<Real,4> Delmean; // Mean flux for [mass, py, pz, E]
                    Array2D<Real, 0, 4, 0, 4> Delvar; // covariance matrix
                    GpuArray<Real,4> rand; // Random normal numbers
                    GpuArray<Real,4> effflux; // effusive flux

                    for (int l=0;l<nspecies;++l) {

                        // Mean effusive fluxes
                        Delmean[0] = alpha[l]*sqrt(rhoL[l]*sqrtTL - rhoR[l]*sqrtTR); // mass

                        Delmean[1] = alpha[l]*(rhoL[l]*vyL*sqrtTL - rhoR[l]*vyR*sqrtTR); // py
                        Delmean[2] = alpha[l]*(rhoL[l]*vzL*sqrtTL - rhoR[l]*vzR*sqrtTR); // pz

                        EL = (2.0*k_B*TL/mass[l]) + 0.5*(vyL*vyL + vzL*vzL);
                        ER = (2.0*k_B*TR/mass[l]) + 0.5*(vyR*vyR + vzR*vzR);
                        Delmean[3] = alpha[l]*(rhoL[l]*sqrtTL*EL - rhoR[l]*sqrtTR*ER); // E
                        
                        // Set the covariance matrix to 0
                        for (int p=0;p<4;++p) {
                            for (int q=0;q<4;++q) {
                                Delvar(p,q) = 0.0;
                            }
                        }

                        // Fill variances
                        Delvar(0,0) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL + rhoR[l]*sqrtTR);
                        
                        pyL = k_B*TL + mass[l]*vyL*vyL;
                        pyR = k_B*TR + mass[l]*vyR*vyR;
                        Delvar(1,1) = alpha[l]*(rhoL[l]*sqrtTL*pyL + rhoR[l]*sqrtTR*pyR);

                        pzL = k_B*TL + mass[l]*vzL*vzL;
                        pzR = k_B*TR + mass[l]*vzR*vzR;
                        Delvar(2,2) = alpha[l]*(rhoL[l]*sqrtTL*pzL + rhoR[l]*sqrtTR*pzR);

                        EL = 24.0*k_B*k_B*TL*TL + 12.0*k_B*mass[l]*TL*(vyL*vyL+vzL*vzL) + mass[l]*mass[l]*(vyL*vyL+vzL*vzL)*(vyL*vyL+vzL*vzL);
                        ER = 24.0*k_B*k_B*TR*TR + 12.0*k_B*mass[l]*TR*(vyR*vyR+vzR*vzR) + mass[l]*mass[l]*(vyR*vyR+vzR*vzR)*(vyR*vyR+vzR*vzR);
                        Delvar(3,3) = alpha[l]*(rhoL[l]*sqrtTL*EL + rhoR[l]*sqrtTR*ER);
                        
                        // Fill covariances
                        Delvar(0,1) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL*vyL + rhoR[l]*sqrtTR*vyR);
                        Delvar(1,0) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL*vyL + rhoR[l]*sqrtTR*vyR);
                        Delvar(0,2) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL*vzL + rhoR[l]*sqrtTR*vzR);
                        Delvar(2,0) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL*vzL + rhoR[l]*sqrtTR*vzR);

                        meL = 4.0*k_B*TL + mass[l]*(vyL*vyL+vzL*vzL);
                        meR = 4.0*k_B*TR + mass[l]*(vyR*vyR+vzR*vzR);
                        Delvar(0,3) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*meL + rhoR[l]*sqrtTR*meR);
                        Delvar(3,0) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*meL + rhoR[l]*sqrtTR*meR);

                        pyeL = vyL*(6.0*k_B*TL + mass[l]*(vyL*vyL+vzL*vzL));
                        pyeR = vyR*(6.0*k_B*TR + mass[l]*(vyR*vyR+vzR*vzR));
                        Delvar(1,3) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*pyeL + rhoR[l]*sqrtTR*pyeR);
                        Delvar(3,1) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*pyeL + rhoR[l]*sqrtTR*pyeR);

                        pzeL = vzL*(6.0*k_B*TL + mass[l]*(vyL*vyL+vzL*vzL));
                        pzeR = vzR*(6.0*k_B*TR + mass[l]*(vyR*vyR+vzR*vzR));
                        Delvar(2,3) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*pzeL + rhoR[l]*sqrtTR*pzeR);
                        Delvar(3,2) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*pzeL + rhoR[l]*sqrtTR*pzeR);

                        // Cholesky factorise the covariance matrix
                        Array2D<Real, 0, 4, 0, 4> sqrtVar;
                        Cholesky(Delvar,4,sqrtVar);

                        // Generate random numbers
                        for (int n=0;n<4;++n) {
                            if (stoch_stress_form) rand[n] = RandomNormal(0.,1.,engine);
                            else rand[n] = 0.0;
                            effflux[n] = 0.0; // initialize effusive flux
                        }
                        
                        // Get effusive fluxes from Langevin integration
                        for (int p=0; p<4; ++p) {
                            for (int q=0; q<4; ++q) {
                                effflux[p] += Delmean[p] + sqrtVar(p,q)*rand[q]; // x = mean + rand*covar
                            }
                        }

                        // Increment total flux
                        xflux(i,j,k,0)  += effflux[0]/vol; // mass flux
                        xflux(i,j,k,5+l) = effflux[0]/vol; // species mass flux
                        xflux(i,j,k,4)  += effflux[3]/vol; // energy flux
                        if (do_1D) {
                        }
                        else if (do_2D) {
                            edgex_v(i,j,k)  += effflux[1]/vol; // y-momentum flux
                        }
                        else {
                            edgex_v(i,j,k)  += effflux[1]/vol; // y-momentum flux
                            edgex_w(i,j,k)  += effflux[2]/vol; // z-momentum flux
                        }

//                        amrex::Print() << i << " " << j << " " << k << " " << l << " " << effflux[0] << " " << effflux[1] << " " << effflux[2] << " " << effflux[3] << "\n";

                    }
                    
                }
            });

        }

    }

}

void applyEffusion(std::array<MultiFab, AMREX_SPACEDIM>& faceflux, 
                   std::array< MultiFab, 2 >& edgeflux_x,
                   MultiFab& cons_in,
                   std::array< MultiFab, AMREX_SPACEDIM >& cumom)
{
    
    BL_PROFILE_VAR("applyEffusion()",applyEffusion);
    
    for ( MFIter mfi(cons_in,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const Array4<Real>& cons = cons_in.array(mfi);
        const Array4<Real>& xflux = faceflux[0].array(mfi);
        const Array4<Real>& ymom = cumom[1].array(mfi);
        const Array4<Real>& zmom = cumom[2].array(mfi);
        const Array4<Real>& edgex_v = edgeflux_x[0].array(mfi);
        const Array4<Real>& edgex_w = edgeflux_x[1].array(mfi);

        if (bx.smallEnd(0) == membrane_cell) { // membrane at the left end (cell to the right of the membrane)

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == membrane_cell) {
                    cons(i,j,k,0) += xflux(i,j,k,0);
                    cons(i,j,k,4) += xflux(i,j,k,4);
                    for (int l=0;l<nspecies;++l) {
                        cons(i,j,k,5+l) += xflux(i,j,k,5+l);
                    }
                    ymom(i,j,k) += edgex_v(i,j,k); 
                    zmom(i,j,k) += edgex_w(i,j,k); 
                }
            });
        }

        else if (bx.bigEnd(0) == membrane_cell - 1) { // membrane at the right end (cell to the left of the membrane)

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == membrane_cell-1) {
                    cons(i,j,k,0) -= xflux(membrane_cell,j,k,0);
                    cons(i,j,k,4) -= xflux(membrane_cell,j,k,4);
                    for (int l=0;l<nspecies;++l) {
                        cons(i,j,k,5+l) -= xflux(membrane_cell,j,k,5+l);
                    }
                    ymom(i,j,k) -= edgex_v(membrane_cell,j,k); 
                    zmom(i,j,k) -= edgex_w(membrane_cell,j,k); 
                }
            });
        
        }

    }

}
