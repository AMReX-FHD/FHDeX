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
    
    for ( MFIter mfi(faceflux[0]); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const Array4<Real>& prim = prim_in.array(mfi);
        const Array4<Real>& cons = cons_in.array(mfi);
        const Array4<Real>& xflux = faceflux[0].array(mfi);
        const Array4<Real>& edgex_v = edgeflux_x[0].array(mfi);
        const Array4<Real>& edgex_w = edgeflux_x[1].array(mfi);
        const Array4<Real>& vely = vel[1].array(mfi);
        const Array4<Real>& velz = vel[2].array(mfi);

        if (bx.smallEnd(0) == membrane_cell && bx.bigEnd(0) == membrane_cell) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == membrane_cell) {

                    Real TL = prim(i-1,j,k,4);
                    Real TR = prim(i,j,k,4); 
                    Real sqrtTL = sqrt(TL);
                    Real sqrtTR = sqrt(TR);

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

                    GpuArray<Real,MAX_SPECIES+3> Delmean; // Mean fluxes for [species, py, pz, E]
                    Delmean[nspecies+0] = 0.0; // py
                    Delmean[nspecies+1] = 0.0; // pz
                    Delmean[nspecies+2] = 0.0; // E
                    Real DelmeanM = 0.0; // total mass

                    for (int l=0;l<nspecies;++l) {
                        Delmean[l] = alpha[l]*sqrt(rhoL[l]*sqrtTL - rhoR[l]*sqrtTR); // mass
                        DelmeanM += Delmean[l];

                        Delmean[nspecies+0] += alpha[l]*(rhoL[l]*vyL*sqrtTL - rhoR[l]*vyR*sqrtTR); // py
                        Delmean[nspecies+1] += alpha[l]*(rhoL[l]*vzL*sqrtTL - rhoR[l]*vzR*sqrtTR); // pz

                        Real EL  = (2.0*k_B*TL/mass[l]) + 0.5*(vyL*vyL + vzL*vzL);
                        Real ER = (2.0*k_B*TR/mass[l]) + 0.5*(vyR*vyR + vzR*vzR);
                        Delmean[nspecies+2] += alpha[l]*(rhoL[l]*sqrtTL*EL - rhoR[l]*sqrtTR*ER); // E
                    }

                    // Construct the covariance matrix
                    Array2D<Real, 1, MAX_SPECIES+3, 1, MAX_SPECIES+3> Delvar;
                    for (int l=0;l<nspecies+3;++l) {
                        for (int n=0;n<nspecies+3;++n) {
                            Delvar(l,n) = 0.0;
                        }
                    }

                    // Fill variances
                    for (int l=0;l<nspecies;++l) {
                        for (int n=0;n<nspecies;++n) {
                            if (l==n) Delvar(l,n) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL + rhoR[l]*sqrtTR);
                            else Delvar(l,n) = 0.0;
                        }
                        
                        Real pyL = k_B*TL + mass[l]*vyL*vyL;
                        Real pyR = k_B*TR + mass[l]*vyR*vyR;
                        Delvar(nspecies+0,nspecies+0) += alpha[l]*(rhoL[l]*sqrtTL*pyL + rhoR[l]*sqrtTR*pyR);
                        Real pzL = k_B*TL + mass[l]*vzL*vzL;
                        Real pzR = k_B*TR + mass[l]*vzR*vzR;
                        Delvar(nspecies+1,nspecies+1) += alpha[l]*(rhoL[l]*sqrtTL*pzL + rhoR[l]*sqrtTR*pzR);

                        Real EL = 24.0*k_B*k_B*TL*TL + 12.0*k_B*mass[l]*TL*(vyL*vyL+vzL*vzL) + mass[l]*mass[l]*(vyL*vyL+vzL*vzL)*(vyL*vyL+vzL*vzL);
                        Real ER = 24.0*k_B*k_B*TR*TR + 12.0*k_B*mass[l]*TR*(vyR*vyR+vzR*vzR) + mass[l]*mass[l]*(vyR*vyR+vzR*vzR)*(vyR*vyR+vzR*vzR);
                        Delvar(nspecies+2,nspecies+2) += alpha[l]*(rhoL[l]*sqrtTL*EL + rhoR[l]*sqrtTR*ER);
                    }

                    // Fill covariances
                    for (int l=0;l<nspecies;++l) {

                        Delvar(l,nspecies+0) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL*vyL + rhoR[l]*sqrtTR*vyR);
                        Delvar(nspecies+0,l) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL*vyL + rhoR[l]*sqrtTR*vyR);
                        Delvar(l,nspecies+1) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL*vzL + rhoR[l]*sqrtTR*vzR);
                        Delvar(nspecies+1,l) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL*vzL + rhoR[l]*sqrtTR*vzR);

                        Real meL = 4.0*k_B*TL + mass[l]*(vyL*vyL+vzL*vzL);
                        Real meR = 4.0*k_B*TR + mass[l]*(vyR*vyR+vzR*vzR);
                        Delvar(l,nspecies+2) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*meL + rhoR[l]*sqrtTR*meR);
                        Delvar(nspecies+2,l) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*meL + rhoR[l]*sqrtTR*meR);

                        Real pyeL = vyL*(6.0*k_B*TL + mass[l]*(vyL*vyL+vzL*vzL));
                        Real pyeR = vyR*(6.0*k_B*TR + mass[l]*(vyR*vyR+vzR*vzR));
                        Delvar(nspecies+2,nspecies+0) += 0.5*alpha[l]*(rhoL[l]*sqrtTL*pyeL + rhoR[l]*sqrtTR*pyeR);
                        Delvar(nspecies+0,nspecies+2) += 0.5*alpha[l]*(rhoL[l]*sqrtTL*pyeL + rhoR[l]*sqrtTR*pyeR);
                        Real pzeL = vzL*(6.0*k_B*TL + mass[l]*(vyL*vyL+vzL*vzL));
                        Real pzeR = vzR*(6.0*k_B*TR + mass[l]*(vyR*vyR+vzR*vzR));
                        Delvar(nspecies+2,nspecies+1) += 0.5*alpha[l]*(rhoL[l]*sqrtTL*pzeL + rhoR[l]*sqrtTR*pzeR);
                        Delvar(nspecies+1,nspecies+2) += 0.5*alpha[l]*(rhoL[l]*sqrtTL*pzeL + rhoR[l]*sqrtTR*pzeR);
                    }
                    
                    // Cholesky factorise the covariance matrix
                    Array2D<Real, 1, MAX_SPECIES+3, 1, MAX_SPECIES+3> sqrtVar;
                    Cholesky(Delvar,nspecies,sqrtVar);

                    // Generate random numbers 
                    GpuArray<Real,MAX_SPECIES+3> rand; // Random normal numbers
                    GpuArray<Real,MAX_SPECIES+3> effflux; // effusive flux
                    for (int n=0;n<nspecies+3;++n) {
                        if (stoch_stress_form == 1) rand[n] = amrex::RandomNormal(0.,1.);
                        else rand[n] = 0.0;
                        effflux[n] = 0.0; // initialize effusive flux
                    }

                    // Get effusive fluxes from Langevin integration
                    for (int ns=0; ns<nspecies+3; ++ns) {
                        for (int ll=0; ll<=ns; ++ll) {
                            effflux[ns] += Delmean[ns] + sqrtVar(ns,ll)*rand[ll]; // x = mean + rand*covar
                        }
                    }

                    // Fill up the fluxes
                    for (int ns=0;ns<nspecies;++ns) { // species and mass flux
                        xflux(i,j,k,ns+5) = effflux[ns]/vol;
                        xflux(i,j,k,0) += effflux[ns]/vol;
                    }
                    xflux(i,j,k,4) = effflux[nspecies+2]/vol; // energy flux
                    edgex_v(i,j,k) = effflux[nspecies+0]/vol; // y-momentum flux
                    edgex_w(i,j,k) = effflux[nspecies+1]/vol; // z-momentum flux
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
    
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {

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
