#include "compressible_functions_stag.H"
#include "compressible_functions.H"
#include "common_functions.H"
#include "rng_functions.H"
#include "reservoirStag_K.H"
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

    AMREX_D_TERM(faceflux[0].setVal(0.0);,
                 faceflux[1].setVal(0.0);,
                 faceflux[2].setVal(0.0););
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
        mass[l] = molmass[l]/(6.023e23);
        alpha[l] = area*transmission[l]*sqrt(k_B/(2.0*3.142*mass[l]));
    }
    
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<Real>& prim = prim_in.array(mfi);
        const Array4<Real>& cons = cons_in.array(mfi);
        const Array4<Real>& xflux = faceflux[0].array(mfi);
        const Array4<Real>& edgex_v = edgeflux_x[0].array(mfi);
        const Array4<Real>& edgex_w = edgeflux_x[1].array(mfi);
        const Array4<Real>& vely = vel[1].array(mfi);
        const Array4<Real>& velz = vel[2].array(mfi);

        if ((bx.smallEnd(0) == membrane_cell) or (bx.bigEnd(0) == membrane_cell - 1)) {

            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {
//            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//            amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, RandomEngine const& engine) noexcept
//            {
//                if (i == membrane_cell) {

                    Real TL = prim(membrane_cell-1,j,k,4);
                    Real TR = prim(membrane_cell,j,k,4); 
                    Real sqrtTL = sqrt(TL);
                    Real sqrtTR = sqrt(TR);

                    Real vyL = vely(membrane_cell-1,j,k);
                    Real vyR = vely(membrane_cell,j,k);
                    Real vzL = velz(membrane_cell-1,j,k);
                    Real vzR = velz(membrane_cell,j,k);

                    GpuArray<Real,MAX_SPECIES> rhoL;
                    GpuArray<Real,MAX_SPECIES> rhoR;
                    for (int l=0;l<nspecies;++l) {
                        rhoL[l] = cons(membrane_cell-1,j,k,5+l);
                        rhoR[l] = cons(membrane_cell,j,k,5+l);
                    }

                    // hard-coded for the case do_1d = 1 and nspecies = 2 (same species) -- copied from the old code for testing purpose
                    // if ((do_1D) and (nspecies==2)) {
                    if (0) {

                        Real tl = TL;
                        Real tr = TR;
                        Real sqrttl = sqrtTL;
                        Real sqrttr = sqrtTR;
                        // species 1
                        Real mm = molmass[0] / 6.02e23;
                        Real fac5 = transmission[0]*(std::pow(k_B,2.5))*6.0/std::sqrt(2*mm*3.142);;
                        Real fac3 = transmission[0]*(std::pow(k_B,1.5))*2.0/std::sqrt(2*mm*3.142);
                        Real fac1 = transmission[0]*(std::pow(k_B,0.5))/std::sqrt(2*mm*3.142);
                        Real rhol = rhoL[0];
                        Real rhor = rhoR[0];
                        Real um = fac3*(sqrttl*tl*rhol-sqrttr*tr*rhor);
                        Real nm = fac1*(sqrttl*rhol-sqrttr*rhor);

                        Real uv = fac5*(sqrttl*tl*tl*rhol+sqrttr*tr*tr*rhor);
                        Real nv = fac1*(sqrttl*rhol+sqrttr*rhor);

                        Real cross = fac3*(sqrttl*tl*rhol+sqrttr*tr*rhor);

                        Real corr = cross/(std::sqrt(uv)*std::sqrt(nv));

                        Real rn1 = amrex::RandomNormal(0.,1.);
                        Real rn2 = amrex::RandomNormal(0.,1.);
                        Real rn3 = rn1*corr + std::sqrt(1-pow(corr,2.))*rn2;
                            
                        xflux(membrane_cell,j,k,5)   = (dt*area*nm + std::sqrt(dt*area*mm*nv)*rn1)/vol; // species mass flux
                        xflux(membrane_cell,j,k,0)  += (dt*area*nm + std::sqrt(dt*area*mm*nv)*rn1)/vol; // mass flux
                        xflux(membrane_cell,j,k,4)  += (dt*area*um + std::sqrt(dt*area*mm*uv)*rn3)/(vol*mm);

                        // species 2
                        mm = molmass[1] / 6.02e23;
                        fac5 = transmission[1]*(std::pow(k_B,2.5))*6.0/std::sqrt(2*mm*3.142);;
                        fac3 = transmission[1]*(std::pow(k_B,1.5))*2.0/std::sqrt(2*mm*3.142);
                        fac1 = transmission[1]*(std::pow(k_B,0.5))/std::sqrt(2*mm*3.142);
                        rhol = rhoL[1];
                        rhor = rhoR[1];
                        um = fac3*(sqrttl*tl*rhol-sqrttr*tr*rhor);
                        nm = fac1*(sqrttl*rhol-sqrttr*rhor);

                        uv = fac5*(sqrttl*tl*tl*rhol+sqrttr*tr*tr*rhor);
                        nv = fac1*(sqrttl*rhol+sqrttr*rhor);

                        cross = fac3*(sqrttl*tl*rhol+sqrttr*tr*rhor);

                        corr = cross/(std::sqrt(uv)*std::sqrt(nv));

                        rn1 = amrex::RandomNormal(0.,1.);
                        rn2 = amrex::RandomNormal(0.,1.);
                        rn3 = rn1*corr + std::sqrt(1-pow(corr,2.))*rn2;
                            
                        xflux(membrane_cell,j,k,6)   = (dt*area*nm + std::sqrt(dt*area*mm*nv)*rn1)/vol; // species mass flux
                        xflux(membrane_cell,j,k,0)  += (dt*area*nm + std::sqrt(dt*area*mm*nv)*rn1)/vol; // mass flux
                        xflux(membrane_cell,j,k,4)  += (dt*area*um + std::sqrt(dt*area*mm*uv)*rn3)/(vol*mm);

                    }
                    
                    // new membrane code for multispecies
                    else {
                    // Compute effusive flux per species
                    for (int l=0;l<nspecies;++l) {
                        
                        if (amrex::Math::abs(transmission[l]) >= 1.e-10) {
                        
                            // Mean effusive accumulation
                            Real efffluxM = 0.0;
                            Real efffluxPy = 0.0;
                            Real efffluxPz = 0.0;
                            Real efffluxE = 0.0;

                            if (do_1D) {
                                efffluxM = dt*alpha[l]*(rhoL[l]*sqrtTL - rhoR[l]*sqrtTR); // mass
                                Real EL = (2.0*k_B*TL/mass[l]);
                                Real ER = (2.0*k_B*TR/mass[l]);
                                efffluxE = dt*alpha[l]*(rhoL[l]*sqrtTL*EL - rhoR[l]*sqrtTR*ER); // E
                            }
                            else if (do_2D) {
                                efffluxM = dt*alpha[l]*(rhoL[l]*sqrtTL - rhoR[l]*sqrtTR); // mass
                                efffluxPy = dt*alpha[l]*(rhoL[l]*vyL*sqrtTL - rhoR[l]*vyR*sqrtTR); // py
                                Real EL = (2.0*k_B*TL/mass[l]) + 0.5*vyL*vyL;
                                Real ER = (2.0*k_B*TR/mass[l]) + 0.5*vyR*vyR;
                                efffluxE = dt*alpha[l]*(rhoL[l]*sqrtTL*EL - rhoR[l]*sqrtTR*ER); // E
                            }
                            else {
                                efffluxM = dt*alpha[l]*(rhoL[l]*sqrtTL - rhoR[l]*sqrtTR); // mass
                                efffluxPy = dt*alpha[l]*(rhoL[l]*vyL*sqrtTL - rhoR[l]*vyR*sqrtTR); // py
                                efffluxPz = dt*alpha[l]*(rhoL[l]*vzL*sqrtTL - rhoR[l]*vzR*sqrtTR); // pz
                                Real EL = (2.0*k_B*TL/mass[l]) + 0.5*(vyL*vyL + vzL*vzL);
                                Real ER = (2.0*k_B*TR/mass[l]) + 0.5*(vyR*vyR + vzR*vzR);
                                efffluxE = dt*alpha[l]*(rhoL[l]*sqrtTL*EL - rhoR[l]*sqrtTR*ER); // E
                            }
                            
                            if (stoch_stress_form == 1) {

                                if (do_1D) {

                                    // Set the covariance matrix to 0
                                    Array2D<Real, 0, 1, 0, 1> Delvar; // covariance matrix [m, E]

                                    // Cholesky factor of the covariance matrix
                                    Array2D<Real, 0, 1, 0, 1> sqrtVar;
                                    for (int p=0;p<2;++p) {
                                        for (int q=0;q<2;++q) {
                                            Delvar(p,q) = 0.0;
                                            sqrtVar(p,q) = 0.0;
                                        }
                                    }

                                    // Fill variances
                                    Delvar(0,0) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL + rhoR[l]*sqrtTR);

                                    Real ELV = 6.0*k_B*k_B*TL*TL;
                                    Real ERV = 6.0*k_B*k_B*TR*TR;
                                    Delvar(1,1) = alpha[l]*(rhoL[l]*sqrtTL*ELV + rhoR[l]*sqrtTR*ERV)/mass[l];
                                    
                                    // Fill covariances
                                    Real meL = 2.0*k_B*TL;
                                    Real meR = 2.0*k_B*TR;
                                    Delvar(0,1) = alpha[l]*(rhoL[l]*sqrtTL*meL + rhoR[l]*sqrtTR*meR);
                                    Delvar(1,0) = alpha[l]*(rhoL[l]*sqrtTL*meL + rhoR[l]*sqrtTR*meR);

                                    Chol1D(Delvar,sqrtVar);

                                    // Generate random numbers
                                    GpuArray<Real,2> rand; // Random normal numbers
                                    for (int n=0;n<2;++n) {
                                        //rand[n] = RandomNormal(0.,1.,engine);
                                        rand[n] = RandomNormal(0.,1.);
                                    }
                                    
                                    // Get effusive accumulation from Langevin integration
                                    for (int q=0; q<2; ++q) {
                                        efffluxM += sqrtVar(0,q)*rand[q]*sqrt(dt); // x = mean + rand*covar
                                    }
                                    for (int q=0; q<2; ++q) {
                                        efffluxE += sqrtVar(1,q)*rand[q]*sqrt(dt); // x = mean + rand*covar
                                    }
                                }

                                else if (do_2D) {
                                    // Set the covariance matrix to 0
                                    Array2D<Real, 0, 2, 0, 2> Delvar; // covariance matrix [m, py, E]
                                    for (int p=0;p<3;++p) {
                                        for (int q=0;q<3;++q) {
                                            Delvar(p,q) = 0.0;
                                        }
                                    }

                                    // Fill variances
                                    Delvar(0,0) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL + rhoR[l]*sqrtTR);
                                    
                                    Real pyL = k_B*TL + mass[l]*vyL*vyL;
                                    Real pyR = k_B*TR + mass[l]*vyR*vyR;
                                    Delvar(1,1) = alpha[l]*(rhoL[l]*sqrtTL*pyL + rhoR[l]*sqrtTR*pyR);

                                    Real ELV = 24.0*k_B*k_B*TL*TL + 12.0*k_B*mass[l]*TL*(vyL*vyL) + mass[l]*mass[l]*(vyL*vyL)*(vyL*vyL);
                                    Real ERV = 24.0*k_B*k_B*TR*TR + 12.0*k_B*mass[l]*TR*(vyR*vyR) + mass[l]*mass[l]*(vyR*vyR)*(vyR*vyR);
                                    Delvar(2,2) = alpha[l]*(rhoL[l]*sqrtTL*ELV + rhoR[l]*sqrtTR*ERV);
                                    
                                    // Fill covariances
                                    Delvar(0,1) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL*vyL + rhoR[l]*sqrtTR*vyR);
                                    Delvar(1,0) = Delvar(0,1);

                                    Real meL = 4.0*k_B*TL + mass[l]*(vyL*vyL);
                                    Real meR = 4.0*k_B*TR + mass[l]*(vyR*vyR);
                                    Delvar(0,2) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*meL + rhoR[l]*sqrtTR*meR);
                                    Delvar(2,0) = Delvar(0,2);

                                    Real pyeL = vyL*(6.0*k_B*TL + mass[l]*(vyL*vyL));
                                    Real pyeR = vyR*(6.0*k_B*TR + mass[l]*(vyR*vyR));
                                    Delvar(1,2) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*pyeL + rhoR[l]*sqrtTR*pyeR);
                                    Delvar(2,1) = Delvar(1,2);

                                    // Cholesky factorise the covariance matrix
                                    Array2D<Real, 0, 2, 0, 2> sqrtVar;
                                    Chol2D(Delvar,sqrtVar);

                                    // Generate random numbers
                                    GpuArray<Real,3> rand; // Random normal numbers
                                    for (int n=0;n<3;++n) {
                                        //rand[n] = RandomNormal(0.,1.,engine);
                                        rand[n] = RandomNormal(0.,1.);
                                    }
                                    
                                    // Get effusive fluxes from Langevin integration
                                    for (int q=0; q<3; ++q) {
                                        efffluxM += sqrtVar(0,q)*rand[q]; // x = mean + rand*covar
                                    }
                                    for (int q=0; q<3; ++q) {
                                        efffluxPy += sqrtVar(1,q)*rand[q]; // x = mean + rand*covar
                                    }
                                    for (int q=0; q<3; ++q) {
                                        efffluxE += sqrtVar(2,q)*rand[q]; // x = mean + rand*covar
                                    }
                                }

                                else {
                                    // Set the covariance matrix to 0
                                    Array2D<Real, 0, 3, 0, 3> Delvar; // covariance matrix
                                    for (int p=0;p<4;++p) {
                                        for (int q=0;q<4;++q) {
                                            Delvar(p,q) = 0.0;
                                        }
                                    }

                                    // Fill variances
                                    Delvar(0,0) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL + rhoR[l]*sqrtTR);
                                    
                                    Real pyL = k_B*TL + mass[l]*vyL*vyL;
                                    Real pyR = k_B*TR + mass[l]*vyR*vyR;
                                    Delvar(1,1) = alpha[l]*(rhoL[l]*sqrtTL*pyL + rhoR[l]*sqrtTR*pyR);

                                    Real pzL = k_B*TL + mass[l]*vzL*vzL;
                                    Real pzR = k_B*TR + mass[l]*vzR*vzR;
                                    Delvar(2,2) = alpha[l]*(rhoL[l]*sqrtTL*pzL + rhoR[l]*sqrtTR*pzR);

                                    Real ELV = 24.0*k_B*k_B*TL*TL + 12.0*k_B*mass[l]*TL*(vyL*vyL+vzL*vzL) + mass[l]*mass[l]*(vyL*vyL+vzL*vzL)*(vyL*vyL+vzL*vzL);
                                    Real ERV = 24.0*k_B*k_B*TR*TR + 12.0*k_B*mass[l]*TR*(vyR*vyR+vzR*vzR) + mass[l]*mass[l]*(vyR*vyR+vzR*vzR)*(vyR*vyR+vzR*vzR);
                                    Delvar(3,3) = alpha[l]*(rhoL[l]*sqrtTL*ELV + rhoR[l]*sqrtTR*ERV);
                                    
                                    // Fill covariances
                                    Delvar(0,1) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL*vyL + rhoR[l]*sqrtTR*vyR);
                                    Delvar(1,0) = Delvar(0,1);
                                    Delvar(0,2) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL*vzL + rhoR[l]*sqrtTR*vzR);
                                    Delvar(2,0) = Delvar(0,2);

                                    Real meL = 4.0*k_B*TL + mass[l]*(vyL*vyL+vzL*vzL);
                                    Real meR = 4.0*k_B*TR + mass[l]*(vyR*vyR+vzR*vzR);
                                    Delvar(0,3) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*meL + rhoR[l]*sqrtTR*meR);
                                    Delvar(3,0) = Delvar(0,3);

                                    Real pyeL = vyL*(6.0*k_B*TL + mass[l]*(vyL*vyL+vzL*vzL));
                                    Real pyeR = vyR*(6.0*k_B*TR + mass[l]*(vyR*vyR+vzR*vzR));
                                    Delvar(1,3) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*pyeL + rhoR[l]*sqrtTR*pyeR);
                                    Delvar(3,1) = Delvar(1,3);

                                    Real pzeL = vzL*(6.0*k_B*TL + mass[l]*(vyL*vyL+vzL*vzL));
                                    Real pzeR = vzR*(6.0*k_B*TR + mass[l]*(vyR*vyR+vzR*vzR));
                                    Delvar(2,3) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*pzeL + rhoR[l]*sqrtTR*pzeR);
                                    Delvar(3,2) = Delvar(2,3);

                                    // Cholesky factorise the covariance matrix
                                    Array2D<Real, 0, 3, 0, 3> sqrtVar;
                                    Chol3D(Delvar,sqrtVar);

                                    // Generate random numbers
                                    GpuArray<Real,4> rand; // Random normal numbers
                                    for (int n=0;n<4;++n) {
                                        //rand[n] = RandomNormal(0.,1.,engine);
                                        rand[n] = RandomNormal(0.,1.);
                                    }
                                    
                                    // Get effusive fluxes from Langevin integration
                                    for (int q=0; q<4; ++q) {
                                        efffluxM += sqrtVar(0,q)*rand[q]; // x = mean + rand*covar
                                    }
                                    for (int q=0; q<4; ++q) {
                                        efffluxPy += sqrtVar(1,q)*rand[q]; // x = mean + rand*covar
                                    }
                                    for (int q=0; q<4; ++q) {
                                        efffluxPz += sqrtVar(2,q)*rand[q]; // x = mean + rand*covar
                                    }
                                    for (int q=0; q<4; ++q) {
                                        efffluxE += sqrtVar(3,q)*rand[q]; // x = mean + rand*covar
                                    }
                                }
                            }

                            // Increment total densities
                            xflux(membrane_cell,j,k,0)  += efffluxM/vol; // mass 
                            xflux(membrane_cell,j,k,5+l) = efffluxM/vol; // species mass
                            xflux(membrane_cell,j,k,4)  += efffluxE/vol; // energy
                            if (do_1D) {
                            }
                            else if (do_2D) {
                                edgex_v(membrane_cell,j,k)  += efffluxPy/vol; // y-momentum
                            }
                            else {
                                edgex_v(membrane_cell,j,k)  += efffluxPy/vol; // y-momentum
                                edgex_w(membrane_cell,j,k)  += efffluxPz/vol; // z-momentum
                            }
                        }

                    }
                    }
                }
                }
                    
//                }
//            });

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

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<Real>& cons = cons_in.array(mfi);
        const Array4<Real>& xflux = faceflux[0].array(mfi);
        const Array4<Real>& ymom = cumom[1].array(mfi);
        const Array4<Real>& zmom = cumom[2].array(mfi);
        const Array4<Real>& edgex_v = edgeflux_x[0].array(mfi);
        const Array4<Real>& edgex_w = edgeflux_x[1].array(mfi);

        if (bx.smallEnd(0) == membrane_cell) { // membrane at the left end (cell to the right of the membrane)

//            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//            {
//                if (i == membrane_cell) {
            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {
                    cons(membrane_cell,j,k,0) += xflux(membrane_cell,j,k,0);
                    cons(membrane_cell,j,k,4) += xflux(membrane_cell,j,k,4);
                    for (int l=0;l<nspecies;++l) {
                        cons(membrane_cell,j,k,5+l) += xflux(membrane_cell,j,k,5+l);
                    }
                    ymom(membrane_cell,j,k) += edgex_v(membrane_cell,j,k); 
                    zmom(membrane_cell,j,k) += edgex_w(membrane_cell,j,k); 
            }
            }
//                }
//            });
        }

        else if (bx.bigEnd(0) == membrane_cell - 1) { // membrane at the right end (cell to the left of the membrane)

//            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//            {
//                if (i == membrane_cell-1) {
            for (auto k = lo.z; k <= hi.z; ++k) {
            for (auto j = lo.y; j <= hi.y; ++j) {
                    cons(membrane_cell-1,j,k,0) -= xflux(membrane_cell,j,k,0);
                    cons(membrane_cell-1,j,k,4) -= xflux(membrane_cell,j,k,4);
                    for (int l=0;l<nspecies;++l) {
                        cons(membrane_cell-1,j,k,5+l) -= xflux(membrane_cell,j,k,5+l);
                    }
                    ymom(membrane_cell-1,j,k) -= edgex_v(membrane_cell,j,k); 
                    zmom(membrane_cell-1,j,k) -= edgex_w(membrane_cell,j,k); 
            }
            }
//                }
//            });
        
        }

    }

}

///////////////////////////////////////////////////////////////////////////
// compute momentum and fluxes at the membrane ////////////////////////////
// for SSA-type membrane; for now membrane is always perpendicular to x ///
///////////////////////////////////////////////////////////////////////////
void 
ComputeFluxMomMembrane(const MultiFab& cons0_in, const MultiFab& prim0_in,
                       const std::array<MultiFab, AMREX_SPACEDIM>& vel0,
                       std::array<MultiFab, AMREX_SPACEDIM>& cumom_mem,
                       std::array<MultiFab, AMREX_SPACEDIM>& faceflux_mem,
                       const amrex::Geometry& geom,
                       const amrex::Real dt)
{
    BL_PROFILE_VAR("ComputeFluxMomMembrane()",ComputeFluxMomMembrane);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    Box dom(geom.Domain());

    Real N_A = Runiv/k_B; // Avagadro's number
    
    GpuArray<Real,MAX_SPECIES> mass;
    for (int l=0;l<nspecies;++l) {
        mass[l] = molmass[l]/(N_A);
    }

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        cumom_mem[d].setVal(0.0);
        faceflux_mem[d].setVal(0.0);
    }

    Real area, vol;
    area = dx[1]*dx[2];
    vol  = dx[0]*dx[1]*dx[2];

    // face-based flux (mass and energy) and normal momentum /////
    //////////////////////////////////////////////////////////////
    for (MFIter mfi(faceflux_mem[0]); mfi.isValid(); ++mfi) {
        
        const Box& bx = mfi.validbox();
        
        const Array4<const Real> cons0   = cons0_in.array(mfi);
        const Array4<const Real> prim0   = prim0_in.array(mfi);
        const Array4<const Real>& xvel0  = (vel0[0]).array(mfi);
        const Array4<const Real>& yvel0  = (vel0[1]).array(mfi);
        const Array4<const Real>& zvel0  = (vel0[2]).array(mfi);
        
        const Array4<Real>& xmom  = (cumom_mem[0]).array(mfi);
        const Array4<Real>& xflux = (faceflux_mem[0]).array(mfi);

        // membrane on hi side of bx
        if ((membrane_cell >= bx.smallEnd(0)) and (membrane_cell <= bx.bigEnd(0))) {
            amrex::ParallelForRNG(bx, 
            [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
            {
                if (i == membrane_cell) {
                    ////////////////// left cell ///////////////////////
                    
                    ////////////////////////////////////////////////////
                    // number of particles in FHD cell (for correction)
                    //Real N = 0.0;
                    //GpuArray<Real,MAX_SPECIES> N_i;
                    //for (int n=0;n<nspecies;++n) {
                    //    N_i[n] = vol*(cons0(i-1,j,k,5+n)/mass[n]);
                    //    N += N_i[n];
                    //}
                    ////////////////////////////////////////////////////
                    
                    Real T = prim0(i-1,j,k,4);
                    Real Vx, Vy, Vz;
                    Vx = 0.5*(xvel0(i-1,j,k) + xvel0(i,j,k)); // normal
//                    Vx = 0.0;
                    Vy = 0.5*(yvel0(i-1,j,k) + yvel0(i-1,j+1,k)); // tangential
                    Vz = 0.5*(zvel0(i-1,j,k) + zvel0(i-1,j,k+1)); // tangential
                    Real mass_cross; // total mass into membrane
                    GpuArray<Real,3> mom_cross; // total momentum into membrane
                    Real en_cross; // total energy into membrane
                    GpuArray<Real,MAX_SPECIES> spec_mass_cross;
                    GpuArray<Real,MAX_SPECIES> rhoYk;
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i-1,j,k,5+n);
                    }
                    if (do_1D) {
                        poisson_process_membrane(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                                 mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vy,Vz,
                                                 transmission,engine);
                    }
                    else if (do_2D) {
                        poisson_process_membrane(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                                 mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vy,Vz,
                                                 transmission,engine);
                    }
                    else {
                        poisson_process_membrane(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                                 mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vy,Vz,
                                                 transmission,engine);
                    }
                    xflux(i,j,k,0) += mass_cross/(dt*area); // update mass flux
                    xflux(i,j,k,4) += en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) += spec_mass_cross[n]/(dt*area); // update species flux
                    }
                    
                    if (do_1D) {
                        xmom(i,j,k) += mom_cross[0]/vol;
                    }
                    else if (do_2D) {
                        xmom(i,j,k) += mom_cross[0]/vol;
//                        xflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
                        xflux(i,j,k,1) += 0.0;
                    }
                    else {
                        xmom(i,j,k) += mom_cross[0]/vol;
//                        xflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
//                        xflux(i,j,k,2) += mom_cross[2]/dt/area; // tangential flux
                        xflux(i,j,k,1) += 0.0;
                        xflux(i,j,k,2) += 0.0;
                    }
                    
                    ////////////////////////////////////////////////////

                    ////////////////// right cell //////////////////////

                    T = prim0(i,j,k,4);
                    Vx = -0.5*(xvel0(i,j,k)+xvel0(i+1,j,k)); // normal
//                    Vx = 0.0;
                    Vy =  0.5*(yvel0(i,j,k)+yvel0(i,j+1,k)); // tangential
                    Vz =  0.5*(zvel0(i,j,k)+zvel0(i,j,k+1)); // tangential
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i,j,k,5+n);
                    }
                    if (do_1D) {
                        poisson_process_membrane(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                                 mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vy,Vz,
                                                 transmission,engine);
                    }
                    else if (do_2D) {
                        poisson_process_membrane(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                                 mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vy,Vz,
                                                 transmission,engine);
                    }
                    else {
                        poisson_process_membrane(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                                 mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vy,Vz,
                                                 transmission,engine);
                    }
                    xflux(i,j,k,0) -= mass_cross/(dt*area); // update mass flux
                    xflux(i,j,k,4) -= en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) -= spec_mass_cross[n]/(dt*area); // update species flux
                    }
                    
                    if (do_1D) {
                        xmom(i,j,k) -= mom_cross[0]/vol;
                    }
                    else if (do_2D) {
                        xmom(i,j,k) -= mom_cross[0]/vol;
//                        xflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
                        xflux(i,j,k,1) -= 0.0;
                    }
                    else {
                        xmom(i,j,k) -= mom_cross[0]/vol;
//                        xflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
//                        xflux(i,j,k,2) -= mom_cross[2]/dt/area; // tangential flux
                        xflux(i,j,k,1) -= 0.0;
                        xflux(i,j,k,2) -= 0.0;
                    }
                    ////////////////////////////////////////////////////
                }
            });
        }
    }

    cumom_mem[0].OverrideSync(geom.periodicity());
    faceflux_mem[0].OverrideSync(geom.periodicity());
}

///////////////////////////////////////////////////////////////////////////
// reset fluxes at the membrane  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void 
ResetMembraneFluxes(const std::array<MultiFab, AMREX_SPACEDIM>& faceflux_mem,
                    std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                    std::array< MultiFab, 2 >& edgeflux_x,
                    std::array< MultiFab, 2 >& edgeflux_y,
                    std::array< MultiFab, 2 >& edgeflux_z,
                    const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("ResetMembraneFluxes()",ResetMembraneFluxes);

    for (MFIter mfi(faceflux[0]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const Array4<Real>& xflux           = (faceflux[0]).array(mfi);
        const Array4<const Real>& xflux_mem = (faceflux_mem[0]).array(mfi);

        if ((membrane_cell >= bx.smallEnd(0)) and (membrane_cell <= bx.bigEnd(0))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == membrane_cell) {
                    xflux(i,j,k,0) = xflux_mem(i,j,k,0); // mass flux
                    xflux(i,j,k,4) = xflux_mem(i,j,k,4); // energy flux
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) = xflux_mem(i,j,k,5+n); // species flux
                    }
                }
            });
        }
    }
}

///////////////////////////////////////////////////////////////////////////
// reset momentum at the membrane /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
void 
ResetMembraneMom(std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                 const std::array<MultiFab, AMREX_SPACEDIM>& cumom_mem,
                 const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("ResetMembraneMom()",ResetMembraneMom);

    for (MFIter mfi(cumom[0]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const Array4<Real>& xmom           = (cumom[0]).array(mfi);
        const Array4<const Real>& xmom_mem = (cumom_mem[0]).array(mfi);

        if ((membrane_cell >= bx.smallEnd(0)) and (membrane_cell <= bx.bigEnd(0))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == membrane_cell) {
                    xmom(i,j,k) = xmom_mem(i,j,k);
//                    xmom(i,j,k) = 0.0;
                }
            });
        }
    }
}

