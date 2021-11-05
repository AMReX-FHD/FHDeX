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
        alpha[l] = area*dt*transmission[l]*sqrt(k_B/(2.0*3.142*mass[l]));
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

                    // Compute effusive flux per species
                    for (int l=0;l<nspecies;++l) {
                        
                        if (amrex::Math::abs(transmission[l]) >= 1.e-10) {
                        
                            //amrex::Print(Print::AllProcs) << i << " " << j << " " << k << " " << l << " " << alpha[l] << " " << rhoL[l] << " " << rhoR[l] << " " << TL << " " << TR << " " << vyL << " " << vyR << " " << vzL << " " << vzR << std::endl;
                            //printf("%d %g %g %g %g %g %g %g %g %g\n",l,alpha[l],rhoL[l],rhoR[l],TL,TR,vyL,vyR,vzL,vzR);
                            
                            // Mean effusive fluxes
                            Real efffluxM = 0.0;
                            Real efffluxPy = 0.0;
                            Real efffluxPz = 0.0;
                            Real efffluxE = 0.0;

                            if (do_1D) {
                                efffluxM = alpha[l]*(rhoL[l]*sqrtTL - rhoR[l]*sqrtTR); // mass
                                Real EL = (2.0*k_B*TL/mass[l]);
                                Real ER = (2.0*k_B*TR/mass[l]);
                                efffluxE = alpha[l]*(rhoL[l]*sqrtTL*EL - rhoR[l]*sqrtTR*ER); // E
                            }
                            else if (do_2D) {
                                efffluxM = alpha[l]*(rhoL[l]*sqrtTL - rhoR[l]*sqrtTR); // mass
                                efffluxPy = alpha[l]*(rhoL[l]*vyL*sqrtTL - rhoR[l]*vyR*sqrtTR); // py
                                Real EL = (2.0*k_B*TL/mass[l]) + 0.5*vyL*vyL;
                                Real ER = (2.0*k_B*TR/mass[l]) + 0.5*vyR*vyR;
                                efffluxE = alpha[l]*(rhoL[l]*sqrtTL*EL - rhoR[l]*sqrtTR*ER); // E
                            }
                            else {
                                efffluxM = alpha[l]*(rhoL[l]*sqrtTL - rhoR[l]*sqrtTR); // mass
                                efffluxPy = alpha[l]*(rhoL[l]*vyL*sqrtTL - rhoR[l]*vyR*sqrtTR); // py
                                efffluxPz = alpha[l]*(rhoL[l]*vzL*sqrtTL - rhoR[l]*vzR*sqrtTR); // pz
                                Real EL = (2.0*k_B*TL/mass[l]) + 0.5*(vyL*vyL + vzL*vzL);
                                Real ER = (2.0*k_B*TR/mass[l]) + 0.5*(vyR*vyR + vzR*vzR);
                                efffluxE = alpha[l]*(rhoL[l]*sqrtTL*EL - rhoR[l]*sqrtTR*ER); // E
                            }
                            
                            if (stoch_stress_form == 1) {

                                if (do_1D) {
                                    // Set the covariance matrix to 0
                                    Array2D<Real, 0, 1, 0, 1> Delvar; // covariance matrix [m, E]
                                    for (int p=0;p<2;++p) {
                                        for (int q=0;q<2;++q) {
                                            Delvar(p,q) = 0.0;
                                        }
                                    }

                                    // Fill variances
                                    Delvar(0,0) = mass[l]*alpha[l]*(rhoL[l]*sqrtTL + rhoR[l]*sqrtTR);

                                    Real ELV = 24.0*k_B*k_B*TL*TL;
                                    Real ERV = 24.0*k_B*k_B*TR*TR;
                                    Delvar(1,1) = alpha[l]*(rhoL[l]*sqrtTL*ELV + rhoR[l]*sqrtTR*ERV);
                                    
                                    // Fill covariances
                                    Real meL = 4.0*k_B*TL;
                                    Real meR = 4.0*k_B*TR;
                                    Delvar(0,1) = 0.5*alpha[l]*(rhoL[l]*sqrtTL*meL + rhoR[l]*sqrtTR*meR);
                                    Delvar(1,0) = Delvar(0,1);

                                    // Cholesky factorise the covariance matrix
                                    Array2D<Real, 0, 1, 0, 1> sqrtVar;
                                    //amrex::Print(Print::AllProcs) << j << " " << k << " " << l << " " << TL << " " << TR << " " << rhoL[l] << " " << rhoR[l] << "\n";
                                    //amrex::Print(Print::AllProcs) << j << " " << k << " " << l << " " << Delvar(0,0) << " " << Delvar(1,1) << " " << Delvar(2,2) << " " << Delvar(3,3) << " " << Delvar(0,1) << " " << Delvar(0,2) << " " << Delvar(0,3) << " " << Delvar(1,2) << " " << Delvar(1,3) << " " << Delvar(2,3) << "\n";
                                    Chol1D(Delvar,sqrtVar);

                                    // Generate random numbers
                                    GpuArray<Real,2> rand; // Random normal numbers
                                    for (int n=0;n<2;++n) {
                                        //rand[n] = RandomNormal(0.,1.,engine);
                                        rand[n] = RandomNormal(0.,1.);
                                    }
                                    
                                    // Get effusive fluxes from Langevin integration
                                    for (int q=0; q<2; ++q) {
                                        efffluxM += sqrtVar(0,q)*rand[q]; // x = mean + rand*covar
                                    }
                                    for (int q=0; q<2; ++q) {
                                        efffluxE += sqrtVar(1,q)*rand[q]; // x = mean + rand*covar
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
                                    //amrex::Print(Print::AllProcs) << j << " " << k << " " << l << " " << TL << " " << TR << " " << rhoL[l] << " " << rhoR[l] << "\n";
                                    //amrex::Print(Print::AllProcs) << j << " " << k << " " << l << " " << Delvar(0,0) << " " << Delvar(1,1) << " " << Delvar(2,2) << " " << Delvar(3,3) << " " << Delvar(0,1) << " " << Delvar(0,2) << " " << Delvar(0,3) << " " << Delvar(1,2) << " " << Delvar(1,3) << " " << Delvar(2,3) << "\n";
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
                                    //amrex::Print(Print::AllProcs) << j << " " << k << " " << l << " " << TL << " " << TR << " " << rhoL[l] << " " << rhoR[l] << "\n";
                                    //amrex::Print(Print::AllProcs) << j << " " << k << " " << l << " " << Delvar(0,0) << " " << Delvar(1,1) << " " << Delvar(2,2) << " " << Delvar(3,3) << " " << Delvar(0,1) << " " << Delvar(0,2) << " " << Delvar(0,3) << " " << Delvar(1,2) << " " << Delvar(1,3) << " " << Delvar(2,3) << "\n";
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

                            // Increment total flux
                            xflux(membrane_cell,j,k,0)  += efffluxM/vol; // mass flux
                            xflux(membrane_cell,j,k,5+l) = efffluxM/vol; // species mass flux
                            xflux(membrane_cell,j,k,4)  += efffluxE/vol; // energy flux
                            if (do_1D) {
                            }
                            else if (do_2D) {
                                edgex_v(membrane_cell,j,k)  += efffluxPy/vol; // y-momentum flux
                            }
                            else {
                                edgex_v(membrane_cell,j,k)  += efffluxPy/vol; // y-momentum flux
                                edgex_w(membrane_cell,j,k)  += efffluxPz/vol; // z-momentum flux
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
