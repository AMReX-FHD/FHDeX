#include "compressible_functions.H"
#include "compressible_functions_stag.H"
#include "common_functions.H"

void calculateFluxStag(const MultiFab& cons_in, const std::array< MultiFab, AMREX_SPACEDIM >& cumom_in, 
                       const MultiFab& prim_in, const std::array< MultiFab, AMREX_SPACEDIM >& facevel_in,
                       const MultiFab& eta_in, const MultiFab& zeta_in, const MultiFab& kappa_in,
                       const MultiFab& chi_in, const MultiFab& D_in,
                       std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in,
                       std::array< MultiFab, 2 >& edgeflux_x_in,
                       std::array< MultiFab, 2 >& edgeflux_y_in,
                       std::array< MultiFab, 2 >& edgeflux_z_in,
                       std::array< MultiFab, 3 >& cenflux_in,
                       std::array<MultiFab, AMREX_SPACEDIM>& stochFlux_in,
                       MultiFab& rancorn_in,
                       const amrex::Geometry geom,
		                   const amrex::Vector< amrex::Real >& stoch_weights,
		                   const amrex::Real* dx, const amrex::Real dt)
{
    BL_PROFILE_VAR("calculateFluxStag()",calculateFluxStag);
    
    // from namelist
    int nspecies_gpu = nspecies;
    int algorithm_type_gpu = algorithm_type;
    int nvars_gpu = nvars;
    int nprimvars_gpu = nprimvars;
    Real Runiv_gpu = Runiv;
    int visc_type_gpu = visc_type;
    int n_cells_z = n_cells[2];
    Real k_B_gpu = k_B;
    
    // from namelist
    GpuArray<Real,MAX_SPECIES> hcv_gpu;
    for (int n=0; n<nspecies; ++n) {
        hcv_gpu[n] = hcv[n];
    }
    
    GpuArray<Real,MAX_SPECIES> hcp_gpu;
    for (int n=0; n<nspecies; ++n) {
        hcp_gpu[n] = hcp[n];
    }
    
    GpuArray<Real,MAX_SPECIES> molmass_gpu;
    for (int n=0; n<nspecies; ++n) {
        molmass_gpu[n] = molmass[n];
    }

    GpuArray<Real,AMREX_SPACEDIM> dx_gpu;
    for (int n=0; n<AMREX_SPACEDIM; ++n) {
        dx_gpu[n] = dx[n];
    }
    
    AMREX_D_TERM(faceflux_in[0].setVal(0.0);,
                 faceflux_in[1].setVal(0.0);,
                 faceflux_in[2].setVal(0.0););

    edgeflux_x_in[0].setVal(0.0);
    edgeflux_x_in[1].setVal(0.0);

    edgeflux_y_in[0].setVal(0.0);
    edgeflux_y_in[1].setVal(0.0);

    edgeflux_z_in[0].setVal(0.0);
    edgeflux_z_in[1].setVal(0.0);

    AMREX_D_TERM(cenflux_in[0].setVal(0.0);,
                 cenflux_in[1].setVal(0.0);,
                 cenflux_in[2].setVal(0.0););

    std::array< MultiFab, AMREX_SPACEDIM > tau_diag; // diagonal stress (defined at cell centers)
    AMREX_D_TERM(tau_diag[0].define(cons_in.boxArray(),cons_in.DistributionMap(),1,ngc);,
                 tau_diag[1].define(cons_in.boxArray(),cons_in.DistributionMap(),1,ngc);,
                 tau_diag[2].define(cons_in.boxArray(),cons_in.DistributionMap(),1,ngc););

    std::array< MultiFab, NUM_EDGE > tau_diagoff; // off diagonal stress at edges
    AMREX_D_TERM(tau_diagoff[0].define(convert(cons_in.boxArray(),nodal_flag_xy),cons_in.DistributionMap(),1,ngc);,
                 tau_diagoff[1].define(convert(cons_in.boxArray(),nodal_flag_yz),cons_in.DistributionMap(),1,ngc);,
                 tau_diagoff[2].define(convert(cons_in.boxArray(),nodal_flag_xz),cons_in.DistributionMap(),1,ngc););

    AMREX_D_TERM(tau_diag[0].setVal(0.0);,
                 tau_diag[1].setVal(0.0);,
                 tau_diag[2].setVal(0.0););

    AMREX_D_TERM(tau_diagoff[0].setVal(0.0);,
                 tau_diagoff[1].setVal(0.0);,
                 tau_diagoff[2].setVal(0.0););

    ////////////////////
    // stochastic fluxes
    ////////////////////
    
    if (stoch_stress_form == 1) {

        Real volinv = 1./(dx[0]*dx[1]*dx[2]);
        Real dtinv = 1./dt;
        
        // Loop over boxes
        for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {

            AMREX_D_TERM(const Array4<Real>& fluxx = faceflux_in[0].array(mfi); ,
                         const Array4<Real>& fluxy = faceflux_in[1].array(mfi); ,
                         const Array4<Real>& fluxz = faceflux_in[2].array(mfi));

            AMREX_D_TERM(const Array4<Real>& ranfluxx = stochFlux_in[0].array(mfi); ,
                         const Array4<Real>& ranfluxy = stochFlux_in[1].array(mfi); ,
                         const Array4<Real>& ranfluxz = stochFlux_in[2].array(mfi));

            const Array4<const Real> prim = prim_in.array(mfi);
            const Array4<const Real> cons = cons_in.array(mfi);

            const Array4<const Real> rancorn = rancorn_in.array(mfi);
        
            const Array4<const Real> eta   = eta_in.array(mfi);
            const Array4<const Real> zeta  = zeta_in.array(mfi);
            const Array4<const Real> kappa = kappa_in.array(mfi);
            const Array4<const Real> chi   = chi_in.array(mfi);
            const Array4<const Real> Dij   = D_in.array(mfi);

            const Box& tbx = mfi.nodaltilebox(0);
            const Box& tby = mfi.nodaltilebox(1);
            const Box& tbz = mfi.nodaltilebox(2);
        
            amrex::ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {


                GpuArray<Real,MAX_SPECIES+5> fweights;
                GpuArray<Real,MAX_SPECIES+5> weiner;
                
                GpuArray<Real,MAX_SPECIES> hk;
                GpuArray<Real,MAX_SPECIES> yy;
                GpuArray<Real,MAX_SPECIES> yyp;
                
                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> DijY_edge;
                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> sqD;
                                   
                Real muxp = (eta(i,j,k)*prim(i,j,k,4) + eta(i-1,j,k)*prim(i-1,j,k,4));
                Real kxp = (kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4) + kappa(i-1,j,k)*prim(i-1,j,k,4)*prim(i-1,j,k,4));

                Real meanT = 0.5*(prim(i,j,k,4)+prim(i-1,j,k,4));

                // Weights for facial fluxes:
                fweights[0] = 0; // No mass flux;
                fweights[1]=sqrt(k_B_gpu*muxp*volinv*dtinv);
                fweights[2]=fweights[1];
                fweights[3]=fweights[1];
                fweights[4]=sqrt(k_B_gpu*kxp*volinv*dtinv);

                // Construct the random increments
                for (int n=0; n<5; ++n) {
                    weiner[n] = fweights[n]*ranfluxx(i,j,k,n);
                }
                
                Real nweight=sqrt(k_B_gpu*volinv*dtinv);
                                
                if (n_cells_z > 1) {

                    // Corner viscosity coefficients in 3D
                    Real muzepp = 0.25*(eta(i,j,k)*prim(i,j,k,4) + eta(i-1,j,k)*prim(i-1,j,k,4) +
                                        eta(i,j+1,k)*prim(i,j+1,k,4) + eta(i-1,j+1,k)*prim(i-1,j+1,k,4) +
                                        eta(i,j,k+1)*prim(i,j,k+1,4) + eta(i-1,j,k+1)*prim(i-1,j,k+1,4) +
                                        eta(i,j+1,k+1)*prim(i,j+1,k+1,4) + eta(i-1,j+1,k+1)*prim(i-1,j+1,k+1,4) )/3.;
                    
                    Real muzemp = 0.25*(eta(i,j-1,k)*prim(i,j-1,k,4) + eta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                        eta(i,j,k)*prim(i,j,k,4) + eta(i-1,j,k)*prim(i-1,j,k,4) +
                                        eta(i,j-1,k+1)*prim(i,j-1,k+1,4) + eta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,4) +
                                        eta(i,j,k+1)*prim(i,j,k+1,4) + eta(i-1,j,k+1)*prim(i-1,j,k+1,4) )/3.;
                    
                    Real muzepm = 0.25*(eta(i,j,k-1)*prim(i,j,k-1,4) + eta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                        eta(i,j+1,k-1)*prim(i,j+1,k-1,4) + eta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,4) +
                                        eta(i,j,k)*prim(i,j,k,4) + eta(i-1,j,k)*prim(i-1,j,k,4) +
                                        eta(i,j+1,k)*prim(i,j+1,k,4) + eta(i-1,j+1,k)*prim(i-1,j+1,k,4) )/3.;
                    
                    Real muzemm = 0.25*(eta(i,j-1,k-1)*prim(i,j-1,k-1,4) + eta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) +
                                        eta(i,j,k-1)*prim(i,j,k-1,4) + eta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                        eta(i,j-1,k)*prim(i,j-1,k,4) + eta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                        eta(i,j,k)*prim(i,j,k,4) + eta(i-1,j,k)*prim(i-1,j,k,4) )/3.;

                    if (amrex::Math::abs(visc_type_gpu) == 3) {

                        muzepp = muzepp + 0.25*(zeta(i,j,k)*prim(i,j,k,4) + zeta(i-1,j,k)*prim(i-1,j,k,4) +
                                                zeta(i,j+1,k)*prim(i,j+1,k,4) + zeta(i-1,j+1,k)*prim(i-1,j+1,k,4) +
                                                zeta(i,j,k+1)*prim(i,j,k+1,4) + zeta(i-1,j,k+1)*prim(i-1,j,k+1,4) +
                                                zeta(i,j+1,k+1)*prim(i,j+1,k+1,4) + zeta(i-1,j+1,k+1)*prim(i-1,j+1,k+1,4) );
                        
                        muzemp = muzemp + 0.25*(zeta(i,j-1,k)*prim(i,j-1,k,4) + zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                                zeta(i,j,k)*prim(i,j,k,4) + zeta(i-1,j,k)*prim(i-1,j,k,4) +
                                                zeta(i,j-1,k+1)*prim(i,j-1,k+1,4) + zeta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,4) +
                                                zeta(i,j,k+1)*prim(i,j,k+1,4) + zeta(i-1,j,k+1)*prim(i-1,j,k+1,4) );
                        
                        muzepm = muzepm + 0.25*(zeta(i,j,k-1)*prim(i,j,k-1,4) + zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                                zeta(i,j+1,k-1)*prim(i,j+1,k-1,4) + zeta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,4) +
                                                zeta(i,j,k)*prim(i,j,k,4) + zeta(i-1,j,k)*prim(i-1,j,k,4) +
                                                zeta(i,j+1,k)*prim(i,j+1,k,4) + zeta(i-1,j+1,k)*prim(i-1,j+1,k,4) );
                        
                        muzemm = muzemm + 0.25*(zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) + zeta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) +
                                                zeta(i,j,k-1)*prim(i,j,k-1,4) + zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                                zeta(i,j-1,k)*prim(i,j-1,k,4) + zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                                zeta(i,j,k)*prim(i,j,k,4) + zeta(i-1,j,k)*prim(i-1,j,k,4) );
                    }

                    weiner[1] = weiner[1] + 0.25*nweight*(sqrt(muzepp)*rancorn(i,j+1,k+1)+
                                                          sqrt(muzemp)*rancorn(i,j,k+1) + sqrt(muzepm)* rancorn(i,j+1,k)+ 
                                                          sqrt(muzemm)*rancorn(i,j,k)); // Random "divergence" stress

                } else if (n_cells_z == 1) {

                    Abort("n_cells_z==1 case for stoch flux not written");
/*                    
          ! Corner viscosity coefficients in 2D
          muzepp = 0.5*(eta(i,j,k)*prim(i,j,k,5) + eta(i-1,j,k)*prim(i-1,j,k,5) + &
               eta(i,j+1,k)*prim(i,j+1,k,5) + eta(i-1,j+1,k)*prim(i-1,j+1,k,5) )/3.
          muzemp = 0.5*(eta(i,j-1,k)*prim(i,j-1,k,5) + eta(i-1,j-1,k)*prim(i-1,j-1,k,5) + &
               eta(i,j,k)*prim(i,j,k,5) + eta(i-1,j,k)*prim(i-1,j,k,5) )/3.

          if (abs(visc_type) .eq. 3) then

             muzepp = muzepp + 0.25*(zeta(i,j,k)*prim(i,j,k,5) + zeta(i-1,j,k)*prim(i-1,j,k,5) + &
                  zeta(i,j+1,k)*prim(i,j+1,k,5) + zeta(i-1,j+1,k)*prim(i-1,j+1,k,5) )
             muzemp = muzemp + 0.25*(zeta(i,j-1,k)*prim(i,j-1,k,5) + zeta(i-1,j-1,k)*prim(i-1,j-1,k,5) + &
                  zeta(i,j,k)*prim(i,j,k,5) + zeta(i-1,j,k)*prim(i-1,j,k,5) )

          endif

          weiner(2) = weiner(2) + 0.5*nweight*(sqrt(muzepp)*rancorn(i,j+1,k)+ &
               sqrt(muzemp)*rancorn(i,j,k)) ! Random "divergence" stress
*/
                }

                for (int n=1; n<5; ++n) {
                    fluxx(i,j,k,n) = fluxx(i,j,k,n) + weiner[n];
                }

                // Viscous heating:
                Real phiflx =  weiner[1]*(prim(i-1,j,k,1)+prim(i,j,k,1)) +
                    weiner[2]*(prim(i-1,j,k,2)+prim(i,j,k,2)) +
                    weiner[3]*(prim(i-1,j,k,3)+prim(i,j,k,3));

                phiflx =  - 0.5*phiflx;

                fluxx(i,j,k,4) = fluxx(i,j,k,4) - phiflx;

                if (algorithm_type_gpu == 2) {

                    for (int n=5; n<5+nspecies_gpu; ++n) {
                        weiner[n] = 0.;
                    }

                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        yy[ns] = amrex::max(0.,amrex::min(1.,prim(i-1,j,k,6+ns)));
                        yyp[ns] = amrex::max(0.,amrex::min(1.,prim(i,j,k,6+ns)));
                    }

                    Real sumy = 0.;
                    Real sumyp = 0.;

                    for (int n=0; n<nspecies_gpu; ++n) {
                        sumy += yy[n];
                        sumyp += yyp[n];
                    }

                    for (int n=0; n<nspecies_gpu; ++n) {
                        yy[n] /= sumy;
                        yyp[n] /= sumyp;
                    }

                    Real MWmix = 0.;

                    for (int ns=0; ns<nspecies_gpu; ++ns) {

                        MWmix = MWmix + 0.5*(yy[ns]+yyp[ns])/molmass_gpu[ns];

                        for (int ll=0; ll<nspecies_gpu; ++ll) {
                            DijY_edge[ns*nspecies_gpu+ll] = 0.5*(Dij(i-1,j,k,ll*nspecies_gpu+ns)*yy[ll] +
                                                                 Dij(i,j,k,ll*nspecies_gpu+ns)*yyp[ll] +
                                                                (Dij(i-1,j,k,ns*nspecies_gpu+ll)*yy[ns] +
                                                                 Dij(i,j,k,ns*nspecies_gpu+ll)*yyp[ns] ));

                        }
                    }
                    
                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        if (amrex::Math::abs(yy[ns]) + amrex::Math::abs(yyp[ns]) <= 1.e-12) {
                            for (int n=0; n<nspecies_gpu; ++n) {
                                DijY_edge[ns*nspecies_gpu+n]=0.;
                                DijY_edge[n*nspecies_gpu+ns]=0.;
                            }
                        }
                    }

                    MWmix = 1. / MWmix;

                    CholeskyDecomp(DijY_edge,nspecies_gpu,sqD);

                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        for (int ll=0; ll<=ns; ++ll) {
                            fweights[5+ll]=sqrt(k_B_gpu*MWmix*volinv/(Runiv_gpu*dt))*sqD[ns*nspecies_gpu+ll];
                            weiner[5+ns] = weiner[5+ns] + fweights[5+ll]*ranfluxx(i,j,k,5+ll);
                        }
                        fluxx(i,j,k,5+ns) = weiner[5+ns];
                    }

                    GetEnthalpies(meanT, hk, hcp_gpu, nspecies_gpu);

                    Real soret = 0.;

                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        soret = soret + (hk[ns] + Runiv_gpu*meanT/molmass_gpu[ns]
                                         *0.5*(chi(i-1,j,k,ns)+chi(i,j,k,ns)))*weiner[5+ns];
                    }
                    fluxx(i,j,k,4) = fluxx(i,j,k,4) + soret;
                }

            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                GpuArray<Real,MAX_SPECIES+5> fweights;
                GpuArray<Real,MAX_SPECIES+5> weiner;
                
                GpuArray<Real,MAX_SPECIES> hk;
                GpuArray<Real,MAX_SPECIES> yy;
                GpuArray<Real,MAX_SPECIES> yyp;
                
                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> DijY_edge;
                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> sqD;
                

                Real muyp = eta(i,j,k)*prim(i,j,k,4) + eta(i,j-1,k)*prim(i,j-1,k,4);
                Real kyp = kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4) + kappa(i,j-1,k)*prim(i,j-1,k,4)*prim(i,j-1,k,4);

                Real meanT = 0.5*(prim(i,j,k,4)+prim(i,j-1,k,4));

                // Weights for facial fluxes:
                fweights[0] = 0; // No mass flux
                fweights[1] = sqrt(k_B_gpu*muyp*volinv*dtinv);
                fweights[2] = sqrt(k_B_gpu*muyp*volinv*dtinv);
                fweights[3] = sqrt(k_B_gpu*muyp*volinv*dtinv);
                fweights[4] = sqrt(k_B_gpu*kyp*volinv*dtinv);

                // Construct the random increments
                for (int n=0; n<5; ++n) {
                    weiner[n] = fweights[n]*ranfluxy(i,j,k,n);
                }

                Real nweight=sqrt(k_B_gpu*volinv*dtinv);
                                
                if (n_cells_z > 1) {

                    // Corner viscosity coefficients 3D
                    Real muzepp = 0.25*(eta(i+1,j-1,k)*prim(i+1,j-1,k,4) + eta(i,j-1,k)*prim(i,j-1,k,4) +
                                        eta(i+1,j,k)*prim(i+1,j,k,4) + eta(i,j,k)*prim(i,j,k,4) +
                                        eta(i+1,j-1,k+1)*prim(i+1,j-1,k+1,4) + eta(i,j-1,k+1)*prim(i,j-1,k+1,4) +
                                        eta(i+1,j,k+1)*prim(i+1,j,k+1,4) + eta(i,j,k+1)*prim(i,j,k+1,4) )/3.;

                    Real muzemp = 0.25*(eta(i-1,j,k)*prim(i-1,j,k,4) + eta(i,j,k)*prim(i,j,k,4) +
                                        eta(i-1,j-1,k)*prim(i-1,j-1,k,4) + eta(i,j-1,k)*prim(i,j-1,k,4) +
                                        eta(i-1,j,k+1)*prim(i-1,j,k+1,4) + eta(i,j,k+1)*prim(i,j,k+1,4) +
                                        eta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,4) + eta(i,j-1,k+1)*prim(i,j-1,k+1,4) )/3.;

                    Real muzepm = 0.25*(eta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,4) + eta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                        eta(i+1,j,k-1)*prim(i+1,j,k-1,4) + eta(i,j,k-1)*prim(i,j,k-1,4) +
                                        eta(i+1,j-1,k)*prim(i+1,j-1,k,4) + eta(i,j-1,k)*prim(i,j-1,k,4) +
                                        eta(i+1,j,k)*prim(i+1,j,k,4) + eta(i,j,k)*prim(i,j,k,4) )/3.;

                    Real muzemm = 0.25*(eta(i-1,j,k-1)*prim(i-1,j,k-1,4) + eta(i,j,k-1)*prim(i,j,k-1,4) +
                                        eta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) + eta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                        eta(i-1,j,k)*prim(i-1,j,k,4) + eta(i,j,k)*prim(i,j,k,4) +
                                        eta(i-1,j-1,k)*prim(i-1,j-1,k,4) + eta(i,j-1,k)*prim(i,j-1,k,4) )/3.;

                    if (amrex::Math::abs(visc_type_gpu) == 3) {

                        muzepp = muzepp + 0.25*(zeta(i+1,j-1,k)*prim(i+1,j-1,k,4) + zeta(i,j-1,k)*prim(i,j-1,k,4) +
                                                zeta(i+1,j,k)*prim(i+1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) +
                                                zeta(i+1,j-1,k+1)*prim(i+1,j-1,k+1,4) + zeta(i,j-1,k+1)*prim(i,j-1,k+1,4) +
                                                zeta(i+1,j,k+1)*prim(i+1,j,k+1,4) + zeta(i,j,k+1)*prim(i,j,k+1,4) );

                        muzemp = muzemp + 0.25*(zeta(i-1,j,k)*prim(i-1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) +
                                                zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) + zeta(i,j-1,k)*prim(i,j-1,k,4) +
                                                zeta(i-1,j,k+1)*prim(i-1,j,k+1,4) + zeta(i,j,k+1)*prim(i,j,k+1,4) +
                                                zeta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,4) + zeta(i,j-1,k+1)*prim(i,j-1,k+1,4) );

                        muzepm =  muzepm +0.25*(zeta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,4) + zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                                zeta(i+1,j,k-1)*prim(i+1,j,k-1,4) + zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                                zeta(i+1,j-1,k)*prim(i+1,j-1,k,4) + zeta(i,j-1,k)*prim(i,j-1,k,4) +
                                                zeta(i+1,j,k)*prim(i+1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) );

                        muzemm = muzemm + 0.25*(zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) + zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                                zeta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) + zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                                zeta(i-1,j,k)*prim(i-1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) +
                                                zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) + zeta(i,j-1,k)*prim(i,j-1,k,4) );
                    }

                    weiner[2] = weiner[2] + 0.25*nweight*(sqrt(muzepp)*rancorn(i+1,j,k+1)+ sqrt(muzemp)*rancorn(i,j,k+1) + 
                                                          sqrt(muzepm)* rancorn(i+1,j,k)+ sqrt(muzemm)*rancorn(i,j,k)); // Random "divergence" stress

                } else if (n_cells_z == 1) {
                    
                    Abort("n_cells_z==1 case for stoch flux not written");
/*
          ! Corner viscosity coefficients 2D
          muzepp = 0.5d0*(eta(i+1,j-1,k)*prim(i+1,j-1,k,5) + eta(i,j-1,k)*prim(i,j-1,k,5) + &
               eta(i+1,j,k)*prim(i+1,j,k,5) + eta(i,j,k)*prim(i,j,k,5) )/3.d0

          muzemp = 0.5d0*(eta(i-1,j,k)*prim(i-1,j,k,5) + eta(i,j,k)*prim(i,j,k,5) + &
               eta(i-1,j-1,k)*prim(i-1,j-1,k,5) + eta(i,j-1,k)*prim(i,j-1,k,5) )/3.d0

          if (abs(visc_type) .eq. 3) then

             muzepp = muzepp + 0.25d0*(zeta(i+1,j-1,k)*prim(i+1,j-1,k,5) + zeta(i,j-1,k)*prim(i,j-1,k,5) + &
                  zeta(i+1,j,k)*prim(i+1,j,k,5) + zeta(i,j,k)*prim(i,j,k,5) )

             muzemp = muzemp + 0.25d0*(zeta(i-1,j,k)*prim(i-1,j,k,5) + zeta(i,j,k)*prim(i,j,k,5) + &
                  zeta(i-1,j-1,k)*prim(i-1,j-1,k,5) + zeta(i,j-1,k)*prim(i,j-1,k,5) )

          endif

          weiner(3) = weiner(3) + 0.5d0*nweight*    &
               (sqrt(muzepp)*rancorn(i+1,j,k) + sqrt(muzemp)*rancorn(i,j,k)) ! Random "divergence" stress
*/
                }

                for (int n=1; n<5; ++n) {
                    fluxy(i,j,k,n) = fluxy(i,j,k,n) + weiner[n];
                }
            
                // Viscous heating:
                Real phiflx =  weiner[1]*(prim(i,j-1,k,1)+prim(i,j,k,1)) +
                    weiner[2]*(prim(i,j-1,k,2)+prim(i,j,k,2)) +
                    weiner[3]*(prim(i,j-1,k,3)+prim(i,j,k,3));
                
                phiflx =  - 0.5*phiflx;

                fluxy(i,j,k,4) = fluxy(i,j,k,4) - phiflx;

                if (algorithm_type_gpu == 2) {

                    for (int n=5; n<5+nspecies_gpu; ++n) {
                        weiner[n] = 0.;
                    }

                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        yy[ns] = amrex::max(0.,amrex::min(1.,prim(i,j-1,k,6+ns)));
                        yyp[ns] = amrex::max(0.,amrex::min(1.,prim(i,j,k,6+ns)));
                    }

                    Real sumy = 0.;
                    Real sumyp = 0.;

                    for (int n=0; n<nspecies_gpu; ++n) {
                        sumy += yy[n];
                        sumyp += yyp[n];
                    }

                    for (int n=0; n<nspecies_gpu; ++n) {
                        yy[n] /= sumy;
                        yyp[n] /= sumyp;
                    }

                    Real MWmix = 0.;

                    for (int ns=0; ns<nspecies_gpu; ++ns) {

                        MWmix = MWmix + 0.5*(yy[ns]+yyp[ns])/molmass_gpu[ns];

                        for (int ll=0; ll<nspecies_gpu; ++ll) {
                            DijY_edge[ns*nspecies_gpu+ll] = 0.5*(Dij(i,j-1,k,ll*nspecies_gpu+ns)*yy[ll] +
                                                                 Dij(i,j,k,ll*nspecies_gpu+ns)*yyp[ll] +
                                                                (Dij(i,j-1,k,ns*nspecies_gpu+ll)*yy[ns] +
                                                                 Dij(i,j,k,ns*nspecies_gpu+ll)*yyp[ns] ));

                        }
                    }
                    
                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        if (amrex::Math::abs(yy[ns]) + amrex::Math::abs(yyp[ns]) <= 1.e-12) {
                            for (int n=0; n<nspecies_gpu; ++n) {
                                DijY_edge[ns*nspecies_gpu+n]=0.;
                                DijY_edge[n*nspecies_gpu+ns]=0.;
                            }
                        }
                    }

                    MWmix = 1. / MWmix;

                    CholeskyDecomp(DijY_edge,nspecies_gpu,sqD);

                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        for (int ll=0; ll<=ns; ++ll) {
                            fweights[5+ll]=sqrt(k_B_gpu*MWmix*volinv/(Runiv_gpu*dt))*sqD[ns*nspecies_gpu+ll];
                            weiner[5+ns] = weiner[5+ns] + fweights[5+ll]*ranfluxy(i,j,k,5+ll);
                        }
                        fluxy(i,j,k,5+ns) = weiner[5+ns];
                    }

                    GetEnthalpies(meanT, hk, hcp_gpu, nspecies_gpu);

                    Real soret = 0.;

                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        soret = soret + (hk[ns] + Runiv_gpu*meanT/molmass_gpu[ns]
                                         *0.5*(chi(i,j-1,k,ns)+chi(i,j,k,ns)))*weiner[5+ns];
                    }
                    fluxy(i,j,k,4) = fluxy(i,j,k,4) + soret;
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                GpuArray<Real,MAX_SPECIES+5> fweights;
                GpuArray<Real,MAX_SPECIES+5> weiner;
                
                GpuArray<Real,MAX_SPECIES> hk;
                GpuArray<Real,MAX_SPECIES> yy;
                GpuArray<Real,MAX_SPECIES> yyp;
                
                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> DijY_edge;
                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> sqD;
                
                if (n_cells_z > 1) {

                    Real muzp = eta(i,j,k)*prim(i,j,k,4) + eta(i,j,k-1)*prim(i,j,k-1,4);
                    Real kzp = kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4) + kappa(i,j,k-1)*prim(i,j,k-1,4)*prim(i,j,k-1,4);

                    Real meanT = 0.5*(prim(i,j,k,4)+prim(i,j,k-1,4));

                    // Weights for facial fluxes:
                    fweights[0] = 0; // No mass flux
                    fweights[1] = sqrt(k_B_gpu*muzp*volinv*dtinv);
                    fweights[2] = sqrt(k_B_gpu*muzp*volinv*dtinv);
                    fweights[3] = sqrt(k_B_gpu*muzp*volinv*dtinv);
                    fweights[4] = sqrt(k_B_gpu*kzp*volinv*dtinv);

                    // Construct the random increments
                    for (int n=0; n<5; ++n) {
                        weiner[n] = fweights[n]*ranfluxz(i,j,k,n);
                    }
                
                    Real nweight=sqrt(k_B_gpu*volinv*dtinv);

                    // Corner viscosity coefficients
                    Real muzepp = 0.25*(eta(i+1,j,k-1)*prim(i+1,j,k-1,4) + eta(i,j,k-1)*prim(i,j,k-1,4) +
                                        eta(i+1,j+1,k-1)*prim(i+1,j+1,k-1,4) + eta(i,j+1,k-1)*prim(i,j+1,k-1,4) +
                                        eta(i+1,j,k)*prim(i+1,j,k,4) + eta(i,j,k)*prim(i,j,k,4) +
                                        eta(i+1,j+1,k)*prim(i+1,j+1,k,4) + eta(i,j+1,k)*prim(i,j+1,k,4) )/3.;

                    Real muzemp = 0.25*(eta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,4) + eta(i,j+1,k-1)*prim(i,j+1,k-1,4) +
                                        eta(i-1,j,k-1)*prim(i-1,j,k-1,4) + eta(i,j,k-1)*prim(i,j,k-1,4) +
                                        eta(i-1,j+1,k)*prim(i-1,j+1,k,4) + eta(i,j+1,k)*prim(i,j+1,k,4) +
                                        eta(i-1,j,k)*prim(i-1,j,k,4) + eta(i,j,k)*prim(i,j,k,4) )/3.;

                    Real muzepm = 0.25*(eta(i+1,j,k)*prim(i+1,j,k,4) + eta(i,j,k-2)*prim(i,j,k,4) +
                                        eta(i+1,j-1,k)*prim(i+1,j-1,k,4) + eta(i,j-1,k-2)*prim(i,j-1,k,4) +
                                        eta(i+1,j,k-1)*prim(i+1,j,k-1,4) + eta(i,j,k-1)*prim(i,j,k-1,4) +
                                        eta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,4) + eta(i,j-1,k-1)*prim(i,j-1,k-1,4) )/3.;

                    Real muzemm = 0.25*(eta(i-1,j-1,k)*prim(i-1,j-1,k,4) + eta(i,j-1,k)*prim(i,j-1,k,4) +
                                        eta(i-1,j,k)*prim(i-1,j,k,4) + eta(i,j,k)*prim(i,j,k,4) +
                                        eta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) + eta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                        eta(i-1,j,k-1)*prim(i-1,j,k-1,4) + eta(i,j,k-1)*prim(i,j,k-1,4) )/3.;

                    if (amrex::Math::abs(visc_type_gpu) == 3) {

                        muzepp = muzepp+ 0.25*(zeta(i+1,j,k-1)*prim(i+1,j,k-1,4) + zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                               zeta(i+1,j+1,k-1)*prim(i+1,j+1,k-1,4) + zeta(i,j+1,k-1)*prim(i,j+1,k-1,4) +
                                               zeta(i+1,j,k)*prim(i+1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) +
                                               zeta(i+1,j+1,k)*prim(i+1,j+1,k,4) + zeta(i,j+1,k)*prim(i,j+1,k,4) );

                        muzemp = muzemp + 0.25*(zeta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,4) + zeta(i,j+1,k-1)*prim(i,j+1,k-1,4) +
                                                zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) + zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                                zeta(i-1,j+1,k)*prim(i-1,j+1,k,4) + zeta(i,j+1,k)*prim(i,j+1,k,4) +
                                                zeta(i-1,j,k)*prim(i-1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) );

                        muzepm = muzepm + 0.25*(zeta(i+1,j,k)*prim(i+1,j,k,4) + zeta(i,j,k-2)*prim(i,j,k,4) +
                                                zeta(i+1,j-1,k)*prim(i+1,j-1,k,4) + zeta(i,j-1,k-2)*prim(i,j-1,k,4) +
                                                zeta(i+1,j,k-1)*prim(i+1,j,k-1,4) + zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                                zeta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,4) + zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) );

                        muzemm = muzemm + 0.25*(zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) + zeta(i,j-1,k)*prim(i,j-1,k,4) +
                                                zeta(i-1,j,k)*prim(i-1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) +
                                                zeta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) + zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                                zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) + zeta(i,j,k-1)*prim(i,j,k-1,4) );

                    }

                    weiner[3] = weiner[3] + 0.25*nweight*   
                        (sqrt(muzepp)*rancorn(i+1,j+1,k)+ sqrt(muzemp)*rancorn(i,j+1,k) + 
                         sqrt(muzepm)* rancorn(i+1,j,k)+ sqrt(muzemm)*rancorn(i,j,k)); // Random "divergence" stress


                    for (int n=1; n<5; ++n) {
                        fluxz(i,j,k,n) = fluxz(i,j,k,n) + weiner[n];
                    }
                    
                    // Viscous heating:
                    Real phiflx =  weiner[1]*(prim(i,j,k-1,1)+prim(i,j,k,1)) +         
                        weiner[2]*(prim(i,j,k-1,2)+prim(i,j,k,2)) +
                        weiner[3]*(prim(i,j,k-1,3)+prim(i,j,k,3));
                    
                    phiflx =  - 0.5*phiflx;

                    fluxz(i,j,k,4) = fluxz(i,j,k,4) - phiflx;

                    if (algorithm_type_gpu == 2) {

                    for (int n=5; n<5+nspecies_gpu; ++n) {
                        weiner[n] = 0.;
                    }

                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        yy[ns] = amrex::max(0.,amrex::min(1.,prim(i,j,k-1,6+ns)));
                        yyp[ns] = amrex::max(0.,amrex::min(1.,prim(i,j,k,6+ns)));
                    }

                    Real sumy = 0.;
                    Real sumyp = 0.;

                    for (int n=0; n<nspecies_gpu; ++n) {
                        sumy += yy[n];
                        sumyp += yyp[n];
                    }

                    for (int n=0; n<nspecies_gpu; ++n) {
                        yy[n] /= sumy;
                        yyp[n] /= sumyp;
                    }

                    Real MWmix = 0.;

                    for (int ns=0; ns<nspecies_gpu; ++ns) {

                        MWmix = MWmix + 0.5*(yy[ns]+yyp[ns])/molmass_gpu[ns];

                        for (int ll=0; ll<nspecies_gpu; ++ll) {
                            DijY_edge[ns*nspecies_gpu+ll] = 0.5*(Dij(i,j,k-1,ll*nspecies_gpu+ns)*yy[ll] +
                                                                 Dij(i,j,k,ll*nspecies_gpu+ns)*yyp[ll] +
                                                                (Dij(i,j,k-1,ns*nspecies_gpu+ll)*yy[ns] +
                                                                 Dij(i,j,k,ns*nspecies_gpu+ll)*yyp[ns] ));

                        }
                    }


                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        if (amrex::Math::abs(yy[ns]) + amrex::Math::abs(yyp[ns]) <= 1.e-12) {
                            for (int n=0; n<nspecies_gpu; ++n) {
                                DijY_edge[ns*nspecies_gpu+n]=0.;
                                DijY_edge[n*nspecies_gpu+ns]=0.;
                            }
                        }
                    }

                    MWmix = 1. / MWmix;

                    CholeskyDecomp(DijY_edge,nspecies_gpu,sqD);

                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        for (int ll=0; ll<=ns; ++ll) {
                            fweights[5+ll]=sqrt(k_B_gpu*MWmix*volinv/(Runiv_gpu*dt))*sqD[ns*nspecies_gpu+ll];
                            weiner[5+ns] = weiner[5+ns] + fweights[5+ll]*ranfluxz(i,j,k,5+ll);
                        }
                        fluxz(i,j,k,5+ns) = weiner[5+ns];
                    }

                    GetEnthalpies(meanT, hk, hcp_gpu, nspecies_gpu);

                    Real soret = 0.;

                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        soret = soret + (hk[ns] + Runiv_gpu*meanT/molmass_gpu[ns]
                                         *0.5*(chi(i,j,k-1,ns)+chi(i,j,k,ns)))*weiner[5+ns];
                    }
                    fluxz(i,j,k,4) = fluxz(i,j,k,4) + soret;
                    
                    }
                }
                
            }); // end lambda function

        } // end MFIter

        // Loop over boxes
        for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {
            
            const Box& bx = mfi.tilebox();

            // call this to enforce flux boundary conditions
            stoch_flux(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                       cons_in[mfi].dataPtr(),  
                       prim_in[mfi].dataPtr(),    
                       faceflux_in[0][mfi].dataPtr(),
                       faceflux_in[1][mfi].dataPtr(),
#if (AMREX_SPACEDIM == 3)
                       faceflux_in[2][mfi].dataPtr(),
#endif
                       stochFlux_in[0][mfi].dataPtr(),
                       stochFlux_in[1][mfi].dataPtr(),
#if (AMREX_SPACEDIM == 3)
                       stochFlux_in[2][mfi].dataPtr(),
#endif
                       rancorn_in[mfi].dataPtr(),
                       eta_in[mfi].dataPtr(),  
                       zeta_in[mfi].dataPtr(),  
                       kappa_in[mfi].dataPtr(),
                       chi_in[mfi].dataPtr(),  
                       D_in[mfi].dataPtr(),  
                       ZFILL(dx), &dt);
        }
    }
        
    ////////////////////
    // diffusive fluxes
    ////////////////////
    
    // Loop over boxes
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {

        AMREX_D_TERM(const Array4<Real>& xflux = faceflux_in[0].array(mfi); ,
                     const Array4<Real>& yflux = faceflux_in[1].array(mfi); ,
                     const Array4<Real>& zflux = faceflux_in[2].array(mfi));

        const Array4<Real>& edgex_v = edgeflux_x_in[0].array(mfi);
        const Array4<Real>& edgex_w = edgeflux_x_in[1].array(mfi);
        const Array4<Real>& edgey_u = edgeflux_y_in[0].array(mfi);
        const Array4<Real>& edgey_w = edgeflux_y_in[1].array(mfi);
        const Array4<Real>& edgez_u = edgeflux_z_in[0].array(mfi);
        const Array4<Real>& edgez_v = edgeflux_z_in[1].array(mfi);

        const Array4<Real>& cenx_u = cenflux_in[0].array(mfi);
        const Array4<Real>& ceny_v = cenflux_in[1].array(mfi);
        const Array4<Real>& cenz_w = cenflux_in[2].array(mfi);

        const Array4<Real> tauxx = tau_diag[0].array(mfi);
        const Array4<Real> tauyy = tau_diag[1].array(mfi);
        const Array4<Real> tauzz = tau_diag[2].array(mfi);

        AMREX_D_TERM(const Array4<Real> tauxy = tau_diagoff[0].array(mfi);,
                     const Array4<Real> tauyz = tau_diagoff[1].array(mfi);,
                     const Array4<Real> tauxz = tau_diagoff[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& momx = cumom_in[0].array(mfi);,
                     Array4<Real const> const& momy = cumom_in[1].array(mfi);,
                     Array4<Real const> const& momz = cumom_in[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& velx = facevel_in[0].array(mfi);,
                     Array4<Real const> const& vely = facevel_in[1].array(mfi);,
                     Array4<Real const> const& velz = facevel_in[2].array(mfi););

        const Array4<const Real> prim = prim_in.array(mfi);
        const Array4<const Real> cons = cons_in.array(mfi);
        
        const Array4<const Real> eta   = eta_in.array(mfi);
        const Array4<const Real> zeta  = zeta_in.array(mfi);
        const Array4<const Real> kappa = kappa_in.array(mfi);
        const Array4<const Real> chi   = chi_in.array(mfi);
        const Array4<const Real> Dij   = D_in.array(mfi);

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        IntVect nd(AMREX_D_DECL(1,1,1));
        const Box& tbn = mfi.tilebox(nd);

        const Box & bx_xy = mfi.tilebox(nodal_flag_xy);
        #if (AMREX_SPACEDIM == 3)
        const Box & bx_xz = mfi.tilebox(nodal_flag_xz);
        const Box & bx_yz = mfi.tilebox(nodal_flag_yz);
        #endif

        const Box& bx = mfi.tilebox();

        Real half = 0.5;

        // Populate diagonal stress
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            Real u_x, v_y, w_z; // velocity gradients
            u_x = (velx(i+1,j,k) - velx(i,j,k))/dx_gpu[0];
            v_y = (vely(i,j+1,k) - vely(i,j,k))/dx_gpu[1];
            w_z = (velz(i,j,k+1) - velz(i,j,k))/dx_gpu[2];

            Real div = u_x + v_y + w_z;
            tauxx(i,j,k) = 2*eta(i,j,k)*u_x + (kappa(i,j,k) - 2*eta(i,j,k)/3.)*div; 
            tauyy(i,j,k) = 2*eta(i,j,k)*v_y + (kappa(i,j,k) - 2*eta(i,j,k)/3.)*div; 
            tauzz(i,j,k) = 2*eta(i,j,k)*w_z + (kappa(i,j,k) - 2*eta(i,j,k)/3.)*div;
            });

        // Populate off-diagonal stress
        amrex::ParallelFor(bx_xy, bx_xz, bx_yz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            Real u_y, v_x; // velocity gradients
            u_y = (velx(i,j,k) - velx(i,j-1,k))/dx_gpu[1];
            v_x = (vely(i,j,k) - vely(i-1,j,k))/dx_gpu[0];
            tauxy(i,j,k) = 0.25*(eta(i-1,j-1,k)+eta(i-1,j,k)+eta(i,j-1,k)+eta(i,j,k))*(u_y+v_x);
            
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            Real u_z, w_x; // velocity gradients
            u_z = (velx(i,j,k) - velx(i,j,k-1))/dx_gpu[2];
            w_x = (velz(i,j,k) - velz(i-1,j,k))/dx_gpu[0];
            tauxz(i,j,k) = 0.25*(eta(i-1,j,k-1)+eta(i-1,j,k)+eta(i,j,k-1)+eta(i,j,k))*(u_z+w_x);
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            Real v_z, w_y; // velocity gradients
            v_z = (vely(i,j,k) - vely(i,j,k-1))/dx_gpu[2];
            w_y = (velz(i,j,k) - velz(i,j-1,k))/dx_gpu[1];
            tauyz(i,j,k) = 0.25*(eta(i,j-1,k-1)+eta(i,j-1,k)+eta(i,j,k-1)+eta(i,j,k))*(v_z+w_y);

        });

        // Loop over faces for flux calculations (4:5+ns)
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            GpuArray<Real,MAX_SPECIES> meanXk;
            GpuArray<Real,MAX_SPECIES> meanYk;
            GpuArray<Real,MAX_SPECIES> dk;
            GpuArray<Real,MAX_SPECIES> Fk;
            GpuArray<Real,MAX_SPECIES> hk;
            GpuArray<Real,MAX_SPECIES> soret;

            xflux(i,j,k,4) += 0.5*velx(i,j,k)*(tauxx(i-1,j,k)+tauxx(i,j,k));
            xflux(i,j,k,4) += 0.25*((vely(i,j+1,k)+vely(i-1,j+1,k))*tauxy(i,j+1,k) + (vely(i,j,k)+vely(i-1,j,k))*tauxy(i,j,k));
            xflux(i,j,k,4) += 0.25*((velz(i,j,k+1)+velz(i-1,j,k+1))*tauxz(i,j,k+1) + (velz(i,j,k)+velz(i-1,j,k))*tauxz(i,j,k));

            Real meanT = 0.5*(prim(i-1,j,k,4)+prim(i,j,k,4));
            Real meanP = 0.5*(prim(i-1,j,k,5)+prim(i,j,k,5));

            if (algorithm_type_gpu == 2) {

                // compute dk
                for (int ns=0; ns<nspecies_gpu; ++ns) {
                    Real term1 = (prim(i,j,k,6+nspecies_gpu+ns)-prim(i-1,j,k,6+nspecies_gpu+ns))/dx_gpu[0];
                    meanXk[ns] = 0.5*(prim(i-1,j,k,6+nspecies_gpu+ns)+prim(i,j,k,6+nspecies_gpu+ns));
                    meanYk[ns] = 0.5*(prim(i-1,j,k,6+ns)+prim(i,j,k,6+ns));
                    Real term2 = (meanXk[ns]-meanYk[ns])*(prim(i,j,k,5)-prim(i-1,j,k,5))/dx_gpu[0]/meanP;
                    dk[ns] = term1 + term2;
                    soret[ns] = 0.5*(chi(i-1,j,k,ns)*prim(i-1,j,k,6+nspecies_gpu+ns)+chi(i,j,k,ns)*prim(i,j,k,6+nspecies_gpu+ns))
                        *(prim(i,j,k,4)-prim(i-1,j,k,4))/dx_gpu[0]/meanT;
                }

                // compute Fk (based on Eqn. 2.5.24, Giovangigli's book)
                for (int kk=0; kk<nspecies_gpu; ++kk) {
                    Fk[kk] = 0.;
                    for (int ll=0; ll<nspecies_gpu; ++ll) {
                        Fk[kk] -= half*(Dij(i-1,j,k,ll*nspecies_gpu+kk)+Dij(i,j,k,ll*nspecies_gpu+kk))*( dk[ll] +soret[ll]);
                    }
                }

                // compute Q (based on Eqn. 2.5.25, Giovangigli's book)
                GetEnthalpies(meanT,hk,hcp_gpu,nspecies_gpu);

                Real Q5 = 0.;
                for (int ns=0; ns<nspecies_gpu; ++ns) {
                    Q5 += (hk[ns] + 0.5 * Runiv_gpu*meanT*(chi(i-1,j,k,ns)+chi(i,j,k,ns))/molmass_gpu[ns])*Fk[ns];
                }
                // heat conduction already included in flux(5)       

                xflux(i,j,k,4) += Q5;

                for (int ns=0; ns<nspecies_gpu; ++ns) {
                    xflux(i,j,k,5+ns) += Fk[ns];
                }
            }
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            
            GpuArray<Real,MAX_SPECIES> meanXk;
            GpuArray<Real,MAX_SPECIES> meanYk;
            GpuArray<Real,MAX_SPECIES> dk;
            GpuArray<Real,MAX_SPECIES> Fk;
            GpuArray<Real,MAX_SPECIES> hk;
            GpuArray<Real,MAX_SPECIES> soret;

            yflux(i,j,k,4) += 0.25*((velx(i+1,j,k)+velx(i+1,j-1,k))*tauxy(i+1,j,k) + (velx(i,j,k)+velx(i,j-1,k))*tauxy(i,j,k));
            yflux(i,j,k,4) += 0.5*vely(i,j,k)*(tauyy(i,j-1,k)+tauyy(i,j,k));
            yflux(i,j,k,4) += 0.25*((velz(i,j,k+1)+velz(i,j-1,k+1))*tauyz(i,j,k+1) + (velz(i,j,k)+velz(i,j-1,k))*tauyz(i,j,k));

            Real meanT = 0.5*(prim(i,j-1,k,4)+prim(i,j,k,4));
            Real meanP = 0.5*(prim(i,j-1,k,5)+prim(i,j,k,5));

            if (algorithm_type_gpu == 2) {
                // compute dk
                for (int ns=0; ns<nspecies_gpu; ++ns) {
                    Real term1 = (prim(i,j,k,6+nspecies_gpu+ns)-prim(i,j-1,k,6+nspecies_gpu+ns))/dx_gpu[1];
                    meanXk[ns] = 0.5*(prim(i,j-1,k,6+nspecies_gpu+ns)+prim(i,j,k,6+nspecies_gpu+ns));
                    meanYk[ns] = 0.5*(prim(i,j-1,k,6+ns)+prim(i,j,k,6+ns));
                    Real term2 = (meanXk[ns]-meanYk[ns])*(prim(i,j,k,5)-prim(i,j-1,k,5))/dx_gpu[1]/meanP;
                    dk[ns] = term1 + term2;
                    soret[ns] = 0.5*(chi(i,j-1,k,ns)*prim(i,j-1,k,6+nspecies_gpu+ns)+chi(i,j,k,ns)*prim(i,j,k,6+nspecies_gpu+ns))
                        *(prim(i,j,k,4)-prim(i,j-1,k,4))/dx_gpu[1]/meanT;
                }

                // compute Fk (based on Eqn. 2.5.24, Giovangigli's book)
                for (int kk=0; kk<nspecies_gpu; ++kk) {
                    Fk[kk] = 0.;
                    for (int ll=0; ll<nspecies_gpu; ++ll) {
                        Fk[kk] -= half*(Dij(i,j-1,k,ll*nspecies_gpu+kk)+Dij(i,j,k,ll*nspecies_gpu+kk))*( dk[ll] +soret[ll]);
                    }
                }

                // compute Q (based on Eqn. 2.5.25, Giovangigli's book)
                GetEnthalpies(meanT,hk,hcp_gpu,nspecies_gpu);

                Real Q5 = 0.0;
                for (int ns=0; ns<nspecies_gpu; ++ns) {
                    Q5 += (hk[ns] + 0.5 * Runiv_gpu*meanT*(chi(i,j-1,k,ns)+chi(i,j,k,ns))/molmass_gpu[ns])*Fk[ns];
                }

                // heat conduction already included in flux(5)

                yflux(i,j,k,4) += Q5;

                for (int ns=0; ns<nspecies_gpu; ++ns) {
                    yflux(i,j,k,5+ns) += Fk[ns];
                }
            }
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            if (n_cells_z > 1) {
            
                GpuArray<Real,MAX_SPECIES> meanXk;
                GpuArray<Real,MAX_SPECIES> meanYk;
                GpuArray<Real,MAX_SPECIES> dk;
                GpuArray<Real,MAX_SPECIES> Fk;
                GpuArray<Real,MAX_SPECIES> hk;
                GpuArray<Real,MAX_SPECIES> soret;
                    
                zflux(i,j,k,4) += 0.25*((velx(i+1,j,k-1)+velx(i+1,j,k))*tauxz(i+1,j,k) + (velx(i,j,k)+velx(i,j,k-1))*tauxz(i,j,k));
                zflux(i,j,k,4) += 0.25*((vely(i,j+1,k-1)+vely(i,j+1,k))*tauyz(i,j+1,k) + (vely(i,j,k)+vely(i,j,k-1))*tauyz(i,j,k));
                zflux(i,j,k,4) += 0.5*velz(i,j,k)*(tauzz(i,j,k-1)+tauzz(i,j,k));

                Real meanT = 0.5*(prim(i,j,k-1,4)+prim(i,j,k,4));
                Real meanP = 0.5*(prim(i,j,k-1,5)+prim(i,j,k,5));

                if (algorithm_type_gpu == 2) {

                    // compute dk
                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        Real term1 = (prim(i,j,k,6+nspecies_gpu+ns)-prim(i,j,k-1,6+nspecies_gpu+ns))/dx_gpu[2];
                        meanXk[ns] = 0.5*(prim(i,j,k-1,6+nspecies_gpu+ns)+prim(i,j,k,6+nspecies_gpu+ns));
                        meanYk[ns] = 0.5*(prim(i,j,k-1,6+ns)+prim(i,j,k,6+ns));
                        Real term2 = (meanXk[ns]-meanYk[ns])*(prim(i,j,k,5)-prim(i,j,k-1,5))/dx_gpu[2]/meanP;
                        dk[ns] = term1 + term2;
                        soret[ns] = 0.5*(chi(i,j,k,ns)*prim(i,j,k-1,6+nspecies_gpu+ns)+chi(i,j,k+1,ns)*prim(i,j,k,6+nspecies_gpu+ns))
                            *(prim(i,j,k,4)-prim(i,j,k-1,4))/dx_gpu[2]/meanT;
                    }

                    // compute Fk (based on Eqn. 2.5.24, Giovangigli's book)
                    for (int kk=0; kk<nspecies_gpu; ++kk) {
                        Fk[kk] = 0.;
                        for (int ll=0; ll<nspecies_gpu; ++ll) {
                            Fk[kk] -= half*(Dij(i,j,k-1,ll*nspecies_gpu+kk)+Dij(i,j,k,ll*nspecies_gpu+kk))*( dk[ll] +soret[ll]);
                        }
                    }

                    // compute Q (based on Eqn. 2.5.25, Giovangigli's book)
                    GetEnthalpies(meanT,hk,hcp_gpu,nspecies_gpu);

                    Real Q5 = 0.0;
                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        Q5 += (hk[ns] + 0.5 * Runiv_gpu*meanT*(chi(i,j,k,ns)+chi(i,j,k,ns))/molmass_gpu[ns])*Fk[ns];
                    }

                    // heat conduction already included in flux(5)
                    zflux(i,j,k,4) += Q5;

                    for (int ns=0; ns<nspecies_gpu; ++ns) {
                        zflux(i,j,k,5+ns) += Fk[ns];
                    }
                }
            
            } // n_cells_z test
        });

        // Loop over edges for momemntum flux calculations [1:3]
        amrex::ParallelFor(bx_xy, bx_xz, bx_yz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            edgey_u(i,j,k) += tauxy(i,j,k);
            edgex_v(i,j,k) += tauxy(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            edgez_u(i,j,k) += tauxz(i,j,k);
            edgex_w(i,j,k) += tauxz(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            edgez_v(i,j,k) += tauyz(i,j,k);
            edgey_w(i,j,k) += tauyz(i,j,k);
        });
        
        // Loop over the center cells and compute fluxes (diagonal momentum terms)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            cenx_u(i,j,k) += tauxx(i,j,k);
            ceny_v(i,j,k) += tauyy(i,j,k);
            cenz_w(i,j,k) += tauzz(i,j,k);
        });
    }

    ////////////////////
    // hyperbolic fluxes
    ////////////////////

    Real wgt2 = 1./12.;
    Real wgt1 = 0.5 + wgt2;

    // Loop over boxes
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {

        AMREX_D_TERM(const Array4<Real>& xflux = faceflux_in[0].array(mfi); ,
                     const Array4<Real>& yflux = faceflux_in[1].array(mfi); ,
                     const Array4<Real>& zflux = faceflux_in[2].array(mfi));

        const Array4<Real>& edgex_v = edgeflux_x_in[0].array(mfi);
        const Array4<Real>& edgex_w = edgeflux_x_in[1].array(mfi);
        const Array4<Real>& edgey_u = edgeflux_y_in[0].array(mfi);
        const Array4<Real>& edgey_w = edgeflux_y_in[1].array(mfi);
        const Array4<Real>& edgez_u = edgeflux_z_in[0].array(mfi);
        const Array4<Real>& edgez_v = edgeflux_z_in[1].array(mfi);

        const Array4<Real>& cenx_u = cenflux_in[0].array(mfi);
        const Array4<Real>& ceny_v = cenflux_in[1].array(mfi);
        const Array4<Real>& cenz_w = cenflux_in[2].array(mfi);

        AMREX_D_TERM(Array4<Real const> const& momx = cumom_in[0].array(mfi);,
                     Array4<Real const> const& momy = cumom_in[1].array(mfi);,
                     Array4<Real const> const& momz = cumom_in[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& velx = facevel_in[0].array(mfi);,
                     Array4<Real const> const& vely = facevel_in[1].array(mfi);,
                     Array4<Real const> const& velz = facevel_in[2].array(mfi););

        const Array4<const Real> prim = prim_in.array(mfi);
        const Array4<const Real> cons = cons_in.array(mfi);

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Box & bx_xy = mfi.tilebox(nodal_flag_xy);
        #if (AMREX_SPACEDIM == 3)
        const Box & bx_xz = mfi.tilebox(nodal_flag_xz);
        const Box & bx_yz = mfi.tilebox(nodal_flag_yz);
        #endif

        const Box& bx = mfi.tilebox();

        if (advection_type == 1) { // interpolate primitive quantities (fix this later for staggered grid -- Ishan)
            
            // Loop over the cells and compute fluxes
            amrex::ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            
                GpuArray<Real,MAX_SPECIES+5> conserved;
                GpuArray<Real,MAX_SPECIES+6> primitive;
                GpuArray<Real,MAX_SPECIES  > Yk;
                    
                for (int l=0; l<nprimvars_gpu; ++l) {
                    primitive[l] = wgt1*(prim(i,j,k,l)+prim(i-1,j,k,l)) - wgt2*(prim(i-2,j,k,l)+prim(i+1,j,k,l));
                }

                Real temp = primitive[4];
                Real rho = primitive[0];
                conserved[0] = rho;

                // want sum of specden == rho
                for (int n=0; n<nspecies_gpu; ++n) {
                    Yk[n] = primitive[6+n];
                }

                Real intenergy;
                GetEnergy(intenergy, Yk, temp, hcv_gpu, nspecies_gpu);

                Real vsqr = primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3];

                conserved[4] = rho*intenergy + 0.5*rho*vsqr;

                xflux(i,j,k,0) += conserved[0]*primitive[1];
                xflux(i,j,k,1) += conserved[0]*(primitive[1]*primitive[1])+primitive[5];
                xflux(i,j,k,2) += conserved[0]*primitive[1]*primitive[2];
                xflux(i,j,k,3) += conserved[0]*primitive[1]*primitive[3];

                xflux(i,j,k,4) += primitive[1]*conserved[4] + primitive[5]*primitive[1];

                if (algorithm_type_gpu == 2) { // Add advection of concentration
                    for (int n=0; n<nspecies_gpu; ++n) {
                        xflux(i,j,k,5+n) += rho*primitive[6+n]*primitive[1];
                    }
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            
                GpuArray<Real,MAX_SPECIES+5> conserved;
                GpuArray<Real,MAX_SPECIES+6> primitive;
                GpuArray<Real,MAX_SPECIES  > Yk;
                    
                for (int l=0; l<nprimvars_gpu; ++l) {
                    primitive[l] = wgt1*(prim(i,j,k,l)+prim(i,j-1,k,l)) - wgt2*(prim(i,j-2,k,l)+prim(i,j+1,k,l));
                }

                Real temp = primitive[4];
                Real rho = primitive[0];
                conserved[0] = rho;

                // want sum of specden == rho
                for (int n=0; n<nspecies_gpu; ++n) {
                    Yk[n] = primitive[6+n];
                }

                Real intenergy;
                GetEnergy(intenergy, Yk, temp, hcv_gpu, nspecies_gpu);

                Real vsqr = primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3];

                conserved[4] = rho*intenergy + 0.5*rho*vsqr;

                yflux(i,j,k,0) += conserved[0]*primitive[2];
                yflux(i,j,k,1) += conserved[0]*primitive[1]*primitive[2];
                yflux(i,j,k,2) += conserved[0]*primitive[2]*primitive[2]+primitive[5];
                yflux(i,j,k,3) += conserved[0]*primitive[3]*primitive[2];

                yflux(i,j,k,4) += primitive[2]*conserved[4] + primitive[5]*primitive[2];

                if (algorithm_type_gpu == 2) { // Add advection of concentration
                    for (int n=0; n<nspecies_gpu; ++n) {
                        yflux(i,j,k,5+n) += rho*primitive[6+n]*primitive[2];
                    }
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            
                GpuArray<Real,MAX_SPECIES+5> conserved;
                GpuArray<Real,MAX_SPECIES+6> primitive;
                GpuArray<Real,MAX_SPECIES  > Yk;
                    
                for (int l=0; l<nprimvars_gpu; ++l) {
                    primitive[l] = wgt1*(prim(i,j,k,l)+prim(i,j,k-1,l)) - wgt2*(prim(i,j,k-2,l)+prim(i,j,k+1,l));
                }

                Real temp = primitive[4];
                Real rho = primitive[0];
                conserved[0] = rho;

                // want sum of specden == rho
                for (int n=0; n<nspecies_gpu; ++n) {
                    Yk[n] = primitive[6+n];
                }

                Real intenergy;
                GetEnergy(intenergy, Yk, temp, hcv_gpu, nspecies_gpu);

                Real vsqr = primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3];

                conserved[4] = rho*intenergy + 0.5*rho*vsqr;

                zflux(i,j,k,0) += conserved[0]*primitive[3];
                zflux(i,j,k,1) += conserved[0]*primitive[1]*primitive[3];
                zflux(i,j,k,2) += conserved[0]*primitive[2]*primitive[3];
                zflux(i,j,k,3) += conserved[0]*primitive[3]*primitive[3]+primitive[5];

                zflux(i,j,k,4) += primitive[3]*conserved[4] + primitive[5]*primitive[3];

                if (algorithm_type_gpu == 2) { // Add advection of concentration
                    for (int n=0; n<nspecies_gpu; ++n) {
                        zflux(i,j,k,5+n) += rho*primitive[6+n]*primitive[3];
                    }
                }

            });
            
        } else if (advection_type == 2) { // interpolate conserved quantitites

            // 1. Loop over the face cells and compute fluxes (all conserved qtys. except momentum; i.e.,[0,4,5-nspecies])
            amrex::ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            
                xflux(i,j,k,0) += momx(i,j,k);
                xflux(i,j,k,4) += 0.5*(cons(i-1,j,k,4)+cons(i,j,k,4))*velx(i,j,k) + 0.25*(prim(i-1,j,k,5)+prim(i,j,k,5))*velx(i,j,k)*(cons(i-1,j,k,0)+cons(i,j,k,0));

                if (algorithm_type_gpu == 2) { // Add advection of concentration
                    for (int n=0; n<nspecies_gpu; ++n) {
                        xflux(i,j,k,5+n) += 0.5*(cons(i-1,j,k,5+n)+cons(i,j,k,5+n))*velx(i,j,k);
                    }
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            
                yflux(i,j,k,0) += momy(i,j,k);
                yflux(i,j,k,4) += 0.5*(cons(i,j-1,k,4)+cons(i,j,k,4))*vely(i,j,k) + 0.25*(prim(i,j-1,k,5)+prim(i,j,k,5))*vely(i,j,k)*(cons(i,j-1,k,0)+cons(i,j,k,0));

                if (algorithm_type_gpu == 2) { // Add advection of concentration
                    for (int n=0; n<nspecies_gpu; ++n) {
                        yflux(i,j,k,5+n) += 0.5*(cons(i,j-1,k,5+n)+cons(i,j,k,5+n))*vely(i,j,k);
                    }
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                zflux(i,j,k,0) += momz(i,j,k);
                zflux(i,j,k,4) += 0.5*(cons(i,j,k-1,4)+cons(i,j,k,4))*velz(i,j,k) + 0.25*(prim(i,j,k-1,5)+prim(i,j,k,5))*velz(i,j,k)*(cons(i,j,k-1,0)+cons(i,j,k,0));

                if (algorithm_type_gpu == 2) { // Add advection of concentration
                    for (int n=0; n<nspecies_gpu; ++n) {
                        zflux(i,j,k,5+n) += 0.5*(cons(i,j,k-1,5+n)+cons(i,j,k,5+n))*velz(i,j,k);
                    }
                }
            });

            // 2. Loop over the edge cells and compute fluxes (off-diagonal momentum terms)
            amrex::ParallelFor(bx_xy, bx_xz, bx_yz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                edgey_u(i,j,k) += 0.25*(momx(i,j-1,k)+momx(i,j,k))*(vely(i-1,j,k)+vely(i,j,k));
                edgex_v(i,j,k) += 0.25*(momy(i-1,j,k)+momy(i,j,k))*(velx(i,j-1,k)+velx(i,j,k));
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                edgez_u(i,j,k) += 0.25*(momx(i,j,k-1)+momx(i,j,k))*(velz(i-1,j,k)+velz(i,j,k));
                edgex_w(i,j,k) += 0.25*(momz(i-1,j,k)+momz(i,j,k))*(velx(i,j,k-1)+velx(i,j,k));
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                edgez_v(i,j,k) += 0.25*(momy(i,j,k-1)+momy(i,j,k))*(velz(i,j-1,k)+velz(i,j,k));
                edgey_w(i,j,k) += 0.25*(momz(i,j-1,k)+momz(i,j,k))*(vely(i,j,k-1)+vely(i,j,k));
            });

            // 3. Loop over the center cells and compute fluxes (diagonal momentum terms)
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                cenx_u(i,j,k) += 0.25*(momx(i,j,k)+momx(i+1,j,k))*(velx(i,j,k)+velx(i+1,j,k)) + prim(i,j,k,5);
                ceny_v(i,j,k) += 0.25*(momy(i,j,k)+momy(i,j+1,k))*(vely(i,j,k)+vely(i,j+1,k)) + prim(i,j,k,5);
                cenz_w(i,j,k) += 0.25*(momz(i,j,k)+momz(i,j,k+1))*(velz(i,j,k)+velz(i,j,k+1)) + prim(i,j,k,5);
            });

        }
    }
}
