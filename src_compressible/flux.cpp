#include "compressible_functions.H"
#include "common_functions.H"

void calculateFlux(const MultiFab& cons_in, const MultiFab& prim_in,
                   const MultiFab& eta_in, const MultiFab& zeta_in, const MultiFab& kappa_in,
                   const MultiFab& chi_in, const MultiFab& D_in,
                   std::array<MultiFab, AMREX_SPACEDIM>& flux_in,
                   std::array<MultiFab, AMREX_SPACEDIM>& stochFlux_in,
                   std::array<MultiFab, AMREX_SPACEDIM>& cornx_in,
                   std::array<MultiFab, AMREX_SPACEDIM>& corny_in,
                   std::array<MultiFab, AMREX_SPACEDIM>& cornz_in,
                   MultiFab& visccorn_in,
                   MultiFab& rancorn_in,
                   const amrex::Geometry& geom,
                   const amrex::Vector< amrex::Real >& /*stoch_weights*/,
                   const amrex::Real dt)
{
    BL_PROFILE_VAR("calculateFlux()",calculateFlux);

    int n_cells_z = n_cells[2];

    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    AMREX_D_TERM(flux_in[0].setVal(0);,
                 flux_in[1].setVal(0);,
                 flux_in[2].setVal(0););

    // ignore for reservoirs and periodic BC
    bool is_lo_x_dirichlet_mass = (bc_mass_lo[0] != 3) and (bc_mass_lo[0] != -1);
    bool is_hi_x_dirichlet_mass = (bc_mass_hi[0] != 3) and (bc_mass_hi[0] != -1);
    bool is_lo_y_dirichlet_mass = (bc_mass_lo[1] != 3) and (bc_mass_lo[1] != -1);
    bool is_hi_y_dirichlet_mass = (bc_mass_hi[1] != 3) and (bc_mass_hi[1] != -1);
    bool is_lo_z_dirichlet_mass = (bc_mass_lo[2] != 3) and (bc_mass_lo[2] != -1);
    bool is_hi_z_dirichlet_mass = (bc_mass_hi[2] != 3) and (bc_mass_hi[2] != -1);

    ////////////////////
    // stochastic fluxes
    ////////////////////

    if (stoch_stress_form == 1) {

        Real volinv = 1./(dx[0]*dx[1]*dx[2]);
        Real dtinv = 1./dt;

        // Loop over boxes
        for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {

            AMREX_D_TERM(const Array4<Real>& fluxx = flux_in[0].array(mfi); ,
                         const Array4<Real>& fluxy = flux_in[1].array(mfi); ,
                         const Array4<Real>& fluxz = flux_in[2].array(mfi));

            AMREX_D_TERM(const Array4<Real>& ranfluxx = stochFlux_in[0].array(mfi); ,
                         const Array4<Real>& ranfluxy = stochFlux_in[1].array(mfi); ,
                         const Array4<Real>& ranfluxz = stochFlux_in[2].array(mfi));

            const Array4<const Real> prim = prim_in.array(mfi);
            // const Array4<const Real> cons = cons_in.array(mfi);

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
                GpuArray<Real,MAX_SPECIES+5> wiener;

                GpuArray<Real,MAX_SPECIES> hk;
                GpuArray<Real,MAX_SPECIES> yy;
                GpuArray<Real,MAX_SPECIES> yyp;

                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> DijY_edge;
                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> sqD;

                Real muxp = (eta(i,j,k)*prim(i,j,k,4) + eta(i-1,j,k)*prim(i-1,j,k,4));
                Real kxp = (kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4) + kappa(i-1,j,k)*prim(i-1,j,k,4)*prim(i-1,j,k,4));

                Real meanT = 0.5*(prim(i,j,k,4)+prim(i-1,j,k,4));

                if ((i == 0) and is_lo_x_dirichlet_mass) {
                    muxp = 2.0*eta(i-1,j,k)*prim(i-1,j,k,4);
                    kxp  = 2.0*kappa(i-1,j,k)*prim(i-1,j,k,4)*prim(i-1,j,k,4);
                    meanT = prim(i-1,j,k,4);
                }
                if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                    muxp = 2.0*eta(i,j,k)*prim(i,j,k,4);
                    kxp  = 2.0*kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4);
                    meanT = prim(i,j,k,4);
                }

                // Weights for facial fluxes:
                fweights[0] = 0; // No mass flux;
                fweights[1]=sqrt(k_B*muxp*volinv*dtinv);
                fweights[2]=fweights[1];
                fweights[3]=fweights[1];
                fweights[4]=sqrt(k_B*kxp*volinv*dtinv);

                // Construct the random increments
                for (int n=0; n<5; ++n) {
                    wiener[n] = fweights[n]*ranfluxx(i,j,k,n);
                }

                Real nweight=sqrt(k_B*volinv*dtinv);

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

                    if ((i == 0) and is_lo_x_dirichlet_mass) {
                        muzepp = 0.5*(eta(i-1,j,k)*prim(i-1,j,k,4) +
                                      eta(i-1,j+1,k)*prim(i-1,j+1,k,4) +
                                      eta(i-1,j,k+1)*prim(i-1,j,k+1,4) +
                                      eta(i-1,j+1,k+1)*prim(i-1,j+1,k+1,4) )/3.;

                        muzemp = 0.5*(eta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                      eta(i-1,j,k)*prim(i-1,j,k,4) +
                                      eta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,4) +
                                      eta(i-1,j,k+1)*prim(i-1,j,k+1,4) )/3.;

                        muzepm = 0.5*(eta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                      eta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,4) +
                                      eta(i-1,j,k)*prim(i-1,j,k,4) +
                                      eta(i-1,j+1,k)*prim(i-1,j+1,k,4) )/3.;

                        muzemm = 0.5*(eta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) +
                                      eta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                      eta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                      eta(i-1,j,k)*prim(i-1,j,k,4) )/3.;
                    }
                    if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                        muzepp = 0.5*(eta(i,j,k)*prim(i,j,k,4) +
                                      eta(i,j+1,k)*prim(i,j+1,k,4) +
                                      eta(i,j,k+1)*prim(i,j,k+1,4) +
                                      eta(i,j+1,k+1)*prim(i,j+1,k+1,4) )/3.;

                        muzemp = 0.5*(eta(i,j-1,k)*prim(i,j-1,k,4) +
                                      eta(i,j,k)*prim(i,j,k,4) +
                                      eta(i,j-1,k+1)*prim(i,j-1,k+1,4) +
                                      eta(i,j,k+1)*prim(i,j,k+1,4) )/3.;

                        muzepm = 0.5*(eta(i,j,k-1)*prim(i,j,k-1,4) +
                                      eta(i,j+1,k-1)*prim(i,j+1,k-1,4) +
                                      eta(i,j,k)*prim(i,j,k,4) +
                                      eta(i,j+1,k)*prim(i,j+1,k,4) )/3.;

                        muzemm = 0.5*(eta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                      eta(i,j,k-1)*prim(i,j,k-1,4) +
                                      eta(i,j-1,k)*prim(i,j-1,k,4) +
                                      eta(i,j,k)*prim(i,j,k,4) )/3.;
                    }



                    if (amrex::Math::abs(visc_type) == 3) {
                        if ((i == 0) and is_lo_x_dirichlet_mass) {
                            muzepp += 0.5*(zeta(i-1,j,k)*prim(i-1,j,k,4) +
                                          zeta(i-1,j+1,k)*prim(i-1,j+1,k,4) +
                                          zeta(i-1,j,k+1)*prim(i-1,j,k+1,4) +
                                          zeta(i-1,j+1,k+1)*prim(i-1,j+1,k+1,4) );

                            muzemp += 0.5*(zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                          zeta(i-1,j,k)*prim(i-1,j,k,4) +
                                          zeta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,4) +
                                          zeta(i-1,j,k+1)*prim(i-1,j,k+1,4) );

                            muzepm += 0.5*(zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                          zeta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,4) +
                                          zeta(i-1,j,k)*prim(i-1,j,k,4) +
                                          zeta(i-1,j+1,k)*prim(i-1,j+1,k,4) );

                            muzemm += 0.5*(zeta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) +
                                          zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                          zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                          zeta(i-1,j,k)*prim(i-1,j,k,4) );
                        }
                        else if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                            muzepp += 0.5*(zeta(i,j,k)*prim(i,j,k,4) +
                                          zeta(i,j+1,k)*prim(i,j+1,k,4) +
                                          zeta(i,j,k+1)*prim(i,j,k+1,4) +
                                          zeta(i,j+1,k+1)*prim(i,j+1,k+1,4) );

                            muzemp += 0.5*(zeta(i,j-1,k)*prim(i,j-1,k,4) +
                                          zeta(i,j,k)*prim(i,j,k,4) +
                                          zeta(i,j-1,k+1)*prim(i,j-1,k+1,4) +
                                          zeta(i,j,k+1)*prim(i,j,k+1,4) );

                            muzepm += 0.5*(zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                          zeta(i,j+1,k-1)*prim(i,j+1,k-1,4) +
                                          zeta(i,j,k)*prim(i,j,k,4) +
                                          zeta(i,j+1,k)*prim(i,j+1,k,4) );

                            muzemm += 0.5*(zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                          zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                          zeta(i,j-1,k)*prim(i,j-1,k,4) +
                                          zeta(i,j,k)*prim(i,j,k,4) );
                        }
                        else {
                            muzepp += 0.25*(zeta(i,j,k)*prim(i,j,k,4) + zeta(i-1,j,k)*prim(i-1,j,k,4) +
                                           zeta(i,j+1,k)*prim(i,j+1,k,4) + zeta(i-1,j+1,k)*prim(i-1,j+1,k,4) +
                                           zeta(i,j,k+1)*prim(i,j,k+1,4) + zeta(i-1,j,k+1)*prim(i-1,j,k+1,4) +
                                           zeta(i,j+1,k+1)*prim(i,j+1,k+1,4) + zeta(i-1,j+1,k+1)*prim(i-1,j+1,k+1,4) );

                            muzemp += 0.25*(zeta(i,j-1,k)*prim(i,j-1,k,4) + zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                           zeta(i,j,k)*prim(i,j,k,4) + zeta(i-1,j,k)*prim(i-1,j,k,4) +
                                           zeta(i,j-1,k+1)*prim(i,j-1,k+1,4) + zeta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,4) +
                                           zeta(i,j,k+1)*prim(i,j,k+1,4) + zeta(i-1,j,k+1)*prim(i-1,j,k+1,4) );

                            muzepm += 0.25*(zeta(i,j,k-1)*prim(i,j,k-1,4) + zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                           zeta(i,j+1,k-1)*prim(i,j+1,k-1,4) + zeta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,4) +
                                           zeta(i,j,k)*prim(i,j,k,4) + zeta(i-1,j,k)*prim(i-1,j,k,4) +
                                           zeta(i,j+1,k)*prim(i,j+1,k,4) + zeta(i-1,j+1,k)*prim(i-1,j+1,k,4) );

                            muzemm += 0.25*(zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) + zeta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) +
                                            zeta(i,j,k-1)*prim(i,j,k-1,4) + zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                            zeta(i,j-1,k)*prim(i,j-1,k,4) + zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                            zeta(i,j,k)*prim(i,j,k,4) + zeta(i-1,j,k)*prim(i-1,j,k,4) );
                        }
                    }

                    Real factor_lo_y = 1.;
                    Real factor_hi_y = 1.;
                    Real factor_lo_z = 1.;
                    Real factor_hi_z = 1.;

                    // 1 = slip
                    // 2 = no-slip
                    if (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) {
                        if (j == 0) {
                            factor_lo_y = (bc_vel_lo[1] == 1) ? std::sqrt(2.0) : 0.;
                        }
                    }
                    if (bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) {
                        if (k == 0) {
                            factor_lo_z = (bc_vel_lo[2] == 1) ? std::sqrt(2.0) : 0.;
                        }
                    }
                    if (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) {
                        if (j == n_cells[1]-1) {
                            factor_hi_y = (bc_vel_hi[1] == 1) ? std::sqrt(2.0) : 0.;
                        }
                    }
                    if (bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) {
                        if (k == n_cells[2]-1) {
                            factor_hi_z = (bc_vel_hi[2] == 1) ? std::sqrt(2.0) : 0.;
                        }
                    }

                    // Random "divergence" stress
                    wiener[1] = wiener[1] + 0.25*nweight*(factor_hi_y*factor_hi_z*sqrt(muzepp)*rancorn(i,j+1,k+1) +
                                                          factor_lo_y*factor_hi_z*sqrt(muzemp)*rancorn(i,j,k+1) +
                                                          factor_hi_y*factor_lo_z*sqrt(muzepm)*rancorn(i,j+1,k) +
                                                          factor_lo_y*factor_lo_z*sqrt(muzemm)*rancorn(i,j,k));

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

          wiener(2) = wiener(2) + 0.5*nweight*(sqrt(muzepp)*rancorn(i,j+1,k)+ &
               sqrt(muzemp)*rancorn(i,j,k)) ! Random "divergence" stress
*/
                }

                // increment momentum fluxes
                for (int n=1; n<4; ++n) {
                    fluxx(i,j,k,n) = fluxx(i,j,k,n) + wiener[n];
                }

                // heat flux
                fluxx(i,j,k,nvars) = fluxx(i,j,k,nvars) + wiener[4];

                // Viscous heating:
                Real phiflxdiag =  wiener[1]*(prim(i-1,j,k,1)+prim(i,j,k,1));

                Real phiflxshear = wiener[2]*(prim(i-1,j,k,2)+prim(i,j,k,2)) +
                                   wiener[3]*(prim(i-1,j,k,3)+prim(i,j,k,3));

                if ((i == 0) and is_lo_x_dirichlet_mass) {
                    phiflxdiag =  2.0*wiener[1]*prim(i-1,j,k,1);
                    phiflxshear = 2.0*wiener[2]*prim(i-1,j,k,2) +
                                  2.0*wiener[3]*prim(i-1,j,k,3);
                }
                if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                    phiflxdiag =  2.0*wiener[1]*prim(i,j,k,1);
                    phiflxshear = 2.0*wiener[2]*prim(i,j,k,2) +
                                  2.0*wiener[3]*prim(i,j,k,3);
                }

                phiflxdiag = -0.5*phiflxdiag;
                phiflxshear = -0.5*phiflxshear;

                fluxx(i,j,k,nvars+1) = fluxx(i,j,k,nvars+1) - phiflxdiag;
                fluxx(i,j,k,nvars+2) = fluxx(i,j,k,nvars+2) - phiflxshear;

                if (algorithm_type == 2) {

                    for (int n=5; n<5+nspecies; ++n) {
                        wiener[n] = 0.;
                    }

                    for (int ns=0; ns<nspecies; ++ns) {
                        yy[ns] = amrex::max(0.,amrex::min(1.,prim(i-1,j,k,6+ns)));
                        yyp[ns] = amrex::max(0.,amrex::min(1.,prim(i,j,k,6+ns)));
                        if ((i == 0) and is_lo_x_dirichlet_mass) {
                            yyp[ns] = amrex::max(0.,amrex::min(1.,prim(i-1,j,k,6+ns)));
                        }
                        if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                            yy[ns] = amrex::max(0.,amrex::min(1.,prim(i,j,k,6+ns)));
                        }
                    }

                    Real sumy = 0.;
                    Real sumyp = 0.;

                    for (int n=0; n<nspecies; ++n) {
                        sumy += yy[n];
                        sumyp += yyp[n];
                    }

                    for (int n=0; n<nspecies; ++n) {
                        yy[n] /= sumy;
                        yyp[n] /= sumyp;
                    }

                    Real MWmix = 0.;

                    for (int ns=0; ns<nspecies; ++ns) {

                        MWmix = MWmix + 0.5*(yy[ns]+yyp[ns])/molmass[ns];

                        for (int ll=0; ll<nspecies; ++ll) {
                            DijY_edge[ns*nspecies+ll] = 0.5*(Dij(i-1,j,k,ll*nspecies+ns)*yy[ll] +
                                                                 Dij(i,j,k,ll*nspecies+ns)*yyp[ll] +
                                                                (Dij(i-1,j,k,ns*nspecies+ll)*yy[ns] +
                                                                 Dij(i,j,k,ns*nspecies+ll)*yyp[ns] ));
                            if ((i == 0) and is_lo_x_dirichlet_mass) {
                                DijY_edge[ns*nspecies+ll] = 0.5*(Dij(i-1,j,k,ll*nspecies+ns)*yy[ll] +
                                                                     Dij(i-1,j,k,ll*nspecies+ns)*yyp[ll] +
                                                                    (Dij(i-1,j,k,ns*nspecies+ll)*yy[ns] +
                                                                     Dij(i-1,j,k,ns*nspecies+ll)*yyp[ns] ));
                            }
                            if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                                DijY_edge[ns*nspecies+ll] = 0.5*(Dij(i,j,k,ll*nspecies+ns)*yy[ll] +
                                                                     Dij(i,j,k,ll*nspecies+ns)*yyp[ll] +
                                                                    (Dij(i,j,k,ns*nspecies+ll)*yy[ns] +
                                                                     Dij(i,j,k,ns*nspecies+ll)*yyp[ns] ));
                            }

                        }
                    }

                    for (int ns=0; ns<nspecies; ++ns) {
                        if (amrex::Math::abs(yy[ns]) + amrex::Math::abs(yyp[ns]) <= 1.e-12) {
                            for (int n=0; n<nspecies; ++n) {
                                DijY_edge[ns*nspecies+n]=0.;
                                DijY_edge[n*nspecies+ns]=0.;
                            }
                        }
                    }

                    MWmix = 1. / MWmix;

                    CholeskyDecomp(DijY_edge,nspecies,sqD);

                    for (int ns=0; ns<nspecies; ++ns) {
                        for (int ll=0; ll<=ns; ++ll) {
                            fweights[5+ll]=sqrt(k_B*MWmix*volinv/(Runiv*dt))*sqD[ns*nspecies+ll];
                            wiener[5+ns] = wiener[5+ns] + fweights[5+ll]*ranfluxx(i,j,k,5+ll);
                        }
                        fluxx(i,j,k,5+ns) = wiener[5+ns];
                    }

                    GetEnthalpies(meanT, hk);

                    Real soret = 0.;

                    for (int ns=0; ns<nspecies; ++ns) {
                        Real soret_s;
                        soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*0.5*(chi(i-1,j,k,ns)+chi(i,j,k,ns)))*wiener[5+ns];
                        if ((i == 0) and is_lo_x_dirichlet_mass) {
                            soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*chi(i-1,j,k,ns))*wiener[5+ns];
                        }
                        if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                            soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*chi(i,j,k,ns))*wiener[5+ns];
                        }
                        soret = soret + soret_s;
                    }
                    fluxx(i,j,k,nvars+3) = fluxx(i,j,k,nvars+3) + soret;
                }

            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                GpuArray<Real,MAX_SPECIES+5> fweights;
                GpuArray<Real,MAX_SPECIES+5> wiener;

                GpuArray<Real,MAX_SPECIES> hk;
                GpuArray<Real,MAX_SPECIES> yy;
                GpuArray<Real,MAX_SPECIES> yyp;

                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> DijY_edge;
                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> sqD;


                Real muyp = eta(i,j,k)*prim(i,j,k,4) + eta(i,j-1,k)*prim(i,j-1,k,4);
                Real kyp = kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4) + kappa(i,j-1,k)*prim(i,j-1,k,4)*prim(i,j-1,k,4);

                Real meanT = 0.5*(prim(i,j,k,4)+prim(i,j-1,k,4));

                if ((j == 0) and is_lo_y_dirichlet_mass) {
                    muyp = 2.0*eta(i,j-1,k)*prim(i,j-1,k,4);
                    kyp  = 2.0*kappa(i,j-1,k)*prim(i,j-1,k,4)*prim(i,j-1,k,4);
                    meanT = prim(i,j-1,k,4);
                }
                if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                    muyp = 2.0*eta(i,j,k)*prim(i,j,k,4);
                    kyp  = 2.0*kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4);
                    meanT = prim(i,j,k,4);
                }

                // Weights for facial fluxes:
                fweights[0] = 0; // No mass flux
                fweights[1] = sqrt(k_B*muyp*volinv*dtinv);
                fweights[2] = sqrt(k_B*muyp*volinv*dtinv);
                fweights[3] = sqrt(k_B*muyp*volinv*dtinv);
                fweights[4] = sqrt(k_B*kyp*volinv*dtinv);

                // Construct the random increments
                for (int n=0; n<5; ++n) {
                    wiener[n] = fweights[n]*ranfluxy(i,j,k,n);
                }

                Real nweight=sqrt(k_B*volinv*dtinv);

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

                    if ((j == 0) and is_lo_y_dirichlet_mass) {
                        muzepp = 0.5*(eta(i+1,j-1,k)*prim(i+1,j-1,k,4) +
                                      eta(i,j-1,k)*prim(i,j-1,k,4) +
                                      eta(i+1,j-1,k+1)*prim(i+1,j-1,k+1,4) +
                                      eta(i,j-1,k+1)*prim(i,j-1,k+1,4) )/3.;

                        muzemp = 0.5*(eta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                      eta(i,j-1,k)*prim(i,j-1,k,4) +
                                      eta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,4) +
                                      eta(i,j-1,k+1)*prim(i,j-1,k+1,4) )/3.;

                        muzepm = 0.5*(eta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,4) +
                                      eta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                      eta(i+1,j-1,k)*prim(i+1,j-1,k,4) +
                                      eta(i,j-1,k)*prim(i,j-1,k,4) )/3.;

                        muzemm = 0.5*(eta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) +
                                      eta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                      eta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                      eta(i,j-1,k)*prim(i,j-1,k,4) )/3.;

                    }
                    if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                        muzepp = 0.5*(eta(i+1,j,k)*prim(i+1,j,k,4) +
                                      eta(i,j,k)*prim(i,j,k,4) +
                                      eta(i+1,j,k+1)*prim(i+1,j,k+1,4) +
                                      eta(i,j,k+1)*prim(i,j,k+1,4) )/3.;

                        muzemp = 0.5*(eta(i-1,j,k)*prim(i-1,j,k,4) +
                                      eta(i,j,k)*prim(i,j,k,4) +
                                      eta(i-1,j,k+1)*prim(i-1,j,k+1,4) +
                                      eta(i,j,k+1)*prim(i,j,k+1,4) )/3.;

                        muzepm = 0.5*(eta(i+1,j,k-1)*prim(i+1,j,k-1,4) +
                                      eta(i,j,k-1)*prim(i,j,k-1,4) +
                                      eta(i+1,j,k)*prim(i+1,j,k,4) +
                                      eta(i,j,k)*prim(i,j,k,4) )/3.;

                        muzemm = 0.5*(eta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                      eta(i,j,k-1)*prim(i,j,k-1,4) +
                                      eta(i-1,j,k)*prim(i-1,j,k,4) +
                                      eta(i,j,k)*prim(i,j,k,4) )/3.;

                    }

                    if (amrex::Math::abs(visc_type) == 3) {

                        if ((j == 0) and is_lo_y_dirichlet_mass) {
                            muzepp += 0.5*(zeta(i+1,j-1,k)*prim(i+1,j-1,k,4) +
                                          zeta(i,j-1,k)*prim(i,j-1,k,4) +
                                          zeta(i+1,j-1,k+1)*prim(i+1,j-1,k+1,4) +
                                          zeta(i,j-1,k+1)*prim(i,j-1,k+1,4) )/3.;

                            muzemp += 0.5*(zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                          zeta(i,j-1,k)*prim(i,j-1,k,4) +
                                          zeta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,4) +
                                          zeta(i,j-1,k+1)*prim(i,j-1,k+1,4) )/3.;

                            muzepm += 0.5*(zeta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,4) +
                                          zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                          zeta(i+1,j-1,k)*prim(i+1,j-1,k,4) +
                                          zeta(i,j-1,k)*prim(i,j-1,k,4) )/3.;

                            muzemm += 0.5*(zeta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) +
                                          zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                          zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                          zeta(i,j-1,k)*prim(i,j-1,k,4) )/3.;

                        }
                        else  if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                            muzepp += 0.5*(zeta(i+1,j,k)*prim(i+1,j,k,4) +
                                          zeta(i,j,k)*prim(i,j,k,4) +
                                          zeta(i+1,j,k+1)*prim(i+1,j,k+1,4) +
                                          zeta(i,j,k+1)*prim(i,j,k+1,4) )/3.;

                            muzemp += 0.5*(zeta(i-1,j,k)*prim(i-1,j,k,4) +
                                          zeta(i,j,k)*prim(i,j,k,4) +
                                          zeta(i-1,j,k+1)*prim(i-1,j,k+1,4) +
                                          zeta(i,j,k+1)*prim(i,j,k+1,4) )/3.;

                            muzepm += 0.5*(zeta(i+1,j,k-1)*prim(i+1,j,k-1,4) +
                                          zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                          zeta(i+1,j,k)*prim(i+1,j,k,4) +
                                          zeta(i,j,k)*prim(i,j,k,4) )/3.;

                            muzemm += 0.5*(zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                          zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                          zeta(i-1,j,k)*prim(i-1,j,k,4) +
                                          zeta(i,j,k)*prim(i,j,k,4) )/3.;

                        }
                        else {
                            muzepp += 0.25*(zeta(i+1,j-1,k)*prim(i+1,j-1,k,4) + zeta(i,j-1,k)*prim(i,j-1,k,4) +
                                           zeta(i+1,j,k)*prim(i+1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) +
                                           zeta(i+1,j-1,k+1)*prim(i+1,j-1,k+1,4) + zeta(i,j-1,k+1)*prim(i,j-1,k+1,4) +
                                           zeta(i+1,j,k+1)*prim(i+1,j,k+1,4) + zeta(i,j,k+1)*prim(i,j,k+1,4) );

                            muzemp += 0.25*(zeta(i-1,j,k)*prim(i-1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) +
                                           zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) + zeta(i,j-1,k)*prim(i,j-1,k,4) +
                                           zeta(i-1,j,k+1)*prim(i-1,j,k+1,4) + zeta(i,j,k+1)*prim(i,j,k+1,4) +
                                           zeta(i-1,j-1,k+1)*prim(i-1,j-1,k+1,4) + zeta(i,j-1,k+1)*prim(i,j-1,k+1,4) );

                            muzepm += 0.25*(zeta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,4) + zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                           zeta(i+1,j,k-1)*prim(i+1,j,k-1,4) + zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                           zeta(i+1,j-1,k)*prim(i+1,j-1,k,4) + zeta(i,j-1,k)*prim(i,j-1,k,4) +
                                           zeta(i+1,j,k)*prim(i+1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) );

                            muzemm += 0.25*(zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) + zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                                    zeta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) + zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                                    zeta(i-1,j,k)*prim(i-1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) +
                                                    zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) + zeta(i,j-1,k)*prim(i,j-1,k,4) );
                        }
                    }

                    Real factor_lo_x = 1.;
                    Real factor_hi_x = 1.;
                    Real factor_lo_z = 1.;
                    Real factor_hi_z = 1.;

                    // 1 = slip
                    // 2 = no-slip
                    if (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) {
                        if (i == 0) {
                            factor_lo_x = (bc_vel_lo[0] == 1) ? std::sqrt(2.0) : 0.;
                        }
                    }
                    if (bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) {
                        if (k == 0) {
                            factor_lo_z = (bc_vel_lo[2] == 1) ? std::sqrt(2.0) : 0.;
                        }
                    }
                    if (bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) {
                        if (i == n_cells[0]-1) {
                            factor_hi_x = (bc_vel_hi[0] == 1) ? std::sqrt(2.0) : 0.;
                        }
                    }
                    if (bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) {
                        if (k == n_cells[2]-1) {
                            factor_hi_z = (bc_vel_hi[2] == 1) ? std::sqrt(2.0) : 0.;
                        }
                    }

                    // Random "divergence" stress
                    wiener[2] = wiener[2] + 0.25*nweight*(factor_hi_x*factor_hi_z*sqrt(muzepp)*rancorn(i+1,j,k+1) +
                                                          factor_lo_x*factor_hi_z*sqrt(muzemp)*rancorn(i,j,k+1) +
                                                          factor_hi_x*factor_lo_z*sqrt(muzepm)*rancorn(i+1,j,k) +
                                                          factor_lo_x*factor_lo_z*sqrt(muzemm)*rancorn(i,j,k));

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

          wiener(3) = wiener(3) + 0.5d0*nweight*    &
               (sqrt(muzepp)*rancorn(i+1,j,k) + sqrt(muzemp)*rancorn(i,j,k)) ! Random "divergence" stress
*/
                }

                for (int n=1; n<4; ++n) {
                    fluxy(i,j,k,n) = fluxy(i,j,k,n) + wiener[n];
                }

                // heat flux
                fluxy(i,j,k,nvars) = fluxy(i,j,k,nvars) + wiener[4];

                // Viscous heating:
                Real phiflxdiag = wiener[2]*(prim(i,j-1,k,2)+prim(i,j,k,2));

                Real phiflxshear = wiener[1]*(prim(i,j-1,k,1)+prim(i,j,k,1)) +
                                   wiener[3]*(prim(i,j-1,k,3)+prim(i,j,k,3));

                if ((j == 0) and is_lo_y_dirichlet_mass) {
                    phiflxdiag =  2.0*wiener[2]*prim(i,j-1,k,2);
                    phiflxshear = 2.0*wiener[1]*prim(i,j-1,k,1) +
                                  2.0*wiener[3]*prim(i,j-1,k,3);
                }
                if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                    phiflxdiag =  2.0*wiener[2]*prim(i,j,k,2);
                    phiflxshear = 2.0*wiener[1]*prim(i,j,k,1) +
                                  2.0*wiener[3]*prim(i,j,k,3);
                }

                phiflxdiag = -0.5*phiflxdiag;
                phiflxshear = -0.5*phiflxshear;

                fluxy(i,j,k,nvars+1) = fluxy(i,j,k,nvars+1) - phiflxdiag;
                fluxy(i,j,k,nvars+2) = fluxy(i,j,k,nvars+2) - phiflxshear;

                if (algorithm_type == 2) {

                    for (int n=5; n<5+nspecies; ++n) {
                        wiener[n] = 0.;
                    }

                    for (int ns=0; ns<nspecies; ++ns) {
                        yy[ns] = amrex::max(0.,amrex::min(1.,prim(i,j-1,k,6+ns)));
                        yyp[ns] = amrex::max(0.,amrex::min(1.,prim(i,j,k,6+ns)));
                        if ((j == 0) and is_lo_y_dirichlet_mass) {
                            yyp[ns] = amrex::max(0.,amrex::min(1.,prim(i,j-1,k,6+ns)));
                        }
                        if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                            yy[ns] = amrex::max(0.,amrex::min(1.,prim(i,j,k,6+ns)));
                        }
                    }

                    Real sumy = 0.;
                    Real sumyp = 0.;

                    for (int n=0; n<nspecies; ++n) {
                        sumy += yy[n];
                        sumyp += yyp[n];
                    }

                    for (int n=0; n<nspecies; ++n) {
                        yy[n] /= sumy;
                        yyp[n] /= sumyp;
                    }

                    Real MWmix = 0.;

                    for (int ns=0; ns<nspecies; ++ns) {

                        MWmix = MWmix + 0.5*(yy[ns]+yyp[ns])/molmass[ns];

                        for (int ll=0; ll<nspecies; ++ll) {
                            DijY_edge[ns*nspecies+ll] = 0.5*(Dij(i,j-1,k,ll*nspecies+ns)*yy[ll] +
                                                                 Dij(i,j,k,ll*nspecies+ns)*yyp[ll] +
                                                                (Dij(i,j-1,k,ns*nspecies+ll)*yy[ns] +
                                                                 Dij(i,j,k,ns*nspecies+ll)*yyp[ns] ));

                            if ((j == 0) and is_lo_y_dirichlet_mass) {
                                DijY_edge[ns*nspecies+ll] = 0.5*(Dij(i,j-1,k,ll*nspecies+ns)*yy[ll] +
                                                                     Dij(i,j-1,k,ll*nspecies+ns)*yyp[ll] +
                                                                    (Dij(i,j-1,k,ns*nspecies+ll)*yy[ns] +
                                                                     Dij(i,j-1,k,ns*nspecies+ll)*yyp[ns] ));
                            }
                            if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                                DijY_edge[ns*nspecies+ll] = 0.5*(Dij(i,j,k,ll*nspecies+ns)*yy[ll] +
                                                                     Dij(i,j,k,ll*nspecies+ns)*yyp[ll] +
                                                                    (Dij(i,j,k,ns*nspecies+ll)*yy[ns] +
                                                                     Dij(i,j,k,ns*nspecies+ll)*yyp[ns] ));
                            }
                        }
                    }

                    for (int ns=0; ns<nspecies; ++ns) {
                        if (amrex::Math::abs(yy[ns]) + amrex::Math::abs(yyp[ns]) <= 1.e-12) {
                            for (int n=0; n<nspecies; ++n) {
                                DijY_edge[ns*nspecies+n]=0.;
                                DijY_edge[n*nspecies+ns]=0.;
                            }
                        }
                    }

                    MWmix = 1. / MWmix;

                    CholeskyDecomp(DijY_edge,nspecies,sqD);

                    for (int ns=0; ns<nspecies; ++ns) {
                        for (int ll=0; ll<=ns; ++ll) {
                            fweights[5+ll]=sqrt(k_B*MWmix*volinv/(Runiv*dt))*sqD[ns*nspecies+ll];
                            wiener[5+ns] = wiener[5+ns] + fweights[5+ll]*ranfluxy(i,j,k,5+ll);
                        }
                        fluxy(i,j,k,5+ns) = wiener[5+ns];
                    }

                    GetEnthalpies(meanT, hk);

                    Real soret = 0.;

                    for (int ns=0; ns<nspecies; ++ns) {
                        Real soret_s;
                        soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*0.5*(chi(i,j-1,k,ns)+chi(i,j,k,ns)))*wiener[5+ns];
                        if ((j == 0) and is_lo_y_dirichlet_mass) {
                            soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*chi(i,j-1,k,ns))*wiener[5+ns];
                        }
                        if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                            soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*chi(i,j,k,ns))*wiener[5+ns];
                        }
                        soret = soret + soret_s;
                    }
                    fluxy(i,j,k,nvars+3) = fluxy(i,j,k,nvars+3) + soret;
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                GpuArray<Real,MAX_SPECIES+5> fweights;
                GpuArray<Real,MAX_SPECIES+5> wiener;

                GpuArray<Real,MAX_SPECIES> hk;
                GpuArray<Real,MAX_SPECIES> yy;
                GpuArray<Real,MAX_SPECIES> yyp;

                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> DijY_edge;
                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> sqD;

                if (n_cells_z > 1) {

                    Real muzp = eta(i,j,k)*prim(i,j,k,4) + eta(i,j,k-1)*prim(i,j,k-1,4);
                    Real kzp = kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4) + kappa(i,j,k-1)*prim(i,j,k-1,4)*prim(i,j,k-1,4);

                    Real meanT = 0.5*(prim(i,j,k,4)+prim(i,j,k-1,4));

                    if ((k == 0) and is_lo_z_dirichlet_mass) {
                        muzp = 2.0*eta(i,j,k-1)*prim(i,j,k-1,4);
                        kzp  = 2.0*kappa(i,j,k-1)*prim(i,j,k-1,4)*prim(i,j,k-1,4);
                        meanT = prim(i,j,k-1,4);
                    }
                    if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                        muzp = 2.0*eta(i,j,k)*prim(i,j,k,4);
                        kzp  = 2.0*kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4);
                        meanT = prim(i,j,k,4);
                    }

                    // Weights for facial fluxes:
                    fweights[0] = 0; // No mass flux
                    fweights[1] = sqrt(k_B*muzp*volinv*dtinv);
                    fweights[2] = sqrt(k_B*muzp*volinv*dtinv);
                    fweights[3] = sqrt(k_B*muzp*volinv*dtinv);
                    fweights[4] = sqrt(k_B*kzp*volinv*dtinv);

                    // Construct the random increments
                    for (int n=0; n<5; ++n) {
                        wiener[n] = fweights[n]*ranfluxz(i,j,k,n);
                    }

                    Real nweight=sqrt(k_B*volinv*dtinv);

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

                    if ((k == 0) and is_lo_z_dirichlet_mass) {
                        muzepp = 0.5*(eta(i+1,j,k-1)*prim(i+1,j,k-1,4) +
                                      eta(i,j,k-1)*prim(i,j,k-1,4) +
                                      eta(i+1,j+1,k-1)*prim(i+1,j+1,k-1,4) +
                                      eta(i,j+1,k-1)*prim(i,j+1,k-1,4) )/3.;

                        muzemp = 0.5*(eta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,4) +
                                      eta(i,j+1,k-1)*prim(i,j+1,k-1,4) +
                                      eta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                      eta(i,j,k-1)*prim(i,j,k-1,4) )/3.;

                        muzepm = 0.5*(eta(i+1,j,k-1)*prim(i+1,j,k-1,4) +
                                      eta(i,j,k-1)*prim(i,j,k-1,4) +
                                      eta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,4) +
                                      eta(i,j-1,k-1)*prim(i,j-1,k-1,4) )/3.;

                        muzemm = 0.5*(eta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) +
                                      eta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                      eta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                      eta(i,j,k-1)*prim(i,j,k-1,4) )/3.;

                    }
                    if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                        muzepp = 0.5*(eta(i+1,j,k)*prim(i+1,j,k,4) +
                                       eta(i,j,k)*prim(i,j,k,4) +
                                       eta(i+1,j+1,k)*prim(i+1,j+1,k,4) +
                                       eta(i,j+1,k)*prim(i,j+1,k,4) )/3.;

                        muzemp = 0.5*(eta(i-1,j+1,k)*prim(i-1,j+1,k,4) +
                                       eta(i,j+1,k)*prim(i,j+1,k,4) +
                                       eta(i-1,j,k)*prim(i-1,j,k,4) +
                                       eta(i,j,k)*prim(i,j,k,4) )/3.;

                        muzepm = 0.5*(eta(i+1,j,k)*prim(i+1,j,k,4) +
                                       eta(i,j,k)*prim(i,j,k,4) +
                                       eta(i+1,j-1,k)*prim(i+1,j-1,k,4) +
                                       eta(i,j-1,k)*prim(i,j-1,k,4) )/3.;

                        muzemm = 0.5*(eta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                       eta(i,j-1,k)*prim(i,j-1,k,4) +
                                       eta(i-1,j,k)*prim(i-1,j,k,4) +
                                       eta(i,j,k)*prim(i,j,k,4) )/3.;

                    }

                    if (amrex::Math::abs(visc_type) == 3) {

                        if ((k == 0) and is_lo_z_dirichlet_mass) {
                            muzepp += 0.5*(zeta(i+1,j,k-1)*prim(i+1,j,k-1,4) +
                                          zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                          zeta(i+1,j+1,k-1)*prim(i+1,j+1,k-1,4) +
                                          zeta(i,j+1,k-1)*prim(i,j+1,k-1,4) )/3.;

                            muzemp += 0.5*(zeta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,4) +
                                          zeta(i,j+1,k-1)*prim(i,j+1,k-1,4) +
                                          zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                          zeta(i,j,k-1)*prim(i,j,k-1,4) )/3.;

                            muzepm += 0.5*(zeta(i+1,j,k-1)*prim(i+1,j,k-1,4) +
                                          zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                          zeta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,4) +
                                          zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) )/3.;

                            muzemm += 0.5*(zeta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) +
                                          zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                          zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) +
                                          zeta(i,j,k-1)*prim(i,j,k-1,4) )/3.;

                        }
                        if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                            muzepp += 0.5*(zeta(i+1,j,k)*prim(i+1,j,k,4) +
                                           zeta(i,j,k)*prim(i,j,k,4) +
                                           zeta(i+1,j+1,k)*prim(i+1,j+1,k,4) +
                                           zeta(i,j+1,k)*prim(i,j+1,k,4) )/3.;

                            muzemp += 0.5*(zeta(i-1,j+1,k)*prim(i-1,j+1,k,4) +
                                           zeta(i,j+1,k)*prim(i,j+1,k,4) +
                                           zeta(i-1,j,k)*prim(i-1,j,k,4) +
                                           zeta(i,j,k)*prim(i,j,k,4) )/3.;

                            muzepm += 0.5*(zeta(i+1,j,k)*prim(i+1,j,k,4) +
                                           zeta(i,j,k)*prim(i,j,k,4) +
                                           zeta(i+1,j-1,k)*prim(i+1,j-1,k,4) +
                                           zeta(i,j-1,k)*prim(i,j-1,k,4) )/3.;

                            muzemm += 0.5*(zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) +
                                           zeta(i,j-1,k)*prim(i,j-1,k,4) +
                                           zeta(i-1,j,k)*prim(i-1,j,k,4) +
                                           zeta(i,j,k)*prim(i,j,k,4) )/3.;

                        }
                        else {
                            muzepp += 0.25*(zeta(i+1,j,k-1)*prim(i+1,j,k-1,4) + zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                           zeta(i+1,j+1,k-1)*prim(i+1,j+1,k-1,4) + zeta(i,j+1,k-1)*prim(i,j+1,k-1,4) +
                                           zeta(i+1,j,k)*prim(i+1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) +
                                           zeta(i+1,j+1,k)*prim(i+1,j+1,k,4) + zeta(i,j+1,k)*prim(i,j+1,k,4) );

                            muzemp += 0.25*(zeta(i-1,j+1,k-1)*prim(i-1,j+1,k-1,4) + zeta(i,j+1,k-1)*prim(i,j+1,k-1,4) +
                                            zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) + zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                            zeta(i-1,j+1,k)*prim(i-1,j+1,k,4) + zeta(i,j+1,k)*prim(i,j+1,k,4) +
                                            zeta(i-1,j,k)*prim(i-1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) );

                            muzepm += 0.25*(zeta(i+1,j,k)*prim(i+1,j,k,4) + zeta(i,j,k-2)*prim(i,j,k,4) +
                                            zeta(i+1,j-1,k)*prim(i+1,j-1,k,4) + zeta(i,j-1,k-2)*prim(i,j-1,k,4) +
                                            zeta(i+1,j,k-1)*prim(i+1,j,k-1,4) + zeta(i,j,k-1)*prim(i,j,k-1,4) +
                                            zeta(i+1,j-1,k-1)*prim(i+1,j-1,k-1,4) + zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) );

                            muzemm += 0.25*(zeta(i-1,j-1,k)*prim(i-1,j-1,k,4) + zeta(i,j-1,k)*prim(i,j-1,k,4) +
                                                    zeta(i-1,j,k)*prim(i-1,j,k,4) + zeta(i,j,k)*prim(i,j,k,4) +
                                                    zeta(i-1,j-1,k-1)*prim(i-1,j-1,k-1,4) + zeta(i,j-1,k-1)*prim(i,j-1,k-1,4) +
                                                    zeta(i-1,j,k-1)*prim(i-1,j,k-1,4) + zeta(i,j,k-1)*prim(i,j,k-1,4) );
                        }

                    }

                    Real factor_lo_x = 1.;
                    Real factor_hi_x = 1.;
                    Real factor_lo_y = 1.;
                    Real factor_hi_y = 1.;

                    // 1 = slip
                    // 2 = no-slip
                    if (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) {
                        if (i == 0) {
                            factor_lo_x = (bc_vel_lo[0] == 1) ? std::sqrt(2.0) : 0.;
                        }
                    }
                    if (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) {
                        if (j == 0) {
                            factor_lo_y = (bc_vel_lo[1] == 1) ? std::sqrt(2.0) : 0.;
                        }
                    }
                    if (bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) {
                        if (i == n_cells[0]-1) {
                            factor_hi_x = (bc_vel_hi[0] == 1) ? std::sqrt(2.0) : 0.;
                        }
                    }
                    if (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) {
                        if (j == n_cells[1]-1) {
                            factor_hi_y = (bc_vel_hi[1] == 1) ? std::sqrt(2.0) : 0.;
                        }
                    }

                    // Random "divergence" stress
                    wiener[3] = wiener[3] + 0.25*nweight*(factor_hi_x*factor_hi_y*sqrt(muzepp)*rancorn(i+1,j+1,k) +
                                                          factor_lo_x*factor_hi_y*sqrt(muzemp)*rancorn(i,j+1,k) +
                                                          factor_hi_x*factor_lo_y*sqrt(muzepm)*rancorn(i+1,j,k) +
                                                          factor_lo_x*factor_lo_y*sqrt(muzemm)*rancorn(i,j,k));

                    for (int n=1; n<4; ++n) {
                        fluxz(i,j,k,n) = fluxz(i,j,k,n) + wiener[n];
                    }

                    // heat flux
                    fluxz(i,j,k,nvars) = fluxz(i,j,k,nvars) + wiener[4];

                    // Viscous heating:
                    Real phiflxdiag = wiener[3]*(prim(i,j,k-1,3)+prim(i,j,k,3));

                    Real phiflxshear = wiener[1]*(prim(i,j,k-1,1)+prim(i,j,k,1)) +
                                       wiener[2]*(prim(i,j,k-1,2)+prim(i,j,k,2));

                    if ((k == 0) and is_lo_z_dirichlet_mass) {
                        phiflxdiag =  2.0*wiener[3]*prim(i,j,k-1,3);
                        phiflxshear = 2.0*wiener[1]*prim(i,j,k-1,1) +
                                      2.0*wiener[2]*prim(i,j,k-1,2);
                    }
                    if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                        phiflxdiag =  2.0*wiener[3]*prim(i,j,k,3);
                        phiflxshear = 2.0*wiener[1]*prim(i,j,k,1) +
                                      2.0*wiener[2]*prim(i,j,k,2);
                    }

                    phiflxdiag = -0.5*phiflxdiag;
                    phiflxshear = -0.5*phiflxshear;

                    fluxz(i,j,k,nvars+1) = fluxz(i,j,k,nvars+1) - phiflxdiag;
                    fluxz(i,j,k,nvars+2) = fluxz(i,j,k,nvars+2) - phiflxshear;

                    if (algorithm_type == 2) {

                    for (int n=5; n<5+nspecies; ++n) {
                        wiener[n] = 0.;
                    }

                    for (int ns=0; ns<nspecies; ++ns) {
                        yy[ns] = amrex::max(0.,amrex::min(1.,prim(i,j,k-1,6+ns)));
                        yyp[ns] = amrex::max(0.,amrex::min(1.,prim(i,j,k,6+ns)));
                        if ((k == 0) and is_lo_z_dirichlet_mass) {
                            yyp[ns] = amrex::max(0.,amrex::min(1.,prim(i,j,k-1,6+ns)));
                        }
                        if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                            yy[ns] = amrex::max(0.,amrex::min(1.,prim(i,j,k,6+ns)));
                        }
                    }

                    Real sumy = 0.;
                    Real sumyp = 0.;

                    for (int n=0; n<nspecies; ++n) {
                        sumy += yy[n];
                        sumyp += yyp[n];
                    }

                    for (int n=0; n<nspecies; ++n) {
                        yy[n] /= sumy;
                        yyp[n] /= sumyp;
                    }

                    Real MWmix = 0.;

                    for (int ns=0; ns<nspecies; ++ns) {

                        MWmix = MWmix + 0.5*(yy[ns]+yyp[ns])/molmass[ns];

                        for (int ll=0; ll<nspecies; ++ll) {
                            DijY_edge[ns*nspecies+ll] = 0.5*(Dij(i,j,k-1,ll*nspecies+ns)*yy[ll] +
                                                                 Dij(i,j,k,ll*nspecies+ns)*yyp[ll] +
                                                                (Dij(i,j,k-1,ns*nspecies+ll)*yy[ns] +
                                                                 Dij(i,j,k,ns*nspecies+ll)*yyp[ns] ));
                            if ((k == 0) and is_lo_z_dirichlet_mass) {
                                DijY_edge[ns*nspecies+ll] = 0.5*(Dij(i,j,k-1,ll*nspecies+ns)*yy[ll] +
                                                                     Dij(i,j,k-1,ll*nspecies+ns)*yyp[ll] +
                                                                    (Dij(i,j,k-1,ns*nspecies+ll)*yy[ns] +
                                                                     Dij(i,j,k-1,ns*nspecies+ll)*yyp[ns] ));
                            }
                            if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                                DijY_edge[ns*nspecies+ll] = 0.5*(Dij(i,j,k,ll*nspecies+ns)*yy[ll] +
                                                                     Dij(i,j,k,ll*nspecies+ns)*yyp[ll] +
                                                                    (Dij(i,j,k,ns*nspecies+ll)*yy[ns] +
                                                                     Dij(i,j,k,ns*nspecies+ll)*yyp[ns] ));
                            }

                        }
                    }


                    for (int ns=0; ns<nspecies; ++ns) {
                        if (amrex::Math::abs(yy[ns]) + amrex::Math::abs(yyp[ns]) <= 1.e-12) {
                            for (int n=0; n<nspecies; ++n) {
                                DijY_edge[ns*nspecies+n]=0.;
                                DijY_edge[n*nspecies+ns]=0.;
                            }
                        }
                    }

                    MWmix = 1. / MWmix;

                    CholeskyDecomp(DijY_edge,nspecies,sqD);

                    for (int ns=0; ns<nspecies; ++ns) {
                        for (int ll=0; ll<=ns; ++ll) {
                            fweights[5+ll]=sqrt(k_B*MWmix*volinv/(Runiv*dt))*sqD[ns*nspecies+ll];
                            wiener[5+ns] = wiener[5+ns] + fweights[5+ll]*ranfluxz(i,j,k,5+ll);
                        }
                        fluxz(i,j,k,5+ns) = wiener[5+ns];
                    }

                    GetEnthalpies(meanT, hk);

                    Real soret = 0.;

                    for (int ns=0; ns<nspecies; ++ns) {
                        Real soret_s;
                        soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*0.5*(chi(i,j,k-1,ns)+chi(i,j,k,ns)))*wiener[5+ns];
                        if ((k == 0) and is_lo_z_dirichlet_mass) {
                            soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*chi(i,j,k-1,ns))*wiener[5+ns];
                        }
                        if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                            soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*chi(i,j,k,ns))*wiener[5+ns];
                        }
                        soret = soret + soret_s;
                    }
                    fluxz(i,j,k,nvars+3) = fluxz(i,j,k,nvars+3) + soret;

                    }
                }

            }); // end lambda function

        } // end MFIter

        StochFlux(flux_in,geom);
        MembraneFlux(flux_in,geom);
    }

    ////////////////////
    // diffusive flxues
    ////////////////////

    // Loop over boxes
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {

        AMREX_D_TERM(const Array4<Real>& fluxx = flux_in[0].array(mfi); ,
                     const Array4<Real>& fluxy = flux_in[1].array(mfi); ,
                     const Array4<Real>& fluxz = flux_in[2].array(mfi));

        const Array4<const Real> prim = prim_in.array(mfi);
        // const Array4<const Real> cons = cons_in.array(mfi);

        const Array4<const Real> eta   = eta_in.array(mfi);
        const Array4<const Real> zeta  = zeta_in.array(mfi);
        const Array4<const Real> kappa = kappa_in.array(mfi);
        const Array4<const Real> chi   = chi_in.array(mfi);
        const Array4<const Real> Dij   = D_in.array(mfi);

        const Array4<Real> cornux = cornx_in[0].array(mfi);
        const Array4<Real> cornvx = cornx_in[1].array(mfi);
        const Array4<Real> cornwx = cornx_in[2].array(mfi);
        const Array4<Real> cornuy = corny_in[0].array(mfi);
        const Array4<Real> cornvy = corny_in[1].array(mfi);
        const Array4<Real> cornwy = corny_in[2].array(mfi);
        const Array4<Real> cornuz = cornz_in[0].array(mfi);
        const Array4<Real> cornvz = cornz_in[1].array(mfi);
        const Array4<Real> cornwz = cornz_in[2].array(mfi);
        const Array4<Real> visccorn = visccorn_in.array(mfi);

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        IntVect nd(AMREX_D_DECL(1,1,1));
        const Box& tbn = mfi.tilebox(nd);

        Real half = 0.5;

        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            GpuArray<Real,MAX_SPECIES> meanXk;
            GpuArray<Real,MAX_SPECIES> meanYk;
            GpuArray<Real,MAX_SPECIES> dk;
            GpuArray<Real,MAX_SPECIES> Fk;
            GpuArray<Real,MAX_SPECIES> hk;
            GpuArray<Real,MAX_SPECIES> soret;

            Real muxp = half*(eta(i,j,k) + eta(i-1,j,k));
            Real kxp = half*(kappa(i,j,k) + kappa(i-1,j,k));

            Real tauxxp = muxp*(prim(i,j,k,1) - prim(i-1,j,k,1))/dx[0];
            Real tauyxp = muxp*(prim(i,j,k,2) - prim(i-1,j,k,2))/dx[0];
            Real tauzxp = muxp*(prim(i,j,k,3) - prim(i-1,j,k,3))/dx[0];

            Real divxp = 0.;

            Real phiflx =  tauxxp*(prim(i-1,j,k,1)+prim(i,j,k,1))
                +  divxp*(prim(i-1,j,k,1)+prim(i,j,k,1))
                +  tauyxp*(prim(i-1,j,k,2)+prim(i,j,k,2))
                +  tauzxp*(prim(i-1,j,k,3)+prim(i,j,k,3));

            Real Qflux = kxp*(prim(i,j,k,4)-prim(i-1,j,k,4))/dx[0];

            if ((i == 0) and is_lo_x_dirichlet_mass) {
                muxp = eta(i-1,j,k);
                kxp  = kappa(i-1,j,k);
                tauxxp = muxp*(prim(i,j,k,1) - prim(i-1,j,k,1))/(0.5*dx[0]);
                tauyxp = muxp*(prim(i,j,k,2) - prim(i-1,j,k,2))/(0.5*dx[0]);
                tauzxp = muxp*(prim(i,j,k,3) - prim(i-1,j,k,3))/(0.5*dx[0]);
                phiflx = 2.0*(tauxxp*(prim(i-1,j,k,1))
                          +  divxp*(prim(i-1,j,k,1))
                          +  tauyxp*(prim(i-1,j,k,2))
                          +  tauzxp*(prim(i-1,j,k,3)));
                Qflux = kxp*(prim(i,j,k,4)-prim(i-1,j,k,4))/(0.5*dx[0]);
            }
            if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                muxp = eta(i,j,k);
                kxp  = kappa(i,j,k);
                tauxxp = muxp*(prim(i,j,k,1) - prim(i-1,j,k,1))/(0.5*dx[0]);
                tauyxp = muxp*(prim(i,j,k,2) - prim(i-1,j,k,2))/(0.5*dx[0]);
                tauzxp = muxp*(prim(i,j,k,3) - prim(i-1,j,k,3))/(0.5*dx[0]);
                phiflx = 2.0*(tauxxp*(prim(i,j,k,1))
                          +  divxp*(prim(i,j,k,1))
                          +  tauyxp*(prim(i,j,k,2))
                          +  tauzxp*(prim(i,j,k,3)));
                Qflux = kxp*(prim(i,j,k,4)-prim(i-1,j,k,4))/(0.5*dx[0]);
            }

            fluxx(i,j,k,1) = fluxx(i,j,k,1) - (tauxxp+divxp);
            fluxx(i,j,k,2) = fluxx(i,j,k,2) - tauyxp;
            fluxx(i,j,k,3) = fluxx(i,j,k,3) - tauzxp;

            // heat flux
            fluxx(i,j,k,nvars) = fluxx(i,j,k,nvars) - Qflux;

            // viscous heating
            fluxx(i,j,k,nvars+1) = fluxx(i,j,k,nvars+1) - (half*phiflx);

            Real meanT = 0.5*(prim(i-1,j,k,4)+prim(i,j,k,4));
            Real meanP = 0.5*(prim(i-1,j,k,5)+prim(i,j,k,5));
            if ((i == 0) and is_lo_x_dirichlet_mass) {
                meanT = prim(i-1,j,k,4);
                meanP = prim(i-1,j,k,5);
            }
            if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                meanT = prim(i,j,k,4);
                meanP = prim(i,j,k,5);
            }

            if (algorithm_type == 2) {

                // compute dk
                for (int ns=0; ns<nspecies; ++ns) {
                    Real term1 = (prim(i,j,k,6+nspecies+ns)-prim(i-1,j,k,6+nspecies+ns))/dx[0];
                    meanXk[ns] = 0.5*(prim(i-1,j,k,6+nspecies+ns)+prim(i,j,k,6+nspecies+ns));
                    meanYk[ns] = 0.5*(prim(i-1,j,k,6+ns)+prim(i,j,k,6+ns));
                    Real term2 = (meanXk[ns]-meanYk[ns])*(prim(i,j,k,5)-prim(i-1,j,k,5))/dx[0]/meanP;
                    dk[ns] = term1 + term2;
                    Real ChiX = 0.5*(chi(i-1,j,k,ns)*prim(i-1,j,k,6+nspecies+ns)+chi(i,j,k,ns)*prim(i,j,k,6+nspecies+ns));
                    soret[ns] = ChiX*(prim(i,j,k,4)-prim(i-1,j,k,4))/dx[0]/meanT;

                    if ((i == 0) and is_lo_x_dirichlet_mass) {
                        term1 = (prim(i,j,k,6+nspecies+ns)-prim(i-1,j,k,6+nspecies+ns))/(0.5*dx[0]);
                        meanXk[ns] = prim(i-1,j,k,6+nspecies+ns);
                        meanYk[ns] = prim(i-1,j,k,6+ns);
                        term2 = (meanXk[ns]-meanYk[ns])*(prim(i,j,k,5)-prim(i-1,j,k,5))/(0.5*dx[0])/meanP;
                        dk[ns] = term1 + term2;
                        ChiX = chi(i-1,j,k,ns)*prim(i-1,j,k,6+nspecies+ns);
                        soret[ns] = ChiX*(prim(i,j,k,4)-prim(i-1,j,k,4))/(0.5*dx[0])/meanT;
                    }
                    if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                        term1 = (prim(i,j,k,6+nspecies+ns)-prim(i-1,j,k,6+nspecies+ns))/(0.5*dx[0]);
                        meanXk[ns] = prim(i,j,k,6+nspecies+ns);
                        meanYk[ns] = prim(i,j,k,6+ns);
                        term2 = (meanXk[ns]-meanYk[ns])*(prim(i,j,k,5)-prim(i-1,j,k,5))/(0.5*dx[0])/meanP;
                        dk[ns] = term1 + term2;
                        ChiX = chi(i,j,k,ns)*prim(i,j,k,6+nspecies+ns);
                        soret[ns] = ChiX*(prim(i,j,k,4)-prim(i-1,j,k,4))/(0.5*dx[0])/meanT;
                    }
                }

                // compute Fk (based on Eqn. 2.5.24, Giovangigli's book)
                for (int kk=0; kk<nspecies; ++kk) {
                    Fk[kk] = 0.;
                    for (int ll=0; ll<nspecies; ++ll) {
                        Real Fks = half*(Dij(i-1,j,k,ll*nspecies+kk)+Dij(i,j,k,ll*nspecies+kk))*( dk[ll] +soret[ll]);
                        if ((i == 0) and is_lo_x_dirichlet_mass) {
                            Fks = Dij(i-1,j,k,ll*nspecies+kk)*( dk[ll] +soret[ll]);
                        }
                        if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                            Fks = Dij(i,j,k,ll*nspecies+kk)*( dk[ll] +soret[ll]);
                        }
                        Fk[kk] = Fk[kk] - Fks;
                    }
                }

                // compute Q (based on Eqn. 2.5.25, Giovangigli's book)
                GetEnthalpies(meanT,hk);

                Real Q5 = 0.;
                for (int ns=0; ns<nspecies; ++ns) {
                    Real Q5s = (hk[ns] + 0.5 * Runiv*meanT*(chi(i-1,j,k,ns)+chi(i,j,k,ns))/molmass[ns])*Fk[ns];
                    if ((i == 0) and is_lo_x_dirichlet_mass) {
                        Q5s = (hk[ns] + Runiv*meanT*chi(i-1,j,k,ns)/molmass[ns])*Fk[ns];
                    }
                    if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                        Q5s = (hk[ns] + Runiv*meanT*chi(i,j,k,ns)/molmass[ns])*Fk[ns];
                    }
                    Q5 = Q5 + Q5s;
                }
                // heat conduction already included in flux(5)

                fluxx(i,j,k,nvars+3) = fluxx(i,j,k,nvars+3) + Q5;

                for (int ns=0; ns<nspecies; ++ns) {
                    fluxx(i,j,k,5+ns) = fluxx(i,j,k,5+ns) + Fk[ns];
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

            Real muyp = half*(eta(i,j,k) + eta(i,j-1,k));
            Real kyp = half*(kappa(i,j,k) + kappa(i,j-1,k));

            Real tauxyp =  muyp*(prim(i,j,k,1) - prim(i,j-1,k,1))/dx[1];
            Real tauyyp =  muyp*(prim(i,j,k,2) - prim(i,j-1,k,2))/dx[1];
            Real tauzyp =  muyp*(prim(i,j,k,3) - prim(i,j-1,k,3))/dx[1];

            Real divyp = 0.;

            Real phiflx = tauxyp*(prim(i,j,k,1)+prim(i,j-1,k,1))
                +  tauyyp*(prim(i,j,k,2)+prim(i,j-1,k,2))
                +  divyp*(prim(i,j,k,2)+prim(i,j-1,k,2))
                +  tauzyp*(prim(i,j,k,3)+prim(i,j-1,k,3));

            Real Qflux = kyp*(prim(i,j,k,4)-prim(i,j-1,k,4))/dx[1];

            if ((j == 0) and is_lo_y_dirichlet_mass) {
                muyp = eta(i,j-1,k);
                kyp  = kappa(i,j-1,k);
                tauxyp = muyp*(prim(i,j,k,1) - prim(i,j-1,k,1))/(0.5*dx[1]);
                tauyyp = muyp*(prim(i,j,k,2) - prim(i,j-1,k,2))/(0.5*dx[1]);
                tauzyp = muyp*(prim(i,j,k,3) - prim(i,j-1,k,3))/(0.5*dx[1]);
                phiflx = 2.0*(tauxyp*(prim(i,j-1,k,1))
                          +  tauyyp*(prim(i,j-1,k,2))
                          +  divyp*(prim(i,j-1,k,2))
                          +  tauzyp*(prim(i,j-1,k,3)));
                Qflux = kyp*(prim(i,j,k,4)-prim(i,j-1,k,4))/(0.5*dx[1]);
            }
            if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                muyp = eta(i,j,k);
                kyp  = kappa(i,j,k);
                tauxyp = muyp*(prim(i,j,k,1) - prim(i,j-1,k,1))/(0.5*dx[1]);
                tauyyp = muyp*(prim(i,j,k,2) - prim(i,j-1,k,2))/(0.5*dx[1]);
                tauzyp = muyp*(prim(i,j,k,3) - prim(i,j-1,k,3))/(0.5*dx[1]);
                phiflx = 2.0*(tauxyp*(prim(i,j,k,1))
                          +  tauyyp*(prim(i,j,k,2))
                          +  divyp*(prim(i,j,k,2))
                          +  tauzyp*(prim(i,j,k,3)));
                Qflux = kyp*(prim(i,j,k,4)-prim(i,j-1,k,4))/(0.5*dx[1]);
            }

            fluxy(i,j,k,1) = fluxy(i,j,k,1) - tauxyp;
            fluxy(i,j,k,2) = fluxy(i,j,k,2) - (tauyyp+divyp);
            fluxy(i,j,k,3) = fluxy(i,j,k,3) - tauzyp;

            // heat flux
            fluxy(i,j,k,nvars) = fluxy(i,j,k,nvars) - Qflux;

            // viscous heating
            fluxy(i,j,k,nvars+1) = fluxy(i,j,k,nvars+1) - (half*phiflx);

            Real meanT = 0.5*(prim(i,j-1,k,4)+prim(i,j,k,4));
            Real meanP = 0.5*(prim(i,j-1,k,5)+prim(i,j,k,5));
            if ((j == 0) and is_lo_y_dirichlet_mass) {
                meanT = prim(i,j-1,k,4);
                meanP = prim(i,j-1,k,5);
            }
            if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                meanT = prim(i,j,k,4);
                meanP = prim(i,j,k,5);
            }

            if (algorithm_type == 2) {
                // compute dk
                for (int ns=0; ns<nspecies; ++ns) {
                    Real term1 = (prim(i,j,k,6+nspecies+ns)-prim(i,j-1,k,6+nspecies+ns))/dx[1];
                    meanXk[ns] = 0.5*(prim(i,j-1,k,6+nspecies+ns)+prim(i,j,k,6+nspecies+ns));
                    meanYk[ns] = 0.5*(prim(i,j-1,k,6+ns)+prim(i,j,k,6+ns));
                    Real term2 = (meanXk[ns]-meanYk[ns])*(prim(i,j,k,5)-prim(i,j-1,k,5))/dx[1]/meanP;
                    dk[ns] = term1 + term2;
                    Real ChiX = 0.5*(chi(i,j-1,k,ns)*prim(i,j-1,k,6+nspecies+ns)+chi(i,j,k,ns)*prim(i,j,k,6+nspecies+ns));
                    soret[ns] = ChiX*(prim(i,j,k,4)-prim(i,j-1,k,4))/dx[1]/meanT;

                    if ((j == 0) and is_lo_y_dirichlet_mass) {
                        term1 = (prim(i,j,k,6+nspecies+ns)-prim(i,j-1,k,6+nspecies+ns))/(0.5*dx[1]);
                        meanXk[ns] = prim(i,j-1,k,6+nspecies+ns);
                        meanYk[ns] = prim(i,j-1,k,6+ns);
                        term2 = (meanXk[ns]-meanYk[ns])*(prim(i,j,k,5)-prim(i,j-1,k,5))/(0.5*dx[1])/meanP;
                        dk[ns] = term1 + term2;
                        ChiX = chi(i,j-1,k,ns)*prim(i,j-1,k,6+nspecies+ns);
                        soret[ns] = ChiX*(prim(i,j,k,4)-prim(i,j-1,k,4))/(0.5*dx[1])/meanT;
                    }
                    if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                        term1 = (prim(i,j,k,6+nspecies+ns)-prim(i,j-1,k,6+nspecies+ns))/(0.5*dx[1]);
                        meanXk[ns] = prim(i,j,k,6+nspecies+ns);
                        meanYk[ns] = prim(i,j,k,6+ns);
                        term2 = (meanXk[ns]-meanYk[ns])*(prim(i,j,k,5)-prim(i,j-1,k,5))/(0.5*dx[1])/meanP;
                        dk[ns] = term1 + term2;
                        ChiX = chi(i,j,k,ns)*prim(i,j,k,6+nspecies+ns);
                        soret[ns] = ChiX*(prim(i,j,k,4)-prim(i,j-1,k,4))/(0.5*dx[1])/meanT;
                    }

                }

                // compute Fk (based on Eqn. 2.5.24, Giovangigli's book)
                for (int kk=0; kk<nspecies; ++kk) {
                    Fk[kk] = 0.;
                    for (int ll=0; ll<nspecies; ++ll) {
                        Real Fks = half*(Dij(i,j-1,k,ll*nspecies+kk)+Dij(i,j,k,ll*nspecies+kk))*( dk[ll] +soret[ll]);
                        if ((j == 0) and is_lo_y_dirichlet_mass) {
                            Fks = Dij(i,j-1,k,ll*nspecies+kk)*( dk[ll] +soret[ll]);
                        }
                        if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                            Fks = Dij(i,j,k,ll*nspecies+kk)*( dk[ll] +soret[ll]);
                        }
                        Fk[kk] = Fk[kk] - Fks;
                    }
                }

                // compute Q (based on Eqn. 2.5.25, Giovangigli's book)
                GetEnthalpies(meanT,hk);

                Real Q5 = 0.0;
                for (int ns=0; ns<nspecies; ++ns) {
                    Real Q5s = (hk[ns] + 0.5 * Runiv*meanT*(chi(i,j-1,k,ns)+chi(i,j,k,ns))/molmass[ns])*Fk[ns];
                    if ((j == 0) and is_lo_y_dirichlet_mass) {
                        Q5s = (hk[ns] + Runiv*meanT*chi(i,j-1,k,ns)/molmass[ns])*Fk[ns];
                    }
                    if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                        Q5s = (hk[ns] + Runiv*meanT*chi(i,j,k,ns)/molmass[ns])*Fk[ns];
                    }
                    Q5 = Q5 + Q5s;
                }

                // heat conduction already included in flux(5)

                fluxy(i,j,k,nvars+3) = fluxy(i,j,k,nvars+3) + Q5;

                for (int ns=0; ns<nspecies; ++ns) {
                    fluxy(i,j,k,5+ns) = fluxy(i,j,k,5+ns) + Fk[ns];
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

            Real muzp = half*(eta(i,j,k) + eta(i,j,k-1));
            Real kzp = half*(kappa(i,j,k) + kappa(i,j,k-1));

            Real tauxzp =  muzp*(prim(i,j,k,1) - prim(i,j,k-1,1))/dx[2];
            Real tauyzp =  muzp*(prim(i,j,k,2) - prim(i,j,k-1,2))/dx[2];
            Real tauzzp =  muzp*(prim(i,j,k,3) - prim(i,j,k-1,3))/dx[2];

            Real divzp = 0.;

            Real phiflx = tauxzp*(prim(i,j,k-1,1)+prim(i,j,k,1))
                +  tauyzp*(prim(i,j,k-1,2)+prim(i,j,k,2))
                +  tauzzp*(prim(i,j,k-1,3)+prim(i,j,k,3))
                +  divzp*(prim(i,j,k-1,3)+prim(i,j,k,3));

            Real Qflux = kzp*(prim(i,j,k,4)-prim(i,j,k-1,4))/dx[2];

            if ((k == 0) and is_lo_z_dirichlet_mass) {
                muzp = eta(i,j,k-1);
                kzp  = kappa(i,j,k-1);
                tauxzp = muzp*(prim(i,j,k,1) - prim(i,j,k-1,1))/(0.5*dx[2]);
                tauyzp = muzp*(prim(i,j,k,2) - prim(i,j,k-1,2))/(0.5*dx[2]);
                tauzzp = muzp*(prim(i,j,k,3) - prim(i,j,k-1,3))/(0.5*dx[2]);
                phiflx = 2.0*(tauxzp*(prim(i,j,k-1,1))
                          +  tauyzp*(prim(i,j,k-1,2))
                          +  tauzzp*(prim(i,j,k-1,3))
                          +  divzp*(prim(i,j,k-1,3)));
                Qflux = kzp*(prim(i,j,k,4)-prim(i,j,k-1,4))/(0.5*dx[2]);
            }
            if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                muzp = eta(i,j,k);
                kzp  = kappa(i,j,k);
                tauxzp = muzp*(prim(i,j,k,1) - prim(i,j,k-1,1))/(0.5*dx[2]);
                tauyzp = muzp*(prim(i,j,k,2) - prim(i,j,k-1,2))/(0.5*dx[2]);
                tauzzp = muzp*(prim(i,j,k,3) - prim(i,j,k-1,3))/(0.5*dx[2]);
                phiflx = 2.0*(tauxzp*(prim(i,j,k,1))
                          +  tauyzp*(prim(i,j,k,2))
                          +  tauzzp*(prim(i,j,k,3))
                          +  divzp*(prim(i,j,k,3)));
                Qflux = kzp*(prim(i,j,k,4)-prim(i,j,k-1,4))/(0.5*dx[2]);
            }

            fluxz(i,j,k,1) = fluxz(i,j,k,1) - tauxzp;
            fluxz(i,j,k,2) = fluxz(i,j,k,2) - tauyzp;
            fluxz(i,j,k,3) = fluxz(i,j,k,3) - (tauzzp+divzp);

            // heat flux
            fluxz(i,j,k,nvars) = fluxz(i,j,k,nvars) - Qflux;

            // viscous heating
            fluxz(i,j,k,nvars+1) = fluxz(i,j,k,nvars+1) - (half*phiflx);

            Real meanT = 0.5*(prim(i,j,k-1,4)+prim(i,j,k,4));
            Real meanP = 0.5*(prim(i,j,k-1,5)+prim(i,j,k,5));
            if ((k == 0) and is_lo_z_dirichlet_mass) {
                meanT = prim(i,j,k-1,4);
                meanP = prim(i,j,k-1,5);
            }
            if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                meanT = prim(i,j,k,4);
                meanP = prim(i,j,k,5);
            }

            if (algorithm_type == 2) {

                // compute dk
                for (int ns=0; ns<nspecies; ++ns) {
                    Real term1 = (prim(i,j,k,6+nspecies+ns)-prim(i,j,k-1,6+nspecies+ns))/dx[2];
                    meanXk[ns] = 0.5*(prim(i,j,k-1,6+nspecies+ns)+prim(i,j,k,6+nspecies+ns));
                    meanYk[ns] = 0.5*(prim(i,j,k-1,6+ns)+prim(i,j,k,6+ns));
                    Real term2 = (meanXk[ns]-meanYk[ns])*(prim(i,j,k,5)-prim(i,j,k-1,5))/dx[2]/meanP;
                    dk[ns] = term1 + term2;
                    Real ChiX = 0.5*(chi(i,j,k,ns)*prim(i,j,k-1,6+nspecies+ns)+chi(i,j,k+1,ns)*prim(i,j,k,6+nspecies+ns));
                    soret[ns] = ChiX*(prim(i,j,k,4)-prim(i,j,k-1,4))/dx[2]/meanT;

                    if ((k == 0) and is_lo_z_dirichlet_mass) {
                        term1 = (prim(i,j,k,6+nspecies+ns)-prim(i,j,k-1,6+nspecies+ns))/(0.5*dx[2]);
                        meanXk[ns] = prim(i,j,k-1,6+nspecies+ns);
                        meanYk[ns] = prim(i,j,k-1,6+ns);
                        term2 = (meanXk[ns]-meanYk[ns])*(prim(i,j,k,5)-prim(i,j,k-1,5))/(0.5*dx[2])/meanP;
                        dk[ns] = term1 + term2;
                        ChiX = chi(i,j,k-1,ns)*prim(i,j,k-1,6+nspecies+ns);
                        soret[ns] = ChiX*(prim(i,j,k,4)-prim(i,j,k-1,4))/(0.5*dx[2])/meanT;
                    }
                    if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                        term1 = (prim(i,j,k,6+nspecies+ns)-prim(i,j,k-1,6+nspecies+ns))/(0.5*dx[2]);
                        meanXk[ns] = prim(i,j,k,6+nspecies+ns);
                        meanYk[ns] = prim(i,j,k,6+ns);
                        term2 = (meanXk[ns]-meanYk[ns])*(prim(i,j,k,5)-prim(i,j,k-1,5))/(0.5*dx[2])/meanP;
                        dk[ns] = term1 + term2;
                        ChiX = chi(i,j,k,ns)*prim(i,j,k,6+nspecies+ns);
                        soret[ns] = ChiX*(prim(i,j,k,4)-prim(i,j,k-1,4))/(0.5*dx[2])/meanT;
                    }
                }

                // compute Fk (based on Eqn. 2.5.24, Giovangigli's book)
                for (int kk=0; kk<nspecies; ++kk) {
                    Fk[kk] = 0.;
                    for (int ll=0; ll<nspecies; ++ll) {
                        Real Fks = half*(Dij(i,j,k-1,ll*nspecies+kk)+Dij(i,j,k,ll*nspecies+kk))*( dk[ll] +soret[ll]);
                        if ((k == 0) and is_lo_z_dirichlet_mass) {
                            Fks = Dij(i,j,k-1,ll*nspecies+kk)*( dk[ll] +soret[ll]);
                        }
                        if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                            Fks = Dij(i,j,k,ll*nspecies+kk)*( dk[ll] +soret[ll]);
                        }
                        Fk[kk] = Fk[kk] - Fks;
                    }
                }

                // compute Q (based on Eqn. 2.5.25, Giovangigli's book)
                GetEnthalpies(meanT,hk);

                Real Q5 = 0.0;
                for (int ns=0; ns<nspecies; ++ns) {
                    Real Q5s = (hk[ns] + 0.5 * Runiv*meanT*(chi(i,j,k,ns)+chi(i,j,k,ns))/molmass[ns])*Fk[ns];
                    if ((k == 0) and is_lo_z_dirichlet_mass) {
                        Q5s = (hk[ns] + Runiv*meanT*chi(i,j,k-1,ns)/molmass[ns])*Fk[ns];
                    }
                    if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                        Q5s = (hk[ns] + Runiv*meanT*chi(i,j,k,ns)/molmass[ns])*Fk[ns];
                    }
                    Q5 = Q5 + Q5s;
                }

                // heat conduction already included in flux(5)
                fluxz(i,j,k,nvars+3) = fluxz(i,j,k,nvars+3) + Q5;

                for (int ns=0; ns<nspecies; ++ns) {
                    fluxz(i,j,k,5+ns) = fluxz(i,j,k,5+ns) + Fk[ns];
                }
            }

            } // n_cells_z test
        });

        if (n_cells_z > 1) {

        amrex::ParallelFor(tbn,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            // Corner viscosity
            Real muxp = 0.125*(eta(i,j-1,k-1) + eta(i-1,j-1,k-1) + eta(i,j,k-1) + eta(i-1,j,k-1)
                               + eta(i,j-1,k) + eta(i-1,j-1,k) + eta(i,j,k) + eta(i-1,j,k));

            Real zetaxp;
            if (amrex::Math::abs(visc_type) == 3) {
                zetaxp = 0.125*(zeta(i,j-1,k-1) + zeta(i-1,j-1,k-1) + zeta(i,j,k-1) + zeta(i-1,j,k-1)+
                                zeta(i,j-1,k) + zeta(i-1,j-1,k) + zeta(i,j,k) + zeta(i-1,j,k));
            } else {
                zetaxp = 0.;
            }

            GpuArray<Real,AMREX_SPACEDIM> DX;
            DX[0] = dx[0];
            DX[1] = dx[1];
            DX[2] = dx[2];

            if ((i == 0) and is_lo_x_dirichlet_mass) {
                DX[0] = 0.5*dx[0];
                muxp = 0.25*(eta(i-1,j-1,k-1) + eta(i-1,j-1,k) + eta(i-1,j,k-1) + eta(i-1,j,k));
                if (amrex::Math::abs(visc_type) == 3) zetaxp = 0.25*(zeta(i-1,j-1,k-1) + zeta(i-1,j-1,k) + zeta(i-1,j,k-1) + zeta(i-1,j,k));
                else zetaxp = 0.;
            }
            if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                DX[0] = 0.5*dx[0];
                muxp = 0.25*(eta(i,j-1,k-1) + eta(i,j-1,k) + eta(i,j,k-1) + eta(i,j,k));
                if (amrex::Math::abs(visc_type) == 3) zetaxp = 0.25*(zeta(i,j-1,k-1) + zeta(i,j-1,k) + zeta(i,j,k-1) + zeta(i,j,k));
                else zetaxp = 0.;
            }
            if ((j == 0) and is_lo_y_dirichlet_mass) {
                DX[1] = 0.5*dx[1];
                muxp = 0.25*(eta(i-1,j-1,k-1) + eta(i-1,j-1,k) + eta(i,j-1,k-1) + eta(i,j-1,k));
                if (amrex::Math::abs(visc_type) == 3) zetaxp = 0.25*(zeta(i-1,j-1,k-1) + zeta(i-1,j-1,k) + zeta(i,j-1,k-1) + zeta(i,j-1,k));
                else zetaxp = 0.;
            }
            if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                DX[1] = 0.5*dx[1];
                muxp = 0.25*(eta(i-1,j,k-1) + eta(i-1,j,k) + eta(i,j,k-1) + eta(i,j,k));
                if (amrex::Math::abs(visc_type) == 3) zetaxp = 0.25*(zeta(i-1,j,k-1) + zeta(i-1,j,k) + zeta(i,j,k-1) + zeta(i,j,k));
                else zetaxp = 0.;
            }
            if ((k == 0) and is_lo_z_dirichlet_mass) {
                DX[2] = 0.5*dx[2];
                muxp = 0.25*(eta(i-1,j-1,k-1) + eta(i-1,j,k-1) + eta(i,j-1,k-1) + eta(i,j,k-1));
                if (amrex::Math::abs(visc_type) == 3) zetaxp = 0.25*(zeta(i-1,j-1,k-1) + zeta(i-1,j,k-1) + zeta(i,j-1,k-1) + zeta(i,j,k-1));
                else zetaxp = 0.;
            }
            if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                DX[2] = 0.5*dx[2];
                muxp = 0.25*(eta(i-1,j-1,k) + eta(i-1,j,k) + eta(i,j-1,k) + eta(i,j,k));
                if (amrex::Math::abs(visc_type) == 3) zetaxp = 0.25*(zeta(i-1,j-1,k) + zeta(i-1,j,k) + zeta(i,j-1,k) + zeta(i,j,k));
                else zetaxp = 0.;
            }

            cornux(i,j,k) = 0.25*muxp*(prim(i,j-1,k-1,1)-prim(i-1,j-1,k-1,1) + prim(i,j,k-1,1)-prim(i-1,j,k-1,1)+
                                         prim(i,j-1,k,1)-prim(i-1,j-1,k,1) + prim(i,j,k,1)-prim(i-1,j,k,1))/DX[0];
            cornvx(i,j,k) = 0.25*muxp*(prim(i,j-1,k-1,2)-prim(i-1,j-1,k-1,2) + prim(i,j,k-1,2)-prim(i-1,j,k-1,2)+
                                         prim(i,j-1,k,2)-prim(i-1,j-1,k,2) + prim(i,j,k,2)-prim(i-1,j,k,2))/DX[0];
            cornwx(i,j,k) = 0.25*muxp*(prim(i,j-1,k-1,3)-prim(i-1,j-1,k-1,3) + prim(i,j,k-1,3)-prim(i-1,j,k-1,3)+
                                         prim(i,j-1,k,3)-prim(i-1,j-1,k,3) + prim(i,j,k,3)-prim(i-1,j,k,3))/DX[0];

            cornuy(i,j,k) = 0.25*muxp* (prim(i-1,j,k-1,1)-prim(i-1,j-1,k-1,1) + prim(i,j,k-1,1)-prim(i,j-1,k-1,1) +
                                          prim(i-1,j,k,1)-prim(i-1,j-1,k,1) + prim(i,j,k,1)-prim(i,j-1,k,1))/DX[1];
            cornvy(i,j,k) = 0.25*muxp* (prim(i-1,j,k-1,2)-prim(i-1,j-1,k-1,2) + prim(i,j,k-1,2)-prim(i,j-1,k-1,2) +
                                          prim(i-1,j,k,2)-prim(i-1,j-1,k,2) + prim(i,j,k,2)-prim(i,j-1,k,2))/DX[1];
            cornwy(i,j,k) = 0.25*muxp* (prim(i-1,j,k-1,3)-prim(i-1,j-1,k-1,3) + prim(i,j,k-1,3)-prim(i,j-1,k-1,3) +
                                          prim(i-1,j,k,3)-prim(i-1,j-1,k,3) + prim(i,j,k,3)-prim(i,j-1,k,3))/DX[1];

            cornuz(i,j,k) = 0.25*muxp*(prim(i-1,j-1,k,1)-prim(i-1,j-1,k-1,1) + prim(i,j-1,k,1)-prim(i,j-1,k-1,1) +
                                         prim(i-1,j,k,1)-prim(i-1,j,k-1,1) + prim(i,j,k,1)-prim(i,j,k-1,1))/DX[2];
            cornvz(i,j,k) = 0.25*muxp*(prim(i-1,j-1,k,2)-prim(i-1,j-1,k-1,2) + prim(i,j-1,k,2)-prim(i,j-1,k-1,2) +
                                         prim(i-1,j,k,2)-prim(i-1,j,k-1,2) + prim(i,j,k,2)-prim(i,j,k-1,2))/DX[2];
            cornwz(i,j,k) = 0.25*muxp*(prim(i-1,j-1,k,3)-prim(i-1,j-1,k-1,3) + prim(i,j-1,k,3)-prim(i,j-1,k-1,3) +
                                         prim(i-1,j,k,3)-prim(i-1,j,k-1,3) + prim(i,j,k,3)-prim(i,j,k-1,3))/DX[2];

            visccorn(i,j,k) =  (muxp/12.+zetaxp/4.)*( // Divergence stress
                (prim(i,  j-1,k-1,1)-prim(i-1,j-1,k-1,1))/DX[0] + (prim(i,j,  k-1,1)-prim(i-1,j  ,k-1,1))/DX[0] +
                (prim(i,  j-1,k  ,1)-prim(i-1,j-1,k,  1))/DX[0] + (prim(i,j,  k,  1)-prim(i-1,j  ,k,  1))/DX[0] +
                (prim(i-1,j  ,k-1,2)-prim(i-1,j-1,k-1,2))/DX[1] + (prim(i,j,  k-1,2)-prim(i  ,j-1,k-1,2))/DX[1] +
                (prim(i-1,j  ,k  ,2)-prim(i-1,j-1,k  ,2))/DX[1] + (prim(i,j,  k,  2)-prim(i  ,j-1,k,  2))/DX[1] +
                (prim(i-1,j-1,k  ,3)-prim(i-1,j-1,k-1,3))/DX[2] + (prim(i,j-1,k,  3)-prim(i  ,j-1,k-1,3))/DX[2] +
                (prim(i-1,j  ,k  ,3)-prim(i-1,j  ,k-1,3))/DX[2] + (prim(i,j,  k,  3)-prim(i  ,j  ,k-1,3))/DX[2]);

        });

        } else if (n_cells_z == 1) {

            Abort("diffusive flux n_cells_z == 1 case not converted yet");

/* OLD FORTRAN CODE TO CONVERT
       do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
       do i = lo(1),hi(1)+1

          ! Corner viscosity
          muxp = 0.25d0*(eta(i,j-1,k) + eta(i-1,j-1,k) + eta(i,j,k) + eta(i-1,j,k))
          if (abs(visc_type) .eq. 3) then
             zetaxp = 0.25d0*(zeta(i,j-1,k) + zeta(i-1,j-1,k) + zeta(i,j,k) + zeta(i-1,j,k))
          else
             zetaxp = 0.0
          endif

          cornux(i,j,k) = 0.5d0*muxp*(prim(i,j-1,k,2)-prim(i-1,j-1,k,2) + prim(i,j,k,2)-prim(i-1,j,k,2))/dx(1)
          cornvx(i,j,k) = 0.5d0*muxp*(prim(i,j-1,k,3)-prim(i-1,j-1,k,3) + prim(i,j,k,3)-prim(i-1,j,k,3))/dx(1)
          cornwx(i,j,k) = 0.5d0*muxp*(prim(i,j-1,k,4)-prim(i-1,j-1,k,4) + prim(i,j,k,4)-prim(i-1,j,k,4))/dx(1)

          cornuy(i,j,k) = 0.5d0*muxp* (prim(i-1,j,k,2)-prim(i-1,j-1,k,2) + prim(i,j,k,2)-prim(i,j-1,k,2))/dx(2)
          cornvy(i,j,k) = 0.5d0*muxp* (prim(i-1,j,k,3)-prim(i-1,j-1,k,3) + prim(i,j,k,3)-prim(i,j-1,k,3))/dx(2)
          cornwy(i,j,k) = 0.5d0*muxp* (prim(i-1,j,k,4)-prim(i-1,j-1,k,4) + prim(i,j,k,4)-prim(i,j-1,k,4))/dx(2)

          cornuz(i,j,k) = 0.d0
          cornvz(i,j,k) = 0.d0
          cornwz(i,j,k) = 0.d0

          visccorn(i,j,k) =  (muxp/6d0+zetaxp/2d0)*( & ! Divergence stress
               (prim(i,j-1,k,2)-prim(i-1,j-1,k,2))/dx(1) + (prim(i,j,k,2)-prim(i-1,j,k,2))/dx(1) + &
               (prim(i-1,j,k,3)-prim(i-1,j-1,k,3))/dx(2) + (prim(i,j,k,3)-prim(i,j-1,k,3))/dx(2))

          ! Copy along z direction
          cornux(i,j,k+1) = cornux(i,j,k)
          cornvx(i,j,k+1) = cornvx(i,j,k)
          cornwx(i,j,k+1) = cornwx(i,j,k)

          cornuy(i,j,k+1) = cornuy(i,j,k)
          cornvy(i,j,k+1) = cornvy(i,j,k)
          cornwy(i,j,k+1) = cornwy(i,j,k)

          cornuz(i,j,k+1) = cornuz(i,j,k)
          cornvz(i,j,k+1) = cornvz(i,j,k)
          cornwz(i,j,k+1) = cornwz(i,j,k)

          visccorn(i,j,k+1) = visccorn(i,j,k)

       end do
       end do
       end do
*/

        } // n_cells_z test

        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            fluxx(i,j,k,1) = fluxx(i,j,k,1) - 0.25*(visccorn(i,j+1,k+1)+visccorn(i,j,k+1) +
                                                      visccorn(i,j+1,k)+visccorn(i,j,k)); // Viscous "divergence" stress

            fluxx(i,j,k,1) = fluxx(i,j,k,1) + .25*
                (cornvy(i,j+1,k+1)+cornvy(i,j,k+1)+cornvy(i,j+1,k)+cornvy(i,j,k)  +
                 cornwz(i,j+1,k+1)+cornwz(i,j,k+1)+cornwz(i,j+1,k)+cornwz(i,j,k));

            fluxx(i,j,k,2) = fluxx(i,j,k,2) - .25*
                (cornuy(i,j+1,k+1)+cornuy(i,j,k+1)+cornuy(i,j+1,k)+cornuy(i,j,k));

            fluxx(i,j,k,3) = fluxx(i,j,k,3) - .25*
                (cornuz(i,j+1,k+1)+cornuz(i,j,k+1)+cornuz(i,j+1,k)+cornuz(i,j,k));

            Real phiflx;

            if ((i == 0) and is_lo_x_dirichlet_mass) {
                phiflx =  0.5*(visccorn(i,j+1,k+1)+visccorn(i,j,k+1) +
                            visccorn(i,j+1,k)+visccorn(i,j,k)
                            -(cornvy(i,j+1,k+1)+cornvy(i,j,k+1)+cornvy(i,j+1,k)+cornvy(i,j,k)  +
                              cornwz(i,j+1,k+1)+cornwz(i,j,k+1)+cornwz(i,j+1,k)+cornwz(i,j,k))) *
                            (prim(i-1,j,k,1));

                phiflx = phiflx + .5*
                            (cornuy(i,j+1,k+1)+cornuy(i,j,k+1)+cornuy(i,j+1,k)+cornuy(i,j,k)) *
                            (prim(i-1,j,k,2));

                phiflx = phiflx + .5*
                            (cornuz(i,j+1,k+1)+cornuz(i,j,k+1)+cornuz(i,j+1,k)+cornuz(i,j,k)) *
                            (prim(i-1,j,k,3));

            }
            else if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                phiflx =  0.5*(visccorn(i,j+1,k+1)+visccorn(i,j,k+1) +
                            visccorn(i,j+1,k)+visccorn(i,j,k)
                            -(cornvy(i,j+1,k+1)+cornvy(i,j,k+1)+cornvy(i,j+1,k)+cornvy(i,j,k)  +
                              cornwz(i,j+1,k+1)+cornwz(i,j,k+1)+cornwz(i,j+1,k)+cornwz(i,j,k))) *
                            (prim(i,j,k,1));

                phiflx = phiflx + .5*
                            (cornuy(i,j+1,k+1)+cornuy(i,j,k+1)+cornuy(i,j+1,k)+cornuy(i,j,k)) *
                            (prim(i,j,k,2));

                phiflx = phiflx + .5*
                            (cornuz(i,j+1,k+1)+cornuz(i,j,k+1)+cornuz(i,j+1,k)+cornuz(i,j,k)) *
                            (prim(i,j,k,3));

            }
            else {
                phiflx =  0.25*(visccorn(i,j+1,k+1)+visccorn(i,j,k+1) +
                            visccorn(i,j+1,k)+visccorn(i,j,k)
                            -(cornvy(i,j+1,k+1)+cornvy(i,j,k+1)+cornvy(i,j+1,k)+cornvy(i,j,k)  +
                              cornwz(i,j+1,k+1)+cornwz(i,j,k+1)+cornwz(i,j+1,k)+cornwz(i,j,k))) *
                            (prim(i-1,j,k,1)+prim(i,j,k,1));

                phiflx = phiflx + .25*
                            (cornuy(i,j+1,k+1)+cornuy(i,j,k+1)+cornuy(i,j+1,k)+cornuy(i,j,k)) *
                            (prim(i-1,j,k,2)+prim(i,j,k,2));

                phiflx = phiflx + .25*
                            (cornuz(i,j+1,k+1)+cornuz(i,j,k+1)+cornuz(i,j+1,k)+cornuz(i,j,k)) *
                            (prim(i-1,j,k,3)+prim(i,j,k,3));

            }

            fluxx(i,j,k,nvars+1) = fluxx(i,j,k,nvars+1)-0.5*phiflx;
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            fluxy(i,j,k,2) = fluxy(i,j,k,2) -
                0.25*(visccorn(i+1,j,k+1)+visccorn(i,j,k+1)+visccorn(i+1,j,k)+visccorn(i,j,k));

            fluxy(i,j,k,2) = fluxy(i,j,k,2) + .25*
                (cornux(i+1,j,k+1)+cornux(i,j,k+1)+cornux(i+1,j,k)+cornux(i,j,k)  +
                 cornwz(i+1,j,k+1)+cornwz(i,j,k+1)+cornwz(i+1,j,k)+cornwz(i,j,k));

            fluxy(i,j,k,1) = fluxy(i,j,k,1) - .25*
                (cornvx(i+1,j,k+1)+cornvx(i,j,k+1)+cornvx(i+1,j,k)+cornvx(i,j,k));

            fluxy(i,j,k,3) = fluxy(i,j,k,3) - .25*
                (cornvz(i+1,j,k+1)+cornvz(i,j,k+1)+cornvz(i+1,j,k)+cornvz(i,j,k));

            Real phiflx;

            if ((j == 0) and is_lo_y_dirichlet_mass) {

                phiflx = 0.5*(visccorn(i+1,j,k+1)+visccorn(i,j,k+1)+visccorn(i+1,j,k)+visccorn(i,j,k)
                               -(cornux(i+1,j,k+1)+cornux(i,j,k+1)+cornux(i+1,j,k)+cornux(i,j,k)  +
                                 cornwz(i+1,j,k+1)+cornwz(i,j,k+1)+cornwz(i+1,j,k)+cornwz(i,j,k))) *
                              (prim(i,j-1,k,2));

                phiflx = phiflx + .5*
                            (cornvx(i+1,j,k+1)+cornvx(i,j,k+1)+cornvx(i+1,j,k)+cornvx(i,j,k)) *
                            (prim(i,j-1,k,1));

                phiflx = phiflx + .5*
                            (cornvz(i+1,j,k+1)+cornvz(i,j,k+1)+cornvz(i+1,j,k)+cornvz(i,j,k)) *
                            (prim(i,j-1,k,3));

            }
            else if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {

                phiflx = 0.5*(visccorn(i+1,j,k+1)+visccorn(i,j,k+1)+visccorn(i+1,j,k)+visccorn(i,j,k)
                               -(cornux(i+1,j,k+1)+cornux(i,j,k+1)+cornux(i+1,j,k)+cornux(i,j,k)  +
                                 cornwz(i+1,j,k+1)+cornwz(i,j,k+1)+cornwz(i+1,j,k)+cornwz(i,j,k))) *
                              (prim(i,j,k,2));

                phiflx = phiflx + .5*
                            (cornvx(i+1,j,k+1)+cornvx(i,j,k+1)+cornvx(i+1,j,k)+cornvx(i,j,k)) *
                            (prim(i,j,k,1));

                phiflx = phiflx + .5*
                            (cornvz(i+1,j,k+1)+cornvz(i,j,k+1)+cornvz(i+1,j,k)+cornvz(i,j,k)) *
                            (prim(i,j,k,3));

            }
            else {
                phiflx = 0.25*(visccorn(i+1,j,k+1)+visccorn(i,j,k+1)+visccorn(i+1,j,k)+visccorn(i,j,k)
                               -(cornux(i+1,j,k+1)+cornux(i,j,k+1)+cornux(i+1,j,k)+cornux(i,j,k)  +
                                 cornwz(i+1,j,k+1)+cornwz(i,j,k+1)+cornwz(i+1,j,k)+cornwz(i,j,k))) *
                              (prim(i,j-1,k,2)+prim(i,j,k,2));

                phiflx = phiflx + .25*
                            (cornvx(i+1,j,k+1)+cornvx(i,j,k+1)+cornvx(i+1,j,k)+cornvx(i,j,k)) *
                            (prim(i,j-1,k,1)+prim(i,j,k,1));

                phiflx = phiflx + .25*
                            (cornvz(i+1,j,k+1)+cornvz(i,j,k+1)+cornvz(i+1,j,k)+cornvz(i,j,k)) *
                            (prim(i,j-1,k,3)+prim(i,j,k,3));

            }

            fluxy(i,j,k,nvars+1) = fluxy(i,j,k,nvars+1)-0.5*phiflx;

        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            if (n_cells_z > 1) {

            fluxz(i,j,k,3) = fluxz(i,j,k,3) -
                0.25*(visccorn(i+1,j+1,k)+visccorn(i,j+1,k)+visccorn(i+1,j,k)+visccorn(i,j,k));

            fluxz(i,j,k,3) = fluxz(i,j,k,3) + .25*
                (cornvy(i+1,j+1,k)+cornvy(i+1,j,k)+cornvy(i,j+1,k)+cornvy(i,j,k)  +
                 cornux(i+1,j+1,k)+cornux(i+1,j,k)+cornux(i,j+1,k)+cornux(i,j,k));

            fluxz(i,j,k,1) = fluxz(i,j,k,1) - .25*
                (cornwx(i+1,j+1,k)+cornwx(i+1,j,k)+cornwx(i,j+1,k)+cornwx(i,j,k));

            fluxz(i,j,k,2) = fluxz(i,j,k,2) - .25*
                (cornwy(i+1,j+1,k)+cornwy(i+1,j,k)+cornwy(i,j+1,k)+cornwy(i,j,k));

            Real phiflx;

            if ((k == 0) and is_lo_z_dirichlet_mass) {

                phiflx = 0.5*(visccorn(i+1,j+1,k)+visccorn(i,j+1,k)+visccorn(i+1,j,k)+visccorn(i,j,k)
                               -(cornvy(i+1,j+1,k)+cornvy(i+1,j,k)+cornvy(i,j+1,k)+cornvy(i,j,k)  +
                                 cornux(i+1,j+1,k)+cornux(i+1,j,k)+cornux(i,j+1,k)+cornux(i,j,k))) *
                                (prim(i,j,k-1,3));

                phiflx = phiflx + .5*
                            (cornwx(i+1,j+1,k)+cornwx(i+1,j,k)+cornwx(i,j+1,k)+cornwx(i,j,k))*
                            (prim(i,j,k-1,1));

                phiflx = phiflx + .5*
                            (cornwy(i+1,j+1,k)+cornwy(i+1,j,k)+cornwy(i,j+1,k)+cornwy(i,j,k)) *
                            (prim(i,j,k-1,2));

            }
            else if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {

                phiflx = 0.5*(visccorn(i+1,j+1,k)+visccorn(i,j+1,k)+visccorn(i+1,j,k)+visccorn(i,j,k)
                               -(cornvy(i+1,j+1,k)+cornvy(i+1,j,k)+cornvy(i,j+1,k)+cornvy(i,j,k)  +
                                 cornux(i+1,j+1,k)+cornux(i+1,j,k)+cornux(i,j+1,k)+cornux(i,j,k))) *
                                (prim(i,j,k,3));

                phiflx = phiflx + .5*
                            (cornwx(i+1,j+1,k)+cornwx(i+1,j,k)+cornwx(i,j+1,k)+cornwx(i,j,k))*
                            (prim(i,j,k,1));

                phiflx = phiflx + .5*
                            (cornwy(i+1,j+1,k)+cornwy(i+1,j,k)+cornwy(i,j+1,k)+cornwy(i,j,k)) *
                            (prim(i,j,k,2));

            }
            else {
                phiflx = 0.25*(visccorn(i+1,j+1,k)+visccorn(i,j+1,k)+visccorn(i+1,j,k)+visccorn(i,j,k)
                               -(cornvy(i+1,j+1,k)+cornvy(i+1,j,k)+cornvy(i,j+1,k)+cornvy(i,j,k)  +
                                 cornux(i+1,j+1,k)+cornux(i+1,j,k)+cornux(i,j+1,k)+cornux(i,j,k))) *
                                (prim(i,j,k-1,3)+prim(i,j,k,3));

                phiflx = phiflx + .25*
                            (cornwx(i+1,j+1,k)+cornwx(i+1,j,k)+cornwx(i,j+1,k)+cornwx(i,j,k))*
                            (prim(i,j,k-1,1)+prim(i,j,k,1));

                phiflx = phiflx + .25*
                            (cornwy(i+1,j+1,k)+cornwy(i+1,j,k)+cornwy(i,j+1,k)+cornwy(i,j,k)) *
                            (prim(i,j,k-1,2)+prim(i,j,k,2));

            }

            fluxz(i,j,k,nvars+1) = fluxz(i,j,k,nvars+1)-0.5*phiflx;

            }

        });

    }

    // Set species flux to zero at the walls (also Dufour)
    BCWallSpeciesFlux(flux_in,geom);

    ////////////////////
    // hyperbolic fluxes
    ////////////////////

    Real wgt2 = 1./12.;
    Real wgt1 = 0.5 + wgt2;

    Real wgta = -0.2;
    Real wgtb = 0.75;
    Real wgtc = 0.5;
    Real wgtd = -0.05;

    // Loop over boxes
    for ( MFIter mfi(cons_in); mfi.isValid(); ++mfi) {

        AMREX_D_TERM(const Array4<Real>& xflux = flux_in[0].array(mfi); ,
                     const Array4<Real>& yflux = flux_in[1].array(mfi); ,
                     const Array4<Real>& zflux = flux_in[2].array(mfi));

        const Array4<const Real> prim = prim_in.array(mfi);
        const Array4<const Real> cons = cons_in.array(mfi);

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        if (advection_type == 1) { // interpolate primitive quantities

            // Loop over the cells and compute fluxes
            amrex::ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                GpuArray<Real,MAX_SPECIES+5> conserved;
                GpuArray<Real,MAX_SPECIES+6> primitive;
                GpuArray<Real,MAX_SPECIES  > Yk;

                for (int l=0; l<nspecies+6; ++l) {
                    primitive[l] = wgt1*(prim(i,j,k,l)+prim(i-1,j,k,l)) - wgt2*(prim(i-2,j,k,l)+prim(i+1,j,k,l));
                }
                if ((i == 0) and is_lo_x_dirichlet_mass) {
                    for (int l=0; l<nspecies+6; ++l) {
                        primitive[l] = prim(i-1,j,k,l);
                    }
                }
                if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                    for (int l=0; l<nspecies+6; ++l) {
                        primitive[l] = prim(i,j,k,l);
                    }
                }

                Real temp = primitive[4];
                Real rho = primitive[0];
                conserved[0] = rho;

                // want sum of specden == rho
                for (int n=0; n<nspecies; ++n) {
                    Yk[n] = primitive[6+n];
                }

                Real intenergy;
                GetEnergy(intenergy, Yk, temp);

                Real vsqr = primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3];

                conserved[4] = rho*intenergy + 0.5*rho*vsqr;

                xflux(i,j,k,0) += conserved[0]*primitive[1];
                xflux(i,j,k,1) += conserved[0]*(primitive[1]*primitive[1])+primitive[5];
                xflux(i,j,k,2) += conserved[0]*primitive[1]*primitive[2];
                xflux(i,j,k,3) += conserved[0]*primitive[1]*primitive[3];

                xflux(i,j,k,4) += primitive[1]*conserved[4] + primitive[5]*primitive[1];

                // also add the diffusive + stochastic contributions from heat flux, viscous heating and Dufour effects
                xflux(i,j,k,4) += xflux(i,j,k,nvars) + xflux(i,j,k,nvars+1) + xflux(i,j,k,nvars+2) + xflux(i,j,k,nvars+3);

                if (algorithm_type == 2) { // Add advection of concentration
                    for (int n=0; n<nspecies; ++n) {
                        xflux(i,j,k,5+n) += rho*primitive[6+n]*primitive[1];
                    }
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                GpuArray<Real,MAX_SPECIES+5> conserved;
                GpuArray<Real,MAX_SPECIES+6> primitive;
                GpuArray<Real,MAX_SPECIES  > Yk;

                for (int l=0; l<nspecies+6; ++l) {
                    primitive[l] = wgt1*(prim(i,j,k,l)+prim(i,j-1,k,l)) - wgt2*(prim(i,j-2,k,l)+prim(i,j+1,k,l));
                }
                if ((j == 0) and is_lo_y_dirichlet_mass) {
                    for (int l=0; l<nspecies+6; ++l) {
                        primitive[l] = prim(i,j-1,k,l);
                    }
                }
                if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                    for (int l=0; l<nspecies+6; ++l) {
                        primitive[l] = prim(i,j,k,l);
                    }
                }

                Real temp = primitive[4];
                Real rho = primitive[0];
                conserved[0] = rho;

                // want sum of specden == rho
                for (int n=0; n<nspecies; ++n) {
                    Yk[n] = primitive[6+n];
                }

                Real intenergy;
                GetEnergy(intenergy, Yk, temp);

                Real vsqr = primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3];

                conserved[4] = rho*intenergy + 0.5*rho*vsqr;

                yflux(i,j,k,0) += conserved[0]*primitive[2];
                yflux(i,j,k,1) += conserved[0]*primitive[1]*primitive[2];
                yflux(i,j,k,2) += conserved[0]*primitive[2]*primitive[2]+primitive[5];
                yflux(i,j,k,3) += conserved[0]*primitive[3]*primitive[2];

                yflux(i,j,k,4) += primitive[2]*conserved[4] + primitive[5]*primitive[2];

                // also add the diffusive + stochastic contributions from heat flux, viscous heating and Dufour effects
                yflux(i,j,k,4) += yflux(i,j,k,nvars) + yflux(i,j,k,nvars+1) + yflux(i,j,k,nvars+2) + yflux(i,j,k,nvars+3);

                if (algorithm_type == 2) { // Add advection of concentration
                    for (int n=0; n<nspecies; ++n) {
                        yflux(i,j,k,5+n) += rho*primitive[6+n]*primitive[2];
                    }
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                GpuArray<Real,MAX_SPECIES+5> conserved;
                GpuArray<Real,MAX_SPECIES+6> primitive;
                GpuArray<Real,MAX_SPECIES  > Yk;

                for (int l=0; l<nspecies+6; ++l) {
                    primitive[l] = wgt1*(prim(i,j,k,l)+prim(i,j,k-1,l)) - wgt2*(prim(i,j,k-2,l)+prim(i,j,k+1,l));
                }
                if ((k == 0) and is_lo_z_dirichlet_mass) {
                    for (int l=0; l<nspecies+6; ++l) {
                        primitive[l] = prim(i,j,k-1,l);
                    }
                }
                if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                    for (int l=0; l<nspecies+6; ++l) {
                        primitive[l] = prim(i,j,k,l);
                    }
                }

                Real temp = primitive[4];
                Real rho = primitive[0];
                conserved[0] = rho;

                // want sum of specden == rho
                for (int n=0; n<nspecies; ++n) {
                    Yk[n] = primitive[6+n];
                }

                Real intenergy;
                GetEnergy(intenergy, Yk, temp);

                Real vsqr = primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3];

                conserved[4] = rho*intenergy + 0.5*rho*vsqr;

                zflux(i,j,k,0) += conserved[0]*primitive[3];
                zflux(i,j,k,1) += conserved[0]*primitive[1]*primitive[3];
                zflux(i,j,k,2) += conserved[0]*primitive[2]*primitive[3];
                zflux(i,j,k,3) += conserved[0]*primitive[3]*primitive[3]+primitive[5];

                zflux(i,j,k,4) += primitive[3]*conserved[4] + primitive[5]*primitive[3];

                // also add the diffusive + stochastic contributions from heat flux, viscous heating and Dufour effects
                zflux(i,j,k,4) += zflux(i,j,k,nvars) + zflux(i,j,k,nvars+1) + zflux(i,j,k,nvars+2) + zflux(i,j,k,nvars+3);

                if (algorithm_type == 2) { // Add advection of concentration
                    for (int n=0; n<nspecies; ++n) {
                        zflux(i,j,k,5+n) += rho*primitive[6+n]*primitive[3];
                    }
                }

            });

        } else if (advection_type == 2) { // interpolate conserved quantitites

            // Loop over the cells and compute fluxes
            amrex::ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                GpuArray<Real,MAX_SPECIES+5> conserved;
                GpuArray<Real,MAX_SPECIES+6> primitive;
                GpuArray<Real,MAX_SPECIES  > Yk;

                // interpolate conserved quantities to faces
                for (int l=0; l<nspecies+5; ++l) {
                    conserved[l] = wgt1*(cons(i,j,k,l)+cons(i-1,j,k,l)) - wgt2*(cons(i-2,j,k,l)+cons(i+1,j,k,l));
                }
                if ((i == 0) and is_lo_x_dirichlet_mass) {
                    for (int l=0; l<nspecies+5; ++l) {
                        conserved[l] = cons(i-1,j,k,l);
                    }
                } else if ((i == 1) and is_lo_x_dirichlet_mass) {
                    for (int l=0; l<nspecies+5; ++l) {
                        conserved[l] = wgta*cons(i-2,j,k,l) + wgtb*cons(i-1,j,k,l) + wgtc*cons(i,j,k,l) + wgtd*cons(i+1,j,k,l);
                    }
                }
                if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                    for (int l=0; l<nspecies+5; ++l) {
                        conserved[l] = cons(i,j,k,l);
                    }
                } else if ((i == n_cells[0]-1) and is_hi_x_dirichlet_mass) {
                    for (int l=0; l<nspecies+5; ++l) {
                        conserved[l] = wgta*cons(i+1,j,k,l) + wgtb*cons(i,j,k,l) + wgtc*cons(i-1,j,k,l) + wgtd*cons(i-2,j,k,l);
                    }
                }

                // compute velocities
                for (int l=1; l<4; ++l) {
                    primitive[l] = conserved[l]/conserved[0];
                }

                // want sum of specden == rho
                for (int n=0; n<nspecies; ++n) {
                    Yk[n] = conserved[5+n]/conserved[0];
                }

                // compute temperature
                Real vsqr = primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3];
                Real intenergy = conserved[4]/conserved[0] - 0.5*vsqr;
                GetTemperature(intenergy, Yk, primitive[4]);

                // compute pressure
                GetPressureGas(primitive[5], Yk, conserved[0], primitive[4]);

                xflux(i,j,k,0) += conserved[0]*primitive[1];
                xflux(i,j,k,1) += conserved[0]*(primitive[1]*primitive[1])+primitive[5];
                xflux(i,j,k,2) += conserved[0]*primitive[1]*primitive[2];
                xflux(i,j,k,3) += conserved[0]*primitive[1]*primitive[3];

                xflux(i,j,k,4) += primitive[1]*conserved[4] + primitive[5]*primitive[1];

                // also add the diffusive + stochastic contributions from heat flux, viscous heating and Dufour effects
                xflux(i,j,k,4) += xflux(i,j,k,nvars) + xflux(i,j,k,nvars+1) + xflux(i,j,k,nvars+2) + xflux(i,j,k,nvars+3);

                if (algorithm_type == 2) { // Add advection of concentration
                    for (int n=0; n<nspecies; ++n) {
                        xflux(i,j,k,5+n) += conserved[5+n]*primitive[1];
                    }
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                GpuArray<Real,MAX_SPECIES+5> conserved;
                GpuArray<Real,MAX_SPECIES+6> primitive;
                GpuArray<Real,MAX_SPECIES  > Yk;

                // interpolate conserved quantities to faces
                for (int l=0; l<nspecies+5; ++l) {
                    conserved[l] = wgt1*(cons(i,j,k,l)+cons(i,j-1,k,l)) - wgt2*(cons(i,j-2,k,l)+cons(i,j+1,k,l));
                }
                if ((j == 0) and is_lo_y_dirichlet_mass) {
                    for (int l=0; l<nspecies+5; ++l) {
                        conserved[l] = cons(i,j-1,k,l);
                    }
                } else if ((j == 1) and is_lo_y_dirichlet_mass) {
                    for (int l=0; l<nspecies+5; ++l) {
                        conserved[l] = wgta*cons(i,j-2,k,l) + wgtb*cons(i,j-1,k,l) + wgtc*cons(i,j,k,l) + wgtd*cons(i,j+1,k,l);
                    }
                }
                if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                    for (int l=0; l<nspecies+5; ++l) {
                        conserved[l] = cons(i,j,k,l);
                    }
                } else if ((j == n_cells[1]-1) and is_hi_y_dirichlet_mass) {
                    for (int l=0; l<nspecies+5; ++l) {
                        conserved[l] = wgta*cons(i,j+1,k,l) + wgtb*cons(i,j,k,l) + wgtc*cons(i,j-1,k,l) + wgtd*cons(i,j-2,k,l);
                    }
                }

                // compute velocities
                for (int l=1; l<4; ++l) {
                    primitive[l] = conserved[l]/conserved[0];
                }

                // want sum of specden == rho
                for (int n=0; n<nspecies; ++n) {
                    Yk[n] = conserved[5+n]/conserved[0];
                }

                // compute temperature
                Real vsqr = primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3];
                Real intenergy = conserved[4]/conserved[0] - 0.5*vsqr;
                GetTemperature(intenergy, Yk, primitive[4]);

                // compute pressure
                GetPressureGas(primitive[5], Yk, conserved[0], primitive[4]);

                yflux(i,j,k,0) += conserved[0]*primitive[2];
                yflux(i,j,k,1) += conserved[0]*primitive[1]*primitive[2];
                yflux(i,j,k,2) += conserved[0]*primitive[2]*primitive[2]+primitive[5];
                yflux(i,j,k,3) += conserved[0]*primitive[3]*primitive[2]  ;

                yflux(i,j,k,4) += primitive[2]*conserved[4] + primitive[5]*primitive[2];

                // also add the diffusive + stochastic contributions from heat flux, viscous heating and Dufour effects
                yflux(i,j,k,4) += yflux(i,j,k,nvars) + yflux(i,j,k,nvars+1) + yflux(i,j,k,nvars+2) + yflux(i,j,k,nvars+3);

                if (algorithm_type == 2) { // Add advection of concentration
                    for (int n=0; n<nspecies; ++n) {
                        yflux(i,j,k,5+n) += conserved[5+n]*primitive[2];
                    }
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                GpuArray<Real,MAX_SPECIES+5> conserved;
                GpuArray<Real,MAX_SPECIES+6> primitive;
                GpuArray<Real,MAX_SPECIES  > Yk;

                // interpolate conserved quantities to faces
                for (int l=0; l<nspecies+5; ++l) {
                    conserved[l] = wgt1*(cons(i,j,k,l)+cons(i,j,k-1,l)) - wgt2*(cons(i,j,k-2,l)+cons(i,j,k+1,l));
                }
                if ((k == 0) and is_lo_z_dirichlet_mass) {
                    for (int l=0; l<nspecies+5; ++l) {
                        conserved[l] = cons(i,j,k-1,l);
                    }
                } else if ((k == 1) and is_lo_z_dirichlet_mass) {
                    for (int l=0; l<nspecies+5; ++l) {
                        conserved[l] = wgta*cons(i,j,k-2,l) + wgtb*cons(i,j,k-1,l) + wgtc*cons(i,j,k,l) + wgtd*cons(i,j,k+1,l);
                    }
                }
                if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                    for (int l=0; l<nspecies+5; ++l) {
                        conserved[l] = cons(i,j,k,l);
                    }
                } else if ((k == n_cells[2]-1) and is_hi_z_dirichlet_mass) {
                    for (int l=0; l<nspecies+5; ++l) {
                        conserved[l] = wgta*cons(i,j,k+1,l) + wgtb*cons(i,j,k,l) + wgtc*cons(i,j,k-1,l) + wgtd*cons(i,j,k-2,l);
                    }
                }

                // compute velocities
                for (int l=1; l<4; ++l) {
                    primitive[l] = conserved[l]/conserved[0];
                }

                // want sum of specden == rho
                for (int n=0; n<nspecies; ++n) {
                    Yk[n] = conserved[5+n]/conserved[0];
                }

                // compute temperature
                Real vsqr = primitive[1]*primitive[1] + primitive[2]*primitive[2] + primitive[3]*primitive[3];
                Real intenergy = conserved[4]/conserved[0] - 0.5*vsqr;
                GetTemperature(intenergy, Yk, primitive[4]);

                // compute pressure
                GetPressureGas(primitive[5], Yk, conserved[0], primitive[4]);


                zflux(i,j,k,0) += conserved[0]*primitive[3];
                zflux(i,j,k,1) += conserved[0]*primitive[1]*primitive[3];
                zflux(i,j,k,2) += conserved[0]*primitive[2]*primitive[3];
                zflux(i,j,k,3) += conserved[0]*primitive[3]*primitive[3]+primitive[5];

                zflux(i,j,k,4) += primitive[3]*conserved[4] + primitive[5]*primitive[3];

                // also add the diffusive + stochastic contributions from heat flux, viscous heating and Dufour effects
                zflux(i,j,k,4) += zflux(i,j,k,nvars) + zflux(i,j,k,nvars+1) + zflux(i,j,k,nvars+2) + zflux(i,j,k,nvars+3);

                if (algorithm_type == 2) { // Add advection of concentration
                    for (int n=0; n<nspecies; ++n) {
                        zflux(i,j,k,5+n) += conserved[5+n]*primitive[3];
                    }
                }
            });

        }
    }
}
