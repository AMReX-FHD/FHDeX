#include "compressible_functions.H"
#include "compressible_functions_stag.H"
#include "common_functions.H"

void calculateFluxStag(const MultiFab& cons_in, const std::array< MultiFab, AMREX_SPACEDIM >& cumom_in,
                       const MultiFab& prim_in, const std::array< MultiFab, AMREX_SPACEDIM >& vel_in,
                       const MultiFab& eta_in, const MultiFab& zeta_in, const MultiFab& kappa_in,
                       const MultiFab& chi_in, const MultiFab& D_in,
                       std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in,
                       std::array< MultiFab, 2 >& edgeflux_x_in,
                       std::array< MultiFab, 2 >& edgeflux_y_in,
                       std::array< MultiFab, 2 >& edgeflux_z_in,
                       std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in,
                       std::array< MultiFab, AMREX_SPACEDIM>& stochface_in,
                       std::array< MultiFab, 2 >& stochedge_x_in,
                       std::array< MultiFab, 2 >& stochedge_y_in,
                       std::array< MultiFab, 2 >& /*stochedge_z_in*/,
                       std::array< MultiFab, AMREX_SPACEDIM>& stochcen_in,
                       const amrex::Geometry& geom,
                       const amrex::Vector< amrex::Real >& /*stoch_weights*/,
                       const amrex::Real dt)
{
    BL_PROFILE_VAR("calculateFluxStag()",calculateFluxStag);

    Box dom(geom.Domain());
    int n_cells_z = n_cells[2];

    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

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

    std::array< MultiFab, AMREX_SPACEDIM > tau_diag_stoch; // diagonal stochastic stress (defined at cell centers)
    AMREX_D_TERM(tau_diag_stoch[0].define(cons_in.boxArray(),cons_in.DistributionMap(),1,ngc);,
                 tau_diag_stoch[1].define(cons_in.boxArray(),cons_in.DistributionMap(),1,ngc);,
                 tau_diag_stoch[2].define(cons_in.boxArray(),cons_in.DistributionMap(),1,ngc););

    std::array< MultiFab, NUM_EDGE > tau_diagoff; // off diagonal stress at edges
    AMREX_D_TERM(tau_diagoff[0].define(convert(cons_in.boxArray(),nodal_flag_xy),cons_in.DistributionMap(),1,ngc);,
                 tau_diagoff[1].define(convert(cons_in.boxArray(),nodal_flag_yz),cons_in.DistributionMap(),1,ngc);,
                 tau_diagoff[2].define(convert(cons_in.boxArray(),nodal_flag_xz),cons_in.DistributionMap(),1,ngc););

    std::array< MultiFab, NUM_EDGE > tau_diagoff_stoch; // off diagonal stochastic stress at edges
    AMREX_D_TERM(tau_diagoff_stoch[0].define(convert(cons_in.boxArray(),nodal_flag_xy),cons_in.DistributionMap(),1,ngc);,
                 tau_diagoff_stoch[1].define(convert(cons_in.boxArray(),nodal_flag_yz),cons_in.DistributionMap(),1,ngc);,
                 tau_diagoff_stoch[2].define(convert(cons_in.boxArray(),nodal_flag_xz),cons_in.DistributionMap(),1,ngc););

    AMREX_D_TERM(tau_diag[0].setVal(0.0);,
                 tau_diag[1].setVal(0.0);,
                 tau_diag[2].setVal(0.0););

    AMREX_D_TERM(tau_diag_stoch[0].setVal(0.0);,
                 tau_diag_stoch[1].setVal(0.0);,
                 tau_diag_stoch[2].setVal(0.0););

    AMREX_D_TERM(tau_diagoff[0].setVal(0.0);,
                 tau_diagoff[1].setVal(0.0);,
                 tau_diagoff[2].setVal(0.0););

    AMREX_D_TERM(tau_diagoff_stoch[0].setVal(0.0);,
                 tau_diagoff_stoch[1].setVal(0.0);,
                 tau_diagoff_stoch[2].setVal(0.0););

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

            const Array4<Real> tauxx_stoch = tau_diag_stoch[0].array(mfi);
            const Array4<Real> tauyy_stoch = tau_diag_stoch[1].array(mfi);
            const Array4<Real> tauzz_stoch = tau_diag_stoch[2].array(mfi);

            AMREX_D_TERM(const Array4<Real> tauxy_stoch = tau_diagoff_stoch[0].array(mfi);,
                         const Array4<Real> tauyz_stoch = tau_diagoff_stoch[1].array(mfi);,
                         const Array4<Real> tauxz_stoch = tau_diagoff_stoch[2].array(mfi););

            AMREX_D_TERM(const Array4<Real>& stochfacex = stochface_in[0].array(mfi); ,
                         const Array4<Real>& stochfacey = stochface_in[1].array(mfi); ,
                         const Array4<Real>& stochfacez = stochface_in[2].array(mfi));

            const Array4<Real>& stochedgex_v = stochedge_x_in[0].array(mfi);
            const Array4<Real>& stochedgex_w = stochedge_x_in[1].array(mfi);
            const Array4<Real>& stochedgey_w = stochedge_y_in[1].array(mfi);

            const Array4<Real>& stochcenx_u = stochcen_in[0].array(mfi);
            const Array4<Real>& stochceny_v = stochcen_in[1].array(mfi);
            const Array4<Real>& stochcenz_w = stochcen_in[2].array(mfi);

            AMREX_D_TERM(Array4<Real const> const& velx = vel_in[0].array(mfi);,
                         Array4<Real const> const& vely = vel_in[1].array(mfi);,
                         Array4<Real const> const& velz = vel_in[2].array(mfi););

            const Array4<const Real> prim = prim_in.array(mfi);

            const Array4<const Real> eta   = eta_in.array(mfi);
            const Array4<const Real> zeta  = zeta_in.array(mfi);
            const Array4<const Real> kappa = kappa_in.array(mfi);
            const Array4<const Real> chi   = chi_in.array(mfi);
            const Array4<const Real> Dij   = D_in.array(mfi);

            const Box& tbx = mfi.nodaltilebox(0);
            const Box& tby = mfi.nodaltilebox(1);
            const Box& tbz = mfi.nodaltilebox(2);

            const Box & bx_xy = mfi.tilebox(nodal_flag_xy);
            #if (AMREX_SPACEDIM == 3)
            const Box & bx_xz = mfi.tilebox(nodal_flag_xz);
            const Box & bx_yz = mfi.tilebox(nodal_flag_yz);
            #endif

            const Box& bx = mfi.growntilebox(1);

            // Populate diagonal stochastic stress
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                Real etaT = eta(i,j,k) * prim(i,j,k,4);
                Real zetaT = zeta(i,j,k) * prim(i,j,k,4);

                Real fac1 = sqrt(2.0 * k_B * etaT * volinv * dtinv);
                Real fac2 =  (-1.0/3.0)*sqrt(2.0 * k_B * etaT * volinv * dtinv);

                if (do_1D) { // 1D

                    fac1 *= sqrt(3.0);
                    fac2 *= sqrt(3.0);

                    Real traceZ = stochcenx_u(i,j,k);

                    tauxx_stoch(i,j,k) = (fac1 * stochcenx_u(i,j,k)) + (fac2 * traceZ);
                    tauyy_stoch(i,j,k) = 0.0;
                    tauzz_stoch(i,j,k) = 0.0;
                }

                else if (do_2D) { // 2D

                    fac2 *= (3.0 + sqrt(3.0))/2.0;

                    Real traceZ = stochcenx_u(i,j,k) + stochceny_v(i,j,k);

                    tauxx_stoch(i,j,k) = (fac1 * stochcenx_u(i,j,k)) + (fac2 * traceZ);
                    tauyy_stoch(i,j,k) = (fac1 * stochceny_v(i,j,k)) + (fac2 * traceZ);
                    tauzz_stoch(i,j,k) = 0.0;
                }

                else { // 3D

                    if (amrex::Math::abs(visc_type) == 3) {
                      fac2 = sqrt(k_B * zetaT * volinv * dtinv / 3.0) - sqrt(2.0 * k_B * etaT * volinv * dtinv)/3.0;
                    }

                    Real traceZ = stochcenx_u(i,j,k) + stochceny_v(i,j,k) + stochcenz_w(i,j,k);

                    tauxx_stoch(i,j,k) = (fac1 * stochcenx_u(i,j,k)) + (fac2 * traceZ);
                    tauyy_stoch(i,j,k) = (fac1 * stochceny_v(i,j,k)) + (fac2 * traceZ);
                    tauzz_stoch(i,j,k) = (fac1 * stochcenz_w(i,j,k)) + (fac2 * traceZ);
                }
            });

            // Populate off-diagonal stress
            amrex::ParallelFor(bx_xy, bx_xz, bx_yz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                if (do_1D) { // 1D
                    tauxy_stoch(i,j,k) = 0.0;
                }
                else { // works for both 2D and 3D
                    Real etaT = 0.25*(eta(i-1,j-1,k)*prim(i-1,j-1,k,4) + eta(i-1,j,k)*prim(i-1,j,k,4) +
                                      eta(i,j-1,k)*prim(i,j-1,k,4) + eta(i,j,k)*prim(i,j,k,4));

                    // Pick boundary values for Dirichlet (stored in ghost)
                    // For corner cases (xy), x wall takes preference
                    if ((j == 0) and is_lo_y_dirichlet_mass) {
                        etaT = 0.5*(eta(i-1,j-1,k)*prim(i-1,j-1,k,4) + eta(i,j-1,k)*prim(i,j-1,k,4));
                    }
                    if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                        etaT = 0.5*(eta(i-1,j,k)*prim(i-1,j,k,4) + eta(i,j,k)*prim(i,j,k,4));
                    }
                    if ((i == 0) and is_lo_x_dirichlet_mass) {
                        etaT = 0.5*(eta(i-1,j-1,k)*prim(i-1,j-1,k,4) + eta(i-1,j,k)*prim(i-1,j,k,4));
                    }
                    if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                        etaT = 0.5*(eta(i,j-1,k)*prim(i,j-1,k,4) + eta(i,j,k)*prim(i,j,k,4));
                    }

                    Real fac = sqrt(2.0 * k_B * etaT * volinv * dtinv);
                    tauxy_stoch(i,j,k) = fac*stochedgex_v(i,j,k);
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                if ((do_1D) or (do_2D)) { // works for 1D and 2D
                    tauxz_stoch(i,j,k) = 0.0;
                }
                else {
                    Real etaT = 0.25*(eta(i-1,j,k-1)*prim(i-1,j,k-1,4) + eta(i-1,j,k)*prim(i-1,j,k,4) +
                                      eta(i,j,k-1)*prim(i,j,k-1,4) + eta(i,j,k)*prim(i,j,k,4));

                    // Pick boundary values for Dirichlet (stored in ghost)
                    // For corner cases (xz), x wall takes preference
                    if ((k == 0) and is_lo_z_dirichlet_mass) {
                        etaT = 0.5*(eta(i-1,j,k-1)*prim(i-1,j,k-1,4) + eta(i,j,k-1)*prim(i,j,k-1,4));
                    }
                    if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                        etaT = 0.5*(eta(i-1,j,k)*prim(i-1,j,k,4) + eta(i,j,k)*prim(i,j,k,4));
                    }
                    if ((i == 0) and is_lo_x_dirichlet_mass) {
                        etaT = 0.5*(eta(i-1,j,k-1)*prim(i-1,j,k-1,4) + eta(i-1,j,k)*prim(i-1,j,k,4));
                    }
                    if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                        etaT = 0.5*(eta(i,j,k-1)*prim(i,j,k-1,4) + eta(i,j,k)*prim(i,j,k,4));
                    }

                    Real fac = sqrt(2.0 * k_B * etaT * volinv * dtinv);
                    tauxz_stoch(i,j,k) = fac*stochedgex_w(i,j,k);
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                if ((do_1D) or (do_2D)) { // works for 1D and 2D
                    tauyz_stoch(i,j,k) = 0.0;
                }
                else {
                    Real etaT = 0.25*(eta(i,j-1,k-1)*prim(i,j-1,k-1,4) + eta(i,j-1,k)*prim(i,j-1,k,4) +
                                      eta(i,j,k-1)*prim(i,j,k-1,4) + eta(i,j,k)*prim(i,j,k,4));

                    // Pick boundary values for Dirichlet (stored in ghost)
                    // For corner cases (yz), y wall takes preference
                    if ((k == 0) and is_lo_z_dirichlet_mass) {
                        etaT = 0.5*(eta(i,j-1,k-1)*prim(i,j-1,k-1,4) + eta(i,j,k-1)*prim(i,j,k-1,4));
                    }
                    if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                        etaT = 0.5*(eta(i,j-1,k)*prim(i,j-1,k,4) + eta(i,j,k)*prim(i,j,k,4));
                    }
                    if ((j == 0) and is_lo_y_dirichlet_mass) {
                        etaT = 0.5*(eta(i,j-1,k-1)*prim(i,j-1,k-1,4) + eta(i,j-1,k)*prim(i,j-1,k,4));
                    }
                    if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                        etaT = 0.5*(eta(i,j,k-1)*prim(i,j,k-1,4) + eta(i,j,k)*prim(i,j,k,4));
                    }

                    Real fac = sqrt(2.0 * k_B * etaT * volinv * dtinv);
                    tauyz_stoch(i,j,k) = fac*stochedgey_w(i,j,k);
                }
            });

            // Loop over faces for flux calculations (4:5+ns)
            amrex::ParallelFor(tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {


                GpuArray<Real,MAX_SPECIES+1> fweights;
                GpuArray<Real,MAX_SPECIES+1> wiener;

                GpuArray<Real,MAX_SPECIES> hk;
                GpuArray<Real,MAX_SPECIES> yy;
                GpuArray<Real,MAX_SPECIES> yyp;

                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> DijY_edge;
                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> sqD;

                Real kxp = (kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4) + kappa(i-1,j,k)*prim(i-1,j,k,4)*prim(i-1,j,k,4));

                Real meanT = 0.5*(prim(i,j,k,4)+prim(i-1,j,k,4));

                if ((i == 0) and is_lo_x_dirichlet_mass) {
                    kxp  = 2.0*kappa(i-1,j,k)*prim(i-1,j,k,4)*prim(i-1,j,k,4);
                    meanT = prim(i-1,j,k,4);
                }
                if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                    kxp  = 2.0*kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4);
                    meanT = prim(i,j,k,4);
                }

                // Weights for facial fluxes:
                fweights[0] = sqrt(k_B*kxp*volinv*dtinv); //energy flux
                wiener[0] = fweights[0]*stochfacex(i,j,k,4);
                // heat flux
                xflux(i,j,k,nvars) = wiener[0];

                // viscous heating
                // diagonal
                xflux(i,j,k,nvars+1) = 0.5*velx(i,j,k)*(tauxx_stoch(i-1,j,k)+tauxx_stoch(i,j,k));
                // shear
                Real visc_shear_heat = 0.0;
                if ((i == 0) and is_lo_x_dirichlet_mass) {
                    visc_shear_heat += 0.5*(vely(i-1,j+1,k)*tauxy_stoch(i,j+1,k)
                                          + vely(i-1,j,k)*tauxy_stoch(i,j,k));
                    visc_shear_heat += 0.5*(velz(i-1,j,k+1)*tauxz_stoch(i,j,k+1)
                                          + velz(i-1,j,k)*tauxz_stoch(i,j,k));
                }
                else if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                    visc_shear_heat += 0.5*(vely(i,j+1,k)*tauxy_stoch(i,j+1,k)
                                          + vely(i,j,k)*tauxy_stoch(i,j,k));
                    visc_shear_heat += 0.5*(velz(i,j,k+1)*tauxz_stoch(i,j,k+1)
                                          + velz(i,j,k)*tauxz_stoch(i,j,k));
                }
                else {
                    visc_shear_heat += 0.25*((vely(i,j+1,k)+vely(i-1,j+1,k))*tauxy_stoch(i,j+1,k)
                                           + (vely(i,j,k)+vely(i-1,j,k))*tauxy_stoch(i,j,k));
                    visc_shear_heat += 0.25*((velz(i,j,k+1)+velz(i-1,j,k+1))*tauxz_stoch(i,j,k+1)
                                           + (velz(i,j,k)+velz(i-1,j,k))*tauxz_stoch(i,j,k));
                }
                xflux(i,j,k,nvars+2) = visc_shear_heat;

                if (algorithm_type == 2) {

                    for (int n=1; n<1+nspecies; ++n) {
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
                            fweights[1+ll] = sqrt(k_B*MWmix*volinv/(Runiv*dt))*sqD[ns*nspecies+ll];
                            wiener[1+ns] = wiener[1+ns] + fweights[1+ll]*stochfacex(i,j,k,5+ll);
                        }
                        xflux(i,j,k,5+ns) = wiener[1+ns];
                    }

                    GetEnthalpies(meanT, hk);

                    Real soret = 0.;

                    for (int ns=0; ns<nspecies; ++ns) {
                        Real soret_s;
                        soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*0.5*(chi(i-1,j,k,ns)+chi(i,j,k,ns)))*wiener[1+ns];
                        if ((i == 0) and is_lo_x_dirichlet_mass) {
                            soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*chi(i-1,j,k,ns))*wiener[1+ns];
                        }
                        if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                            soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*chi(i,j,k,ns))*wiener[1+ns];
                        }
                        soret += soret_s;
                    }
                    xflux(i,j,k,nvars+3) = soret;
                }

            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                GpuArray<Real,MAX_SPECIES+1> fweights;
                GpuArray<Real,MAX_SPECIES+1> wiener;

                GpuArray<Real,MAX_SPECIES> hk;
                GpuArray<Real,MAX_SPECIES> yy;
                GpuArray<Real,MAX_SPECIES> yyp;

                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> DijY_edge;
                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> sqD;

                Real kyp = kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4) + kappa(i,j-1,k)*prim(i,j-1,k,4)*prim(i,j-1,k,4);

                Real meanT = 0.5*(prim(i,j,k,4)+prim(i,j-1,k,4));

                if ((j == 0) and is_lo_y_dirichlet_mass) {
                    kyp  = 2.0*kappa(i,j-1,k)*prim(i,j-1,k,4)*prim(i,j-1,k,4);
                    meanT = prim(i,j-1,k,4);
                }
                if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                    kyp  = 2.0*kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4);
                    meanT = prim(i,j,k,4);
                }

                // viscous heating
                // diagonal
                yflux(i,j,k,nvars+1) = 0.5*vely(i,j,k)*(tauyy_stoch(i,j-1,k)+tauyy_stoch(i,j,k));
                // shear
                Real visc_shear_heat = 0.0;
                if ((j == 0) and is_lo_y_dirichlet_mass) {
                    visc_shear_heat += 0.5*(velx(i+1,j-1,k)*tauxy_stoch(i+1,j,k)
                                           + velx(i,j-1,k)*tauxy_stoch(i,j,k));
                    visc_shear_heat += 0.5*(velz(i,j-1,k+1)*tauyz_stoch(i,j,k+1)
                                          + velz(i,j-1,k)*tauyz_stoch(i,j,k));
                }
                else if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                    visc_shear_heat += 0.5*(velx(i+1,j,k)*tauxy_stoch(i+1,j,k)
                                         +  velx(i,j,k)*tauxy_stoch(i,j,k));
                    visc_shear_heat += 0.5*(velz(i,j,k+1)*tauyz_stoch(i,j,k+1)
                                         +  velz(i,j,k)*tauyz_stoch(i,j,k));
                }
                else {
                    visc_shear_heat += 0.25*((velx(i+1,j,k)+velx(i+1,j-1,k))*tauxy_stoch(i+1,j,k)
                                           + (velx(i,j,k)+velx(i,j-1,k))*tauxy_stoch(i,j,k));
                    visc_shear_heat += 0.25*((velz(i,j,k+1)+velz(i,j-1,k+1))*tauyz_stoch(i,j,k+1)
                                           + (velz(i,j,k)+velz(i,j-1,k))*tauyz_stoch(i,j,k));
                }
                yflux(i,j,k,nvars+2) = visc_shear_heat;

                if (do_1D) { // 1D
                    yflux(i,j,k,nvars) = 0.0;
                    yflux(i,j,k,nvars+3) = 0.0;
                    for (int ns=0; ns<nspecies; ++ns) {
                        yflux(i,j,k,5+ns) = 0.0;
                    }
                }
                else { // works for 2D and 3D

                    // Weights for facial fluxes:
                    fweights[0] = sqrt(k_B*kyp*volinv*dtinv);
                    wiener[0] = fweights[0]*stochfacey(i,j,k,4);
                    // heat flux
                    yflux(i,j,k,nvars) = wiener[0];

                    if (algorithm_type == 2) {

                        for (int n=1; n<1+nspecies; ++n) {
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
                                fweights[1+ll] = sqrt(k_B*MWmix*volinv/(Runiv*dt))*sqD[ns*nspecies+ll];
                                wiener[1+ns] = wiener[1+ns] + fweights[1+ll]*stochfacey(i,j,k,5+ll);
                            }
                            yflux(i,j,k,5+ns) = wiener[1+ns];
                        }

                        GetEnthalpies(meanT, hk);

                        Real soret = 0.;

                        for (int ns=0; ns<nspecies; ++ns) {
                            Real soret_s;
                            soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*0.5*(chi(i,j-1,k,ns)+chi(i,j,k,ns)))*wiener[1+ns];
                            if ((j == 0) and is_lo_y_dirichlet_mass) {
                                soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*chi(i,j-1,k,ns))*wiener[1+ns];
                            }
                            if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                                soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*chi(i,j,k,ns))*wiener[1+ns];
                            }
                            soret += soret_s;
                        }
                        yflux(i,j,k,nvars+3) = soret;
                    }
                }
            },

            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                GpuArray<Real,MAX_SPECIES+1> fweights;
                GpuArray<Real,MAX_SPECIES+1> wiener;

                GpuArray<Real,MAX_SPECIES> hk;
                GpuArray<Real,MAX_SPECIES> yy;
                GpuArray<Real,MAX_SPECIES> yyp;

                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> DijY_edge;
                GpuArray<Real,MAX_SPECIES*MAX_SPECIES> sqD;

                Real kzp = kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4) + kappa(i,j,k-1)*prim(i,j,k-1,4)*prim(i,j,k-1,4);

                Real meanT = 0.5*(prim(i,j,k,4)+prim(i,j,k-1,4));

                if ((k == 0) and is_lo_z_dirichlet_mass) {
                    kzp  = 2.0*kappa(i,j,k-1)*prim(i,j,k-1,4)*prim(i,j,k-1,4);
                    meanT = prim(i,j,k-1,4);
                }
                if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                    kzp  = 2.0*kappa(i,j,k)*prim(i,j,k,4)*prim(i,j,k,4);
                    meanT = prim(i,j,k,4);
                }

                // viscous heating
                // diagonal
                zflux(i,j,k,nvars+1) = 0.5*velz(i,j,k)*(tauzz_stoch(i,j,k-1)+tauzz_stoch(i,j,k));
                // shear
                Real visc_shear_heat = 0.0;
                if ((k == 0) and is_lo_z_dirichlet_mass) {
                    visc_shear_heat += 0.5*(velx(i+1,j,k-1)*tauxz_stoch(i+1,j,k)
                                          + velx(i,j,k-1)*tauxz_stoch(i,j,k));
                    visc_shear_heat += 0.5*(vely(i,j+1,k-1)*tauyz_stoch(i,j+1,k)
                                          + vely(i,j,k-1)*tauyz_stoch(i,j,k));
                }
                else if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                    visc_shear_heat += 0.5*(velx(i+1,j,k)*tauxz_stoch(i+1,j,k)
                                          + velx(i,j,k)*tauxz_stoch(i,j,k));
                    visc_shear_heat += 0.5*(vely(i,j+1,k)*tauyz_stoch(i,j+1,k)
                                          + vely(i,j,k)*tauyz_stoch(i,j,k));
                }
                else {
                    visc_shear_heat += 0.25*((velx(i+1,j,k-1)+velx(i+1,j,k))*tauxz_stoch(i+1,j,k)
                                           + (velx(i,j,k)+velx(i,j,k-1))*tauxz_stoch(i,j,k));
                    visc_shear_heat += 0.25*((vely(i,j+1,k-1)+vely(i,j+1,k))*tauyz_stoch(i,j+1,k)
                                           + (vely(i,j,k)+vely(i,j,k-1))*tauyz_stoch(i,j,k));
                }
                zflux(i,j,k,nvars+2) = visc_shear_heat;

                if ((do_1D) or (do_2D)) { // works for 1D and 2D
                    zflux(i,j,k,nvars) = 0.0;
                    zflux(i,j,k,nvars+3) = 0.0;
                    for (int ns=0; ns<nspecies; ++ns) {
                        zflux(i,j,k,5+ns) = 0.0;
                    }
                }
                else { // 3D

                    // Weights for facial fluxes:
                    fweights[0] = sqrt(k_B*kzp*volinv*dtinv);
                    wiener[0] = fweights[0]*stochfacez(i,j,k,4);
                    // heat flux
                    zflux(i,j,k,nvars) = wiener[0];


                    if (algorithm_type == 2) {

                        for (int n=1; n<1+nspecies; ++n) {
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
                                fweights[1+ll] = sqrt(k_B*MWmix*volinv/(Runiv*dt))*sqD[ns*nspecies+ll];
                                wiener[1+ns] = wiener[1+ns] + fweights[1+ll]*stochfacez(i,j,k,5+ll);
                            }
                            zflux(i,j,k,5+ns) = wiener[1+ns];
                        }

                        GetEnthalpies(meanT, hk);

                        Real soret = 0.;

                        for (int ns=0; ns<nspecies; ++ns) {
                            Real soret_s;
                            soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*0.5*(chi(i,j,k-1,ns)+chi(i,j,k,ns)))*wiener[1+ns];
                            if ((k == 0) and is_lo_z_dirichlet_mass) {
                                soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*chi(i,j,k-1,ns))*wiener[1+ns];
                            }
                            if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                                soret_s = (hk[ns] + Runiv*meanT/molmass[ns]*chi(i,j,k,ns))*wiener[1+ns];
                            }
                            soret += soret_s;
                        }
                        zflux(i,j,k,nvars+3) = soret;

                    }
                }

            }); // end lambda function

            // Loop over edges for momemntum flux calculations [1:3]
            amrex::ParallelFor(bx_xy, bx_xz, bx_yz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                edgey_u(i,j,k) = tauxy_stoch(i,j,k);
                edgex_v(i,j,k) = tauxy_stoch(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                edgez_u(i,j,k) = tauxz_stoch(i,j,k);
                edgex_w(i,j,k) = tauxz_stoch(i,j,k);
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                edgez_v(i,j,k) = tauyz_stoch(i,j,k);
                edgey_w(i,j,k) = tauyz_stoch(i,j,k);
            });

            // Loop over the center cells and compute fluxes (diagonal momentum terms)
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                cenx_u(i,j,k) = tauxx_stoch(i,j,k);
                ceny_v(i,j,k) = tauyy_stoch(i,j,k);
                cenz_w(i,j,k) = tauzz_stoch(i,j,k);
            });
        } // end MFIter


        // Enforce flux boundary conditions
        StochFluxStag(faceflux_in,cenflux_in,edgeflux_x_in,edgeflux_y_in,edgeflux_z_in,geom);
        if (membrane_cell >= 0) {
            StochFluxMem(faceflux_in,edgeflux_x_in,edgeflux_y_in,edgeflux_z_in);
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

        AMREX_D_TERM(Array4<Real const> const& velx = vel_in[0].array(mfi);,
                     Array4<Real const> const& vely = vel_in[1].array(mfi);,
                     Array4<Real const> const& velz = vel_in[2].array(mfi););

        const Array4<const Real> prim = prim_in.array(mfi);

        const Array4<const Real> eta   = eta_in.array(mfi);
        const Array4<const Real> zeta  = zeta_in.array(mfi);
        const Array4<const Real> kappa = kappa_in.array(mfi);
        const Array4<const Real> chi   = chi_in.array(mfi);
        const Array4<const Real> Dij   = D_in.array(mfi);

        const Box& tbx = mfi.nodaltilebox(0);
        const Box& tby = mfi.nodaltilebox(1);
        const Box& tbz = mfi.nodaltilebox(2);

        const Box & bx_xy = mfi.tilebox(nodal_flag_xy);
        #if (AMREX_SPACEDIM == 3)
        const Box & bx_xz = mfi.tilebox(nodal_flag_xz);
        const Box & bx_yz = mfi.tilebox(nodal_flag_yz);
        #endif

        const Box& bx = mfi.growntilebox(1);

        Real half = 0.5;

        // Populate diagonal stress
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            Real u_x, v_y, w_z; // velocity gradients
            if (do_1D) { // 1D
                u_x = (velx(i+1,j,k) - velx(i,j,k))/dx[0];

                Real div = u_x; // divergence
                if (amrex::Math::abs(visc_type) == 3) {
                  tauxx(i,j,k) = 2*eta(i,j,k)*u_x + (zeta(i,j,k) - 2*eta(i,j,k)/3.)*div;
                  tauyy(i,j,k) = 0.0;
                  tauzz(i,j,k) = 0.0;
                }
                else {
                  tauxx(i,j,k) = 2*eta(i,j,k)*u_x + (0.0 - 2*eta(i,j,k)/3.)*div;
                  tauyy(i,j,k) = 0.0;
                  tauzz(i,j,k) = 0.0;
                }
            }
            else if (do_2D) { // 2D
                u_x = (velx(i+1,j,k) - velx(i,j,k))/dx[0];
                v_y = (vely(i,j+1,k) - vely(i,j,k))/dx[1];

                Real div = u_x + v_y; // divergence
                if (amrex::Math::abs(visc_type) == 3) {
                  tauxx(i,j,k) = 2*eta(i,j,k)*u_x + (zeta(i,j,k) - 2*eta(i,j,k)/3.)*div;
                  tauyy(i,j,k) = 2*eta(i,j,k)*v_y + (zeta(i,j,k) - 2*eta(i,j,k)/3.)*div;
                  tauzz(i,j,k) = 0.0;
                }
                else {
                  tauxx(i,j,k) = 2*eta(i,j,k)*u_x + (0.0 - 2*eta(i,j,k)/3.)*div;
                  tauyy(i,j,k) = 2*eta(i,j,k)*v_y + (0.0 - 2*eta(i,j,k)/3.)*div;
                  tauzz(i,j,k) = 0.0;
                }
            }
            else { // 3D
                u_x = (velx(i+1,j,k) - velx(i,j,k))/dx[0];
                v_y = (vely(i,j+1,k) - vely(i,j,k))/dx[1];
                w_z = (velz(i,j,k+1) - velz(i,j,k))/dx[2];

                Real div = u_x + v_y + w_z; // divergence
                if (amrex::Math::abs(visc_type) == 3) {
                  tauxx(i,j,k) = 2*eta(i,j,k)*u_x + (zeta(i,j,k) - 2*eta(i,j,k)/3.)*div;
                  tauyy(i,j,k) = 2*eta(i,j,k)*v_y + (zeta(i,j,k) - 2*eta(i,j,k)/3.)*div;
                  tauzz(i,j,k) = 2*eta(i,j,k)*w_z + (zeta(i,j,k) - 2*eta(i,j,k)/3.)*div;
                }
                else {
                  tauxx(i,j,k) = 2*eta(i,j,k)*u_x + (0.0 - 2*eta(i,j,k)/3.)*div;
                  tauyy(i,j,k) = 2*eta(i,j,k)*v_y + (0.0 - 2*eta(i,j,k)/3.)*div;
                  tauzz(i,j,k) = 2*eta(i,j,k)*w_z + (0.0 - 2*eta(i,j,k)/3.)*div;
                }
            }

        });

        // Populate off-diagonal stress
        amrex::ParallelFor(bx_xy, bx_xz, bx_yz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            if (do_1D) { // 1D
                tauxy(i,j,k) = 0.0;
            }
            else { // works for both 2D and 3D
                Real u_y, v_x, eta_interp; // velocity gradients
                u_y = (velx(i,j,k) - velx(i,j-1,k))/dx[1];
                v_x = (vely(i,j,k) - vely(i-1,j,k))/dx[0];
                eta_interp = 0.25*(eta(i-1,j-1,k)+eta(i-1,j,k)+eta(i,j-1,k)+eta(i,j,k));
                // Pick boundary values for Dirichlet (stored in ghost)
                // For corner cases (xy), x wall takes preference
                if ((j == 0) and is_lo_y_dirichlet_mass) {
                    u_y = (velx(i,j,k) - velx(i,j-1,k))/(0.5*dx[1]);
                    eta_interp = 0.5*(eta(i-1,j-1,k)+eta(i,j-1,k));
                }
                if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                    u_y = (velx(i,j,k) - velx(i,j-1,k))/(0.5*dx[1]);
                    eta_interp = 0.5*(eta(i-1,j,k)+eta(i,j,k));
                }
                if ((i == 0) and is_lo_x_dirichlet_mass) {
                    v_x = (vely(i,j,k) - vely(i-1,j,k))/(0.5*dx[0]);
                    eta_interp = 0.5*(eta(i-1,j-1,k)+eta(i-1,j,k));
                }
                if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                    v_x = (vely(i,j,k) - vely(i-1,j,k))/(0.5*dx[0]);
                    eta_interp = 0.5*(eta(i,j-1,k)+eta(i,j,k));
                }
                tauxy(i,j,k) = eta_interp*(u_y+v_x);
            }
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            if ((do_1D) or (do_2D)) { // works for 1D and 2D
                tauxz(i,j,k) = 0.0;
            }
            else {
                Real u_z, w_x, eta_interp; // velocity gradients
                u_z = (velx(i,j,k) - velx(i,j,k-1))/dx[2];
                w_x = (velz(i,j,k) - velz(i-1,j,k))/dx[0];
                eta_interp = 0.25*(eta(i-1,j,k-1)+eta(i-1,j,k)+eta(i,j,k-1)+eta(i,j,k));
                // Pick boundary values for Dirichlet (stored in ghost)
                // For corner cases (xz), x wall takes preference
                if ((k == 0) and is_lo_z_dirichlet_mass) {
                    u_z = (velx(i,j,k) - velx(i,j,k-1))/(0.5*dx[2]);
                    eta_interp = 0.5*(eta(i-1,j,k-1)+eta(i,j,k-1));
                }
                if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                    u_z = (velx(i,j,k) - velx(i,j,k-1))/(0.5*dx[2]);
                    eta_interp = 0.5*(eta(i-1,j,k)+eta(i,j,k));
                }
                if ((i == 0) and is_lo_x_dirichlet_mass) {
                    w_x = (velz(i,j,k) - velz(i-1,j,k))/(0.5*dx[0]);
                    eta_interp = 0.5*(eta(i-1,j,k-1)+eta(i-1,j,k));
                }
                if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                    w_x = (velz(i,j,k) - velz(i-1,j,k))/(0.5*dx[0]);
                    eta_interp = 0.5*(eta(i,j,k-1)+eta(i,j,k));
                }
                tauxz(i,j,k) = eta_interp*(u_z+w_x);

            }
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            if ((do_1D) or (do_2D)) { // works for 1D and 2D
                tauyz(i,j,k) = 0.0;
            }
            else {
                Real v_z, w_y, eta_interp; // velocity gradients
                v_z = (vely(i,j,k) - vely(i,j,k-1))/dx[2];
                w_y = (velz(i,j,k) - velz(i,j-1,k))/dx[1];
                eta_interp = 0.25*(eta(i,j-1,k-1)+eta(i,j-1,k)+eta(i,j,k-1)+eta(i,j,k));
                // Pick boundary values for Dirichlet (stored in ghost)
                // For corner cases (yz), y wall takes preference
                if ((k == 0) and is_lo_z_dirichlet_mass) {
                    v_z = (vely(i,j,k) - vely(i,j,k-1))/(0.5*dx[2]);
                    eta_interp = 0.5*(eta(i,j-1,k-1)+eta(i,j,k-1));
                }
                if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                    v_z = (vely(i,j,k) - vely(i,j,k-1))/(0.5*dx[2]);
                    eta_interp = 0.5*(eta(i,j-1,k)+eta(i,j,k));
                }
                if ((j == 0) and is_lo_y_dirichlet_mass) {
                    w_y = (velz(i,j,k) - velz(i,j-1,k))/(0.5*dx[1]);
                    eta_interp = 0.5*(eta(i,j-1,k-1)+eta(i,j-1,k));
                }
                if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                    w_y = (velz(i,j,k) - velz(i,j-1,k))/(0.5*dx[1]);
                    eta_interp = 0.5*(eta(i,j,k-1)+eta(i,j,k));
                }
                tauyz(i,j,k) = eta_interp*(v_z+w_y);
            }
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

            Real kxp   = 0.5*(kappa(i-1,j,k)+kappa(i,j,k));
            Real meanT = 0.5*(prim(i-1,j,k,4)+prim(i,j,k,4));
            Real meanP = 0.5*(prim(i-1,j,k,5)+prim(i,j,k,5));
            if ((i == 0) and is_lo_x_dirichlet_mass) {
                kxp   = kappa(i-1,j,k);
                meanT = prim(i-1,j,k,4);
                meanP = prim(i-1,j,k,5);
            }
            if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                kxp   = kappa(i,j,k);
                meanT = prim(i,j,k,4);
                meanP = prim(i,j,k,5);
            }

            // viscous heating (automatically taken care of setting shear stress to zero above for 1D and 2D)
            // diagonal
            xflux(i,j,k,nvars+1) -= 0.5*velx(i,j,k)*(tauxx(i-1,j,k)+tauxx(i,j,k));
            // shear
            Real visc_shear_heat = 0.0;
            if ((i == 0) and is_lo_x_dirichlet_mass) {
                visc_shear_heat -= 0.5*(vely(i-1,j+1,k)*tauxy(i,j+1,k)
                                      + vely(i-1,j,k)*tauxy(i,j,k));
                visc_shear_heat -= 0.5*(velz(i-1,j,k+1)*tauxz(i,j,k+1)
                                      + velz(i-1,j,k)*tauxz(i,j,k));
                // heat flux
                xflux(i,j,k,nvars) -= kxp*(prim(i,j,k,4)-prim(i-1,j,k,4))/(0.5*dx[0]);
            }
            else if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                visc_shear_heat -= 0.5*(vely(i,j+1,k)*tauxy(i,j+1,k)
                                      + vely(i,j,k)*tauxy(i,j,k));
                visc_shear_heat -= 0.5*(velz(i,j,k+1)*tauxz(i,j,k+1)
                                      + velz(i,j,k)*tauxz(i,j,k));
                // heat flux
                xflux(i,j,k,nvars) -= kxp*(prim(i,j,k,4)-prim(i-1,j,k,4))/(0.5*dx[0]);
            }
            else {
                visc_shear_heat -= 0.25*((vely(i,j+1,k)+vely(i-1,j+1,k))*tauxy(i,j+1,k)
                                       + (vely(i,j,k)+vely(i-1,j,k))*tauxy(i,j,k));
                visc_shear_heat -= 0.25*((velz(i,j,k+1)+velz(i-1,j,k+1))*tauxz(i,j,k+1)
                                       + (velz(i,j,k)+velz(i-1,j,k))*tauxz(i,j,k));
                // heat flux
                xflux(i,j,k,nvars) -= kxp*(prim(i,j,k,4)-prim(i-1,j,k,4))/dx[0];
            }
            xflux(i,j,k,nvars+2) += visc_shear_heat;

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
                        Fk[kk] -= Fks;
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
                    Q5 += Q5s;
                }
                // heat conduction already included in flux(5)
                xflux(i,j,k,nvars+3) += Q5;

                for (int ns=0; ns<nspecies; ++ns) {
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

            // viscous heating (automatically taken care of setting shear stress to zero above for 1D and 2D)
            // diagonal
            yflux(i,j,k,nvars+1) -= 0.5*vely(i,j,k)*(tauyy(i,j-1,k)+tauyy(i,j,k));
            // shear
            Real visc_shear_heat = 0.0;
            if ((j == 0) and is_lo_y_dirichlet_mass) {
                visc_shear_heat -= 0.5*(velx(i+1,j-1,k)*tauxy(i+1,j,k)
                                      + velx(i,j-1,k)*tauxy(i,j,k));
                visc_shear_heat -= 0.5*(velz(i,j-1,k+1)*tauyz(i,j,k+1)
                                      + velz(i,j-1,k)*tauyz(i,j,k));
            }
            else if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                visc_shear_heat -= 0.5*(velx(i+1,j,k)*tauxy(i+1,j,k)
                                      + velx(i,j,k)*tauxy(i,j,k));
                visc_shear_heat -= 0.5*(velz(i,j,k+1)*tauyz(i,j,k+1)
                                      + velz(i,j,k)*tauyz(i,j,k));
            }
            else {
                visc_shear_heat -= 0.25*((velx(i+1,j,k)+velx(i+1,j-1,k))*tauxy(i+1,j,k)
                                        + (velx(i,j,k)+velx(i,j-1,k))*tauxy(i,j,k));
                visc_shear_heat -= 0.25*((velz(i,j,k+1)+velz(i,j-1,k+1))*tauyz(i,j,k+1)
                                        + (velz(i,j,k)+velz(i,j-1,k))*tauyz(i,j,k));
            }
            yflux(i,j,k,nvars+2) += visc_shear_heat;

            if (do_1D) { // 1D
                yflux(i,j,k,nvars) -= 0.0;
                yflux(i,j,k,nvars+3) += 0.0;
                for (int ns=0; ns<nspecies; ++ns) {
                    yflux(i,j,k,5+ns) += 0.0;
                }
            }
            else { // works for 2D and 3D
                Real kyp, meanT, meanP;
                if ((j == 0) and is_lo_y_dirichlet_mass) {
                    kyp   = kappa(i,j-1,k);
                    meanT = prim(i,j-1,k,4);
                    meanP = prim(i,j-1,k,5);
                    // heat flux
                    yflux(i,j,k,nvars) -= kyp*(prim(i,j,k,4)-prim(i,j-1,k,4))/(0.5*dx[1]);
                }
                else if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                    kyp   = kappa(i,j,k);
                    meanT = prim(i,j,k,4);
                    meanP = prim(i,j,k,5);
                    // heat flux
                    yflux(i,j,k,nvars) -= kyp*(prim(i,j,k,4)-prim(i,j-1,k,4))/(0.5*dx[1]);
                }
                else {
                    kyp   = 0.5*(kappa(i,j-1,k)+kappa(i,j,k));
                    meanT = 0.5*(prim(i,j-1,k,4)+prim(i,j,k,4));
                    meanP = 0.5*(prim(i,j-1,k,5)+prim(i,j,k,5));
                    // heat flux
                    yflux(i,j,k,nvars) -= kyp*(prim(i,j,k,4)-prim(i,j-1,k,4))/dx[1];
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
                            Fk[kk] -= Fks;
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
                        Q5 += Q5s;
                    }

                    // heat conduction already included in flux(5)

                    yflux(i,j,k,nvars+3) += Q5;

                    for (int ns=0; ns<nspecies; ++ns) {
                        yflux(i,j,k,5+ns) += Fk[ns];
                    }
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

                // viscous heating (automatically taken care of setting shear stress to zero above for 1D and 2D)
                // diagonal
                zflux(i,j,k,nvars+1) -= 0.5*velz(i,j,k)*(tauzz(i,j,k-1)+tauzz(i,j,k));
                // shear
                Real visc_shear_heat = 0.0;
                if ((k == 0) and is_lo_z_dirichlet_mass) {
                    visc_shear_heat -= 0.5*(velx(i+1,j,k-1)*tauxz(i+1,j,k)
                                           + velx(i,j,k-1)*tauxz(i,j,k));
                    visc_shear_heat -= 0.5*(vely(i,j+1,k-1)*tauyz(i,j+1,k)
                                           + vely(i,j,k-1)*tauyz(i,j,k));
                }
                else if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                    visc_shear_heat -= 0.5*(velx(i+1,j,k)*tauxz(i+1,j,k)
                                           + velx(i,j,k)*tauxz(i,j,k));
                    visc_shear_heat -= 0.5*(vely(i,j+1,k)*tauyz(i,j+1,k)
                                           + vely(i,j,k)*tauyz(i,j,k));
                }
                else {
                    visc_shear_heat -= 0.25*((velx(i+1,j,k-1)+velx(i+1,j,k))*tauxz(i+1,j,k)
                                           + (velx(i,j,k)+velx(i,j,k-1))*tauxz(i,j,k));
                    visc_shear_heat -= 0.25*((vely(i,j+1,k-1)+vely(i,j+1,k))*tauyz(i,j+1,k)
                                       + (vely(i,j,k)+vely(i,j,k-1))*tauyz(i,j,k));
                }
                zflux(i,j,k,nvars+2) += visc_shear_heat;

                if ((do_1D) or (do_2D)) { // works for 1D and 2D
                    zflux(i,j,k,nvars) -= 0.0;
                    zflux(i,j,k,nvars+3) += 0.0;
                    for (int ns=0; ns<nspecies; ++ns) {
                        zflux(i,j,k,5+ns) += 0.0;
                    }
                }
                else { // 3D
                    Real kzp, meanT, meanP;
                    if ((k == 0) and is_lo_z_dirichlet_mass) {
                        kzp   = kappa(i,j,k-1);
                        meanT = prim(i,j,k-1,4);
                        meanP = prim(i,j,k-1,5);
                        // heat flux
                        zflux(i,j,k,nvars) -= kzp*(prim(i,j,k,4)-prim(i,j,k-1,4))/(0.5*dx[2]);
                    }
                    else if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                        kzp   = kappa(i,j,k);
                        meanT = prim(i,j,k,4);
                        meanP = prim(i,j,k,5);
                        // heat flux
                        zflux(i,j,k,nvars) -= kzp*(prim(i,j,k,4)-prim(i,j,k-1,4))/(0.5*dx[2]);
                    }
                    else {
                        kzp   = 0.5*(kappa(i,j,k-1)+kappa(i,j,k));
                        meanT = 0.5*(prim(i,j,k-1,4)+prim(i,j,k,4));
                        meanP = 0.5*(prim(i,j,k-1,5)+prim(i,j,k,5));
                        // heat flux
                        zflux(i,j,k,nvars) -= kzp*(prim(i,j,k,4)-prim(i,j,k-1,4))/dx[2];
                    }

                    if (algorithm_type == 2) {

                        // compute dk
                        for (int ns=0; ns<nspecies; ++ns) {
                            Real term1 = (prim(i,j,k,6+nspecies+ns)-prim(i,j,k-1,6+nspecies+ns))/dx[2];
                            meanXk[ns] = 0.5*(prim(i,j,k-1,6+nspecies+ns)+prim(i,j,k,6+nspecies+ns));
                            meanYk[ns] = 0.5*(prim(i,j,k-1,6+ns)+prim(i,j,k,6+ns));
                            Real term2 = (meanXk[ns]-meanYk[ns])*(prim(i,j,k,5)-prim(i,j,k-1,5))/dx[2]/meanP;
                            dk[ns] = term1 + term2;
                            Real ChiX = 0.5*(chi(i,j,k-1,ns)*prim(i,j,k-1,6+nspecies+ns)+chi(i,j,k,ns)*prim(i,j,k,6+nspecies+ns));
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
                                Fk[kk] -= Fks;
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
                            Q5 += Q5s;
                        }

                        // heat conduction already included in flux(5)
                        zflux(i,j,k,nvars+3) += Q5;

                        for (int ns=0; ns<nspecies; ++ns) {
                            zflux(i,j,k,5+ns) += Fk[ns];
                        }
                    }
                }
            } // n_cells_z test
        });

        // Loop over edges for momemntum flux calculations [1:3]
        amrex::ParallelFor(bx_xy, bx_xz, bx_yz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            edgey_u(i,j,k) -= tauxy(i,j,k);
            edgex_v(i,j,k) -= tauxy(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            edgez_u(i,j,k) -= tauxz(i,j,k);
            edgex_w(i,j,k) -= tauxz(i,j,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            edgez_v(i,j,k) -= tauyz(i,j,k);
            edgey_w(i,j,k) -= tauyz(i,j,k);
        });

        // Loop over the center cells and compute fluxes (diagonal momentum terms)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            cenx_u(i,j,k) -= tauxx(i,j,k);
            ceny_v(i,j,k) -= tauyy(i,j,k);
            cenz_w(i,j,k) -= tauzz(i,j,k);
        });
    }

    // Set species flux to zero at the walls
    BCWallReservoirFluxStag(faceflux_in,cenflux_in,geom);

    ////////////////////
    // hyperbolic fluxes
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

        AMREX_D_TERM(Array4<Real const> const& momx = cumom_in[0].array(mfi);,
                     Array4<Real const> const& momy = cumom_in[1].array(mfi);,
                     Array4<Real const> const& momz = cumom_in[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& velx = vel_in[0].array(mfi);,
                     Array4<Real const> const& vely = vel_in[1].array(mfi);,
                     Array4<Real const> const& velz = vel_in[2].array(mfi););

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

        const Box& bx = mfi.growntilebox(1);

        /**********************************************/
        // Advection Models (advection_type)
        // Type -1: turn off advective fluxes
        // Type 1: species flux = j*Y; energy flux = j*(E + p/rho)
        // Type 0: species flux = j*Y; energy flux = j*(sum_{k}h_k*Y_k + 0.5*u.u)
        // Type 2: species flux = v*(\rho Y); energy flux = v*(\rhoE + p)
        // this will work directly for 1D and 2D as all the velocities in the y- and z-directions are always zero

        // 1. Loop over the face cells and compute fluxes of rho, rhoY, rhoE
        amrex::ParallelFor(tbx, tby, tbz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            if ((i == 0) and is_lo_x_dirichlet_mass) {
                xflux(i,j,k,0) += 0.0;
                xflux(i,j,k,4) += 0.0;
                if (algorithm_type == 2) {
                    for (int n=0; n<nspecies; ++n) {
                        xflux(i,j,k,5+n) += 0.0;
                    }
                }
            }
            else if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                xflux(i,j,k,0) += 0.0;
                xflux(i,j,k,4) += 0.0;
                if (algorithm_type == 2) {
                    for (int n=0; n<nspecies; ++n) {
                        xflux(i,j,k,5+n) += 0.0;
                    }
                }
            }
            else {
                if (advection_type >= 0) {
                    xflux(i,j,k,0) += momx(i,j,k);

                    Real meanT, meanRho, meanP, meanE;
                    GpuArray<Real,MAX_SPECIES> Yk;
                    GpuArray<Real,MAX_SPECIES> hk;

                    // temperature, density and pressure at the face
                    meanT   = 0.5*(prim(i-1,j,k,4) + prim(i,j,k,4));
                    meanRho = 0.5*(prim(i-1,j,k,0) + prim(i,j,k,0));
                    meanP   = 0.5*(prim(i-1,j,k,5) + prim(i,j,k,5));

                    // enthalpy and energy at the face
                    GetEnthalpies(meanT, hk);
                    meanE = 0.5*(cons(i-1,j,k,4) + cons(i,j,k,4))/meanRho;

                    // add energy flux
                    if (advection_type == 1) {
                        xflux(i,j,k,4) += momx(i,j,k)*(meanE + (meanP/meanRho));
                    }
                    else if (advection_type == 2) {
                        xflux(i,j,k,4) += 0.5*(cons(i-1,j,k,4)+cons(i,j,k,4))*velx(i,j,k) +
                                          0.5*(prim(i-1,j,k,5)+prim(i,j,k,5))*velx(i,j,k);
                    }

                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            // concentration advection
                            Yk[n] = 0.5*(prim(i-1,j,k,6+n)+prim(i,j,k,6+n));
                            if ((advection_type == 0) or (advection_type == 1)) {
                                xflux(i,j,k,5+n) += Yk[n]*momx(i,j,k);
                            }
                            else if (advection_type == 2) {
                                xflux(i,j,k,5+n) += 0.5*(cons(i-1,j,k,5+n)+cons(i,j,k,5+n))*velx(i,j,k);
                            }

                            // enthalpy advection (advection_type == 0)
                            if (advection_type == 0) {
                                xflux(i,j,k,4) += momx(i,j,k)*Yk[n]*hk[n];
                            }
                        }
                    }

                    if (advection_type == 0) {
                        // Evaluate KE/rho = 1/2(v.v) on neighboring cells of this face
                        Real ke_rho_P = 0.; // i
                        Real ke_rho_M = 0.; // i-1
                        ke_rho_P += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                        ke_rho_P += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                        ke_rho_P += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                        ke_rho_P *= (0.125/cons(i,j,k,0)/cons(i,j,k,0));
                        ke_rho_M += (momx(i,j,k) + momx(i-1,j,k))*(momx(i,j,k) + momx(i-1,j,k));
                        ke_rho_M += (momy(i-1,j+1,k) + momy(i-1,j,k))*(momy(i-1,j+1,k) + momy(i-1,j,k));
                        ke_rho_M += (momz(i-1,j,k+1) + momz(i-1,j,k))*(momz(i-1,j,k+1) + momz(i-1,j,k));
                        ke_rho_M *= (0.125/cons(i-1,j,k,0)/cons(i-1,j,k,0));

                        // add mom*KE/rho to energy flux
                        xflux(i,j,k,4) += momx(i,j,k)*0.5*(ke_rho_P+ke_rho_M);
                    }
                }
            }
            // add the diffusive + stochastic contributions from heat flux, viscous heating and Dufour effects
            xflux(i,j,k,4) += xflux(i,j,k,nvars) + xflux(i,j,k,nvars+1) + xflux(i,j,k,nvars+2) + xflux(i,j,k,nvars+3);
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            if ((j == 0) and is_lo_y_dirichlet_mass) {
                yflux(i,j,k,0) += 0.0;
                yflux(i,j,k,4) += 0.0;
                if (algorithm_type == 2) {
                    for (int n=0; n<nspecies; ++n) {
                        yflux(i,j,k,5+n) += 0.0;
                    }
                }
            }
            else if  ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                yflux(i,j,k,0) += 0.0;
                yflux(i,j,k,4) += 0.0;
                if (algorithm_type == 2) {
                    for (int n=0; n<nspecies; ++n) {
                        yflux(i,j,k,5+n) += 0.0;
                    }
                }
            }
            else {
                if (advection_type >= 0) {
                    yflux(i,j,k,0) += momy(i,j,k);

                    Real meanT, meanRho, meanP, meanE;
                    GpuArray<Real,MAX_SPECIES> Yk;
                    GpuArray<Real,MAX_SPECIES> hk;

                    // temperature, density and pressure at the face
                    meanT   = 0.5*(prim(i,j-1,k,4) + prim(i,j,k,4));
                    meanRho = 0.5*(prim(i,j-1,k,0) + prim(i,j,k,0));
                    meanP   = 0.5*(prim(i,j-1,k,5) + prim(i,j,k,5));

                    // enthalpy and energy at the face
                    GetEnthalpies(meanT, hk);
                    meanE = 0.5*(cons(i,j-1,k,4) + cons(i,j,k,4))/meanRho;

                    // add energy flux
                    if (advection_type == 1) {
                        yflux(i,j,k,4) += momy(i,j,k)*(meanE + (meanP/meanRho));
                    }
                    else if (advection_type == 2) {
                        yflux(i,j,k,4) += 0.5*(cons(i,j-1,k,4)+cons(i,j,k,4))*vely(i,j,k) +
                                          0.5*(prim(i,j-1,k,5)+prim(i,j,k,5))*vely(i,j,k);
                    }

                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            // concentration advection
                            Yk[n] = 0.5*(prim(i,j-1,k,6+n) + prim(i,j,k,6+n));
                            if ((advection_type == 0) or (advection_type == 1)) {
                                yflux(i,j,k,5+n) += Yk[n]*momy(i,j,k);
                            }
                            else if (advection_type == 2) {
                                yflux(i,j,k,5+n) += 0.5*(cons(i,j-1,k,5+n)+cons(i,j,k,5+n))*vely(i,j,k);
                            }

                            // enthalpy advection (advection_type == 0)
                            if (advection_type == 0) {
                                yflux(i,j,k,4) += momy(i,j,k)*Yk[n]*hk[n];
                            }
                        }
                    }

                    if (advection_type == 0) {
                        // Evaluate KE/rho = 1/2(v.v) on neighboring cells of this face
                        Real ke_rho_P = 0.; // i
                        Real ke_rho_M = 0.; // i-1
                        ke_rho_P += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                        ke_rho_P += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                        ke_rho_P += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                        ke_rho_P *= (0.125/cons(i,j,k,0)/cons(i,j,k,0));
                        ke_rho_M += (momx(i+1,j-1,k) + momx(i,j-1,k))*(momx(i+1,j-1,k) + momx(i,j-1,k));
                        ke_rho_M += (momy(i,j,k) + momy(i,j-1,k))*(momy(i,j,k) + momy(i,j-1,k));
                        ke_rho_M += (momz(i,j-1,k+1) + momz(i,j-1,k))*(momz(i,j-1,k+1) + momz(i,j-1,k));
                        ke_rho_M *= (0.125/cons(i,j-1,k,0)/cons(i,j-1,k,0));

                        // add mom*KE/rho to energy flux
                        yflux(i,j,k,4) += momy(i,j,k)*0.5*(ke_rho_P+ke_rho_M);
                    }
                }
            }
            // add the diffusive + stochastic contributions from heat flux, viscous heating and Dufour effects
            yflux(i,j,k,4) += yflux(i,j,k,nvars) + yflux(i,j,k,nvars+1) + yflux(i,j,k,nvars+2) + yflux(i,j,k,nvars+3);
        },

        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            if ((k == 0) and is_lo_z_dirichlet_mass) {
                zflux(i,j,k,0) += 0.0;
                zflux(i,j,k,4) += 0.0;
                if (algorithm_type == 2) {
                    for (int n=0; n<nspecies; ++n) {
                        zflux(i,j,k,5+n) += 0.0;
                    }
                }
            }
            else if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                zflux(i,j,k,0) += 0.0;
                zflux(i,j,k,4) += 0.0;
                if (algorithm_type == 2) {
                    for (int n=0; n<nspecies; ++n) {
                        zflux(i,j,k,5+n) += 0.0;
                    }
                }
            }
            else {
                if (advection_type >= 0) {
                    zflux(i,j,k,0) += momz(i,j,k);

                    Real meanT, meanRho, meanP, meanE;
                    GpuArray<Real,MAX_SPECIES> Yk;
                    GpuArray<Real,MAX_SPECIES> hk;

                    // temperature, density and pressure at the face
                    meanT   = 0.5*(prim(i,j,k-1,4) + prim(i,j,k,4));
                    meanRho = 0.5*(prim(i,j,k-1,0) + prim(i,j,k,0));
                    meanP   = 0.5*(prim(i,j,k-1,5) + prim(i,j,k,5));

                    // enthalpy and energy at the face
                    GetEnthalpies(meanT, hk);
                    meanE = 0.5*(cons(i,j,k-1,4) + cons(i,j,k,4))/meanRho;

                    // add energy flux
                    if (advection_type == 1) {
                        zflux(i,j,k,4) += momz(i,j,k)*(meanE + (meanP/meanRho));
                    }
                    else if (advection_type == 2) {
                        zflux(i,j,k,4) += 0.5*(cons(i,j,k-1,4)+cons(i,j,k,4))*velz(i,j,k) +
                                          0.5*(prim(i,j,k-1,5)+prim(i,j,k,5))*velz(i,j,k);
                    }

                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            // concentration advection
                            Yk[n] = 0.5*(prim(i,j,k-1,6+n) + prim(i,j,k,6+n));
                            if ((advection_type == 0) or (advection_type == 1)) {
                                zflux(i,j,k,5+n) += Yk[n]*momz(i,j,k);
                            }
                            else if (advection_type == 2) {
                                zflux(i,j,k,5+n) += 0.5*(cons(i,j,k-1,5+n)+cons(i,j,k,5+n))*velz(i,j,k);
                            }

                            // enthalpy advection (advection_type == 0)
                            if (advection_type == 0) {
                                zflux(i,j,k,4) += momz(i,j,k)*Yk[n]*hk[n];
                            }
                        }
                    }

                    if (advection_type == 0) {
                        // Evaluate KE/rho = 1/2(v.v) on neighboring cells of this face
                        Real ke_rho_P = 0.; // i
                        Real ke_rho_M = 0.; // i-1
                        ke_rho_P += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                        ke_rho_P += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                        ke_rho_P += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                        ke_rho_P *= (0.125/cons(i,j,k,0)/cons(i,j,k,0));
                        ke_rho_M += (momx(i+1,j,k-1) + momx(i,j,k-1))*(momx(i+1,j,k-1) + momx(i,j,k-1));
                        ke_rho_M += (momy(i,j+1,k-1) + momy(i,j,k-1))*(momy(i,j+1,k-1) + momy(i,j,k-1));
                        ke_rho_M += (momz(i,j,k) + momz(i,j,k-1))*(momz(i,j,k) + momz(i,j,k-1));
                        ke_rho_M *= (0.125/cons(i,j-1,k,0)/cons(i,j-1,k,0));

                        // add mom*KE/rho to energy flux
                        zflux(i,j,k,4) += momz(i,j,k)*0.5*(ke_rho_P+ke_rho_M);
                    }
                }
            }
            // add the diffusive + stochastic contributions from heat flux, viscous heating and Dufour effects
            zflux(i,j,k,4) += zflux(i,j,k,nvars) + zflux(i,j,k,nvars+1) + zflux(i,j,k,nvars+2) + zflux(i,j,k,nvars+3);
        });

        // 2. Loop over the edge cells and compute fluxes (off-diagonal momentum terms)
        amrex::ParallelFor(bx_xy, bx_xz, bx_yz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            if (advection_type >= 0) {
                // Pick boundary values for Dirichlet (stored in ghost)
                // For corner cases (xy), x wall takes preference
                Real y_u = 0.25*(momx(i,j-1,k)+momx(i,j,k))*(vely(i-1,j,k)+vely(i,j,k));
                Real x_v = 0.25*(momy(i-1,j,k)+momy(i,j,k))*(velx(i,j-1,k)+velx(i,j,k));
                if ((j == 0) and is_lo_y_dirichlet_mass) {
                    y_u = 0.5*(momx(i,j-1,k))*(vely(i-1,j,k)+vely(i,j,k));
                    x_v = 0.5*(momy(i-1,j,k)+momy(i,j,k))*(velx(i,j-1,k));
                }
                if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                    y_u = 0.5*(momx(i,j,k))*(vely(i-1,j,k)+vely(i,j,k));
                    x_v = 0.5*(momy(i-1,j,k)+momy(i,j,k))*(velx(i,j,k));
                }
                if ((i == 0) and is_lo_x_dirichlet_mass) {
                    y_u = 0.5*(momx(i,j-1,k)+momx(i,j,k))*(vely(i-1,j,k));
                    x_v = 0.5*(momy(i-1,j,k))*(velx(i,j-1,k)+velx(i,j,k));
                }
                if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                    y_u = 0.5*(momx(i,j-1,k)+momx(i,j,k))*(vely(i,j,k));
                    x_v = 0.5*(momy(i,j,k))*(velx(i,j-1,k)+velx(i,j,k));
                }
                edgey_u(i,j,k) += y_u;
                edgex_v(i,j,k) += x_v;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            if (advection_type >= 0) {
                // Pick boundary values for Dirichlet (stored in ghost)
                // For corner cases (xz), x wall takes preference
                Real z_u = 0.25*(momx(i,j,k-1)+momx(i,j,k))*(velz(i-1,j,k)+velz(i,j,k));
                Real x_w = 0.25*(momz(i-1,j,k)+momz(i,j,k))*(velx(i,j,k-1)+velx(i,j,k));
                if ((k == 0) and is_lo_z_dirichlet_mass) {
                    z_u = 0.5*(momx(i,j,k-1))*(velz(i-1,j,k)+velz(i,j,k));
                    x_w = 0.5*(momz(i-1,j,k)+momz(i,j,k))*(velx(i,j,k-1));
                }
                if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                    z_u = 0.5*(momx(i,j,k))*(velz(i-1,j,k)+velz(i,j,k));
                    x_w = 0.5*(momz(i-1,j,k)+momz(i,j,k))*(velx(i,j,k));
                }
                if ((i == 0) and is_lo_x_dirichlet_mass) {
                    z_u = 0.5*(momx(i,j,k-1)+momx(i,j,k))*(velz(i-1,j,k));
                    x_w = 0.5*(momz(i-1,j,k))*(velx(i,j,k-1)+velx(i,j,k));
                }
                if ((i == n_cells[0]) and is_hi_x_dirichlet_mass) {
                    z_u = 0.5*(momx(i,j,k-1)+momx(i,j,k))*(velz(i,j,k));
                    x_w = 0.5*(momz(i,j,k))*(velx(i,j,k-1)+velx(i,j,k));
                }
                edgez_u(i,j,k) += z_u;
                edgex_w(i,j,k) += x_w;
            }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            if (advection_type >= 0) {
                // Pick boundary values for Dirichlet (stored in ghost)
                // For corner cases (yz), y wall takes preference
                Real z_v = 0.25*(momy(i,j,k-1)+momy(i,j,k))*(velz(i,j-1,k)+velz(i,j,k));
                Real y_w = 0.25*(momz(i,j-1,k)+momz(i,j,k))*(vely(i,j,k-1)+vely(i,j,k));
                if ((k == 0) and is_lo_z_dirichlet_mass) {
                    z_v = 0.5*(momy(i,j,k-1))*(velz(i,j-1,k)+velz(i,j,k));
                    y_w = 0.5*(momz(i,j-1,k)+momz(i,j,k))*(vely(i,j,k-1));
                }
                if ((k == n_cells[2]) and is_hi_z_dirichlet_mass) {
                    z_v = 0.5*(momy(i,j,k))*(velz(i,j-1,k)+velz(i,j,k));
                    y_w = 0.5*(momz(i,j-1,k)+momz(i,j,k))*(vely(i,j,k));
                }
                if ((j == 0) and is_lo_y_dirichlet_mass) {
                    z_v = 0.5*(momy(i,j,k-1)+momy(i,j,k))*(velz(i,j-1,k));
                    y_w = 0.5*(momz(i,j-1,k))*(vely(i,j,k-1)+vely(i,j,k));
                }
                if ((j == n_cells[1]) and is_hi_y_dirichlet_mass) {
                    z_v = 0.5*(momy(i,j,k-1)+momy(i,j,k))*(velz(i,j,k));
                    y_w = 0.5*(momz(i,j,k))*(vely(i,j,k-1)+vely(i,j,k));
                }
                edgez_v(i,j,k) += z_v;
                edgey_w(i,j,k) += y_w;
            }
        });

        // 3. Loop over the center cells and compute fluxes (diagonal momentum terms)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            if (advection_type >= 0) {
                if (do_1D) { // 1D
                    cenx_u(i,j,k) += 0.25*(momx(i,j,k)+momx(i+1,j,k))*(velx(i,j,k)+velx(i+1,j,k)) + prim(i,j,k,5);
                    ceny_v(i,j,k) += 0.0;
                    cenz_w(i,j,k) += 0.0;
                }
                else if (do_2D) { // 2D
                    cenx_u(i,j,k) += 0.25*(momx(i,j,k)+momx(i+1,j,k))*(velx(i,j,k)+velx(i+1,j,k)) + prim(i,j,k,5);
                    ceny_v(i,j,k) += 0.25*(momy(i,j,k)+momy(i,j+1,k))*(vely(i,j,k)+vely(i,j+1,k)) + prim(i,j,k,5);
                    cenz_w(i,j,k) += 0.0;
                }
                else { // 3D
                    cenx_u(i,j,k) += 0.25*(momx(i,j,k)+momx(i+1,j,k))*(velx(i,j,k)+velx(i+1,j,k)) + prim(i,j,k,5);
                    ceny_v(i,j,k) += 0.25*(momy(i,j,k)+momy(i,j+1,k))*(vely(i,j,k)+vely(i,j+1,k)) + prim(i,j,k,5);
                    cenz_w(i,j,k) += 0.25*(momz(i,j,k)+momz(i,j,k+1))*(velz(i,j,k)+velz(i,j,k+1)) + prim(i,j,k,5);
                }
            }
        });
    }
}
