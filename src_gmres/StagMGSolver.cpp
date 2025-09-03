#include "gmres_functions.H"
#include "common_functions.H"

using namespace amrex;

StagMGSolver::StagMGSolver() {}

void StagMGSolver::Define(const BoxArray& ba_in,
                          const DistributionMapping& dmap_in,
                          const Geometry& geom_in) {

    BL_PROFILE_VAR("StagMGSolver::Define()",StagMGSolver_Define);

    // get the problem domain and boxarray at level 0
    pd_base = geom_in.Domain();
    ba_base = ba_in;

    dmap = dmap_in;

    // compute the number of multigrid levels assuming stag_mg_minwidth is the length of the
    // smallest dimension of the smallest grid at the coarsest multigrid level
    nlevs_mg = ComputeNlevsMG(ba_base);
    if (stag_mg_verbosity >= 3) {
        Print() << "Total number of multigrid levels: " << nlevs_mg << std::endl;
    }

    alpha_fc_mg.resize(nlevs_mg);
    rhs_fc_mg.resize(nlevs_mg);
    phi_fc_mg.resize(nlevs_mg);
    Lphi_fc_mg.resize(nlevs_mg);
    resid_fc_mg.resize(nlevs_mg);
    beta_ed_mg.resize(nlevs_mg);

    beta_cc_mg.resize(nlevs_mg);
    gamma_cc_mg.resize(nlevs_mg);

    dx_mg.resize(nlevs_mg);

    geom_mg.resize(nlevs_mg);

    const Real* dx = geom_in.CellSize();

    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    Vector<int> is_periodic(AMREX_SPACEDIM);
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        is_periodic[i] = geom_in.isPeriodic(i);
    }

    for (int n=0; n<nlevs_mg; ++n) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            // compute dx at this level of multigrid
            dx_mg[n][d] = dx[d] * pow(2,n);
        }

        // create the problem domain for this multigrid level
        Box pd(pd_base);
        pd.coarsen(pow(2,n));

        geom_mg[n].define(pd,&real_box,CoordSys::cartesian,is_periodic.data());

        // create the boxarray for this multigrid level
        BoxArray ba(ba_base);
        ba.coarsen(pow(2,n));

        if ( n == 0 && !(ba == ba_base) ) {
            Abort("Finest multigrid level boxarray and coarsest problem boxarrays do not match");
        }

        // build multifabs used in multigrid coarsening
        beta_cc_mg[n].define(ba, dmap, 1, 1);
        gamma_cc_mg[n].define(ba, dmap, 1, 1);

        for (int d=0; d<AMREX_SPACEDIM; d++) {
            alpha_fc_mg[n][d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);
            rhs_fc_mg[n][d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
            phi_fc_mg[n][d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
            Lphi_fc_mg[n][d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
            resid_fc_mg[n][d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);

            // Put in to fix FPE traps
            Lphi_fc_mg[n][d].setVal(0);
            rhs_fc_mg[n][d].setVal(0);
            phi_fc_mg[n][d].setVal(0);
        }

        // build beta_ed_mg
        if (AMREX_SPACEDIM == 2) {
            beta_ed_mg[n][0].define(convert(ba, nodal_flag), dmap, 1, 0);
        }
        else if (AMREX_SPACEDIM == 3) {
            for (int d=0; d<AMREX_SPACEDIM; d++)
                beta_ed_mg[n][d].define(convert(ba, nodal_flag_edge[d]), dmap, 1, 0);
        }
    } // end loop over multigrid levels

}


// solve "(theta*alpha*I - L) phi = rhs" using multigrid with Gauss-Seidel relaxation
// if amrex::Math::abs(visc_type) = 1, L = div beta grad
// if amrex::Math::abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
// if amrex::Math::abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
// if visc_type > 1 we assume constant coefficients
// if visc_type < 1 we assume variable coefficients
// beta_cc, and gamma_cc are cell-centered
// alpha_fc, phi_fc, and rhs_fc are face-centered
// beta_ed is nodal (2d) or edge-centered (3d)
// phi_fc must come in initialized to some value, preferably a reasonable guess
void StagMGSolver::Solve(const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                         const MultiFab & beta_cc,
                         const std::array<MultiFab, NUM_EDGE> & beta_ed,
                         const MultiFab & gamma_cc,
                         std::array<MultiFab, AMREX_SPACEDIM> & phi_fc,
                         const std::array<MultiFab, AMREX_SPACEDIM> & rhs_fc,
                         const Real & theta_alpha)
{
    BL_PROFILE_VAR("StagMGSolver::Solve()",StagMGSolver_Solve);

    if (stag_mg_verbosity >= 1) {
        Print() << "Begin call to stag_mg_solver\n";
    }

    // initial and current residuals
    Vector<Real> resid0(AMREX_SPACEDIM);
    Vector<Real> resid0_l2(AMREX_SPACEDIM);
    Vector<Real> resid(AMREX_SPACEDIM);
    Vector<Real> resid_l2(AMREX_SPACEDIM);
    Real resid_temp;

    int n, color_start, color_end;

    // copy level 1 coefficients into mg array of coefficients
    MultiFab::Copy(beta_cc_mg[0],  beta_cc,  0, 0, 1, 1);
    MultiFab::Copy(gamma_cc_mg[0], gamma_cc, 0, 0, 1, 1);

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(alpha_fc_mg[0][d], alpha_fc[d], 0, 0, 1, 0);
        // multiply alpha_fc_mg by theta_alpha
        alpha_fc_mg[0][d].mult(theta_alpha,0,1,0);
    }

    MultiFab::Copy(beta_ed_mg[0][0], beta_ed[0], 0, 0, 1, 0);
    if (AMREX_SPACEDIM == 3) {
        MultiFab::Copy(beta_ed_mg[0][1], beta_ed[1], 0, 0, 1, 0);
        MultiFab::Copy(beta_ed_mg[0][2], beta_ed[2], 0, 0, 1, 0);
    }

    // coarsen coefficients
    for (n=1; n<nlevs_mg; ++n) {
        // need ghost cells set to zero to prevent intermediate NaN states
        // that cause some compilers to fail
        beta_cc_mg[n].setVal(0.);
        gamma_cc_mg[n].setVal(0.);

        // cc_restriction on beta_cc_mg and gamma_cc_mg
        // NOTE: CCRestriction calls FillBoundary

        CCRestriction(beta_cc_mg[n],  beta_cc_mg[n-1], geom_mg[n]);
        CCRestriction(gamma_cc_mg[n], gamma_cc_mg[n-1], geom_mg[n]);

        // stag_restriction on alpha_fc_mg
        StagRestriction(alpha_fc_mg[n], alpha_fc_mg[n-1], 1);

        // NOTE: StagRestriction, NodalRestriction, and EdgeRestriction do not
        // call FillBoundary => Do them here for now

        for (int d=0; d<AMREX_SPACEDIM; d++) {
            alpha_fc_mg[n][d].FillBoundary(geom_mg[n].periodicity());
        }

#if (AMREX_SPACEDIM == 2)
        // nodal_restriction on beta_ed_mg
        NodalRestriction(beta_ed_mg[n][0],beta_ed_mg[n-1][0]);
#elif (AMREX_SPACEDIM == 3)
        // edge_restriction on beta_ed_mg
        EdgeRestriction(beta_ed_mg[n],beta_ed_mg[n-1]);
#endif
    }

    /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Now we solve the homogeneous problem
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        // initialize phi_fc_mg = phi_fc as an initial guess
        MultiFab::Copy(phi_fc_mg[0][d],phi_fc[d],0,0,1,0);

        // set values on physical boundaries
        MultiFabPhysBCDomainVel(phi_fc_mg[0][d], geom_mg[0],d);

        // fill periodic ghost cells
        phi_fc_mg[0][d].FillBoundary(geom_mg[0].periodicity());

        // fill physical ghost cells
        MultiFabPhysBCMacVel(phi_fc_mg[0][d], geom_mg[0], d);

        // set rhs_fc_mg at level 1 by copying in passed-in rhs_fc
        MultiFab::Copy(rhs_fc_mg[0][d], rhs_fc[d], 0, 0, 1, 0);
    }

    // compute norm of initial residual
    // first compute viscous part of Lphi
    StagApplyOp(geom_mg[0],beta_cc_mg[0],gamma_cc_mg[0],beta_ed_mg[0],
                phi_fc_mg[0],Lphi_fc_mg[0],alpha_fc_mg[0],dx_mg[0].data(),1.);

    // now subtract the rest of the RHS from Lphi.
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        // compute Lphi - rhs
        MultiFab::Subtract(Lphi_fc_mg[0][d],rhs_fc_mg[0][d],0,0,1,1);

        // compute L0 norm of Lphi - rhs
        resid0[d] = Lphi_fc_mg[0][d].norm0();
// FIXME - need to write an L2 norm for staggered fields
//        resid0_l2[d] = Lphi_fc_mg[0][d].norm2();
        if (stag_mg_verbosity >= 2) {
            Print() << "Initial residual " << d << " " << resid0[d] << std::endl;
        }
    }

    if ( std::all_of(resid0.begin(), resid0.end(), [](Real x){return x==0.;}) ) {
        if (stag_mg_verbosity >= 1) {
            Print() << "Initial residual is zero; exiting staggered multigrid solver" << std::endl;
        }
        return;
    }

    // if some (but not all) of the residuals are zero
    // set the zero residuals to the maximum so the multigrid will begin work
    if ( std::any_of(resid0.begin(), resid0.end(), [](Real x){return x==0.;}) ) {
        std::fill(resid0.begin(),    resid0.end(),    *max_element(resid0.begin(),    resid0.end()));
// FIXME - need to write an L2 norm for staggered fields
//        std::fill(resid0_l2.begin(), resid0_l2.end(), *max_element(resid0_l2.begin(), resid0_l2.end()));
    }

    if (stag_mg_smoother == 0) {
        color_start = 0;
        color_end = 0;
    }
    else {
        color_start = 1;
        color_end = 2*AMREX_SPACEDIM;
    }

    for (int vcycle=1; vcycle<=stag_mg_max_vcycles; ++vcycle) {

        if (stag_mg_verbosity >= 2) {
            Print() << "Begin V-Cycle " << vcycle << std::endl;
        }

        // set phi to zero at coarser levels as initial guess for residual equation
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            for (n=1; n<nlevs_mg; ++n) {
                phi_fc_mg[n][d].setVal(0.);
            }
        }

        // down the V-cycle
        for (n=0; n<=nlevs_mg-2; ++n) {

            // print out residual
            if (stag_mg_verbosity >= 3) {

                // compute Lphi
                StagApplyOp(geom_mg[n],beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                            phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),1.);

                // now subtract the rest of the RHS from Lphi.
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    // compute Lphi - rhs, and report residual
                    MultiFab::Subtract(Lphi_fc_mg[n][d],rhs_fc_mg[n][d],0,0,1,0);
                    resid_temp = Lphi_fc_mg[n][d].norm0();
                    Print() << "Residual for comp " << d << " before    smooths at level "
                            << n << " " << resid_temp << std::endl;
                }
            }

            for (int m=1; m<=stag_mg_nsmooths_down; ++m) {

                // do the smooths
                for (int color=color_start; color<=color_end; ++color) {

                    // the form of weighted Jacobi we are using is
                    // phi^{k+1} = phi^k + omega*D^{-1}*(rhs-Lphi)
                    // where D is the diagonal matrix containing the diagonal elements of L

                    // compute Lphi
                    StagApplyOp(geom_mg[n],beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                                phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),1.,color);

                    // update phi = phi + omega*D^{-1}*(rhs-Lphi)
                    StagMGUpdate(phi_fc_mg[n],rhs_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],
                                 beta_cc_mg[n],beta_ed_mg[n],gamma_cc_mg[n],dx_mg[n].data(),color);

                    for (int d=0; d<AMREX_SPACEDIM; ++d) {

                        // set values on physical boundaries
                        MultiFabPhysBCDomainVel(phi_fc_mg[n][d], geom_mg[n],d);

                        // fill periodic ghost cells
                        phi_fc_mg[n][d].FillBoundary(geom_mg[n].periodicity());

                        // fill physical ghost cells
                        MultiFabPhysBCMacVel(phi_fc_mg[n][d], geom_mg[n],d);
                    }

                } // end loop over colors

                // print out residual
                if (stag_mg_verbosity >= 4) {

                    // compute Lphi
                    StagApplyOp(geom_mg[n],beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                                phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),1.);

                    // now subtract the rest of the RHS from Lphi.
                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        // compute Lphi - rhs, and report residual
                        MultiFab::Subtract(Lphi_fc_mg[n][d],rhs_fc_mg[n][d],0,0,1,0);
                        resid_temp = Lphi_fc_mg[n][d].norm0();
                        Print() << "Residual for comp " << d << " after    smooth " << m << " at level "
                                << n << " " << resid_temp << std::endl;
                    }
                }

            } // end loop over nsmooths

            /////////////////
            // compute residual

            // compute Lphi
            StagApplyOp(geom_mg[n],beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                        phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),1.);

            // now subtract the rest of the RHS from Lphi.
            for (int d=0; d<AMREX_SPACEDIM; ++d) {

                // compute Lphi - rhs, and then multiply by -1
                MultiFab::Subtract(Lphi_fc_mg[n][d],rhs_fc_mg[n][d],0,0,1,0);
                Lphi_fc_mg[n][d].mult(-1.,0,1,0);
                if (stag_mg_verbosity >= 3) {
                    resid_temp = Lphi_fc_mg[n][d].norm0();
                    Print() << "Residual for comp " << d << " after all smooths at level "
                            << n << " " << resid_temp << std::endl;
                }

                // set values on physical boundaries
                MultiFabPhysBCDomainVel(Lphi_fc_mg[n][d], geom_mg[n],d);

                // fill periodic ghost cells
                Lphi_fc_mg[n][d].FillBoundary(geom_mg[n].periodicity());

                // fill physical ghost cells
                MultiFabPhysBCMacVel(Lphi_fc_mg[n][d], geom_mg[n],d);
            }

            // restrict/coarsen residual and put it in rhs_fc
            StagRestriction(rhs_fc_mg[n+1],Lphi_fc_mg[n]);

            for (int d=0; d<AMREX_SPACEDIM; d++) {
                // set residual to zero on physical boundaries
                MultiFabPhysBCDomainVel(rhs_fc_mg[n+1][d], geom_mg[n+1], d);
            }

        }  // end loop over nlevs_mg (bottom of V-cycle)

        // bottom solve
        n = nlevs_mg-1;

        if (stag_mg_verbosity >= 3) {
            Print() << "Begin bottom solve" << std::endl;
        }

        ////////////////////////////
        // just do smooths at the current level as the bottom solve

        // print out residual
        if (stag_mg_verbosity >= 3) {

            // compute Lphi
            StagApplyOp(geom_mg[n],beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                        phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),1.);

            // now subtract the rest of the RHS from Lphi.
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                // compute Lphi - rhs, and report residual

                MultiFab::Subtract(Lphi_fc_mg[n][d],rhs_fc_mg[n][d],0,0,1,0);
                resid_temp = Lphi_fc_mg[n][d].norm0();
                Print() << "Residual for comp " << d << " before    smooths at level "
                        << n << " " << resid_temp << std::endl;
            }
        }

        for (int m=1; m<=stag_mg_nsmooths_bottom; ++m) {

            // do the smooths
            for (int color=color_start; color<=color_end; ++color) {

                // the form of weighted Jacobi we are using is
                // phi^{k+1} = phi^k + omega*D^{-1}*(rhs-Lphi)
                // where D is the diagonal matrix containing the diagonal elements of L

                // compute Lphi
                StagApplyOp(geom_mg[n],beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                            phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),1.,color);

                // update phi = phi + omega*D^{-1}*(rhs-Lphi)
                StagMGUpdate(phi_fc_mg[n],rhs_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],
                             beta_cc_mg[n],beta_ed_mg[n],gamma_cc_mg[n],dx_mg[n].data(),color);

                for (int d=0; d<AMREX_SPACEDIM; d++) {

                    // set values on physical boundaires
                    MultiFabPhysBCDomainVel(phi_fc_mg[n][d], geom_mg[n],d);

                    // fill periodic ghost cells
                    phi_fc_mg[n][d].FillBoundary(geom_mg[n].periodicity());

                    // fill physical ghost cells
                    MultiFabPhysBCMacVel(phi_fc_mg[n][d], geom_mg[n],d);
                }

            } // end loop over colors

        } // end loop over nsmooths

        ////////////////////
        // compute residual

        // compute Lphi
        StagApplyOp(geom_mg[n],beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                    phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),1.);

        // now subtract the rest of the RHS from Lphi.
        for (int d=0; d<AMREX_SPACEDIM; ++d) {

            // compute Lphi - rhs, and then multiply by -1
            MultiFab::Subtract(Lphi_fc_mg[n][d],rhs_fc_mg[n][d],0,0,1,0);
            Lphi_fc_mg[n][d].mult(-1.,0,1,0);
            if (stag_mg_verbosity >= 3) {
                resid_temp = Lphi_fc_mg[n][d].norm0();
                Print() << "Residual for comp " << d << " after all smooths at level "
                        << n << " " << resid_temp << std::endl;
            }

            // set values on physical boundaries
            MultiFabPhysBCDomainVel(Lphi_fc_mg[n][d], geom_mg[n],d);

            // fill periodic ghost cells
            Lphi_fc_mg[n][d].FillBoundary(geom_mg[n].periodicity());

            // fill physical ghost cells
            MultiFabPhysBCMacVel(Lphi_fc_mg[n][d], geom_mg[n],d);
        }

        if (stag_mg_verbosity >= 3) {
            Print() << "End bottom solve" << std::endl;
        }

        // up the V-cycle
        for (n=nlevs_mg-2; n>=0; --n) {

            // prolongate/interpolate correction to update phi
            StagProlongation(phi_fc_mg[n+1],phi_fc_mg[n]);

            for (int d=0; d<AMREX_SPACEDIM; ++d) {

                // set values on physical boundaries
                MultiFabPhysBCDomainVel(phi_fc_mg[n][d], geom_mg[n],d);

                // fill periodic ghost cells
                phi_fc_mg[n][d].FillBoundary(geom_mg[n].periodicity());

                // fill physical ghost cells
                MultiFabPhysBCMacVel(phi_fc_mg[n][d], geom_mg[n],d);
            }

            // print out residual
            if (stag_mg_verbosity >= 3) {

                // compute Lphi
                StagApplyOp(geom_mg[n],beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                            phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),1.);

                // now subtract the rest of the RHS from Lphi.
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    // compute Lphi - rhs, and report residual
                    MultiFab::Subtract(Lphi_fc_mg[n][d],rhs_fc_mg[n][d],0,0,1,0);
                    resid_temp = Lphi_fc_mg[n][d].norm0();
                    Print() << "Residual for comp " << d << " before    smooths at level "
                            << n << " " << resid_temp << std::endl;
                }
            }

            for (int m=1; m<=stag_mg_nsmooths_up; ++m) {

                // do the smooths
                for (int color=color_start; color<=color_end; ++color) {

                    // the form of weighted Jacobi we are using is
                    // phi^{k+1} = phi^k + omega*D^{-1}*(rhs-Lphi)
                    // where D is the diagonal matrix containing the diagonal elements of L

                    // compute Lphi
                    StagApplyOp(geom_mg[n],beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                                phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),1.,color);

                    // update phi = phi + omega*D^{-1}*(rhs-Lphi)
                    StagMGUpdate(phi_fc_mg[n],rhs_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],
                                 beta_cc_mg[n],beta_ed_mg[n],gamma_cc_mg[n],dx_mg[n].data(),color);

                    for (int d=0; d<AMREX_SPACEDIM; d++) {

                        // set values on physical boundaries
                        MultiFabPhysBCDomainVel(phi_fc_mg[n][d], geom_mg[n],d);

                        // fill periodic ghost cells
                        phi_fc_mg[n][d].FillBoundary(geom_mg[n].periodicity());

                        // fill physical ghost cells
                        MultiFabPhysBCMacVel(phi_fc_mg[n][d], geom_mg[n],d);
                    }

                } // end loop over colors

                // print out residual
                if (stag_mg_verbosity >= 4) {

                    // compute Lphi
                    StagApplyOp(geom_mg[n],beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                                phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),1.);

                    // now subtract the rest of the RHS from Lphi.
                    for (int d=0; d<AMREX_SPACEDIM; ++d) {
                        // compute Lphi - rhs, and report residual
                        MultiFab::Subtract(Lphi_fc_mg[n][d],rhs_fc_mg[n][d],0,0,1,0);
                        resid_temp = Lphi_fc_mg[n][d].norm0();
                        Print() << "Residual for comp " << d << " after    smooth " << m << " at level "
                                << n << " " << resid_temp << std::endl;
                    }
                }

            } // end loop over stag_mg_nsmooths_up

            if (stag_mg_verbosity >= 3) {

                // compute Lphi
                StagApplyOp(geom_mg[n],beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                            phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),1.);

                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    // compute Lphi - rhs, and report residual
                    MultiFab::Subtract(Lphi_fc_mg[n][d],rhs_fc_mg[n][d],0,0,1,0);
                    resid_temp = Lphi_fc_mg[n][d].norm0();
                    Print() << "Residual for comp " << d << " after all smooths at level "
                            << n << " " << resid_temp << std::endl;
                }
            }

        } // end loop over nlevs_mg (top of V-cycle)

        if (stag_mg_verbosity >= 2) {
            Print() << "End   V-Cycle " << vcycle << std::endl;
        }

        // compute norm of residual

        // compute Lphi
        StagApplyOp(geom_mg[0],beta_cc_mg[0],gamma_cc_mg[0],beta_ed_mg[0],
                    phi_fc_mg[0],Lphi_fc_mg[0],alpha_fc_mg[0],dx_mg[0].data(),1.);

        // compute Lphi - rhs
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            // compute Lphi - rhs
            MultiFab::Subtract(Lphi_fc_mg[0][d],rhs_fc_mg[0][d],0,0,1,0);
        }

        // compute L0 norm of Lphi - rhs and determine if the problem is solved
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            resid[d] = Lphi_fc_mg[0][d].norm0();
// FIXME - need to write an L2 norm for staggered fields
//            resid_l2[d] = Lphi_fc_mg[0][d].norm2();
            if (stag_mg_verbosity >= 2) {
                Print() << "Residual     " << d << " " << resid[d] << std::endl;
                Print() << "resid/resid0 " << d << " " << resid[d]/resid0[d] << std::endl;
            }
        }
        if (stag_mg_verbosity >= 1) {
// FIXME - need to write an L2 norm for staggered fields
/*
            Real sum1=0., sum2=0.;
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                sum1 += pow(resid_l2[d],2.0);
                sum2 += pow(resid0_l2[d],2.0);
            }
            Print() << "StagMG: L2 |r|/|r0|: " << vcycle
                    << sqrt(sum1)/sqrt(sum2);
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                Print() << " " << resid_l2[d]/resid0_l2[d];
            }
            Print() << std::endl;
*/
        }

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            resid[d] /= resid0[d];
        }

        if ( std::all_of(resid.begin(), resid.end(), [](Real x){return x <= stag_mg_rel_tol;}) ) {
            if (stag_mg_verbosity >= 1) {
                Print() << "Solved in " << vcycle << " staggered V-cycles" << std::endl;
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    Print() << "resid/resid0 " << d << " " << resid[d] << std::endl;
                }
            }
            break;
        }

        if (vcycle == stag_mg_max_vcycles) {
            if (stag_mg_verbosity >= 1) {
                Print() << "Exiting staggered multigrid; maximum number of V-Cycles reached" << std::endl;
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    Print() << "resid/resid0 " << d << " " << resid[d] << std::endl;
                }
            }
        }

    } // end loop over stag_mg_max_vcycles

    //////////////////////////////////
    // Done with multigrid
    //////////////////////////////////

    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        // copy solution back into phi_fc
        MultiFab::Copy(phi_fc[d],phi_fc_mg[0][d],0,0,1,0);

        // set values on physical boundaries
        MultiFabPhysBCDomainVel(phi_fc[d], geom_mg[0],d);

        // fill periodic ghost cells
        phi_fc[d].FillBoundary(geom_mg[0].periodicity());

        // fill physical ghost cells
        MultiFabPhysBCMacVel(phi_fc[d], geom_mg[0],d);
    }

    // vcycle_counter += AMREX_SPACEDIM*stag_mg_max_vcycles;

    if (stag_mg_verbosity >= 1) {
        Print() << "End call to stag_mg_solver\n";
    }
}

// compute the number of multigrid levels assuming minwidth is the length of the
// smallest dimension of the smallest grid at the coarsest multigrid level
int StagMGSolver::ComputeNlevsMG(const BoxArray& ba) {

    BL_PROFILE_VAR("ComputeNlevsMG()",ComputeNlevsMG);

    int nlevs_mg = -1;

    for (int i=0; i<ba.size(); ++i) {
        Box bx = ba.get(i);
        IntVect iv = bx.bigEnd() - bx.smallEnd() + IntVect(1);

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            int temp = iv[d];
            int rdir = 1;
            while (temp%2 == 0 && temp/2 >= stag_mg_minwidth) {
                temp /= 2;
                ++rdir;
            }

            if (nlevs_mg == -1) {
                nlevs_mg = rdir;
            }
            else {
                nlevs_mg = amrex::min(rdir,nlevs_mg);
            }
        }
    }

    return nlevs_mg;
}

void StagMGSolver::CCRestriction(MultiFab& phi_c, const MultiFab& phi_f, const Geometry& geom_c)
{
    BL_PROFILE_VAR("CCRestriction()",CCRestriction);

    // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for ( MFIter mfi(phi_c,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // Get the index space of the valid region
        const Box& bx = mfi.tilebox();

        Array4<Real      > const& phi_c_fab = phi_c.array(mfi);
        Array4<Real const> const& phi_f_fab = phi_f.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
#if (AMREX_SPACEDIM==2)
        {
            phi_c_fab(i,j,k) = 0.25*(  phi_f_fab(2*i,2*j  ,k) + phi_f_fab(2*i+1,2*j  ,k)
                                     + phi_f_fab(2*i,2*j+1,k) + phi_f_fab(2*i+1,2*j+1,k) );
        });
#elif (AMREX_SPACEDIM == 3)
        {
            phi_c_fab(i,j,k) = 0.125*(  phi_f_fab(2*i,2*j  ,2*k  ) + phi_f_fab(2*i+1,2*j  ,2*k  )
                                      + phi_f_fab(2*i,2*j+1,2*k  ) + phi_f_fab(2*i+1,2*j+1,2*k  )
                                      + phi_f_fab(2*i,2*j  ,2*k+1) + phi_f_fab(2*i+1,2*j  ,2*k+1)
                                      + phi_f_fab(2*i,2*j+1,2*k+1) + phi_f_fab(2*i+1,2*j+1,2*k+1) );
        });
#endif
    }

    phi_c.FillBoundary(geom_c.periodicity());
}

void StagMGSolver::StagRestriction(std::array< MultiFab, AMREX_SPACEDIM >& phi_c,
                     const std::array< MultiFab, AMREX_SPACEDIM >& phi_f,
                     int simple_stencil)
{

    BL_PROFILE_VAR("StagRestriction()",StagRestriction);

    // loop over boxes (note we are not passing in a cell-centered MultiFab)
    for ( MFIter mfi(phi_c[0],TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // since the MFIter is built on a nodal MultiFab we need to build the
        // nodal tileboxes for each direction in this way
        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x);,
                     Box bx_y = mfi.tilebox(nodal_flag_y);,
                     Box bx_z = mfi.tilebox(nodal_flag_z););

        const Box& index_bounds = amrex::getIndexBounds(AMREX_D_DECL(bx_x, bx_y, bx_z));

        AMREX_D_TERM(Array4<Real> const& phix_c_fab = phi_c[0].array(mfi);,
                     Array4<Real> const& phiy_c_fab = phi_c[1].array(mfi);,
                     Array4<Real> const& phiz_c_fab = phi_c[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& phix_f_fab = phi_f[0].array(mfi);,
                     Array4<Real const> const& phiy_f_fab = phi_f[1].array(mfi);,
                     Array4<Real const> const& phiz_f_fab = phi_f[2].array(mfi););

        if (simple_stencil == 0) {

#if (AMREX_SPACEDIM == 2)
            amrex::ParallelFor(bx_x, bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                phix_c_fab(i,j,k) = 0.25*( phix_f_fab(2*i  ,2*j,k) + phix_f_fab(2*i  ,2*j+1,k))
                                 + 0.125*( phix_f_fab(2*i+1,2*j,k) + phix_f_fab(2*i+1,2*j+1,k)
                                          +phix_f_fab(2*i-1,2*j,k) + phix_f_fab(2*i-1,2*j+1,k));
            },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                phiy_c_fab(i,j,k) = 0.25*( phiy_f_fab(2*i,2*j  ,k) + phiy_f_fab(2*i+1,2*j  ,k))
                                 + 0.125*( phiy_f_fab(2*i,2*j+1,k) + phiy_f_fab(2*i+1,2*j+1,k)
                                          +phiy_f_fab(2*i,2*j-1,k) + phiy_f_fab(2*i+1,2*j-1,k));
            });
#elif (AMREX_SPACEDIM == 3)
            amrex::ParallelFor(bx_x, bx_y, bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                phix_c_fab(i,j,k) = 0.125* ( phix_f_fab(2*i  ,2*j,2*k  ) + phix_f_fab(2*i  ,2*j+1,2*k  )
                                            +phix_f_fab(2*i  ,2*j,2*k+1) + phix_f_fab(2*i  ,2*j+1,2*k+1) )
                                 + 0.0625* ( phix_f_fab(2*i+1,2*j,2*k  ) + phix_f_fab(2*i+1,2*j+1,2*k  )
                                            +phix_f_fab(2*i+1,2*j,2*k+1) + phix_f_fab(2*i+1,2*j+1,2*k+1) )
                                 + 0.0625* ( phix_f_fab(2*i-1,2*j,2*k  ) + phix_f_fab(2*i-1,2*j+1,2*k  )
                                            +phix_f_fab(2*i-1,2*j,2*k+1) + phix_f_fab(2*i-1,2*j+1,2*k+1) );
            },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                phiy_c_fab(i,j,k) = 0.125* ( phiy_f_fab(2*i,2*j  ,2*k  ) + phiy_f_fab(2*i+1,2*j  ,2*k  )
                                            +phiy_f_fab(2*i,2*j  ,2*k+1) + phiy_f_fab(2*i+1,2*j  ,2*k+1) )
                                 + 0.0625* ( phiy_f_fab(2*i,2*j+1,2*k  ) + phiy_f_fab(2*i+1,2*j+1,2*k  )
                                            +phiy_f_fab(2*i,2*j+1,2*k+1) + phiy_f_fab(2*i+1,2*j+1,2*k+1) )
                                 + 0.0625* ( phiy_f_fab(2*i,2*j-1,2*k  ) + phiy_f_fab(2*i+1,2*j-1,2*k  )
                                            +phiy_f_fab(2*i,2*j-1,2*k+1) + phiy_f_fab(2*i+1,2*j-1,2*k+1) );
            },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                phiz_c_fab(i,j,k) = 0.125* ( phiz_f_fab(2*i,2*j  ,2*k  ) + phiz_f_fab(2*i+1,2*j  ,2*k  )
                                            +phiz_f_fab(2*i,2*j+1,2*k  ) + phiz_f_fab(2*i+1,2*j+1,2*k  ) )
                                 + 0.0625* ( phiz_f_fab(2*i,2*j  ,2*k+1) + phiz_f_fab(2*i+1,2*j  ,2*k+1)
                                            +phiz_f_fab(2*i,2*j+1,2*k+1) + phiz_f_fab(2*i+1,2*j+1,2*k+1) )
                                 + 0.0625* ( phiz_f_fab(2*i,2*j  ,2*k-1) + phiz_f_fab(2*i+1,2*j  ,2*k-1)
                                            +phiz_f_fab(2*i,2*j+1,2*k-1) + phiz_f_fab(2*i+1,2*j+1,2*k-1) );
            });
#endif
        }
        else if (simple_stencil == 1) {

#if (AMREX_SPACEDIM == 2)
            amrex::ParallelFor(bx_x, bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                phix_c_fab(i,j,k) = 0.5*(phix_f_fab(2*i,2*j,k) + phix_f_fab(2*i,2*j+1,k));
            },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                phiy_c_fab(i,j,k) = 0.5*(phiy_f_fab(2*i,2*j,k) + phiy_f_fab(2*i+1,2*j,k));
            });
#elif (AMREX_SPACEDIM == 3)
            amrex::ParallelFor(bx_x, bx_y, bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                phix_c_fab(i,j,k) = 0.25* ( phix_f_fab(2*i,2*j,2*k  ) + phix_f_fab(2*i,2*j+1,2*k  )
                                           +phix_f_fab(2*i,2*j,2*k+1) + phix_f_fab(2*i,2*j+1,2*k+1) );

            },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                phiy_c_fab(i,j,k) = 0.25* ( phiy_f_fab(2*i,2*j,2*k  ) + phiy_f_fab(2*i+1,2*j,2*k  )
                                           +phiy_f_fab(2*i,2*j,2*k+1) + phiy_f_fab(2*i+1,2*j,2*k+1) );

            },
                               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                phiz_c_fab(i,j,k) = 0.25* ( phiz_f_fab(2*i,2*j  ,2*k) + phiz_f_fab(2*i+1,2*j  ,2*k)
                                           +phiz_f_fab(2*i,2*j+1,2*k) + phiz_f_fab(2*i+1,2*j+1,2*k) );

            });
#endif
        }
    }
}

void StagMGSolver::NodalRestriction(MultiFab& phi_c, const MultiFab& phi_f)
{
    BL_PROFILE_VAR("NodalRestriction()",NodalRestriction);

    IntVect nodal(AMREX_D_DECL(1,1,1));

    // loop over boxes (note we are not passing in a cell-centered MultiFab)
    for ( MFIter mfi(phi_c,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // note this is NODAL
        const Box& bx = mfi.tilebox();

        Array4<Real      > const& phi_c_fab = phi_c.array(mfi);
        Array4<Real const> const& phi_f_fab = phi_f.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            phi_c_fab(i,j,k) = phi_f_fab(2*i,2*j,2*k);
        });
    }
}

void StagMGSolver::EdgeRestriction(std::array< MultiFab, NUM_EDGE >& phi_c,
                     const std::array< MultiFab, NUM_EDGE >& phi_f)
{
    BL_PROFILE_VAR("EdgeRestriction()",EdgeRestriction);

    if (AMREX_SPACEDIM != 3) {
        Abort("Edge restriction can only be called for 3D!");
    }

    // loop over boxes (note we are not passing in a cell-centered MultiFab)
    for ( MFIter mfi(phi_c[0],TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // since the MFIter is built on a nodal MultiFab we need to build the
        // nodal tileboxes for each edge in this way
        Box bx_xy = mfi.tilebox(nodal_flag_xy);
        Box bx_xz = mfi.tilebox(nodal_flag_xz);
        Box bx_yz = mfi.tilebox(nodal_flag_yz);

        const Box& index_bounds = amrex::getIndexBounds(bx_xy, bx_xz, bx_yz);

        Array4<Real> const& phixy_c_fab = phi_c[0].array(mfi);
        Array4<Real> const& phixz_c_fab = phi_c[1].array(mfi);
        Array4<Real> const& phiyz_c_fab = phi_c[2].array(mfi);

        Array4<Real const> const& phixy_f_fab = phi_f[0].array(mfi);
        Array4<Real const> const& phixz_f_fab = phi_f[1].array(mfi);
        Array4<Real const> const& phiyz_f_fab = phi_f[2].array(mfi);

        amrex::ParallelFor(bx_xy, bx_xz, bx_yz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            phixy_c_fab(i,j,k) = 0.5*(phixy_f_fab(2*i,2*j,2*k)+phixy_f_fab(2*i,2*j,2*k+1));
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            phixz_c_fab(i,j,k) =  0.5*(phixz_f_fab(2*i,2*j,2*k)+phixz_f_fab(2*i,2*j+1,2*k));
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            phiyz_c_fab(i,j,k) =  0.5*(phiyz_f_fab(2*i,2*j,2*k)+phiyz_f_fab(2*i+1,2*j,2*k));
        });
    }
}

void StagMGSolver::StagProlongation(const std::array< MultiFab, AMREX_SPACEDIM >& phi_c_in,
                                    std::array< MultiFab, AMREX_SPACEDIM >& phi_f_in)
{

    BL_PROFILE_VAR("StagProlongation()",StagProlongation);

    // loop over boxes (note we are not passing in a cell-centered MultiFab)
    for ( MFIter mfi(phi_f_in[0],TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // since the MFIter is built on a nodal MultiFab we need to build the
        // nodal tileboxes for each direction in this way
        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x);,
                     Box bx_y = mfi.tilebox(nodal_flag_y);,
                     Box bx_z = mfi.tilebox(nodal_flag_z););

        AMREX_D_TERM(Array4<Real const> const& phix_c = phi_c_in[0].array(mfi);,
                     Array4<Real const> const& phiy_c = phi_c_in[1].array(mfi);,
                     Array4<Real const> const& phiz_c = phi_c_in[2].array(mfi););

        AMREX_D_TERM(Array4<Real> const& phix_f = phi_f_in[0].array(mfi);,
                     Array4<Real> const& phiy_f = phi_f_in[1].array(mfi);,
                     Array4<Real> const& phiz_f = phi_f_in[2].array(mfi););

#if (AMREX_SPACEDIM == 2)
        amrex::ParallelFor(bx_x, bx_y, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int joff = pow(-1,j%2+1);

            if (i%2 == 0) {
                // linear interpolation
                phix_f(i,j,k) = phix_f(i,j,k)
                    + 0.75*phix_c(i/2,j/2     ,k)
                    + 0.25*phix_c(i/2,j/2+joff,k);
            } else {
                // bilinear interpolation
                phix_f(i,j,k) = phix_f(i,j,k)
                    + 0.375*phix_c(i/2  ,j/2     ,k)
                    + 0.125*phix_c(i/2  ,j/2+joff,k)
                    + 0.375*phix_c(i/2+1,j/2     ,k)
                    + 0.125*phix_c(i/2+1,j/2+joff,k);
            }
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ioff = pow(-1,i%2+1);

            if (j%2 == 0) {
                // linear interpolation
                phiy_f(i,j,k) = phiy_f(i,j,k)
                    + 0.75*phiy_c(i/2     ,j/2,k)
                    + 0.25*phiy_c(i/2+ioff,j/2,k);
            } else {
                // bilinear interpolation
                phiy_f(i,j,k) = phiy_f(i,j,k)
                    + 0.375*phiy_c(i/2     ,j/2  ,k)
                    + 0.125*phiy_c(i/2+ioff,j/2  ,k)
                    + 0.375*phiy_c(i/2     ,j/2+1,k)
                    + 0.125*phiy_c(i/2+ioff,j/2+1,k);
            }
        });

#elif (AMREX_SPACEDIM == 3)

        Real nine16 = 9./16.;
        Real three16 = 3./16.;
        Real one16 = 1./16.;
        Real nine32 = 9./32.;
        Real three32 = 3./32.;
        Real one32 = 1./32.;

        amrex::ParallelFor(bx_x, bx_y, bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            int joff = pow(-1,j%2+1);
            int koff = pow(-1,k%2+1);

            if (i%2 == 0) {
                // bilinear in the yz plane
                phix_f(i,j,k) = phix_f(i,j,k)
                    + nine16 *phix_c(i/2,j/2     ,k/2     )
                    + three16*phix_c(i/2,j/2+joff,k/2     )
                    + three16*phix_c(i/2,j/2     ,k/2+koff)
                    + one16  *phix_c(i/2,j/2+joff,k/2+koff);
            } else {
                // bilinear in the yz plane, linear in x
                phix_f(i,j,k) = phix_f(i,j,k)
                    + nine32 *phix_c(i/2  ,j/2     ,k/2     )
                    + three32*phix_c(i/2  ,j/2+joff,k/2     )
                    + three32*phix_c(i/2  ,j/2     ,k/2+koff)
                    + one32  *phix_c(i/2  ,j/2+joff,k/2+koff)
                    + nine32 *phix_c(i/2+1,j/2     ,k/2     )
                    + three32*phix_c(i/2+1,j/2+joff,k/2     )
                    + three32*phix_c(i/2+1,j/2     ,k/2+koff)
                    + one32  *phix_c(i/2+1,j/2+joff,k/2+koff);
            }
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            int ioff = pow(-1,i%2+1);
            int koff = pow(-1,k%2+1);

            if (j%2 == 0) {
                // bilinear in the xz plane
                phiy_f(i,j,k) = phiy_f(i,j,k)
                    + nine16* phiy_c(i/2     ,j/2,k/2     )
                    + three16*phiy_c(i/2+ioff,j/2,k/2     )
                    + three16*phiy_c(i/2     ,j/2,k/2+koff)
                    + one16*  phiy_c(i/2+ioff,j/2,k/2+koff);
            } else {
                // bilinear in the yz plane, linear in y
                phiy_f(i,j,k) = phiy_f(i,j,k)
                    + nine32* phiy_c(i/2     ,j/2  ,k/2     )
                    + three32*phiy_c(i/2+ioff,j/2  ,k/2     )
                    + three32*phiy_c(i/2     ,j/2  ,k/2+koff)
                    + one32*  phiy_c(i/2+ioff,j/2  ,k/2+koff)
                    + nine32* phiy_c(i/2     ,j/2+1,k/2     )
                    + three32*phiy_c(i/2+ioff,j/2+1,k/2     )
                    + three32*phiy_c(i/2     ,j/2+1,k/2+koff)
                    + one32*  phiy_c(i/2+ioff,j/2+1,k/2+koff);
            }
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            int ioff = pow(-1,i%2+1);
            int joff = pow(-1,j%2+1);

            if (k%2 == 0) {
                // bilinear in the xy plane
                phiz_f(i,j,k) = phiz_f(i,j,k)
                    + nine16* phiz_c(i/2     ,j/2     ,k/2)
                    + three16*phiz_c(i/2+ioff,j/2     ,k/2)
                    + three16*phiz_c(i/2     ,j/2+joff,k/2)
                    + one16*  phiz_c(i/2+ioff,j/2+joff,k/2);
            } else {
                // bilinear in the xy plane, linear in z
                phiz_f(i,j,k) = phiz_f(i,j,k)
                    + nine32* phiz_c(i/2     ,j/2     ,k/2  )
                    + three32*phiz_c(i/2+ioff,j/2     ,k/2  )
                    + three32*phiz_c(i/2     ,j/2+joff,k/2  )
                    + one32*  phiz_c(i/2+ioff,j/2+joff,k/2  )
                    + nine32* phiz_c(i/2     ,j/2     ,k/2+1)
                    + three32*phiz_c(i/2+ioff,j/2     ,k/2+1)
                    + three32*phiz_c(i/2     ,j/2+joff,k/2+1)
                    + one32*  phiz_c(i/2+ioff,j/2+joff,k/2+1);
            }
        });
#endif
    }
}

AMREX_GPU_HOST_DEVICE
inline
void stag_mg_update_visc_p1 (Box const& tbx,
                             AMREX_D_DECL(Box const& xbx,
                                          Box const& ybx,
                                          Box const& zbx),
                             AMREX_D_DECL(Array4<Real> const& phix,
                                          Array4<Real> const& phiy,
                                          Array4<Real> const& phiz),
                             AMREX_D_DECL(Array4<Real const> const& rhsx,
                                          Array4<Real const> const& rhsy,
                                          Array4<Real const> const& rhsz),
                             AMREX_D_DECL(Array4<Real const> const& Lpx,
                                          Array4<Real const> const& Lpy,
                                          Array4<Real const> const& Lpz),
                             AMREX_D_DECL(Array4<Real const> const& alphax,
                                          Array4<Real const> const& alphay,
                                          Array4<Real const> const& alphaz),
                             Array4<Real const> const& beta,
                             Array4<Real const> const& gamma,
                             AMREX_D_DECL(bool do_x,
                                          bool do_y,
                                          bool do_z),
                             int offset,  int color, Real stag_mg_omega,
                             const GpuArray<Real, AMREX_SPACEDIM> & dx) noexcept
{
    // xbx, ybx, and zbx are the face-centered boxes

    // if running on the host
    // tlo is the minimal box containins the union of the face-centered grid boxes

    // if running on the gpu, tlo is a box with a single point that comes
    // from the union of the face-centered grid boxes

    const auto tlo = lbound(tbx);
    const auto thi = ubound(tbx);

    // if running on the host, x/y/zlo and x/y/zhi are set to
    // the lower/uppser bounds of x/y/zbx

    // if running on the gpu, x/y/zlo and x/y/zhi are set to
    // the single point defined by tlo, unless tlo is outside of the union
    // of the face-centered grid boxes, in which case they are set to
    // values that make sure the loop is not entered

    AMREX_D_TERM(const auto xlo = amrex::elemwiseMax(tlo, lbound(xbx));,
                 const auto ylo = amrex::elemwiseMax(tlo, lbound(ybx));,
                 const auto zlo = amrex::elemwiseMax(tlo, lbound(zbx)););

    AMREX_D_TERM(const auto xhi = amrex::elemwiseMin(thi, ubound(xbx));,
                 const auto yhi = amrex::elemwiseMin(thi, ubound(ybx));,
                 const auto zhi = amrex::elemwiseMin(thi, ubound(zbx)););

    int ioff;
    Real fac;
    Real dxsqinv = 1./(dx[0]*dx[0]);

    if (do_x) {

        for (int k = xlo.z; k <= xhi.z; ++k) {
        for (int j = xlo.y; j <= xhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (xlo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = xlo.x+ioff; i <= xhi.x; i+=offset) {
            fac = alphax(i,j,k) + 2.*AMREX_SPACEDIM*beta(i,j,k) * dxsqinv;
            phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac;
        }
        }
        }
    }

    if (do_y) {

        for (int k = ylo.z; k <= yhi.z; ++k) {
        for (int j = ylo.y; j <= yhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (ylo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = ylo.x+ioff; i <= yhi.x; i+=offset) {
            fac = alphay(i,j,k) + 2.*AMREX_SPACEDIM*beta(i,j,k) * dxsqinv;
            phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac;
        }
        }
        }
    }

#if (AMREX_SPACEDIM == 3)
    if (do_z) {

        for (int k = zlo.z; k <= zhi.z; ++k) {
        for (int j = zlo.y; j <= zhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (zlo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = zlo.x+ioff; i <= zhi.x; i+=offset) {
            fac = alphaz(i,j,k) + 2.*AMREX_SPACEDIM*beta(i,j,k) * dxsqinv;
            phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac;
        }
        }
        }
    }
#endif

}

AMREX_GPU_HOST_DEVICE
inline
void stag_mg_update_visc_m1 (Box const& tbx,
                             AMREX_D_DECL(Box const& xbx,
                                          Box const& ybx,
                                          Box const& zbx),
                             AMREX_D_DECL(Array4<Real> const& phix,
                                          Array4<Real> const& phiy,
                                          Array4<Real> const& phiz),
                             AMREX_D_DECL(Array4<Real const> const& rhsx,
                                          Array4<Real const> const& rhsy,
                                          Array4<Real const> const& rhsz),
                             AMREX_D_DECL(Array4<Real const> const& Lpx,
                                          Array4<Real const> const& Lpy,
                                          Array4<Real const> const& Lpz),
                             AMREX_D_DECL(Array4<Real const> const& alphax,
                                          Array4<Real const> const& alphay,
                                          Array4<Real const> const& alphaz),
                             Array4<Real const> const& beta,
                             Array4<Real const> const& beta_xy,
#if (AMREX_SPACEDIM == 3)
                             Array4<Real const> const& beta_xz,
                             Array4<Real const> const& beta_yz,
#endif
                             AMREX_D_DECL(bool do_x,
                                          bool do_y,
                                          bool do_z),
                             int offset,  int color, Real stag_mg_omega,
                             const GpuArray<Real, AMREX_SPACEDIM> & dx) noexcept
{
    // xbx, ybx, and zbx are the face-centered boxes

    // if running on the host
    // tlo is the minimal box containins the union of the face-centered grid boxes

    // if running on the gpu, tlo is a box with a single point that comes
    // from the union of the face-centered grid boxes

    const auto tlo = lbound(tbx);
    const auto thi = ubound(tbx);

    // if running on the host, x/y/zlo and x/y/zhi are set to
    // the lower/uppser bounds of x/y/zbx

    // if running on the gpu, x/y/zlo and x/y/zhi are set to
    // the single point defined by tlo, unless tlo is outside of the union
    // of the face-centered grid boxes, in which case they are set to
    // values that make sure the loop is not entered

    AMREX_D_TERM(const auto xlo = amrex::elemwiseMax(tlo, lbound(xbx));,
                 const auto ylo = amrex::elemwiseMax(tlo, lbound(ybx));,
                 const auto zlo = amrex::elemwiseMax(tlo, lbound(zbx)););

    AMREX_D_TERM(const auto xhi = amrex::elemwiseMin(thi, ubound(xbx));,
                 const auto yhi = amrex::elemwiseMin(thi, ubound(ybx));,
                 const auto zhi = amrex::elemwiseMin(thi, ubound(zbx)););

    int ioff;
    Real fac;
    Real dxsqinv = 1./(dx[0]*dx[0]);

    if (do_x) {

        for (int k = xlo.z; k <= xhi.z; ++k) {
        for (int j = xlo.y; j <= xhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (xlo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = xlo.x+ioff; i <= xhi.x; i+=offset) {

            fac = alphax(i,j,k) +
                ( beta(i,j,k)+beta(i-1,j,k)
                  +beta_xy(i,j,k)+beta_xy(i,j+1,k)
#if (AMREX_SPACEDIM == 3)
                  +beta_xz(i,j,k)+beta_xz(i,j,k+1)
#endif
                    ) * dxsqinv;

                phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac;
        }
        }
        }
    }

    if (do_y) {

        for (int k = ylo.z; k <= yhi.z; ++k) {
        for (int j = ylo.y; j <= yhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (ylo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = ylo.x+ioff; i <= yhi.x; i+=offset) {

            fac = alphay(i,j,k) +
                ( beta(i,j,k)+beta(i,j-1,k)
                  +beta_xy(i,j,k)+beta_xy(i+1,j,k)
#if (AMREX_SPACEDIM == 3)
                  +beta_yz(i,j,k)+beta_yz(i,j,k+1)
#endif
                    ) * dxsqinv;

            phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac;

        }
        }
        }
    }

#if (AMREX_SPACEDIM == 3)
    if (do_z) {

        for (int k = zlo.z; k <= zhi.z; ++k) {
        for (int j = zlo.y; j <= zhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (zlo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = zlo.x+ioff; i <= zhi.x; i+=offset) {

            fac = alphaz(i,j,k) +
                ( beta(i,j,k)+beta(i,j,k-1)
                  +beta_xz(i,j,k)+beta_xz(i+1,j,k)
                  +beta_yz(i,j,k)+beta_yz(i,j+1,k) ) * dxsqinv;

            phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac;

        }
        }
        }
    }
#endif

}

AMREX_GPU_HOST_DEVICE
inline
void stag_mg_update_visc_p2 (Box const& tbx,
                             AMREX_D_DECL(Box const& xbx,
                                          Box const& ybx,
                                          Box const& zbx),
                             AMREX_D_DECL(Array4<Real> const& phix,
                                          Array4<Real> const& phiy,
                                          Array4<Real> const& phiz),
                             AMREX_D_DECL(Array4<Real const> const& rhsx,
                                          Array4<Real const> const& rhsy,
                                          Array4<Real const> const& rhsz),
                             AMREX_D_DECL(Array4<Real const> const& Lpx,
                                          Array4<Real const> const& Lpy,
                                          Array4<Real const> const& Lpz),
                             AMREX_D_DECL(Array4<Real const> const& alphax,
                                          Array4<Real const> const& alphay,
                                          Array4<Real const> const& alphaz),
                             Array4<Real const> const& beta,
                             Array4<Real const> const& gamma,
                             AMREX_D_DECL(bool do_x,
                                          bool do_y,
                                          bool do_z),
                             int offset,  int color, Real stag_mg_omega,
                             const GpuArray<Real, AMREX_SPACEDIM> & dx) noexcept
{
    // xbx, ybx, and zbx are the face-centered boxes

    // if running on the host
    // tlo is the minimal box containins the union of the face-centered grid boxes

    // if running on the gpu, tlo is a box with a single point that comes
    // from the union of the face-centered grid boxes

    const auto tlo = lbound(tbx);
    const auto thi = ubound(tbx);

    // if running on the host, x/y/zlo and x/y/zhi are set to
    // the lower/uppser bounds of x/y/zbx

    // if running on the gpu, x/y/zlo and x/y/zhi are set to
    // the single point defined by tlo, unless tlo is outside of the union
    // of the face-centered grid boxes, in which case they are set to
    // values that make sure the loop is not entered

    AMREX_D_TERM(const auto xlo = amrex::elemwiseMax(tlo, lbound(xbx));,
                 const auto ylo = amrex::elemwiseMax(tlo, lbound(ybx));,
                 const auto zlo = amrex::elemwiseMax(tlo, lbound(zbx)););

    AMREX_D_TERM(const auto xhi = amrex::elemwiseMin(thi, ubound(xbx));,
                 const auto yhi = amrex::elemwiseMin(thi, ubound(ybx));,
                 const auto zhi = amrex::elemwiseMin(thi, ubound(zbx)););

    int ioff;
    Real fac;
    Real dxsqinv = 1./(dx[0]*dx[0]);

    if (do_x) {

        for (int k = xlo.z; k <= xhi.z; ++k) {
        for (int j = xlo.y; j <= xhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (xlo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = xlo.x+ioff; i <= xhi.x; i+=offset) {
            fac = alphax(i,j,k) + 2.*(1.+AMREX_SPACEDIM)*beta(i,j,k) * dxsqinv;
            phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac;
        }
        }
        }
    }

    if (do_y) {

        for (int k = ylo.z; k <= yhi.z; ++k) {
        for (int j = ylo.y; j <= yhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (ylo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = ylo.x+ioff; i <= yhi.x; i+=offset) {
            fac = alphay(i,j,k) + 2.*(1.+AMREX_SPACEDIM)*beta(i,j,k) * dxsqinv;
            phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac  ;
        }
        }
        }
    }

#if (AMREX_SPACEDIM == 3)
    if (do_z) {

        for (int k = zlo.z; k <= zhi.z; ++k) {
        for (int j = zlo.y; j <= zhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (zlo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = zlo.x+ioff; i <= zhi.x; i+=offset) {
            fac = alphaz(i,j,k) + 2.*(1.+AMREX_SPACEDIM)*beta(i,j,k) * dxsqinv;
            phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac;
        }
        }
        }
    }
#endif

}

AMREX_GPU_HOST_DEVICE
inline
void stag_mg_update_visc_m2 (Box const& tbx,
                             AMREX_D_DECL(Box const& xbx,
                                          Box const& ybx,
                                          Box const& zbx),
                             AMREX_D_DECL(Array4<Real> const& phix,
                                          Array4<Real> const& phiy,
                                          Array4<Real> const& phiz),
                             AMREX_D_DECL(Array4<Real const> const& rhsx,
                                          Array4<Real const> const& rhsy,
                                          Array4<Real const> const& rhsz),
                             AMREX_D_DECL(Array4<Real const> const& Lpx,
                                          Array4<Real const> const& Lpy,
                                          Array4<Real const> const& Lpz),
                             AMREX_D_DECL(Array4<Real const> const& alphax,
                                          Array4<Real const> const& alphay,
                                          Array4<Real const> const& alphaz),
                             Array4<Real const> const& beta,
                             Array4<Real const> const& beta_xy,
#if (AMREX_SPACEDIM == 3)
                             Array4<Real const> const& beta_xz,
                             Array4<Real const> const& beta_yz,
#endif
                             AMREX_D_DECL(bool do_x,
                                          bool do_y,
                                          bool do_z),
                             int offset,  int color, Real stag_mg_omega,
                             const GpuArray<Real, AMREX_SPACEDIM> & dx) noexcept
{
    // xbx, ybx, and zbx are the face-centered boxes

    // if running on the host
    // tlo is the minimal box containins the union of the face-centered grid boxes

    // if running on the gpu, tlo is a box with a single point that comes
    // from the union of the face-centered grid boxes

    const auto tlo = lbound(tbx);
    const auto thi = ubound(tbx);

    // if running on the host, x/y/zlo and x/y/zhi are set to
    // the lower/uppser bounds of x/y/zbx

    // if running on the gpu, x/y/zlo and x/y/zhi are set to
    // the single point defined by tlo, unless tlo is outside of the union
    // of the face-centered grid boxes, in which case they are set to
    // values that make sure the loop is not entered

    AMREX_D_TERM(const auto xlo = amrex::elemwiseMax(tlo, lbound(xbx));,
                 const auto ylo = amrex::elemwiseMax(tlo, lbound(ybx));,
                 const auto zlo = amrex::elemwiseMax(tlo, lbound(zbx)););

    AMREX_D_TERM(const auto xhi = amrex::elemwiseMin(thi, ubound(xbx));,
                 const auto yhi = amrex::elemwiseMin(thi, ubound(ybx));,
                 const auto zhi = amrex::elemwiseMin(thi, ubound(zbx)););

    int ioff;
    Real fac;
    Real dxsqinv = 1./(dx[0]*dx[0]);

    if (do_x) {

        for (int k = xlo.z; k <= xhi.z; ++k) {
        for (int j = xlo.y; j <= xhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (xlo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = xlo.x+ioff; i <= xhi.x; i+=offset) {
            fac = alphax(i,j,k) +
                ( 2.*beta(i,j,k)+2.*beta(i-1,j,k)
                  +beta_xy(i,j,k)+beta_xy(i,j+1,k)
#if (AMREX_SPACEDIM == 3)
                  +beta_xz(i,j,k)+beta_xz(i,j,k+1)
#endif
                    ) * dxsqinv;
            phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac;
        }
        }
        }
    }

    if (do_y) {

        for (int k = ylo.z; k <= yhi.z; ++k) {
        for (int j = ylo.y; j <= yhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (ylo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = ylo.x+ioff; i <= yhi.x; i+=offset) {
            fac = alphay(i,j,k) +
                ( 2.*beta(i,j,k)+2.*beta(i,j-1,k)
                  +beta_xy(i,j,k)+beta_xy(i+1,j,k)
#if (AMREX_SPACEDIM == 3)
                  +beta_yz(i,j,k)+beta_yz(i,j,k+1)
#endif
                    ) * dxsqinv;
            phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac;
        }
        }
        }
    }

#if (AMREX_SPACEDIM == 3)
    if (do_z) {

        for (int k = zlo.z; k <= zhi.z; ++k) {
        for (int j = zlo.y; j <= zhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (zlo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = zlo.x+ioff; i <= zhi.x; i+=offset) {
            fac = alphaz(i,j,k) +
                ( 2.*beta(i,j,k)+2.*beta(i,j,k-1)
                  +beta_xz(i,j,k)+beta_xz(i+1,j,k)
                  +beta_yz(i,j,k)+beta_yz(i,j+1,k) ) * dxsqinv;
            phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac;
        }
        }
        }
    }
#endif

}

AMREX_GPU_HOST_DEVICE
inline
void stag_mg_update_visc_p3 (Box const& tbx,
                             AMREX_D_DECL(Box const& xbx,
                                          Box const& ybx,
                                          Box const& zbx),
                             AMREX_D_DECL(Array4<Real> const& phix,
                                          Array4<Real> const& phiy,
                                          Array4<Real> const& phiz),
                             AMREX_D_DECL(Array4<Real const> const& rhsx,
                                          Array4<Real const> const& rhsy,
                                          Array4<Real const> const& rhsz),
                             AMREX_D_DECL(Array4<Real const> const& Lpx,
                                          Array4<Real const> const& Lpy,
                                          Array4<Real const> const& Lpz),
                             AMREX_D_DECL(Array4<Real const> const& alphax,
                                          Array4<Real const> const& alphay,
                                          Array4<Real const> const& alphaz),
                             Array4<Real const> const& beta,
                             Array4<Real const> const& gamma,
                             AMREX_D_DECL(bool do_x,
                                          bool do_y,
                                          bool do_z),
                             int offset,  int color, Real stag_mg_omega,
                             const GpuArray<Real, AMREX_SPACEDIM> & dx) noexcept
{
    // xbx, ybx, and zbx are the face-centered boxes

    // if running on the host
    // tlo is the minimal box containins the union of the face-centered grid boxes

    // if running on the gpu, tlo is a box with a single point that comes
    // from the union of the face-centered grid boxes

    const auto tlo = lbound(tbx);
    const auto thi = ubound(tbx);

    // if running on the host, x/y/zlo and x/y/zhi are set to
    // the lower/uppser bounds of x/y/zbx

    // if running on the gpu, x/y/zlo and x/y/zhi are set to
    // the single point defined by tlo, unless tlo is outside of the union
    // of the face-centered grid boxes, in which case they are set to
    // values that make sure the loop is not entered

    AMREX_D_TERM(const auto xlo = amrex::elemwiseMax(tlo, lbound(xbx));,
                 const auto ylo = amrex::elemwiseMax(tlo, lbound(ybx));,
                 const auto zlo = amrex::elemwiseMax(tlo, lbound(zbx)););

    AMREX_D_TERM(const auto xhi = amrex::elemwiseMin(thi, ubound(xbx));,
                 const auto yhi = amrex::elemwiseMin(thi, ubound(ybx));,
                 const auto zhi = amrex::elemwiseMin(thi, ubound(zbx)););

    int ioff;
    Real fac;
    Real dxsqinv = 1./(dx[0]*dx[0]);
    Real fac2 = (AMREX_SPACEDIM == 2) ? 14./3. : 20./3.;

    if (do_x) {

        for (int k = xlo.z; k <= xhi.z; ++k) {
        for (int j = xlo.y; j <= xhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (xlo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = xlo.x+ioff; i <= xhi.x; i+=offset) {
            fac = alphax(i,j,k)+(fac2*beta(i,j,k)+2.*gamma(i,j,k)) * dxsqinv;
            phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac;
        }
        }
        }
    }

    if (do_y) {

        for (int k = ylo.z; k <= yhi.z; ++k) {
        for (int j = ylo.y; j <= yhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (ylo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = ylo.x+ioff; i <= yhi.x; i+=offset) {
            fac = alphay(i,j,k)+(fac2*beta(i,j,k)+2.*gamma(i,j,k)) * dxsqinv;
            phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac;
        }
        }
        }
    }

#if (AMREX_SPACEDIM == 3)
    if (do_z) {

        for (int k = zlo.z; k <= zhi.z; ++k) {
        for (int j = zlo.y; j <= zhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (zlo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = zlo.x+ioff; i <= zhi.x; i+=offset) {
            fac = alphaz(i,j,k)+(fac2*beta(i,j,k)+2.*gamma(i,j,k)) * dxsqinv;
            phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac;
        }
        }
        }
    }
#endif

}

AMREX_GPU_HOST_DEVICE
inline
void stag_mg_update_visc_m3 (Box const& tbx,
                             AMREX_D_DECL(Box const& xbx,
                                          Box const& ybx,
                                          Box const& zbx),
                             AMREX_D_DECL(Array4<Real> const& phix,
                                          Array4<Real> const& phiy,
                                          Array4<Real> const& phiz),
                             AMREX_D_DECL(Array4<Real const> const& rhsx,
                                          Array4<Real const> const& rhsy,
                                          Array4<Real const> const& rhsz),
                             AMREX_D_DECL(Array4<Real const> const& Lpx,
                                          Array4<Real const> const& Lpy,
                                          Array4<Real const> const& Lpz),
                             AMREX_D_DECL(Array4<Real const> const& alphax,
                                          Array4<Real const> const& alphay,
                                          Array4<Real const> const& alphaz),
                             Array4<Real const> const& beta,
                             Array4<Real const> const& beta_xy,
#if (AMREX_SPACEDIM == 3)
                             Array4<Real const> const& beta_xz,
                             Array4<Real const> const& beta_yz,
#endif
                             Array4<Real const> const& gamma,
                             AMREX_D_DECL(bool do_x,
                                          bool do_y,
                                          bool do_z),
                             int offset,  int color, Real stag_mg_omega,
                             const GpuArray<Real, AMREX_SPACEDIM> & dx) noexcept
{
    // xbx, ybx, and zbx are the face-centered boxes

    // if running on the host
    // tlo is the minimal box containins the union of the face-centered grid boxes

    // if running on the gpu, tlo is a box with a single point that comes
    // from the union of the face-centered grid boxes

    const auto tlo = lbound(tbx);
    const auto thi = ubound(tbx);

    // if running on the host, x/y/zlo and x/y/zhi are set to
    // the lower/uppser bounds of x/y/zbx

    // if running on the gpu, x/y/zlo and x/y/zhi are set to
    // the single point defined by tlo, unless tlo is outside of the union
    // of the face-centered grid boxes, in which case they are set to
    // values that make sure the loop is not entered

    AMREX_D_TERM(const auto xlo = amrex::elemwiseMax(tlo, lbound(xbx));,
                 const auto ylo = amrex::elemwiseMax(tlo, lbound(ybx));,
                 const auto zlo = amrex::elemwiseMax(tlo, lbound(zbx)););

    AMREX_D_TERM(const auto xhi = amrex::elemwiseMin(thi, ubound(xbx));,
                 const auto yhi = amrex::elemwiseMin(thi, ubound(ybx));,
                 const auto zhi = amrex::elemwiseMin(thi, ubound(zbx)););

    int ioff;
    Real fac;
    Real dxsqinv = 1./(dx[0]*dx[0]);
    Real fourthirds = 4./3.;

    if (do_x) {

        for (int k = xlo.z; k <= xhi.z; ++k) {
        for (int j = xlo.y; j <= xhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (xlo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = xlo.x+ioff; i <= xhi.x; i+=offset) {

                   fac = alphax(i,j,k) +
                        ( fourthirds*beta(i,j,k)+gamma(i,j,k)
                        +fourthirds*beta(i-1,j,k)+gamma(i-1,j,k)
                        +beta_xy(i,j,k)+beta_xy(i,j+1,k)
#if (AMREX_SPACEDIM == 3)
                        +beta_xz(i,j,k)+beta_xz(i,j,k+1)
#endif
                            ) * dxsqinv;

                   phix(i,j,k) = phix(i,j,k) + stag_mg_omega*(rhsx(i,j,k)-Lpx(i,j,k)) / fac;
        }
        }
        }
    }

    if (do_y) {

        for (int k = ylo.z; k <= yhi.z; ++k) {
        for (int j = ylo.y; j <= yhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (ylo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = ylo.x+ioff; i <= yhi.x; i+=offset) {

                   fac = alphay(i,j,k) +
                        ( fourthirds*beta(i,j,k)+gamma(i,j,k)
                        +fourthirds*beta(i,j-1,k)+gamma(i,j-1,k)
                        +beta_xy(i,j,k)+beta_xy(i+1,j,k)
#if (AMREX_SPACEDIM == 3)
                        +beta_yz(i,j,k)+beta_yz(i,j,k+1)
#endif
                            ) * dxsqinv;

                   phiy(i,j,k) = phiy(i,j,k) + stag_mg_omega*(rhsy(i,j,k)-Lpy(i,j,k)) / fac;
        }
        }
        }
    }

#if (AMREX_SPACEDIM == 3)
    if (do_z) {

        for (int k = zlo.z; k <= zhi.z; ++k) {
        for (int j = zlo.y; j <= zhi.y; ++j) {
        ioff = 0;
        if (offset == 2 && (zlo.x+j+k)%2 != (color+1)%2 ) {
          ioff = 1;
        }
        AMREX_PRAGMA_SIMD
        for (int i = zlo.x+ioff; i <= zhi.x; i+=offset) {

                   fac = alphaz(i,j,k) +
                        ( fourthirds*beta(i,j,k)+gamma(i,j,k)
                        +fourthirds*beta(i,j,k-1)+gamma(i,j,k-1)
                        +beta_xz(i,j,k)+beta_xz(i+1,j,k)
                        +beta_yz(i,j,k)+beta_yz(i,j+1,k) ) * dxsqinv;

                   phiz(i,j,k) = phiz(i,j,k) + stag_mg_omega*(rhsz(i,j,k)-Lpz(i,j,k)) / fac;
        }
        }
        }
    }
#endif

}

void StagMGSolver::StagMGUpdate (std::array< MultiFab, AMREX_SPACEDIM >& phi_fc,
                                 const std::array< MultiFab, AMREX_SPACEDIM >& rhs_fc,
                                 const std::array< MultiFab, AMREX_SPACEDIM >& Lphi_fc,
                                 const std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
                                 const MultiFab& beta_cc,
                                 const std::array< MultiFab, NUM_EDGE >& beta_ed,
                                 const MultiFab& gamma_cc,
                                 const Real* dx,
                                 const int& color)
{
    BL_PROFILE_VAR("StagMGUpdate()",StagMGUpdate);

    AMREX_D_DECL(bool do_x,do_y,do_z);

    int offset = 1;

    if (color == 0) {
        AMREX_D_TERM(do_x = true;,
                     do_y = true;,
                     do_z = true;);

    }
    else if (color == 1 || color == 2) {
        AMREX_D_TERM(do_x = true;,
                     do_y = false;,
                     do_z = false;);
        offset = 2;
    }
    else if (color == 3 || color == 4) {
        AMREX_D_TERM(do_x = false;,
                     do_y = true;,
                     do_z = false;);
        offset = 2;
    }
#if (AMREX_SPACEDIM == 3)
    else if (color == 5 || color == 6) {
        AMREX_D_TERM(do_x = false;,
                     do_y = false;,
                     do_z = true;);
        offset = 2;
    }
#endif
    else {
        Abort("StagMGUpdate: Invalid Color");
    }

    GpuArray<Real,AMREX_SPACEDIM> dx_gpu{AMREX_D_DECL(dx[0], dx[1], dx[2])};

    // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for ( MFIter mfi(beta_cc,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        // Get the index space of the valid region
        const Box& bx = mfi.tilebox();

        AMREX_D_TERM(Array4<Real> const& phix_fab = phi_fc[0].array(mfi);,
                     Array4<Real> const& phiy_fab = phi_fc[1].array(mfi);,
                     Array4<Real> const& phiz_fab = phi_fc[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& rhsx_fab = rhs_fc[0].array(mfi);,
                     Array4<Real const> const& rhsy_fab = rhs_fc[1].array(mfi);,
                     Array4<Real const> const& rhsz_fab = rhs_fc[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& Lphix_fab = Lphi_fc[0].array(mfi);,
                     Array4<Real const> const& Lphiy_fab = Lphi_fc[1].array(mfi);,
                     Array4<Real const> const& Lphiz_fab = Lphi_fc[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& alphax_fab = alpha_fc[0].array(mfi);,
                     Array4<Real const> const& alphay_fab = alpha_fc[1].array(mfi);,
                     Array4<Real const> const& alphaz_fab = alpha_fc[2].array(mfi););

        Array4<Real const> const& beta_cc_fab = beta_cc.array(mfi);
        Array4<Real const> const& gamma_cc_fab = gamma_cc.array(mfi);

        Array4<Real const> const& beta_xy_fab = beta_ed[0].array(mfi);
#if (AMREX_SPACEDIM == 3)
        Array4<Real const> const& beta_xz_fab = beta_ed[1].array(mfi);
        Array4<Real const> const& beta_yz_fab = beta_ed[2].array(mfi);
#endif

        AMREX_D_TERM(const Box& bx_x = mfi.nodaltilebox(0);,
                     const Box& bx_y = mfi.nodaltilebox(1);,
                     const Box& bx_z = mfi.nodaltilebox(2););

        const Box& index_bounds = amrex::getIndexBounds(AMREX_D_DECL(bx_x,bx_y,bx_z));

        Real omega = stag_mg_omega;

        if (visc_type == 1) {

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(index_bounds, tbx,
            {
                stag_mg_update_visc_p1(tbx, AMREX_D_DECL(bx_x,bx_y,bx_z),
                                       AMREX_D_DECL(phix_fab,phiy_fab,phiz_fab),
                                       AMREX_D_DECL(rhsx_fab,rhsy_fab,rhsz_fab),
                                       AMREX_D_DECL(Lphix_fab,Lphiy_fab,Lphiz_fab),
                                       AMREX_D_DECL(alphax_fab,alphay_fab,alphaz_fab),
                                       beta_cc_fab, gamma_cc_fab,
                                       AMREX_D_DECL(do_x,do_y,do_z),
                                       offset, color, omega, dx_gpu);
            });

        }
        else if (visc_type == -1) {

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(index_bounds, tbx,
            {
                stag_mg_update_visc_m1(tbx, AMREX_D_DECL(bx_x,bx_y,bx_z),
                                       AMREX_D_DECL(phix_fab,phiy_fab,phiz_fab),
                                       AMREX_D_DECL(rhsx_fab,rhsy_fab,rhsz_fab),
                                       AMREX_D_DECL(Lphix_fab,Lphiy_fab,Lphiz_fab),
                                       AMREX_D_DECL(alphax_fab,alphay_fab,alphaz_fab),
                                       beta_cc_fab, beta_xy_fab,
#if (AMREX_SPACEDIM == 3)
                                       beta_xz_fab, beta_yz_fab,
#endif
                                       AMREX_D_DECL(do_x,do_y,do_z),
                                       offset, color, omega, dx_gpu);
            });

        }
        else if (visc_type == 2) {

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(index_bounds, tbx,
            {
                stag_mg_update_visc_p2(tbx, AMREX_D_DECL(bx_x,bx_y,bx_z),
                                       AMREX_D_DECL(phix_fab,phiy_fab,phiz_fab),
                                       AMREX_D_DECL(rhsx_fab,rhsy_fab,rhsz_fab),
                                       AMREX_D_DECL(Lphix_fab,Lphiy_fab,Lphiz_fab),
                                       AMREX_D_DECL(alphax_fab,alphay_fab,alphaz_fab),
                                       beta_cc_fab, gamma_cc_fab,
                                       AMREX_D_DECL(do_x,do_y,do_z),
                                       offset, color, omega, dx_gpu);
            });

        }
        else if (visc_type == -2) {

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(index_bounds, tbx,
            {
                stag_mg_update_visc_m2(tbx, AMREX_D_DECL(bx_x,bx_y,bx_z),
                                       AMREX_D_DECL(phix_fab,phiy_fab,phiz_fab),
                                       AMREX_D_DECL(rhsx_fab,rhsy_fab,rhsz_fab),
                                       AMREX_D_DECL(Lphix_fab,Lphiy_fab,Lphiz_fab),
                                       AMREX_D_DECL(alphax_fab,alphay_fab,alphaz_fab),
                                       beta_cc_fab, beta_xy_fab,
#if (AMREX_SPACEDIM == 3)
                                       beta_xz_fab, beta_yz_fab,
#endif
                                       AMREX_D_DECL(do_x,do_y,do_z),
                                       offset, color, omega, dx_gpu);
            });
        }
        else if (visc_type == 3) {

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(index_bounds, tbx,
            {
                stag_mg_update_visc_p3(tbx, AMREX_D_DECL(bx_x,bx_y,bx_z),
                                       AMREX_D_DECL(phix_fab,phiy_fab,phiz_fab),
                                       AMREX_D_DECL(rhsx_fab,rhsy_fab,rhsz_fab),
                                       AMREX_D_DECL(Lphix_fab,Lphiy_fab,Lphiz_fab),
                                       AMREX_D_DECL(alphax_fab,alphay_fab,alphaz_fab),
                                       beta_cc_fab, gamma_cc_fab,
                                       AMREX_D_DECL(do_x,do_y,do_z),
                                       offset, color, omega, dx_gpu);
            });

        }
        else if (visc_type == -3) {

            AMREX_LAUNCH_HOST_DEVICE_LAMBDA(index_bounds, tbx,
            {
                stag_mg_update_visc_m3(tbx, AMREX_D_DECL(bx_x,bx_y,bx_z),
                                       AMREX_D_DECL(phix_fab,phiy_fab,phiz_fab),
                                       AMREX_D_DECL(rhsx_fab,rhsy_fab,rhsz_fab),
                                       AMREX_D_DECL(Lphix_fab,Lphiy_fab,Lphiz_fab),
                                       AMREX_D_DECL(alphax_fab,alphay_fab,alphaz_fab),
                                       beta_cc_fab, beta_xy_fab,
#if (AMREX_SPACEDIM == 3)
                                       beta_xz_fab, beta_yz_fab,
#endif
                                       gamma_cc_fab,
                                       AMREX_D_DECL(do_x,do_y,do_z),
                                       offset, color, omega, dx_gpu);
            });
        }
    }
}
