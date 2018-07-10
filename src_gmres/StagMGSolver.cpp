#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"

#include "common_functions.H"
#include "common_namespace.H"
#include <AMReX_VisMF.H>

using namespace amrex;
using namespace gmres;
using namespace common;

// solve "(theta*alpha*I - L) phi = rhs" using multigrid with Gauss-Seidel relaxation
// if abs(visc_type) = 1, L = div beta grad
// if abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
// if abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
// if visc_type > 1 we assume constant coefficients
// if visc_type < 1 we assume variable coefficients
// beta_cc, and gamma_cc are cell-centered
// alpha_fc, phi_fc, and rhs_fc are face-centered
// beta_ed is nodal (2d) or edge-centered (3d)
// phi_fc must come in initialized to some value, preferably a reasonable guess
void StagMGSolver(const std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
                  const MultiFab& beta_cc,
                  const std::array< MultiFab, NUM_EDGE >& beta_ed,
                  const MultiFab& gamma_cc,
                  std::array< MultiFab, AMREX_SPACEDIM >& phi_fc,
                  const std::array< MultiFab, AMREX_SPACEDIM >& rhs_fc,
                  const Real& theta,
                  const Geometry& geom)
{
    Vector<int> is_periodic(AMREX_SPACEDIM);
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        is_periodic[i] = geom.isPeriodic(i);
    }

    // get the problem domain and boxarray at level 0
    Box pd_base = geom.Domain();
    BoxArray ba_base = beta_cc.boxArray();

    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    // compute the number of multigrid levels assuming stag_mg_minwidth is the length of the
    // smallest dimension of the smallest grid at the coarsest multigrid level
    int nlevs_mg = ComputeNlevsMG(ba_base);
    if (stag_mg_verbosity >= 3) {
        Print() << "Total number of multigrid levels: " << nlevs_mg << std::endl;
    }

    Vector<Geometry> geom_mg(nlevs_mg);
    int n;

    //////////////////////////////////
    // allocate multifabs used in multigrid coarsening

    // cell-centered
    Vector<MultiFab>  beta_cc_mg(nlevs_mg);
    Vector<MultiFab> gamma_cc_mg(nlevs_mg);

    // face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > alpha_fc_mg(nlevs_mg);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >   rhs_fc_mg(nlevs_mg);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >   phi_fc_mg(nlevs_mg);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >  Lphi_fc_mg(nlevs_mg);
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > resid_fc_mg(nlevs_mg);
    Vector<std::array< MultiFab, NUM_EDGE       > >  beta_ed_mg(nlevs_mg); // nodal in 2D, edge-based in 3D
    //////////////////////////////////

    // initial and current residuals
    Vector<Real> resid0(AMREX_SPACEDIM);
    Vector<Real> resid0_l2(AMREX_SPACEDIM);
    Vector<Real> resid(AMREX_SPACEDIM);
    Vector<Real> resid_l2(AMREX_SPACEDIM);
    Real resid_temp;

    int color_start, color_end;

    DistributionMapping dmap = beta_cc.DistributionMap();

    const Real* dx = geom.CellSize();
    Vector<std::array< Real, AMREX_SPACEDIM > > dx_mg(nlevs_mg);

    for (n=0; n<nlevs_mg; ++n) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            // compute dx at this level of multigrid
            dx_mg[n][d] = dx[d] * pow(2,n);
        }

        // create the problem domain for this multigrid level
        Box pd = pd_base.coarsen(pow(2,n));

        geom_mg[n].define(pd,&real_box,CoordSys::cartesian,is_periodic.data());

        // create the boxarray for this multigrid level
        BoxArray ba(ba_base); 
        ba.coarsen(pow(2,n));
     
        if ( n == 0 && !(ba == ba_base) ) {
            Abort("Finest multigrid level boxarray and coarsest problem boxarrays do not match");
        }

        // build multifabs used in multigrid coarsening
         beta_cc_mg[n].define(ba,dmap,1,1);
        gamma_cc_mg[n].define(ba,dmap,1,1);

        AMREX_D_TERM(alpha_fc_mg[n][0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                     alpha_fc_mg[n][1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                     alpha_fc_mg[n][2].define(convert(ba,nodal_flag_z), dmap, 1, 0););
        AMREX_D_TERM(  rhs_fc_mg[n][0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                       rhs_fc_mg[n][1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                       rhs_fc_mg[n][2].define(convert(ba,nodal_flag_z), dmap, 1, 0););
        AMREX_D_TERM(  phi_fc_mg[n][0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                       phi_fc_mg[n][1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                       phi_fc_mg[n][2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
        AMREX_D_TERM( Lphi_fc_mg[n][0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                      Lphi_fc_mg[n][1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                      Lphi_fc_mg[n][2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
        AMREX_D_TERM(resid_fc_mg[n][0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                     resid_fc_mg[n][1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                     resid_fc_mg[n][2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

        // build beta_ed_mg
        if (AMREX_SPACEDIM == 2) {
            beta_ed_mg[n][0].define(convert(ba,nodal_flag), dmap, 1, 0);
        }
        else if (AMREX_SPACEDIM == 3) {
            beta_ed_mg[n][0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
            beta_ed_mg[n][1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
            beta_ed_mg[n][2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
        }
    } // end loop over multigrid levels

    // copy level 1 coefficients into mg array of coefficients
    MultiFab::Copy(beta_cc_mg[0],beta_cc,0,0,1,1);
    MultiFab::Copy(gamma_cc_mg[0],gamma_cc,0,0,1,1);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(alpha_fc_mg[0][d],alpha_fc[d],0,0,1,0);
        // multiply alpha_fc_mg by theta
        alpha_fc_mg[0][d].mult(theta,0,1,0);
    }
    MultiFab::Copy(beta_ed_mg[0][0],beta_ed[0],0,0,1,0);
    if (AMREX_SPACEDIM == 3) {
        MultiFab::Copy(beta_ed_mg[0][1],beta_ed[1],0,0,1,0);
        MultiFab::Copy(beta_ed_mg[0][2],beta_ed[2],0,0,1,0);
    }
    
    // coarsen coefficients
    for (n=1; n<nlevs_mg; ++n) {
        // need ghost cells set to zero to prevent intermediate NaN states
        // that cause some compilers to fail
        beta_cc_mg[n].setVal(0.);
        gamma_cc_mg[n].setVal(0.);

        // cc_restriction on beta_cc_mg and gamma_cc_mg
        CCRestriction( beta_cc_mg[n], beta_cc_mg[n-1],geom_mg[n]);
        CCRestriction(gamma_cc_mg[n],gamma_cc_mg[n-1],geom_mg[n]);

        // stag_restriction on alpha_fc_mg
        StagRestriction(alpha_fc_mg[n],alpha_fc_mg[n-1],1);

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

        // fill periodic ghost cells
        phi_fc_mg[0][d].FillBoundary(geom_mg[0].periodicity());

    }

    // set rhs_fc_mg at level 1 by copying in passed-in rhs_fc
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(rhs_fc_mg[0][d],rhs_fc[d],0,0,1,0);
    }

    // compute norm of initial residual
    // first compute Lphi
    StagApplyOp(beta_cc_mg[0],gamma_cc_mg[0],beta_ed_mg[0],
                phi_fc_mg[0],Lphi_fc_mg[0],alpha_fc_mg[0],dx_mg[0].data());
    
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        // compute Lphi - rhs
        MultiFab::Subtract(Lphi_fc_mg[0][d],rhs_fc_mg[0][d],0,0,1,0);

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
                StagApplyOp(beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                            phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data());
    
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
                    StagApplyOp(beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                                phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),color);

                    // update phi = phi + omega*D^{-1}*(rhs-Lphi)
                    StagMGUpdate(phi_fc_mg[n],rhs_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],
                                 beta_cc_mg[n],beta_ed_mg[n],gamma_cc_mg[n],dx_mg[n].data(),color);

                    for (int d=0; d<AMREX_SPACEDIM; ++d) {

                        // fill periodic ghost cells
                        phi_fc_mg[n][d].FillBoundary(geom_mg[n].periodicity());

                    }

                } // end loop over colors

                // print out residual
                if (stag_mg_verbosity >= 4) {

                    // compute Lphi
                    StagApplyOp(beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                                phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data());
    
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
            StagApplyOp(beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                        phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data());

            for (int d=0; d<AMREX_SPACEDIM; ++d) {

                // compute Lphi - rhs, and then multiply by -1
                MultiFab::Subtract(Lphi_fc_mg[n][d],rhs_fc_mg[n][d],0,0,1,0);
                Lphi_fc_mg[n][d].mult(-1.,0,1,0);
                if (stag_mg_verbosity >= 3) {
                    resid_temp = Lphi_fc_mg[n][d].norm0();
                    Print() << "Residual for comp " << d << " after all smooths at level " 
                            << n << " " << resid_temp << std::endl;
                }

                // fill periodic ghost cells
                Lphi_fc_mg[n][d].FillBoundary(geom_mg[n].periodicity());

            }

            // restrict/coarsen residual and put it in rhs_fc
            StagRestriction(rhs_fc_mg[n+1],Lphi_fc_mg[n]);

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
            StagApplyOp(beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                        phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data());

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
                StagApplyOp(beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                            phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),color);

                // update phi = phi + omega*D^{-1}*(rhs-Lphi)
                StagMGUpdate(phi_fc_mg[n],rhs_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],
                             beta_cc_mg[n],beta_ed_mg[n],gamma_cc_mg[n],dx_mg[n].data(),color);

                for (int d=0; d<AMREX_SPACEDIM; ++d) {

                    // fill periodic ghost cells
                    phi_fc_mg[n][d].FillBoundary(geom_mg[n].periodicity());

                }
            } // end loop over colors

            // print out residual
            if (stag_mg_verbosity >= 4) {
                
                // compute Lphi
                StagApplyOp(beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                            phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data());
                
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    // compute Lphi - rhs, and report residual
                    MultiFab::Subtract(Lphi_fc_mg[n][d],rhs_fc_mg[n][d],0,0,1,0);
                    resid_temp = Lphi_fc_mg[n][d].norm0();
                    Print() << "Residual for comp " << d << " after    smooth " << m << " at level " 
                            << n << " " << resid_temp << std::endl;
                }
            }

        } // end loop over nsmooths

        ////////////////////
        // compute residual

        // compute Lphi
        StagApplyOp(beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                    phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data());
        
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            
            // compute Lphi - rhs, and then multiply by -1
            MultiFab::Subtract(Lphi_fc_mg[n][d],rhs_fc_mg[n][d],0,0,1,0);
            Lphi_fc_mg[n][d].mult(-1.,0,1,0);
            if (stag_mg_verbosity >= 3) {
                resid_temp = Lphi_fc_mg[n][d].norm0();
                Print() << "Residual for comp " << d << " after all smooths at level " 
                        << n << " " << resid_temp << std::endl;
            }

            // fill periodic ghost cells
            Lphi_fc_mg[n][d].FillBoundary(geom_mg[n].periodicity());

        }

        if (stag_mg_verbosity >= 3) {
            Print() << "End bottom solve" << std::endl;
        }

        // up the V-cycle
        for (n=nlevs_mg-2; n>=0; --n) {

            // prolongate/interpolate correction to update phi
            StagProlongation(phi_fc_mg[n+1],phi_fc_mg[n]);

            for (int d=0; d<AMREX_SPACEDIM; ++d) {

                // fill periodic ghost cells
                phi_fc_mg[n][d].FillBoundary(geom_mg[n].periodicity());

            }

            // print out residual
            if (stag_mg_verbosity >= 3) {

                // compute Lphi
                StagApplyOp(beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                            phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data());
    
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
                    StagApplyOp(beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                                phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data(),color);

                    // update phi = phi + omega*D^{-1}*(rhs-Lphi)
                    StagMGUpdate(phi_fc_mg[n],rhs_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],
                                 beta_cc_mg[n],beta_ed_mg[n],gamma_cc_mg[n],dx_mg[n].data(),color);

                    for (int d=0; d<AMREX_SPACEDIM; ++d) {

                        // fill periodic ghost cells
                        phi_fc_mg[n][d].FillBoundary(geom_mg[n].periodicity());

                    }

                } // end loop over colors

                // print out residual
                if (stag_mg_verbosity >= 4) {

                    // compute Lphi
                    StagApplyOp(beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                                phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data());
    
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
                StagApplyOp(beta_cc_mg[n],gamma_cc_mg[n],beta_ed_mg[n],
                            phi_fc_mg[n],Lphi_fc_mg[n],alpha_fc_mg[n],dx_mg[n].data());
        
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
        StagApplyOp(beta_cc_mg[0],gamma_cc_mg[0],beta_ed_mg[0],
                    phi_fc_mg[0],Lphi_fc_mg[0],alpha_fc_mg[0],dx_mg[0].data());

        // compute Lphi - rhs
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            // compute Lphi - rhs
            MultiFab::Subtract(Lphi_fc_mg[0][d],rhs_fc_mg[0][d],0,0,1,0);
        }
	
	//////// TEST L2 INNER PROD /////////////////////
	// for (int d=0; d<AMREX_SPACEDIM; ++d) {
	//   phi_fc_mg[0][d].setVal(1.);
	//   Lphi_fc_mg[0][d].setVal(1.);
	// }
	// StagInnerProd(phi_fc_mg[0],0,Lphi_fc_mg[0],0,resid_l2);
	// for (int d=0; d<AMREX_SPACEDIM; ++d) {
	//   Print() << resid_l2[d] << "\n";
	// }
	// // amrex::Abort("Exit L2 Test");
	// exit(0);
	//////////////////////////////////////////////////

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

        // fill periodic ghost cells
        phi_fc[d].FillBoundary(geom.periodicity());

    }

    // vcycle_counter += AMREX_SPACEDIM*stag_mg_max_vcycles;

    if (stag_mg_verbosity >= 1) {
        Print() << "\nEnd call to stag_mg_solver\n";
    }

}

// compute the number of multigrid levels assuming minwidth is the length of the
// smallest dimension of the smallest grid at the coarsest multigrid level
int ComputeNlevsMG(const BoxArray& ba) {

    int nlevs_mg = -1;

    for (int i=0; i<ba.size(); ++i) {
        Box bx = ba.get(i);
        IntVect iv = bx.bigEnd() - bx.smallEnd() + IntVect(1);

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            int temp = iv[d];
            int rdir = 1;
            while (temp%2 == 0 && temp/stag_mg_minwidth != 1) {
                temp /= 2;
                ++rdir;
            }

            if (nlevs_mg == -1) {
                nlevs_mg = rdir;
            }
            else {
                nlevs_mg = std::min(rdir,nlevs_mg);
            }
        }
    }

    return nlevs_mg;
}

void CCRestriction(MultiFab& phi_c, const MultiFab& phi_f, const Geometry& geom_c)
{
    // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for ( MFIter mfi(phi_c); mfi.isValid(); ++mfi ) {

        // Get the index space of the valid region
        const Box& validBox = mfi.validbox();

        cc_restriction(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                       BL_TO_FORTRAN_3D(phi_c[mfi]),
                       BL_TO_FORTRAN_3D(phi_f[mfi]));
    }

    phi_c.FillBoundary(geom_c.periodicity());

}

void StagRestriction(std::array< MultiFab, AMREX_SPACEDIM >& phi_c, 
                     const std::array< MultiFab, AMREX_SPACEDIM >& phi_f,
                     int simple_stencil)
{

    // loop over boxes (note we are not passing in a cell-centered MultiFab)
    for ( MFIter mfi(phi_c[0]); mfi.isValid(); ++mfi ) {

        // Get the index space of the valid region
        // there are no cell-centered MultiFabs so use this to get
        // a cell-centered box
        const Box& validBox = amrex::enclosedCells(mfi.validbox());

        stag_restriction(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                         BL_TO_FORTRAN_3D(phi_c[0][mfi]),
                         BL_TO_FORTRAN_3D(phi_f[0][mfi]),
                         BL_TO_FORTRAN_3D(phi_c[1][mfi]),
                         BL_TO_FORTRAN_3D(phi_f[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                         BL_TO_FORTRAN_3D(phi_c[2][mfi]),
                         BL_TO_FORTRAN_3D(phi_f[2][mfi]),
#endif
                         &simple_stencil);
    }
}

void NodalRestriction(MultiFab& phi_c, const MultiFab& phi_f)
{
    // loop over boxes (note we are not passing in a cell-centered MultiFab)
    for ( MFIter mfi(phi_c); mfi.isValid(); ++mfi ) {

        // Get the index space of the valid region
        // there are no cell-centered MultiFabs so use this to get
        // a cell-centered box
        const Box& validBox = amrex::enclosedCells(mfi.validbox());

        nodal_restriction(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                          BL_TO_FORTRAN_3D(phi_c[mfi]),
                          BL_TO_FORTRAN_3D(phi_f[mfi]));

    }
}

void EdgeRestriction(std::array< MultiFab, NUM_EDGE >& phi_c, 
                     const std::array< MultiFab, NUM_EDGE >& phi_f)
{
    if (AMREX_SPACEDIM != 3) {
        Abort("Edge restriction can only be called for 3D!");
    }

    // loop over boxes (note we are not passing in a cell-centered MultiFab)
    for ( MFIter mfi(phi_c[0]); mfi.isValid(); ++mfi ) {

        // Get the index space of the valid region
        // there are no cell-centered MultiFabs so use this to get
        // a cell-centered box
        const Box& validBox = amrex::enclosedCells(mfi.validbox());

        edge_restriction(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                         BL_TO_FORTRAN_3D(phi_c[0][mfi]),
                         BL_TO_FORTRAN_3D(phi_f[0][mfi]),
                         BL_TO_FORTRAN_3D(phi_c[1][mfi]),
                         BL_TO_FORTRAN_3D(phi_f[1][mfi]),
                         BL_TO_FORTRAN_3D(phi_c[2][mfi]),
                         BL_TO_FORTRAN_3D(phi_f[2][mfi]));
    }
}

void StagProlongation(const std::array< MultiFab, AMREX_SPACEDIM >& phi_c, 
                      std::array< MultiFab, AMREX_SPACEDIM >& phi_f)
{

    // loop over boxes (note we are not passing in a cell-centered MultiFab)
    for ( MFIter mfi(phi_f[0]); mfi.isValid(); ++mfi ) {

        // Get the index space of the valid region
        // there are no cell-centered MultiFabs so use this to get
        // a cell-centered box
        const Box& validBox = amrex::enclosedCells(mfi.validbox());

        stag_prolongation(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                          BL_TO_FORTRAN_3D(phi_c[0][mfi]),
                          BL_TO_FORTRAN_3D(phi_f[0][mfi]),
                          BL_TO_FORTRAN_3D(phi_c[1][mfi]),
                          BL_TO_FORTRAN_3D(phi_f[1][mfi])
#if (AMREX_SPACEDIM == 3)
                        , BL_TO_FORTRAN_3D(phi_c[2][mfi]),
                          BL_TO_FORTRAN_3D(phi_f[2][mfi])
#endif
                          );
    }

}

void StagMGUpdate(std::array< MultiFab, AMREX_SPACEDIM >& phi_fc,
                  const std::array< MultiFab, AMREX_SPACEDIM >& rhs_fc,
                  const std::array< MultiFab, AMREX_SPACEDIM >& Lphi_fc,
                  const std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
                  const MultiFab& beta_cc,
                  const std::array< MultiFab, NUM_EDGE >& beta_ed,
                  const MultiFab& gamma_cc,
                  const Real* dx,
                  const int& color)
{

    // loop over boxes (make sure mfi takes a cell-centered multifab as an argument)
    for ( MFIter mfi(beta_cc); mfi.isValid(); ++mfi ) {

        // Get the index space of the valid region
        const Box& validBox = mfi.validbox();


        stag_mg_update(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
                       BL_TO_FORTRAN_3D(phi_fc[0][mfi]),
                       BL_TO_FORTRAN_3D(phi_fc[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                       BL_TO_FORTRAN_3D(phi_fc[2][mfi]),
#endif
                       BL_TO_FORTRAN_3D(rhs_fc[0][mfi]),
                       BL_TO_FORTRAN_3D(rhs_fc[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                       BL_TO_FORTRAN_3D(rhs_fc[2][mfi]),
#endif
                       BL_TO_FORTRAN_3D(Lphi_fc[0][mfi]),
                       BL_TO_FORTRAN_3D(Lphi_fc[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                       BL_TO_FORTRAN_3D(Lphi_fc[2][mfi]),
#endif
                       BL_TO_FORTRAN_3D(alpha_fc[0][mfi]),
                       BL_TO_FORTRAN_3D(alpha_fc[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                       BL_TO_FORTRAN_3D(alpha_fc[2][mfi]),
#endif
                       BL_TO_FORTRAN_3D(beta_cc[mfi]),
#if (AMREX_SPACEDIM == 2)
                       BL_TO_FORTRAN_3D(beta_ed[0][mfi]),
#elif (AMREX_SPACEDIM == 3)
                       BL_TO_FORTRAN_3D(beta_ed[0][mfi]),
                       BL_TO_FORTRAN_3D(beta_ed[1][mfi]),
                       BL_TO_FORTRAN_3D(beta_ed[2][mfi]),
#endif
                       BL_TO_FORTRAN_3D(gamma_cc[mfi]),
                       dx, &color);



    }
}
