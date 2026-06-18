#include <AMReX_VisMF.H>

#include "common_functions.H"

#include "gmres_functions.H"


#include <ib_functions.H>


#include <IBParticleContainer.H>




void IBGMRES(std::array<MultiFab, AMREX_SPACEDIM> & b_u, const MultiFab & b_p,
             std::array<MultiFab, AMREX_SPACEDIM> & x_u, MultiFab & x_p,
             std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
             MultiFab & beta, std::array<MultiFab, NUM_EDGE> & beta_ed,
             MultiFab & gamma, Real theta_alpha,
             const IBParticleContainer & ib_pc,
             const Geometry & geom,
             Real & norm_pre_rhs) {

    BL_PROFILE_VAR("GMRES()", GMRES);

    if (gmres_verbose >= 1) {
        Print() << "Begin call to GMRES" << std::endl;
    }

    Vector<Real> cs(gmres_max_inner);
    Vector<Real> sn(gmres_max_inner);
    Vector<Real>  y(gmres_max_inner);
    Vector<Real>  s(gmres_max_inner + 1);

    Vector<Vector<Real>> H(gmres_max_inner + 1, Vector<Real>(gmres_max_inner));

    int outer_iter, total_iter, i_copy; // for looping iteration
    int i=0;

    Real norm_b;            // |b|;           computed once at beginning
    Real norm_pre_b;        // |M^-1 b|;      computed once at beginning
    Real norm_resid;        // |M^-1 (b-Ax)|; computed at beginning of each outer iteration
    Real norm_init_resid;   // |M^-1 (b-Ax)|; computed once at beginning
    Real norm_resid_Stokes; // |b-Ax|;        computed at beginning of each outer iteration
    Real norm_init_Stokes;  // |b-Ax|;        computed once at beginning
    Real norm_u_noprecon;   // u component of norm_resid_Stokes
    Real norm_p_noprecon;   // p component of norm_resid_Stokes
    Real norm_lambda_noprecon;
    Real norm_resid_est;

    Real norm_u; // temporary norms used to build full-state norm
    Real norm_p; // temporary norms used to build full-state norm
    Real norm_lambda;

    Vector<Real> inner_prod_vel(AMREX_SPACEDIM);
    Real inner_prod_pres;
    Real inner_prod_lambda;

    BoxArray ba              = b_p.boxArray();
    DistributionMapping dmap = b_p.DistributionMap();


    // # of ghost cells must match x_u so higher-order stencils can work
    std::array< MultiFab, AMREX_SPACEDIM > r_u;
    std::array< MultiFab, AMREX_SPACEDIM > w_u;
    std::array< MultiFab, AMREX_SPACEDIM > tmp_u;
    std::array< MultiFab, AMREX_SPACEDIM > V_u;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
          r_u[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, x_u[d].nGrow());
          w_u[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);
        tmp_u[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);
          V_u[d].define(convert(ba, nodal_flag_dir[d]), dmap, gmres_max_inner + 1, 0);
    }


    // # of ghost cells must match x_p so higher-order stencils can work
    MultiFab r_p  (ba, dmap,                  1, x_p.nGrow());
    MultiFab w_p  (ba, dmap,                  1, 0);
    MultiFab tmp_p(ba, dmap,                  1, 0);
    MultiFab V_p  (ba, dmap,gmres_max_inner + 1, 0); // Krylov vectors


    //___________________________________________________________________________
    // Get all the immersed-boudary particle indices (used to iterate below)
    Vector<IBP_info> ibp_info;
    // TODO: make `x_lambda` a referenced argument
    std::map<std::pair<int, int>, Vector<RealVect>>         b_lambda;
    std::map<std::pair<int, int>, Vector<RealVect>>         x_lambda;
    std::map<std::pair<int, int>, Vector<RealVect>>         r_lambda;
    std::map<std::pair<int, int>, Vector<RealVect>>         w_lambda;
    std::map<std::pair<int, int>, Vector<RealVect>>       tmp_lambda;
    std::map<std::pair<int, int>, Vector<Vector<RealVect>>> V_lambda;

    std::map<std::pair<int, int>, Vector<RealVect>> marker_pos;

    int ibpc_lev = 0; // assume single level for now
    int ib_grow  = 6; // using the 6-point stencil

    // NOTE: use `ib_pc` BoxArray to collect IB particle data
    MultiFab dummy(ib_pc.ParticleBoxArray(ibpc_lev),
                   ib_pc.ParticleDistributionMap(ibpc_lev), 1, 1);
    for (MFIter mfi(dummy, ib_pc.tile_size); mfi.isValid(); ++mfi){
        IBParticleContainer::PairIndex index(mfi.index(), mfi.LocalTileIndex());
        ib_pc.IBParticleInfo(ibp_info, ibpc_lev, index, true);
    }

    Vector<std::pair<int, int>> part_indices(ibp_info.size());
    for (int i=0; i<ibp_info.size(); ++i) {
        part_indices[i] = ibp_info[i].asPairIndex();

        // Pre-allocate particle arrays
        const Vector<RealVect> marker_positions = ib_pc.MarkerPositions(0, part_indices[i]);
        // ... initialized to (0..0) by default constructor
          b_lambda[part_indices[i]].resize(marker_positions.size());
          x_lambda[part_indices[i]].resize(marker_positions.size());
          r_lambda[part_indices[i]].resize(marker_positions.size());
          w_lambda[part_indices[i]].resize(marker_positions.size());
        tmp_lambda[part_indices[i]].resize(marker_positions.size());

          V_lambda[part_indices[i]].resize(gmres_max_inner + 1);
        for (int j=0; j<gmres_max_inner+1; ++j)
            V_lambda[part_indices[i]][j].resize(marker_positions.size());

        // Fill these with initial values
        marker_pos[part_indices[i]] = marker_positions;
    }


    // DEBUG:
    Print() << "Found " << part_indices.size() << " many IB particles in rank 0:"
            << std::endl;
    for (const auto & pid : part_indices)
        Print() << "[" << pid.first << ", " << pid.second << "]";
    Print() << std::endl << std::endl;



    /****************************************************************************
     *                                                                          *
     * Preflight work: apply scaling and compute perconditioned norms_b         *
     *                                                                          *
     ***************************************************************************/


    // TODO: Assumptions: dx=dy=dz (approx) as well as constant viscosity
    const Real * dx = geom.CellSize();
    const Real c1 = dx[0]/visc_coef;
    const Real c2 = dx[0]*dx[0];
    const Real c3 = 1./dx[0];
    const Real c4 = 1./(visc_coef*dx[0]);

    scale_factor = c1;

    // apply scaling factor
    if (scale_factor != 1.) {
        theta_alpha = theta_alpha*scale_factor;

        // we will solve for scale*x_p so we need to scale the initial guess
        x_p.mult(scale_factor, 0, 1, x_p.nGrow());

        // scale the rhs:
        for (int d=0; d<AMREX_SPACEDIM; ++d)
            b_u[d].mult(scale_factor, 0, 1, b_u[d].nGrow());

        // scale the viscosities:
        beta.mult(scale_factor, 0, 1, beta.nGrow());
        gamma.mult(scale_factor, 0, 1, gamma.nGrow());
        for (int d=0; d<NUM_EDGE; ++d)
            beta_ed[d].mult(scale_factor, 0, 1, beta_ed[d].nGrow());
    }


    //___________________________________________________________________________
    // First application of preconditioner


    // Dummy cell-centered MultiFab (used to construct MFIters)
    MultiFab dummy_iter(ba, dmap, 1, 0);

    // ApplyPrecon(b_u, b_p, tmp_u, tmp_p,
    //             alpha_fc, beta, beta_ed, gamma, theta_alpha,
    //             geom);
    IBMPrecon(b_u, b_p, tmp_u, tmp_p, alpha_fc, beta, beta_ed, gamma, theta_alpha,
              ib_pc, part_indices, tmp_lambda, b_lambda, geom);


    // preconditioned norm_b: norm_pre_b
    StagL2Norm(tmp_u, 0, norm_u);
    CCL2Norm(tmp_p, 0, norm_p);
    MarkerL2Norm(part_indices, dummy_iter, geom, marker_pos, tmp_lambda, norm_lambda);
    norm_p       = p_norm_weight*norm_p;
    norm_lambda  = p_norm_weight*norm_lambda; // TODO: use p_norm_weight for now
    norm_pre_b   = sqrt(norm_u*norm_u + norm_p*norm_p + norm_lambda*norm_lambda);
    norm_pre_rhs = norm_pre_b;


    // calculate the l2 norm of rhs
    StagL2Norm(b_u, 0, norm_u);
    CCL2Norm(b_p, 0, norm_p);
    MarkerL2Norm(part_indices, dummy_iter, geom, marker_pos, b_lambda, norm_lambda);
    norm_p      = p_norm_weight*norm_p;
    norm_lambda = p_norm_weight*norm_lambda; // TODO: use p_norm_weight for now
    norm_b      = sqrt(norm_u*norm_u + norm_p*norm_p + norm_lambda*norm_lambda);

    //! If norm_b=0 we should return zero as the solution and "return" from this routine
    // It is important to use gmres_abs_tol and not 0 since sometimes due to roundoff we
    // get a nonzero number that should really be zero
    if (gmres_verbose >= 1) {
        // Useful to print out to give expected scale for gmres_abs_tol
        Print() << "GMRES.cpp: GMRES called with ||rhs||=" << norm_b << std::endl;
    }

    if (norm_b <= gmres_abs_tol) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            x_u[d].setVal(0.);
        }
        x_p.setVal(0.);

        if (gmres_verbose >= 1) {
            Print() << "GMRES.cpp: converged in 0 iterations since rhs=0" << std::endl;
        }
        return;
    }



    ///////////////////
    // begin outer iteration
    ///////////////////



    total_iter = 0;
    outer_iter = 0;

    do {

        //_______________________________________________________________________
        // Compute tmp = b - Ax

        // Calculate tmp = Ax
        // Fluid part: ........................................... (v, p) = (Av - Gp, -Dv)
        ApplyMatrix(tmp_u, tmp_p, x_u, x_p,
                    alpha_fc, beta, beta_ed, gamma, theta_alpha,
                    geom);

        // IBM part: ............................. (v, lambda) = (Av - Gp - S lambda, -Jv)
        ApplyIBM(tmp_u, tmp_lambda, x_u, ib_pc, part_indices, x_lambda,
                 ib_grow, ibpc_lev, geom);

        // tmp = b - Ax
        // Fluid part: ........................... tmp_u = b_u - (Ax_u - Gx_p - Sx_lambda)
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Subtract(tmp_u[d], b_u[d], 0, 0, 1, 0);
            tmp_u[d].mult(-1., 0, 1, 0);
        }

        // Pressure part: ............................................ tmp_p = b_p - (-Dv)
        MultiFab::Subtract(tmp_p, b_p, 0, 0, 1, 0);
        tmp_p.mult(-1., 0, 1, 0);

        // IBM part: ....................................... tmp_lambda = b_lambda - (-Jv)
        MarkerInvSub(part_indices, tmp_lambda, b_lambda);


        //_______________________________________________________________________
        // un-preconditioned residuals
        StagL2Norm(tmp_u, 0, norm_u_noprecon);
        CCL2Norm(tmp_p, 0, norm_p_noprecon);
        MarkerL2Norm(part_indices, dummy_iter, geom, marker_pos,
                     tmp_lambda, norm_lambda_noprecon);

        norm_p_noprecon      = p_norm_weight*norm_p_noprecon;
        norm_lambda_noprecon = p_norm_weight*norm_lambda_noprecon; // TODO: use p_norm_weight for now
        norm_resid_Stokes    = sqrt(norm_u_noprecon*norm_u_noprecon
                                    + norm_p_noprecon*norm_p_noprecon
                                    + norm_lambda_noprecon*norm_lambda_noprecon);

        if (outer_iter == 0)
            norm_init_Stokes = norm_resid_Stokes;


        //_______________________________________________________________________
        // Print verbose output
        if (gmres_verbose >= 2) {
            Print() << "total Iters = " << total_iter << std::endl;
            Print() << "r/(r_0,b) = " << norm_resid_Stokes/norm_init_Stokes << "  "
                    << norm_resid_Stokes/norm_b << std::endl;
        }

        if (gmres_verbose >= 3) {
            Print() << "un-Precond. rel. resid. (u,v,p,lambda) = "
                    << norm_resid_Stokes/norm_init_Stokes    << "  "
                    << norm_u_noprecon/norm_init_Stokes      << "  "
                    << norm_p_noprecon/norm_init_Stokes      << "  "
                    << norm_lambda_noprecon/norm_init_Stokes << std::endl;
        }


        //_______________________________________________________________________
        // solve for r = M^{-1} tmp
        // We should not be counting these toward the number of mg cycles performed
        // ApplyPrecon(tmp_u, tmp_p, r_u, r_p,
        //             alpha_fc, beta, beta_ed, gamma, theta_alpha,
        //             geom);
        IBMPrecon(tmp_u, tmp_p, r_u, r_p, alpha_fc, beta, beta_ed, gamma, theta_alpha,
                  ib_pc, part_indices, r_lambda, tmp_lambda, geom);


        // resid = sqrt(dot_product(r, r))
        StagL2Norm(r_u, 0, norm_u);
        CCL2Norm(r_p, 0, norm_p);
        MarkerL2Norm(part_indices, dummy_iter, geom, marker_pos, r_lambda, norm_lambda);
        norm_p      = p_norm_weight*norm_p;
        norm_lambda = p_norm_weight*norm_lambda; // TODO: use p_norm_weight for now
        norm_resid  = sqrt(norm_u*norm_u + norm_p*norm_p + norm_lambda*norm_lambda);


        // If first iteration, save the initial preconditioned residual
        if (outer_iter==0) {
            norm_init_resid = norm_resid;
            norm_resid_est  = norm_resid;
        }


        //_______________________________________________________________________
        // Print verbose output
        if (gmres_verbose >= 3) {
            Print() << "Precond. rel. res. (u,v,p,lambda) = "
                    << norm_resid/norm_init_resid  << "  "
                    << norm_u/norm_init_resid      << "  "
                    << norm_p/norm_init_resid      << "  "
                    << norm_lambda/norm_init_resid << std::endl;
        }


        //_______________________________________________________________________
        // We need to test the residual now and exit OuterLoop if converged
        if (total_iter >= gmres_max_iter) {
            if (gmres_verbose >= 1) {
                Print() << "GMRES did not converge in max number of total inner iterations: Exiting"
                        << std::endl;
            }

            break;

        } else if (total_iter >= gmres_min_iter) {
            // other options
            if(norm_resid <= gmres_rel_tol*amrex::min(norm_pre_b, norm_init_resid)) {
                if (gmres_verbose >= 2) {
                    Print() << "GMRES converged: Outer = " << outer_iter << ",  Inner = " << i
                            << " Total=" << total_iter << std::endl;
                }

                if (norm_resid_Stokes >= 10*gmres_rel_tol*amrex::min(norm_b, norm_init_Stokes)) {
                    Print() << "GMRES.cpp: Warning: gmres may not have converged: |r|/|b|= "
                            << norm_resid_Stokes/norm_b << " |r|/|r0|="
                            << norm_resid_Stokes/norm_init_Stokes << std::endl;
                }

                // Only exit if the *true* preconditioned residual is less than tolerance:
                // Do not trust the gmres estimate

                break; // exit OuterLoop

            } else if (norm_resid <= gmres_abs_tol) {

                if (gmres_verbose >= 2) {
                    Print() << "GMRES converged: Outer = " << outer_iter << ",  Inner = " << i
                            << " Total=" << total_iter << std::endl;
                }

                break; // exit OuterLoop
            }
        }


        if (outer_iter >= gmres_max_outer) {
            Print() << "GMRES did not converge in max number of outer iterations: Exiting" << std::endl;
            break; // exit OuterLoop
        }

        outer_iter = outer_iter + 1;


        //_______________________________________________________________________
        // Print verbose output
        if (gmres_verbose >= 3) {
            Print() << "Begin outer iteration " << outer_iter << std::endl;
        }


        //_______________________________________________________________________
        // Create the first basis in Krylov space: V(1) = r / norm(r)
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Copy(V_u[d], r_u[d], 0, 0, 1, 0);
            V_u[d].mult(1./norm_resid, 0, 1, 0);
        }

        MultiFab::Copy(V_p, r_p, 0, 0, 1, 0);
        V_p.mult(1./norm_resid, 0, 1, 0);

        MarkerCopy(part_indices, 0, V_lambda, r_lambda);
        MarkerMult(part_indices, 0, 1./norm_resid, V_lambda);

        // s = norm(r) * e_0
        std::fill(s.begin(), s.end(), 0.);
        s[0] = norm_resid;



        ///////////////////////
        // begin inner iteration
        ///////////////////////



        // i is the inner iteration loop index
        for (i=0; i<gmres_max_inner; ++i) {


            //___________________________________________________________________
            // Print verbose output
            if (gmres_verbose >= 3) {
                Print() << "Begin inner iteration " << i+1 << std::endl;
            }


            total_iter = total_iter + 1;
            i_copy     = i;


            //___________________________________________________________________
            // tmp=A*V(i)

            // use r_p and r_u as temporaries to hold ith component of V
            for (int d=0; d<AMREX_SPACEDIM; ++d)
                MultiFab::Copy(r_u[d], V_u[d], i, 0, 1, 0);

            MultiFab::Copy(r_p, V_p, i, 0, 1, 0);

            MarkerCopy(part_indices, i, r_lambda, V_lambda);

            ApplyMatrix(tmp_u, tmp_p, r_u, r_p,
                        alpha_fc, beta, beta_ed, gamma, theta_alpha,
                        geom);

            ApplyIBM(tmp_u, tmp_lambda, r_u, ib_pc, part_indices, r_lambda,
                     ib_grow, ibpc_lev, geom);


            //___________________________________________________________________
            // w = M^{-1} A*V(i)
            // ApplyPrecon(tmp_u, tmp_p, w_u, w_p,
            //             alpha_fc, beta, beta_ed, gamma, theta_alpha,
            //             geom);
            IBMPrecon(tmp_u, tmp_p, w_u, w_p,
                      alpha_fc, beta, beta_ed, gamma, theta_alpha,
                      ib_pc, part_indices, w_lambda, tmp_lambda, geom);


            //___________________________________________________________________
            // Form Hessenberg matrix H
            for (int k=0; k<=i; ++k) {
                // H(k,i) = dot_product(w, V(k))
                //        = dot_product(w_u, V_u(k))+dot_product(w_p, V_p(k))
                StagInnerProd(w_u, 0, V_u, k, inner_prod_vel);
                CCInnerProd(w_p, 0, V_p, k, inner_prod_pres);
                MarkerInnerProd(part_indices, k, dummy_iter, geom, marker_pos,
                                w_lambda, V_lambda, inner_prod_lambda);
                H[k][i] = std::accumulate(inner_prod_vel.begin(), inner_prod_vel.end(), 0.)
                          + pow(p_norm_weight, 2.0)*inner_prod_pres
                          + pow(p_norm_weight, 2.0)*inner_prod_lambda; // TODO: use p_norm_weight for now


                // w = w - H(k,i) * V(k)
                // use tmp_u and tmp_p as temporaries to hold kth component of V(k)
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    MultiFab::Copy(tmp_u[d], V_u[d], k, 0, 1, 0);
                    tmp_u[d].mult(H[k][i], 0, 1, 0);
                    MultiFab::Subtract(w_u[d], tmp_u[d], 0, 0, 1, 0);
                }

                MultiFab::Copy(tmp_p, V_p, k, 0, 1, 0);
                tmp_p.mult(H[k][i], 0, 1, 0);
                MultiFab::Subtract(w_p, tmp_p, 0, 0, 1, 0);

                MarkerCopy(part_indices, k, tmp_lambda, V_lambda);
                MarkerMult(part_indices, H[k][i], tmp_lambda);
                MarkerSub(part_indices, w_lambda, tmp_lambda);
            }

            // H(i+1,i) = norm(w)
            StagL2Norm(w_u, 0, norm_u);
            CCL2Norm(w_p, 0, norm_p);
            MarkerL2Norm(part_indices, dummy_iter, geom, marker_pos, w_lambda, norm_lambda);
            norm_p      = p_norm_weight*norm_p;
            norm_lambda = p_norm_weight*norm_lambda; // TODO: use p_norm_weight for now
            H[i+1][i]   = sqrt(norm_u*norm_u + norm_p*norm_p + norm_lambda*norm_lambda);


            //___________________________________________________________________
            // V(i+1) = w / H(i+1,i)
            if (H[i+1][i] != 0.) {
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    MultiFab::Copy(V_u[d], w_u[d], 0, i+1, 1, 0);
                    V_u[d].mult(1./H[i+1][i], i+1, 1, 0);
                }

                MultiFab::Copy(V_p, w_p, 0, i+1, 1, 0);
                V_p.mult(1./H[i+1][i], i+1, 1, 0);

                MarkerCopy(part_indices, i+1, V_lambda, w_lambda);
                MarkerMult(part_indices, i+1, 1./H[i+1][i], V_lambda);

            } else { Abort("GMRES.cpp: error in orthogonalization"); }


            //___________________________________________________________________
            // solve least square problem
            LeastSquares(i, H, cs, sn, s);
            norm_resid_est = amrex::Math::abs(s[i+1]);


            //___________________________________________________________________
            // Print verbose output
            if (gmres_verbose >= 2) {
                Print() << "Total iter " << total_iter << ", est. rel. resid. |Pr|/(Pr0, b)= "
                        << norm_resid_est/norm_init_resid << "  "
                        << norm_resid_est/norm_pre_b      << std::endl;
            }

            if (gmres_verbose >= 3) {
                Print() << "Inner. rel. res. (u,v,p,lambda) = "
                        << norm_resid/norm_init_resid  << "  "
                        << norm_u/norm_init_resid      << "  "
                        << norm_p/norm_init_resid      << "  "
                        << norm_lambda/norm_init_resid << std::endl;
            }



            //___________________________________________________________________
            // Inner loop termination condition
            if (total_iter >= gmres_max_iter) {
                break; // exit InnerLoop
            } else if (total_iter >= gmres_min_iter) {
                if ((norm_resid_est <= gmres_rel_tol*amrex::min(norm_pre_b, norm_init_resid))
                    || (norm_resid_est <= gmres_abs_tol)) {
                    break; // exit InnerLoop
                }
            }

        } // end of inner loop



        //_______________________________________________________________________
        // Update the solution
        // first, solve for y
        SolveUTriangular(i_copy-1, H, s, y);

        // then, x = x + dot(V(1:i),y(1:i))
        // UpdateSol(x_u, x_p, V_u, V_p, y, i_copy);
        UpdateSolIBM(part_indices, x_u, x_p, x_lambda, V_u, V_p, V_lambda, y, i_copy);

    } while (true); // end of outer loop (do iter=1,gmres_max_outer)

    // AJN - this is here since I notice epsilon roundoff errors building up
    //       just enough to destroy the asymmetry in time-advancement codes that
    //       ultimately causes lack of convergence in subsequent gmres calls
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        x_u[d].OverrideSync(geom.periodicity());


    // apply scaling factor
    if (scale_factor != 1.) {
        theta_alpha = theta_alpha/scale_factor;
        // the solution we got is scale*x_p

        x_p.mult(1./scale_factor, 0, 1, x_p.nGrow());
        // unscale the rhs

        for (int d=0; d<AMREX_SPACEDIM; ++d)
            b_u[d].mult(1./scale_factor, 0, 1, b_u[d].nGrow());

        // unscale the viscosities
        beta.mult(1./scale_factor, 0, 1, beta.nGrow());
        gamma.mult(1./scale_factor, 0, 1, gamma.nGrow());
        for (int d=0; d<NUM_EDGE; ++d)
            beta_ed[d].mult(1./scale_factor, 0, 1, beta_ed[d].nGrow());
    }

    if (gmres_verbose >= 1) {
        Print() << "Done with GMRES:" << std::endl;
        Print() << "  total ITERs = " << total_iter << std::endl;
        Print() << "  residual/(norm_b,initial) = " << norm_resid/norm_b << "  "
                << norm_resid/norm_init_resid << std::endl;
    }

}



void IBMPrecon(const std::array<MultiFab, AMREX_SPACEDIM> & b_u, const MultiFab & b_p,
               std::array<MultiFab, AMREX_SPACEDIM> & x_u, MultiFab & x_p,
               std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
               const MultiFab & beta, const std::array<MultiFab, NUM_EDGE> & beta_ed,
               const MultiFab & gamma, const Real & theta_alpha,
               const IBParticleContainer & ib_pc,
               const Vector<std::pair<int, int>> & pindex_list,
               std::map<std::pair<int, int>, Vector<RealVect>> & marker_forces,
               const std::map<std::pair<int, int>, Vector<RealVect>> & marker_W,
               const Geometry & geom)
{

    BL_PROFILE_VAR("IBMPrecon()", IBMPrecon);

    BoxArray ba              = b_p.boxArray();
    DistributionMapping dmap = b_p.DistributionMap();

    Real         mean_val_pres;
    Vector<Real> mean_val_umac(AMREX_SPACEDIM);



    // TODO: Assumptions: dx=dy=dz (approx) as well as constant viscosity
    const Real * dx = geom.CellSize();
    const Real c1 = dx[0]/visc_coef;
    const Real c2 = dx[0]*dx[0];
    const Real c3 = 1./dx[0];
    const Real c4 = 1./(visc_coef*dx[0]);




    MultiFab phi     (ba, dmap, 1, 1);
    MultiFab mac_rhs (ba, dmap, 1, 0);
    MultiFab zero_fab(ba, dmap, 1, 0);
    MultiFab x_p_tmp (ba, dmap, 1, 1);

    // set zero_fab_fc to 0
    zero_fab.setVal(0.);

    // build alphainv_fc, one_fab_fc, and zero_fab_fc
    std::array< MultiFab, AMREX_SPACEDIM > alphainv_fc;
    std::array< MultiFab, AMREX_SPACEDIM > one_fab_fc;
    std::array< MultiFab, AMREX_SPACEDIM > zero_fab_fc;
    std::array< MultiFab, AMREX_SPACEDIM > gradp;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        alphainv_fc[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);
         one_fab_fc[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);
        zero_fab_fc[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);
              gradp[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);

        // set alphainv_fc to 1/alpha_fc
        // set one_fab_fc to 1
        // set zero_fab_fc to 0
        alphainv_fc[d].setVal(1.);
        alphainv_fc[d].divide(alpha_fc[d],0,1,0);
         one_fab_fc[d].setVal(1.);
        zero_fab_fc[d].setVal(0.);
    }


    int ib_grow  = 6; // using the 6-point stencil
    int ib_level = 0; // assume single level for now


    // set the initial guess for Phi in the Poisson solve to 0
    // set x_u = 0 as initial guess
    phi.setVal(0.);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        x_u[d].setVal(0.);

    x_p.setVal(0.);


    std::array<MultiFab, AMREX_SPACEDIM> JLS_V;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        JLS_V[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, ib_grow);
        JLS_V[d].setVal(0.);
    }

    MultiFab JLS_P(ba, dmap, 1, ib_grow);
    JLS_P.setVal(0.);

    std::map<std::pair<int, int>, Vector<RealVect>> JLS;
    for (const auto & pindex : pindex_list){
        // initialized to (0..0)
        JLS[pindex].resize(marker_W.at(pindex).size());
    }


    // 1 = projection preconditioner
    // 2 = lower triangular preconditioner
    // 3 = upper triangular preconditioner
    // 4 = block diagonal preconditioner
    // 5 = Uzawa-type approximation (see paper)
    // 6 = upper triangular + viscosity-based BFBt Schur complement (from Georg Stadler)

    // projection preconditioner
    if (amrex::Math::abs(precon_type) == 1) {

        //_______________________________________________________________________
        // Temporary arrays

        std::array<MultiFab, AMREX_SPACEDIM> Ag;
        std::array<MultiFab, AMREX_SPACEDIM> AGphi;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            Ag[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, ib_grow);
            Ag[d].setVal(0.);

            AGphi[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, ib_grow);
            AGphi[d].setVal(0.);
        }



        /************************************************************************
         *                                                                      *
         * Fluid velocity part (also used later)                                *
         * Calculates: x_u = A^{-1}g - GLp^{-1}(DA^{-1}g + h)                   *
         *                                                                      *
         ***********************************************************************/

        //_______________________________________________________________________
        // STEP 1: Solve for an intermediate state, x_u^star, using an implicit
        // viscous solve
        //         x_u^star = A^{-1} b_u

        // x_u^star = A^{-1} b_u ......................................... x_u = A^{-1}g
        StagMGSolver(alpha_fc, beta, beta_ed, gamma, Ag, b_u, theta_alpha, geom);

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            Ag[d].FillBoundary(geom.periodicity());

            MultiFab::Copy(x_u[d], Ag[d], 0, 0, 1, x_u[d].nGrow());
            x_u[d].FillBoundary(geom.periodicity()); // Just in case
        }


        //_______________________________________________________________________
        // STEP 2: Construct RHS for pressure Poisson problem

        // set mac_rhs = D(x_u^star) ................................ mac_rhs = DA^{-1}g
        ComputeDiv(mac_rhs, x_u, 0, 0, 1, geom, 0);

        // add b_p to mac_rhs ................................... mac_rhs = DA^{-1}g + h
        MultiFab::Add(mac_rhs, b_p, 0, 0, 1, 0);

        //_______________________________________________________________________
        // STEP 3: Compute x_u

        // use multigrid to solve for Phi ......................... phi = Lp^{-1}mac_rhs
        // x_u^star is only passed in to get a norm for absolute residual criteria
        MacProj macproj;
        macproj.Define(ba,dmap,geom);
        macproj.solve(alphainv_fc, mac_rhs, phi, geom);

        // x_u = x_u^star - (alpha I)^-1 grad Phi ...... x_u = A^{-1}g - GLp^{-1}mac_rhs
        SubtractWeightedGradP(x_u, alphainv_fc, phi, gradp, geom);


        /************************************************************************
         *                                                                      *
         * Immersed boundary part                                               *
         *                                                                      *
         ***********************************************************************/

        //_______________________________________________________________________
        // Temporary arrays
        std::map<std::pair<int, int>, Vector<RealVect>> JAgW;
        std::map<std::pair<int, int>, Vector<RealVect>> JAGphi;
        std::map<std::pair<int, int>, Vector<RealVect>> JLS_rhs;
        for (const auto & pindex : pindex_list){
            // initialized to (0..0)
            JAgW[pindex].resize(marker_W.at(pindex).size());
            JAGphi[pindex].resize(marker_W.at(pindex).size());
            JLS_rhs[pindex].resize(marker_W.at(pindex).size());
        }



        //_______________________________________________________________________
        // Pure pressure gradient complement: G\phi and A^{-1}G\phi
        std::array<MultiFab, AMREX_SPACEDIM> Gphi;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            Gphi[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, ib_grow);
            Gphi[d].setVal(0.);
        }

        SubtractWeightedGradP(Gphi, alphainv_fc, phi, gradp, geom);

        for (int d=0; d<AMREX_SPACEDIM; ++d)
            Gphi[d].FillBoundary(geom.periodicity());

        // StagMGSolver(alpha_fc, beta, beta_ed, gamma, AGphi, Gphi, theta_alpha, geom);

        // for (int d=0; d<AMREX_SPACEDIM; ++d) {
        //     AGphi[d].mult(-1., 0, 1, ib_grow);
        //     AGphi[d].FillBoundary(geom.periodicity());
        // }


        //_______________________________________________________________________
        // JAgW, and JAGphi preconditioner terms

        // J-interpolated terms: A^{-1}g, A^{-1}G\phi, sourced above
        for (const auto & pindex : pindex_list) {
                  auto & jagw = JAgW.at(pindex);
            const auto & W    = marker_W.at(pindex);

            ib_pc.InterpolateMarkers(ib_level, pindex, jagw, Ag);
            for (auto & elt : jagw) elt = c3*elt;

            for (int i=0; i<jagw.size(); ++i)
                jagw[i] = jagw[i] + W[i]; // ........................ JAgW = JA^{-1}g + W


            auto & jagphi = JAGphi.at(pindex); // ................. JAGphi = JA^{-1}G\phi
            // ib_pc.InterpolateMarkers(ib_level, pindex, jagphi, AGphi);
            ib_pc.InterpolateMarkers(ib_level, pindex, jagphi, Gphi);
            for (auto & elt : jagphi) elt = c3*elt;
        }


        //_______________________________________________________________________
        // JLS preconditioner term

        // RHS term for preconditioner
        for (const auto & pindex : pindex_list){
                  auto & jls    = JLS_rhs.at(pindex);
            const auto & jagw   = JAgW.at(pindex);
            const auto & jagphi = JAGphi.at(pindex);

            for (int i=0; i<jls.size(); ++i)
                jls[i] = jagphi[i] + jagw[i]; // .. JLS_rhs = JA^{-1}G\phi + JA^{-1}g + W
        }

        // Preconditioner guess: L^{-1} ~ A^{-1} => (JL^{-1}S)^{-1} = JAS

        std::array<MultiFab, AMREX_SPACEDIM> spread_rhs;
        std::array<MultiFab, AMREX_SPACEDIM> spread_weights;
        std::array<MultiFab, AMREX_SPACEDIM> AS_rhs;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            spread_rhs[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, ib_grow);
            spread_rhs[d].setVal(0.);

            spread_weights[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, ib_grow);
            spread_weights[d].setVal(0.);

            AS_rhs[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, ib_grow);
            AS_rhs[d].setVal(0.);
        }


        // const Real * dx = geom.CellSize();
        // Real invvol     = 1.;
        // for (int d=0; d<AMREX_SPACEDIM; ++d)
        //     invvol = invvol/dx[d];

        Real c  = 0.326;
        Real D1 = 1.e2; // SJ
        Real D2 = 1.e2; // JS

        for (const auto & pindex : pindex_list) {
            const auto & jls = JLS_rhs.at(pindex);

            // ib_pc.SpreadMarkers(ib_level, pindex, jls, spread_rhs, spread_weights);
            ib_pc.InvInterpolateMarkers(ib_level, pindex, jls, spread_rhs);
        }

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            // spread_rhs[d].mult(c2*invvol, 0, 1, ib_grow);
            spread_rhs[d].mult(c2/D1, 0, 1, ib_grow);
            // spread_rhs[d].mult(1./5, 0, 1, ib_grow);
            spread_rhs[d].FillBoundary(geom.periodicity());
            spread_weights[d].FillBoundary(geom.periodicity());
        }


        // Apply A (Helmhotz) operator
        StagApplyOp(geom, beta, gamma, beta_ed, spread_rhs, AS_rhs, alpha_fc, dx, theta_alpha);
        // for (int d=0; d<AMREX_SPACEDIM; ++d)
        //     MultiFab::Copy(AS_rhs[d], spread_rhs[d], 0, 0, 1, ib_grow);

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            AS_rhs[d].FillBoundary(geom.periodicity());
        }


        // Precon: ....................... JLS = - JAS (JA^{-1}G\phi + JA^{-1}g + W )
        for (const auto & pindex : pindex_list) {
            auto & jls = JLS.at(pindex);

            ib_pc.InterpolateMarkers(ib_level, pindex, jls, AS_rhs);
            // ib_pc.InterpolateMarkers(ib_level, pindex, jls, AS_rhs, spread_weights);
            // ib_pc.InvSpreadMarkers(ib_level, pindex, jls, AS_rhs);

            // for (auto & marker : jls) marker = c3*invvol * marker;
            for (auto & marker : jls) marker = c3/D2 * marker;
        }


        // for (const auto & pindex : pindex_list) {
        //           auto & jls     = JLS.at(pindex);
        //     const auto & jls_rhs = JLS_rhs.at(pindex);

        //     for (int i=0; i<jls.size(); ++i) {
        //         jls[i] = jls_rhs[i];
        //     }
        // }


        // Real spring_constant = 1.;
        // for (const auto & pindex : pindex_list) {
        //     auto & jls = JLS.at(pindex);

        //     ib_pc.InterpolateMarkers(ib_level, pindex, jls, x_u);
        //     // ib_pc.InterpolateMarkers(ib_level, pindex, jls, AS_rhs, spread_weights);

        //     for (auto & marker : jls)
        //         marker = -c3*spring_constant * marker;

        //     const auto & W = marker_W.at(pindex);

        //     for (int i=0; i<jls.size(); ++i)
        //         jls[i] = jls[i] - W[i];
        // }



        //_______________________________________________________________________
        // Compute JLS preconditioner contributions for the velocity and pressure
        // JLS_V = A^{-1}S JLS
        // JLS_P = (\theta\rho_0\Lp^{-1}-\mu_0 1) DA^{-1}S JLS

        std::array<MultiFab, AMREX_SPACEDIM> JLS_V_rhs;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            JLS_V_rhs[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, ib_grow);
            JLS_V_rhs[d].setVal(0.);

            spread_weights[d].setVal(0.);
        }

        MultiFab JLS_P_rhs(ba, dmap, 1, ib_grow);
        JLS_P_rhs.setVal(0.);


        // Velocity Part:

        for (const auto & pindex : pindex_list) {
            const auto & jls = JLS.at(pindex);

            ib_pc.SpreadMarkers(ib_level, pindex, jls, JLS_V_rhs, spread_weights);
        }

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            // Scale the Spreading operator by c2
            JLS_V_rhs[d].mult(c2, 0, 1, ib_grow);

            JLS_V_rhs[d].FillBoundary(geom.periodicity());
            spread_weights[d].FillBoundary(geom.periodicity());
        }

        StagMGSolver(alpha_fc, beta, beta_ed, gamma, JLS_V, JLS_V_rhs, theta_alpha, geom);

        for (int d=0; d<AMREX_SPACEDIM; ++d)
            JLS_V[d].FillBoundary(geom.periodicity());


        // Pressure Part:
        // Note that `JLS_V` (above) is already equal to the A^{-1}S part in `DA^{-1}S`
        // This also means the IBM preconditioner has the property that `div(u) = h`

        // ................................................ JLS_P_rhs = DA^{-1}S JLS
        ComputeDiv(JLS_P_rhs, JLS_V, 0, 0, 1, geom, 0);
        JLS_P_rhs.FillBoundary(geom.periodicity());

        // use multigrid to solve for Phi ............. JLS_P = Lp^{-1} DA^{-1}S JLS

        macproj.Solve(alphainv_fc, JLS_P_rhs, JLS_P, geom);

        // x_u = x_u^star - (alpha I)^-1 grad Phi ...... x_u = A^{-1}g - GLp^{-1}mac_rhs
        SubtractWeightedGradP(JLS_V, alphainv_fc, JLS_P, gradp, geom);

        for (int d=0; d<AMREX_SPACEDIM; ++d)
            JLS_V[d].FillBoundary(geom.periodicity());


        JLS_P.mult(theta_alpha, 0, 1, 0);
        MultiFab::Subtract(JLS_P, JLS_P_rhs, 0, 0, 1, 0);
        JLS_P.FillBoundary(geom.periodicity());




        /************************************************************************
         *                                                                      *
         * Pressure part                                                        *
         * Calculates: x_p = { -(DA^{-1}g + h) + Lp^{-1}(DA^{-1}g + h)          *
         *                   { L_alpha Lp^{-1}(DA^{-1}g + h)                    *
         *                                                                      *
         ***********************************************************************/

        ////////////////////
        // STEP 4: Compute x_p by applying the Schur complement approximation
        ////////////////////

        if (visc_schur_approx == 0) {
            // if precon_type = +1, or theta_alpha=0 then x_p = theta_alpha*Phi - c*beta*(mac_rhs)
            // if precon_type = -1                   then x_p = theta_alpha*Phi - c*beta*L_alpha Phi

            if (precon_type == 1 || theta_alpha == 0) {
                // first set x_p = -mac_rhs ...................... x_p = -(DA^{-1}g + h)
                MultiFab::Copy(x_p, mac_rhs, 0, 0, 1, 0);
                x_p.mult(-1., 0, 1, 0);
            } else {
                // first set x_p = -L_alpha Phi .... x_p = L_alpha Lp^{-1}(DA^{-1}g + h)
                CCApplyNegLap(phi, x_p, alphainv_fc, geom);
            }

            if ( amrex::Math::abs(visc_type) == 1 || amrex::Math::abs(visc_type) == 2) {
                // multiply x_p by beta; x_p = -beta L_alpha Phi
                MultiFab::Multiply(x_p, beta, 0, 0, 1, 0);

                if (amrex::Math::abs(visc_type) == 2) {
                    // multiply by c=2; x_p = -2*beta L_alpha Phi
                    x_p.mult(2., 0, 1, 0);
                }
            } else if (amrex::Math::abs(visc_type) == 3) {

                // multiply x_p by gamma, use mac_rhs a temparary to save x_p
                MultiFab::Copy(mac_rhs, x_p, 0, 0, 1, 0);
                MultiFab::Multiply(mac_rhs, gamma, 0, 0, 1, 0);
                // multiply x_p by beta; x_p = -beta L_alpha Phi
                MultiFab::Multiply(x_p, beta, 0, 0, 1, 0);
                // multiply by c=4/3; x_p = -(4/3) beta L_alpha Phi
                x_p.mult(4./3., 0, 1, 0);
                // x_p = -(4/3) beta L_alpha Phi - gamma L_alpha Phi
                MultiFab::Add(x_p, mac_rhs, 0, 0, 1, 0);
            }

            // multiply Phi by theta_alpha
            phi.mult(theta_alpha, 0, 1, 0);

            // add theta_alpha*Phi to x_p
            MultiFab::Add(x_p, phi, 0, 0, 1, 0);

        } else { Abort("StagApplyOp: visc_schur_approx != 0 not supported"); }

    } else { Abort("StagApplyOp: unsupposed precon_type"); }



    /****************************************************************************
     *                                                                          *
     * Add Immersed boundary part                                               *
     *                                                                          *
     ***************************************************************************/

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Add(x_u[d], JLS_V[d], 0, 0, 1, 0);
        x_u[d].FillBoundary(geom.periodicity());
    }

    MultiFab::Add(x_p, JLS_P, 0, 0, 1, 0);
    x_p.FillBoundary(geom.periodicity());


    for (const auto & pindex : pindex_list) {
        auto & force = marker_forces.at(pindex);
        for (auto & elt : force)
            elt = RealVect{0,0,0};
    }


    for (const auto & pindex : pindex_list) {
              auto & force = marker_forces.at(pindex);
        const auto & jls   = JLS.at(pindex);

        for (int i=0; i<jls.size(); ++i)
            force[i] = force[i] + jls[i];
    }



    ////////////////////
    // STEP 5: Handle null-space issues in MG solvers
    ////////////////////

    // subtract off mean value: Single level only! No need for ghost cells
    SumStag(geom, x_u, mean_val_umac, true);
    SumCC(x_p, 0, mean_val_pres, true);

    // The pressure Poisson problem is always singular:
    x_p.plus(-mean_val_pres, 0, 1, 0);

    // The velocity problem is also singular under these cases
    if (theta_alpha == 0.) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            if (geom.isPeriodic(d)) {
                x_u[d].plus(-mean_val_umac[d], 0, 1, 0);
            }
        }
    }
}



void ApplyIBM(      std::array<MultiFab, AMREX_SPACEDIM>            & b_u,
                    std::map<std::pair<int, int>, Vector<RealVect>> & b_lambda,
              const std::array<MultiFab, AMREX_SPACEDIM>            & x_u,
              const IBParticleContainer                             & ib_pc,
              const Vector<std::pair<int, int>>                     & part_indices,
              const std::map<std::pair<int, int>, Vector<RealVect>> & x_lambda,
              int ib_grow, int ibpc_lev, const Geometry & geom ) {



    // TODO: Assumptions: dx=dy=dz (approx) as well as constant viscosity
    const Real * dx = geom.CellSize();
    const Real c1 = dx[0]/visc_coef;
    const Real c2 = dx[0]*dx[0];
    const Real c3 = 1./dx[0];
    const Real c4 = 1./(visc_coef*dx[0]);



    // Temporary containers:
    std::array<MultiFab, AMREX_SPACEDIM> SLambda;
    std::array<MultiFab, AMREX_SPACEDIM> spread_weights;
    std::array<MultiFab, AMREX_SPACEDIM> x_u_buffer;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        SLambda[d].define(x_u[d].boxArray(), x_u[d].DistributionMap(), 1, ib_grow);
        SLambda[d].setVal(0.);

        x_u_buffer[d].define(x_u[d].boxArray(), x_u[d].DistributionMap(), 1, ib_grow);
        x_u_buffer[d].setVal(0.);

        spread_weights[d].define(x_u[d].boxArray(), x_u[d].DistributionMap(), 1, ib_grow);
        spread_weights[d].setVal(0.);
    }


    // ............................................................ SLambda = S lambda
    for (const auto & pid : part_indices) {
        const auto & lambda = x_lambda.at(pid);

        ib_pc.SpreadMarkers(ibpc_lev, pid, lambda, SLambda, spread_weights);
    }

    for (int d=0; d<AMREX_SPACEDIM; ++d) {

        // Scale the S lambda operator by c2
        SLambda[d].mult(c2, 0, 1, ib_grow);

        SLambda[d].FillBoundary(geom.periodicity());
        spread_weights[d].FillBoundary(geom.periodicity());
    }

    Real norm_sl = 0;
    Real norm_bu = 0;
    StagL2Norm(SLambda, 0, norm_sl);
    StagL2Norm(b_u, 0, norm_bu);
    Print() << "norm SLambda = " << norm_sl <<  " norm bu = " << norm_bu << std::endl;


    // ........................................................ v = Av - Gp - S lambda
    // ...... (computed elsewhere, before this function call) ------^^^^^^^

    for (int d=0; d<AMREX_SPACEDIM; ++d)
        MultiFab::Subtract(b_u[d], SLambda[d], 0, 0, 1, 0);

    // Buffer x_u so that it has enough ghost cells (to cover the kernel)
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(x_u_buffer[d], x_u[d], 0, 0, 1, 0);
        x_u_buffer[d].FillBoundary(geom.periodicity());
    }

    // .................................................................. lambda = -Jv
    for (const auto & pid : part_indices) {
        auto & bl = b_lambda.at(pid);

        ib_pc.InterpolateMarkers(ibpc_lev, pid, bl, x_u_buffer);
        for (auto & elt : bl) elt = -c3*elt;
    }
}



void MarkerAdd(Vector<RealVect> & a, const Vector<RealVect> & b)  {

    for (int i=0; i<a.size(); ++i)
        a[i] = a[i] + b[i];
}



void MarkerAdd(const Vector<std::pair<int, int>> & part_indices,
                     std::map<std::pair<int, int>, Vector<RealVect>> & a,
               const std::map<std::pair<int, int>, Vector<RealVect>> & b) {

    for (const auto & pid : part_indices) {
              auto & a_markers = a.at(pid);
        const auto & b_markers = b.at(pid);

        MarkerAdd(a_markers, b_markers);
    }
}



void MarkerAdd(const Vector<std::pair<int, int>> & part_indices, int comp,
                     std::map<std::pair<int, int>,        Vector<RealVect>>  & a,
               const std::map<std::pair<int, int>, Vector<Vector<RealVect>>> & b) {

    for (const auto & pid : part_indices) {
              auto & a_markers = a.at(pid);
        const auto & b_markers = b.at(pid);

        MarkerAdd(a_markers, b_markers[comp]);
    }
}



void MarkerSub(Vector<RealVect> & a, const Vector<RealVect> & b) {

    for (int i=0; i<a.size(); ++i)
        a[i] = a[i] - b[i];
}



void MarkerSub(const Vector<std::pair<int, int>> & part_indices,
                     std::map<std::pair<int, int>, Vector<RealVect>> & a,
               const std::map<std::pair<int, int>, Vector<RealVect>> & b) {

    for (const auto & pid : part_indices) {
              auto & a_markers = a.at(pid);
        const auto & b_markers = b.at(pid);

        MarkerSub(a_markers, b_markers);
    }
}



void MarkerInvSub(Vector<RealVect> & a, const Vector<RealVect> & b) {

    for (int i=0; i<a.size(); ++i)
        a[i] = b[i] - a[i];
}



void MarkerInvSub(const Vector<std::pair<int, int>> & part_indices,
                        std::map<std::pair<int, int>, Vector<RealVect>> & a,
                  const std::map<std::pair<int, int>, Vector<RealVect>> & b) {

    for (const auto & pid : part_indices) {
              auto & a_markers = a.at(pid);
        const auto & b_markers = b.at(pid);

        MarkerInvSub(a_markers, b_markers);
    }
}



void MarkerMult(Real factor, Vector<RealVect> & a) {

    for (auto & elt : a)
        elt = elt*factor;
}



void MarkerMult(const Vector<std::pair<int, int>> & part_indices, Real factor,
                std::map<std::pair<int, int>, Vector<RealVect>> & a ) {

    for (const auto & pid : part_indices) {
        auto & a_markers = a.at(pid);

        MarkerMult(factor, a_markers);
    }
}



void MarkerMult(const Vector<std::pair<int, int>> & part_indices, int comp, Real factor,
                std::map<std::pair<int, int>, Vector<Vector<RealVect>>> & a ) {

    for (const auto & pid : part_indices) {
        auto & a_markers = a.at(pid);

        MarkerMult(factor, a_markers[comp]);
    }
}



void MarkerInnerProd(const Box & bx, const Geometry & geom,
                     const Vector<RealVect> & marker_pos,
                     const Vector<RealVect> & a, const Vector<RealVect> & b,
                     Real & v) {

    v = 0;

    for (int i=0; i<a.size(); ++i) {
        IntVect cc_index = geom.CellIndex(marker_pos[i].dataPtr());

        if (bx.contains(cc_index)) {
            Real vi = a[i].dotProduct(b[i]);
            v       = v + vi;
        }
    }
}



void MarkerInnerProd(const Vector<std::pair<int, int>> & part_indices,
                     const MultiFab & cc_iter, const Geometry & geom,
                     const std::map<std::pair<int, int>, Vector<RealVect>> & marker_pos,
                     const std::map<std::pair<int, int>, Vector<RealVect>> & a,
                     const std::map<std::pair<int, int>, Vector<RealVect>> & b,
                     Real & v) {

    v = 0.;

    for (MFIter mfi(cc_iter); mfi.isValid(); ++ mfi) {
        for (const auto & pid : part_indices) {

            const Box & bx = mfi.tilebox();

            Real l2_norm = 0.;
            MarkerInnerProd(bx, geom, marker_pos.at(pid),
                            a.at(pid), b.at(pid),
                            l2_norm);
            v = v + l2_norm;
        }
    }

    ParallelDescriptor::ReduceRealSum(v);
}



void MarkerInnerProd(const Vector<std::pair<int, int>> & part_indices, int comp,
                     const MultiFab & cc_iter, const Geometry & geom,
                     const std::map<std::pair<int, int>, Vector<RealVect>> & marker_pos,
                     const std::map<std::pair<int, int>,        Vector<RealVect>>  & a,
                     const std::map<std::pair<int, int>, Vector<Vector<RealVect>>> & b,
                     Real & v) {

    v = 0.;

    for (MFIter mfi(cc_iter); mfi.isValid(); ++ mfi) {
        for (const auto & pid : part_indices) {

            const Box & bx = mfi.tilebox();

            Real l2_norm = 0.;
            MarkerInnerProd(bx, geom, marker_pos.at(pid),
                            a.at(pid), b.at(pid)[comp],
                            l2_norm);
            v = v + l2_norm;
        }
    }

    ParallelDescriptor::ReduceRealSum(v);
}




void MarkerL2Norm(const Box & bx, const Geometry & geom, const Vector<RealVect> & marker_pos,
                  const Vector<RealVect> & markers, Real & norm_l2) {

    norm_l2 = 0.;
    MarkerInnerProd(bx, geom, marker_pos, markers, markers, norm_l2);
    norm_l2 = sqrt(norm_l2);

}



void MarkerL2Norm(const Vector<std::pair<int, int>> & part_indices,
                  const MultiFab & cc_iter, const Geometry & geom,
                  const std::map<std::pair<int, int>, Vector<RealVect>> & marker_pos,
                  const std::map<std::pair<int, int>, Vector<RealVect>> & b, Real & v) {

    v = 0.;

    for (MFIter mfi(cc_iter); mfi.isValid(); ++ mfi) {
        for (const auto & pid : part_indices) {

            const Box & bx = mfi.tilebox();

            Real l2_norm = 0.;
            MarkerL2Norm(bx, geom, marker_pos.at(pid), b.at(pid), l2_norm);
            v = v + l2_norm;
        }
    }

    ParallelDescriptor::ReduceRealSum(v);
}



void MarkerCopy(Vector<RealVect> & a, const Vector<RealVect> & b) {

    for (int i=0; i<a.size(); ++i)
        a[i] = b[i];
}



void MarkerCopy(const Vector<std::pair<int, int>> & part_indices,
                      std::map<std::pair<int, int>, Vector<RealVect>> & a,
                const std::map<std::pair<int, int>, Vector<RealVect>> & b) {

    for (const auto & pid : part_indices) {
              auto & a_markers = a.at(pid);
        const auto & b_markers = b.at(pid);

        MarkerCopy(a_markers, b_markers);
    }
}



void MarkerCopy(const Vector<std::pair<int, int>> & part_indices, int comp,
                      std::map<std::pair<int, int>, Vector<Vector<RealVect>>> & a,
                const std::map<std::pair<int, int>,        Vector<RealVect>>  & b) {

    for (const auto & pid : part_indices) {
              auto & a_markers = a.at(pid);
        const auto & b_markers = b.at(pid);

        MarkerCopy(a_markers[comp], b_markers);
    }
}



void MarkerCopy(const Vector<std::pair<int, int>> & part_indices, int comp,
                      std::map<std::pair<int, int>,        Vector<RealVect>>  & a,
                const std::map<std::pair<int, int>, Vector<Vector<RealVect>>> & b) {

    for (const auto & pid : part_indices) {
              auto & a_markers = a.at(pid);
        const auto & b_markers = b.at(pid);

        MarkerCopy(a_markers, b_markers[comp]);
    }
}



void UpdateSolIBM(const Vector<std::pair<int, int>>                       & part_indices,
                  std::array<MultiFab, AMREX_SPACEDIM>                    & x_u,
                  MultiFab                                                & x_p,
                  std::map<std::pair<int, int>,        Vector<RealVect>>  & x_lambda,
                  std::array<MultiFab, AMREX_SPACEDIM>                    & V_u,
                  MultiFab                                                & V_p,
                  std::map<std::pair<int, int>, Vector<Vector<RealVect>>> & V_lambda,
                  const Vector<Real>                                      & y,
                  int i ) {

    // set V(i) = V(i)*y(i)
    // set x = x + V(i)

    for (int iter=0; iter<=i; ++iter) {

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            V_u[d].mult(y[iter], iter, 1, 0);
            MultiFab::Add(x_u[d], V_u[d], iter, 0, 1, 0);
        }

        V_p.mult(y[iter], iter, 1, 0);
        MultiFab::Add(x_p, V_p, iter, 0, 1, 0);

        MarkerMult(part_indices, iter, y[iter], V_lambda);
        MarkerAdd(part_indices, iter, x_lambda, V_lambda);
    }
}
