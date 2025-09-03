#include "GMRES.H"

GMRES::GMRES (const BoxArray& ba_in,
              const DistributionMapping& dmap_in,
              const Geometry& geom_in) {

    BL_PROFILE_VAR("GMRES::GMRES()", GMRES);

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        r_u[d]        .define(convert(ba_in, nodal_flag_dir[d]), dmap_in, 1,                 1);
        w_u[d]        .define(convert(ba_in, nodal_flag_dir[d]), dmap_in, 1,                 0);
        tmp_u[d]      .define(convert(ba_in, nodal_flag_dir[d]), dmap_in, 1,                 0);
        scr_u[d]      .define(convert(ba_in, nodal_flag_dir[d]), dmap_in, 1,                 0);
        V_u[d]        .define(convert(ba_in, nodal_flag_dir[d]), dmap_in, gmres_max_inner+1, 0);
        alphainv_fc[d].define(convert(ba_in, nodal_flag_dir[d]), dmap_in, 1, 0);
    }

    r_p.define  (ba_in, dmap_in,                  1, 1);
    w_p.define  (ba_in, dmap_in,                  1, 0);
    tmp_p.define(ba_in, dmap_in,                  1, 0);
    scr_p.define(ba_in, dmap_in,                  1, 0);
    V_p.define  (ba_in, dmap_in,gmres_max_inner + 1, 0); // Krylov vectors

    StagSolver.Define(ba_in,dmap_in,geom_in);
    Pcon.Define(ba_in,dmap_in,geom_in);
}


void GMRES::Solve (std::array<MultiFab, AMREX_SPACEDIM> & b_u, MultiFab & b_p,
                   std::array<MultiFab, AMREX_SPACEDIM> & x_u, MultiFab & x_p,
                   std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                   MultiFab & beta, std::array<MultiFab, NUM_EDGE> & beta_ed,
                   MultiFab & gamma,
                   Real theta_alpha,
                   const Geometry & geom,
                   Real & norm_pre_rhs)
{

    BL_PROFILE_VAR("GMRES::Solve()", GMRES_Solve);

    if (gmres_verbose >= 1) {
        Print() << "Begin call to GMRES" << std::endl;
    }

    Vector<Real> cs(gmres_max_inner);
    Vector<Real> sn(gmres_max_inner);
    Vector<Real>  y(gmres_max_inner);
    Vector<Real>  s(gmres_max_inner+1);

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
    Real norm_resid_est;

    Real norm_u; // temporary norms used to build full-state norm
    Real norm_p; // temporary norms used to build full-state norm

    Vector<Real> inner_prod_vel(AMREX_SPACEDIM);
    Real inner_prod_pres;

    //////////////////////////////////////
    // account for inhomogeneous boundary conditions, e.g., moving walls
    // use r_u, tmp_u, r_p, tmp_p as temporary storage

    bool inhomogeneous_fix = false;
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (wallspeed_x_lo[i] != 0.) inhomogeneous_fix = true;
        if (wallspeed_x_hi[i] != 0.) inhomogeneous_fix = true;
        if (wallspeed_y_lo[i] != 0.) inhomogeneous_fix = true;
        if (wallspeed_y_hi[i] != 0.) inhomogeneous_fix = true;
#if (AMREX_SPACEDIM == 3)
        if (wallspeed_z_lo[i] != 0.) inhomogeneous_fix = true;
        if (wallspeed_z_hi[i] != 0.) inhomogeneous_fix = true;
#endif
    }

    if (inhomogeneous_fix) {

        // fill MultiFabs for velocity and pressure to zero
        // then fill ghost cells
        // then apply the operator
        // subtract the result from the rhs

        int is_inhomogeneous = 1;

        r_p.setVal(0.);
        MultiFabPhysBC(r_p, geom, 0, 1, PRES_BC_COMP);

        for (int i=0; i<AMREX_SPACEDIM; ++i ) {
            r_u[i].setVal(0.);
            MultiFabPhysBCMacVel(r_u[i], geom, i, is_inhomogeneous);
        }

        ApplyMatrix(tmp_u, tmp_p, r_u, r_p, alpha_fc, beta, beta_ed, gamma, theta_alpha, geom,
                    is_inhomogeneous);

        MultiFab::Subtract(b_p, tmp_p, 0, 0, 1, 0);
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            MultiFab::Subtract(b_u[i], tmp_u[i], 0, 0, 1, 0);
        }

    }
    //////////////////////////////////////

    //////////////////////////////////////

    /****************************************************************************
     *                                                                          *
     * Preflight work: apply scaling and compute perconditioned norms_b         *
     *                                                                          *
     ***************************************************************************/

    // set alphainv_fc to 1/alpha_fc
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        alphainv_fc[d].setVal(1.);
        alphainv_fc[d].divide(alpha_fc[d],0,1,0);
    }

    // apply scaling factor
    if (scale_factor != 1.) {
        theta_alpha = theta_alpha*scale_factor;

        // we will solve for scale*x_p so we need to scale the initial guess
        x_p.mult(scale_factor, 0, 1, x_p.nGrow());

        // scale the rhs:
        for (int d=0; d<AMREX_SPACEDIM; ++d)
            b_u[d].mult(scale_factor,0,1,b_u[d].nGrow());

        // scale the viscosities:
        beta.mult(scale_factor, 0, 1, beta.nGrow());
        gamma.mult(scale_factor, 0, 1, gamma.nGrow());
        for (int d=0; d<NUM_EDGE; ++d)
            beta_ed[d].mult(scale_factor, 0, 1, beta_ed[d].nGrow());
    }


    // First application of preconditioner
    Pcon.Apply(b_u, b_p, tmp_u, tmp_p, alpha_fc, alphainv_fc,
               beta, beta_ed, gamma, theta_alpha, geom, StagSolver);


    // preconditioned norm_b: norm_pre_b
    StagL2Norm(tmp_u, 0, scr_u, norm_u);
    CCL2Norm(tmp_p, 0, scr_p, norm_p);
    norm_p       = p_norm_weight*norm_p;
    norm_pre_b   = sqrt(norm_u*norm_u + norm_p*norm_p);
    norm_pre_rhs = norm_pre_b;


    // calculate the l2 norm of rhs
    StagL2Norm(b_u, 0, scr_u, norm_u);
    CCL2Norm(b_p, 0, scr_p, norm_p);
    norm_p = p_norm_weight*norm_p;
    norm_b = sqrt(norm_u*norm_u + norm_p*norm_p);


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
        ApplyMatrix(tmp_u, tmp_p, x_u, x_p, alpha_fc, beta, beta_ed, gamma, theta_alpha, geom);

        // tmp = b - Ax
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Subtract(tmp_u[d], b_u[d], 0, 0, 1, 0);
            tmp_u[d].mult(-1., 0, 1, 0);
        }
        MultiFab::Subtract(tmp_p, b_p, 0, 0, 1, 0);
        tmp_p.mult(-1., 0, 1, 0);


        //_______________________________________________________________________
        // un-preconditioned residuals
        StagL2Norm(tmp_u, 0, scr_u, norm_u_noprecon);
        CCL2Norm(tmp_p, 0, scr_p, norm_p_noprecon);
        norm_p_noprecon   = p_norm_weight*norm_p_noprecon;
        norm_resid_Stokes = sqrt(norm_u_noprecon*norm_u_noprecon + norm_p_noprecon*norm_p_noprecon);

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
            Print() << "un-Precond. rel. resid. (u,v,p) = " << norm_resid_Stokes/norm_init_Stokes
                    << "  " << norm_u_noprecon/norm_init_Stokes
                    << "  " << norm_p_noprecon/norm_init_Stokes << std::endl;
        }


        //_______________________________________________________________________
        // solve for r = M^{-1} tmp
        // We should not be counting these toward the number of mg cycles performed
        Pcon.Apply(tmp_u, tmp_p, r_u, r_p, alpha_fc, alphainv_fc,
                   beta, beta_ed, gamma, theta_alpha, geom, StagSolver);


        // resid = sqrt(dot_product(r, r))
        StagL2Norm(r_u, 0, scr_u, norm_u);
        CCL2Norm(r_p, 0, scr_p, norm_p);
        norm_p     = p_norm_weight*norm_p;
        norm_resid = sqrt(norm_u*norm_u + norm_p*norm_p);


        // If first iteration, save the initial preconditioned residual
        if (outer_iter==0) {
            norm_init_resid = norm_resid;
            norm_resid_est = norm_resid;
        }


        //_______________________________________________________________________
        // Print verbose output
        if (gmres_verbose >= 3) {
            Print() << "Precond. rel. res. (u,v,p) = " << norm_resid/norm_init_resid << "  "
                    << norm_u/norm_init_resid << "  " << norm_p/norm_init_resid << std::endl;
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

            ApplyMatrix(tmp_u, tmp_p, r_u, r_p, alpha_fc, beta, beta_ed, gamma, theta_alpha, geom);


            //___________________________________________________________________
            // w = M^{-1} A*V(i)
            Pcon.Apply(tmp_u, tmp_p, w_u, w_p, alpha_fc, alphainv_fc,
                       beta, beta_ed, gamma, theta_alpha, geom, StagSolver);


            //___________________________________________________________________
            // Form Hessenberg matrix H
            for (int k=0; k<=i; ++k) {
                // H(k,i) = dot_product(w, V(k))
                //        = dot_product(w_u, V_u(k))+dot_product(w_p, V_p(k))
                StagInnerProd(w_u, 0, V_u, k, scr_u, inner_prod_vel);
                CCInnerProd(w_p, 0, V_p, k, scr_p, inner_prod_pres);
                H[k][i] = std::accumulate(inner_prod_vel.begin(), inner_prod_vel.end(), 0.)
                          + pow(p_norm_weight, 2.0)*inner_prod_pres;


                // w = w - H(k,i) * V(k)
                // use tmp_u and tmp_p as temporaries to hold kth component of V(k)
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    MultiFab::Copy(tmp_u[d], V_u[d], k, 0, 1, 0);
                    tmp_u[d].mult(H[k][i], 0, 1, 0);
                    MultiFab::Subtract(w_u[d], tmp_u[d], 0, 0, 1, 0);
                }
                MultiFab::Copy(tmp_p, V_p, k, 0, 1, 0);
                tmp_p.mult(H[k][i], 0, 1, 0);
                MultiFab::Subtract(w_p,tmp_p, 0, 0, 1, 0);
            }

            // H(i+1,i) = norm(w)
            StagL2Norm(w_u, 0, scr_u, norm_u);
            CCL2Norm(w_p, 0, scr_p, norm_p);
            norm_p    = p_norm_weight*norm_p;
            H[i+1][i] = sqrt(norm_u*norm_u + norm_p*norm_p);


            //___________________________________________________________________
            // V(i+1) = w / H(i+1,i)
            if (H[i+1][i] != 0.) {
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    MultiFab::Copy(V_u[d], w_u[d], 0, i + 1, 1, 0);
                    V_u[d].mult(1./H[i+1][i], i + 1, 1, 0);
                }
                MultiFab::Copy(V_p, w_p, 0, i + 1, 1, 0);
                V_p.mult(1./H[i+1][i], i + 1, 1, 0);
            } else {
                Abort("GMRES.cpp: error in orthogonalization");
            }


            //___________________________________________________________________
            // solve least square problem
            LeastSquares(i, H, cs, sn, s);
            norm_resid_est = amrex::Math::abs(s[i+1]);


            //___________________________________________________________________
            // Print verbose output
            if (gmres_verbose >= 2) {
                Print() << "Total iter " << total_iter << ",  est. rel. resid. |Pr|/(Pr0,b)= "
                        << norm_resid_est/norm_init_resid << "  "
                        << norm_resid_est/norm_pre_b << std::endl;
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
        UpdateSol(x_u,x_p,V_u,V_p,y,i_copy);

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
            b_u[d].mult(1./scale_factor,0,1,b_u[d].nGrow());

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

void UpdateSol(std::array<MultiFab, AMREX_SPACEDIM>& x_u,
               MultiFab& x_p,
               std::array<MultiFab, AMREX_SPACEDIM>& V_u,
               MultiFab& V_p,
               Vector<Real>& y,
               int i)
{
    // set V(i) = V(i)*y(i)
    // set x = x + V(i)
    for (int iter=0; iter<=i; ++iter) {
        V_p.mult(y[iter],iter,1,0);
        MultiFab::Add(x_p,V_p,iter,0,1,0);
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            V_u[d].mult(y[iter],iter,1,0);
            MultiFab::Add(x_u[d],V_u[d],iter,0,1,0);
        }
    }
}

void LeastSquares(int i,
                  Vector<Vector<Real>>& H,
                  Vector<Real>& cs,
                  Vector<Real>& sn,
                  Vector<Real>& s)
{
    Real temp;

    // apply Givens rotation
    for (int k=0; k<=i-1; ++k) {
        temp      =  cs[k]*H[k][i] + sn[k]*H[k+1][i];
        H[k+1][i] = -sn[k]*H[k][i] + cs[k]*H[k+1][i];
        H[k][i] = temp;
    }

    // form i-th rotation matrix
    RotMat(H[i][i], H[i+1][i], cs[i], sn[i]);

    // approximate residual norm
    temp = cs[i]*s[i];
    s[i+1] = -sn[i]*s[i];
    s[i] = temp;
    H[i][i] = cs[i]*H[i][i] + sn[i]*H[i+1][i];
    H[i+1][i] = 0.;
}

void RotMat(Real a, Real b,
            Real& cs, Real& sn)
{
    Real temp;

    if (b == 0.) {
        cs = 1.;
        sn = 0.;
    }
    else if (amrex::Math::abs(b) > amrex::Math::abs(a)) {
        temp = a/b;
        sn = 1./sqrt(1.+temp*temp);
        cs = temp*sn;
    }
    else {
        temp = b/a;
        cs = 1./sqrt(1.+temp*temp);
        sn = temp*cs;
    }
}

void SolveUTriangular(int k, Vector<Vector<Real>>& H, Vector<Real>& s, Vector<Real>& y)
{
    Real dot;

    y[k+1] = s[k+1]/H[k+1][k+1];
    for (int i=k; i>=0; --i) {

        dot = 0.;
        for (int j=i+1; j<= k+1; ++j) {
            dot += H[i][j]*y[j];
        }

        y[i] = (s[i] - dot) / H[i][i];
    }


}
