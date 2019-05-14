#include <AMReX_VisMF.H>

#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"


#include <IBParticleContainer.H>


using namespace common;
using namespace gmres;


void IBGMRES(std::array<MultiFab, AMREX_SPACEDIM> & b_u, const MultiFab & b_p,
             std::array<MultiFab, AMREX_SPACEDIM> & x_u, MultiFab & x_p,
             const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
             MultiFab & beta, std::array<MultiFab, NUM_EDGE> & beta_ed,
             MultiFab & gamma, Real theta_alpha,
             IBParticleContainer & ib_pc,
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
    Real norm_resid_est;

    Real norm_u; // temporary norms used to build full-state norm
    Real norm_p; // temporary norms used to build full-state norm

    Vector<Real> inner_prod_vel(AMREX_SPACEDIM);
    Real inner_prod_pres;

    BoxArray ba = b_p.boxArray();
    DistributionMapping dmap = b_p.DistributionMap();


    // # of ghost cells must match x_u so higher-order stencils can work
    std::array< MultiFab, AMREX_SPACEDIM > r_u;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        r_u[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, x_u[d].nGrow());

    std::array< MultiFab, AMREX_SPACEDIM > w_u;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        w_u[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);

    std::array< MultiFab, AMREX_SPACEDIM > tmp_u;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        tmp_u[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);

    std::array< MultiFab, AMREX_SPACEDIM > V_u;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        V_u[d].define(convert(ba, nodal_flag_dir[d]), dmap, gmres_max_inner + 1, 0);


    // # of ghost cells must match x_p so higher-order stencils can work
    MultiFab r_p  (ba, dmap,                  1, x_p.nGrow());
    MultiFab w_p  (ba, dmap,                  1, 0);
    MultiFab tmp_p(ba, dmap,                  1, 0);
    MultiFab V_p  (ba, dmap,gmres_max_inner + 1, 0); // Krylov vectors


    //___________________________________________________________________________
    // Get all the immersed-boudary particle indices (used to iterate below)
    Vector<IBP_info> ibp_info;
    // Vector<RealVect> marker_positions;
    // TODO: make this a referenced argument
    std::map<std::pair<int, int>, Vector<RealVect>> marker_forces;
    std::map<std::pair<int, int>, Vector<RealVect>> marker_W;



    // TODO: assuming only 1 level for now
    int ibpc_lev = 0;

    MultiFab dummy(ib_pc.ParticleBoxArray(ibpc_lev),
                   ib_pc.ParticleDistributionMap(ibpc_lev), 1, 1);
    for (MFIter mfi(dummy, ib_pc.tile_size); mfi.isValid(); ++mfi){
        IBParticleContainer::PairIndex index(mfi.index(), mfi.LocalTileIndex());
        ib_pc.IBParticleInfo(ibp_info, ibpc_lev, index, true);
    }

    Vector<std::pair<int, int>> part_indices(ibp_info.size());
    for (int i=0; i<ibp_info.size(); ++i) {
        part_indices[i] = ibp_info[i].asPairIndex();

        // Pre-allocate force arrays
        const Vector<RealVect> marker_positions = ib_pc.MarkerPositions(0, part_indices[i]);
        // ... initialized to (0..0)
        marker_forces[part_indices[i]].resize(marker_positions.size());
        // for (auto & elt : marker_forces[part_indices[i]])
        //     elt = RealVect{1, 1, 1};
        marker_W[part_indices[i]].resize(marker_positions.size());
    }


    Print() << "Found " << part_indices.size() << " many IB particles in rank 0:"
            << std::endl;
    for (const auto & pid : part_indices)
        Print() << "[" << pid.first << ", " << pid.second << "]";
    Print() << std::endl << std::endl;


    // We don't need the marker positions here at all :P
    // for (const auto & pindex : part_indices) {
    //     const Vector<RealVect> pmarkers = ib_pc.MarkerPositions(0, pindex);
    //     for (const auto & marker : pmarkers)
    //         marker_positions.push_back(marker);
    // }


    int ib_grow =  ib_pc.get_nghost() + 6;
    //  using the 6-point stencil ----- ^

    // Buffer velocity to allow for enough ghost cells
    std::array<MultiFab, AMREX_SPACEDIM> u_buffer;
    std::array<MultiFab, AMREX_SPACEDIM> spread_f;
    std::array<MultiFab, AMREX_SPACEDIM> u_precon_f;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        u_buffer[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, ib_grow);
        spread_f[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, ib_grow);
        u_precon_f[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, ib_grow);

        MultiFab::Copy(u_buffer[d], x_u[d], 0, 0, 1, 0);
        u_buffer[d].FillBoundary(geom.periodicity());
    }
    MultiFab p_buffer(ba, dmap, 1, ib_grow);
    MultiFab::Copy(p_buffer, x_p, 0, 0, 1, 0);
    p_buffer.FillBoundary(geom.periodicity());



    /****************************************************************************
     *                                                                          *
     * Preflight work: apply scaling and compute perconditioned norms_b         *
     *                                                                          *
     ***************************************************************************/


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

    // 1. Fluid Precon
    ApplyPrecon(b_u, b_p, tmp_u, tmp_p, alpha_fc, beta, beta_ed, gamma, theta_alpha, geom);

    // 2. IB Precon

    // 2.a Lambda = J*u
    for (const auto & pid : part_indices)
        ib_pc.InterpolateMarkers(0, pid, marker_forces[pid], u_buffer);

    // 2.b Lambda = J*u - J A^{-1} G p
    std::array<MultiFab, AMREX_SPACEDIM> Gp;
    std::array<MultiFab, AMREX_SPACEDIM> AGp;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        Gp[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, ib_grow);
        AGp[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, ib_grow);
    }

    ComputeGrad(p_buffer, Gp, 0, 0, 1, geom);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        Gp[d].FillBoundary(geom.periodicity());

    StagMGSolver(alpha_fc, beta, beta_ed, gamma, AGp, Gp, theta_alpha, geom);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        AGp[d].FillBoundary(geom.periodicity());

    for (const auto & pid : part_indices) {
        Vector<RealVect> mf_tmp(marker_forces[pid].size());
        ib_pc.InterpolateMarkers(0, pid, mf_tmp, AGp);
        for (int i=0; i<marker_forces[pid].size(); ++i) 
            marker_forces[pid][i] = marker_forces[pid][i] + mf_tmp[i];
    }

    // 2.c spread_f = S*Lambda
    for (const auto & pid : part_indices) {
        ib_pc.SpreadMarkers(0, pid, marker_forces[pid], spread_f);
        for (int d=0; d<AMREX_SPACEDIM; ++d)
            spread_f[d].FillBoundary(geom.periodicity());
    }

    // 2.d u_precon_f = A^{-1} spread_f = A^{-1} S Lambda
    StagMGSolver(alpha_fc, beta, beta_ed, gamma, u_precon_f, spread_f, theta_alpha, geom);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        u_precon_f[d].FillBoundary(geom.periodicity());

    for (int d=0; d<AMREX_SPACEDIM; ++d)
        VisMF::Write(u_precon_f[d], "u_precon_f_"+std::to_string(d));




    // preconditioned norm_b: norm_pre_b
    StagL2Norm(tmp_u, 0, norm_u);
    CCL2Norm(tmp_p, 0, norm_p);
    norm_p       = p_norm_weight*norm_p;
    norm_pre_b   = sqrt(norm_u*norm_u + norm_p*norm_p);
    norm_pre_rhs = norm_pre_b;


    // calculate the l2 norm of rhs
    StagL2Norm(b_u, 0, norm_u);
    CCL2Norm(b_p, 0, norm_p);
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
        StagL2Norm(tmp_u, 0, norm_u_noprecon);
        CCL2Norm(tmp_p, 0, norm_p_noprecon);
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
        ApplyPrecon(tmp_u, tmp_p, r_u, r_p, alpha_fc, beta, beta_ed, gamma, theta_alpha, geom);


        // resid = sqrt(dot_product(r, r))
        StagL2Norm(r_u, 0, norm_u);
        CCL2Norm(r_p, 0, norm_p);
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
            if(norm_resid <= gmres_rel_tol*std::min(norm_pre_b, norm_init_resid)) {
                if (gmres_verbose >= 2) {
                    Print() << "GMRES converged: Outer = " << outer_iter << ",  Inner = " << i
                            << " Total=" << total_iter << std::endl;
                }

                if (norm_resid_Stokes >= 10*gmres_rel_tol*std::min(norm_b, norm_init_Stokes)) {
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
            ApplyPrecon(tmp_u, tmp_p, w_u, w_p, alpha_fc, beta, beta_ed, gamma, theta_alpha, geom);


            //___________________________________________________________________
            // Form Hessenberg matrix H
            for (int k=0; k<=i; ++k) {
                // H(k,i) = dot_product(w, V(k))
                //        = dot_product(w_u, V_u(k))+dot_product(w_p, V_p(k))
                StagInnerProd(w_u, 0, V_u, k, inner_prod_vel);
                CCInnerProd(w_p, 0, V_p, k, inner_prod_pres);
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
            StagL2Norm(w_u, 0, norm_u);
            CCL2Norm(w_p, 0, norm_p);
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
            norm_resid_est = std::abs(s[i+1]);


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
                if ((norm_resid_est <= gmres_rel_tol*std::min(norm_pre_b, norm_init_resid))
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



void IBMPrecon(const std::array<MultiFab, AMREX_SPACEDIM> & b_u, const MultiFab & b_p,
               std::array<MultiFab, AMREX_SPACEDIM> & x_u, MultiFab & x_p,
               const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
               const MultiFab & beta, const std::array<MultiFab, NUM_EDGE> & beta_ed,
               const MultiFab & gamma, const Real & theta_alpha,
               IBParticleContainer & ib_pc,
               const Geometry & geom)
{

    BL_PROFILE_VAR("ApplyPrecon()",ApplyPrecon);

    BoxArray ba = b_p.boxArray();
    DistributionMapping dmap = b_p.DistributionMap();

    Real         mean_val_pres;
    Vector<Real> mean_val_umac(AMREX_SPACEDIM);


    MultiFab phi     (ba,dmap,1,1);
    MultiFab mac_rhs (ba,dmap,1,0);
    MultiFab zero_fab(ba,dmap,1,0);
    MultiFab x_p_tmp (ba,dmap,1,1);

    // set zero_fab_fc to 0
    zero_fab.setVal(0.);

    // build alphainv_fc, one_fab_fc, zero_fab_fc, and b_u_tmp
    std::array< MultiFab, AMREX_SPACEDIM > alphainv_fc;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        alphainv_fc[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);

    std::array< MultiFab, AMREX_SPACEDIM > one_fab_fc;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        one_fab_fc[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);

    std::array< MultiFab, AMREX_SPACEDIM > zero_fab_fc;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        zero_fab_fc[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);

    std::array< MultiFab, AMREX_SPACEDIM > b_u_tmp;
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        b_u_tmp[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 0);


    // set alphainv_fc to 1/alpha_fc
    // set one_fab_fc to 1
    // set zero_fab_fc to 0
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        alphainv_fc[d].setVal(1.);
        alphainv_fc[d].divide(alpha_fc[d],0,1,0);
        one_fab_fc[d].setVal(1.);
        zero_fab_fc[d].setVal(0.);
    }

    // set the initial guess for Phi in the Poisson solve to 0
    // set x_u = 0 as initial guess
    phi.setVal(0.);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        x_u[d].setVal(0.);
    }

    // 1 = projection preconditioner
    // 2 = lower triangular preconditioner
    // 3 = upper triangular preconditioner
    // 4 = block diagonal preconditioner
    // 5 = Uzawa-type approximation (see paper)
    // 6 = upper triangular + viscosity-based BFBt Schur complement (from Georg Stadler)

    // projection preconditioner
    if (abs(precon_type) == 1) {

        ////////////////////
        // STEP 1: Solve for an intermediate state, x_u^star, using an implicit viscous solve
        //         x_u^star = A^{-1} b_u
        ////////////////////

        // x_u^star = A^{-1} b_u
        StagMGSolver(alpha_fc,beta,beta_ed,gamma,x_u,b_u,theta_alpha,geom);

        ////////////////////
        // STEP 2: Construct RHS for pressure Poisson problem
        ////////////////////

        // set mac_rhs = D(x_u^star)
        ComputeDiv(mac_rhs,x_u,0,0,1,geom,0);

        // add b_p to mac_rhs
        MultiFab::Add(mac_rhs,b_p,0,0,1,0);

        ////////////////////
        // STEP 3: Compute x_u
        ////////////////////

        // use multigrid to solve for Phi
        // x_u^star is only passed in to get a norm for absolute residual criteria
        MacProj(alphainv_fc,mac_rhs,phi,geom);

        // x_u = x_u^star - (alpha I)^-1 grad Phi
        SubtractWeightedGradP(x_u,alphainv_fc,phi,geom);

        ////////////////////
        // STEP 4: Compute x_p by applying the Schur complement approximation
        ////////////////////

        if (visc_schur_approx == 0) {
            // if precon_type = +1, or theta_alpha=0 then x_p = theta_alpha*Phi - c*beta*(mac_rhs)
            // if precon_type = -1                   then x_p = theta_alpha*Phi - c*beta*L_alpha Phi

            if (precon_type == 1 || theta_alpha == 0) {
                // first set x_p = -mac_rhs
                MultiFab::Copy(x_p,mac_rhs,0,0,1,0);
                x_p.mult(-1.,0,1,0);
            }
            else {
                // first set x_p = -L_alpha Phi
                CCApplyOp(phi,x_p,zero_fab,alphainv_fc,geom);
            }

            if ( abs(visc_type) == 1 || abs(visc_type) == 2) {
                // multiply x_p by beta; x_p = -beta L_alpha Phi
                MultiFab::Multiply(x_p,beta,0,0,1,0);

                if (abs(visc_type) == 2) {
                    // multiply by c=2; x_p = -2*beta L_alpha Phi
                    x_p.mult(2.,0,1,0);
                }
            }
            else if (abs(visc_type) == 3) {

                // multiply x_p by gamma, use mac_rhs a temparary to save x_p
                MultiFab::Copy(mac_rhs,x_p,0,0,1,0);
                MultiFab::Multiply(mac_rhs,gamma,0,0,1,0);
                // multiply x_p by beta; x_p = -beta L_alpha Phi
                MultiFab::Multiply(x_p,beta,0,0,1,0);
                // multiply by c=4/3; x_p = -(4/3) beta L_alpha Phi
                x_p.mult(4./3.,0,1,0);
                // x_p = -(4/3) beta L_alpha Phi - gamma L_alpha Phi
                MultiFab::Add(x_p,mac_rhs,0,0,1,0);
            }

            // multiply Phi by theta_alpha
            phi.mult(theta_alpha,0,1,0);

            // add theta_alpha*Phi to x_p
            MultiFab::Add(x_p,phi,0,0,1,0);
        }
        else {
            Abort("StagApplyOp: visc_schur_approx != 0 not supported");
        }
    }
    else {
        Abort("StagApplyOp: unsupposed precon_type");
    }

    ////////////////////
    // STEP 5: Handle null-space issues in MG solvers
    ////////////////////

    // subtract off mean value: Single level only! No need for ghost cells
    SumStag(x_u,0,mean_val_umac,true);
    SumCC(x_p,0,mean_val_pres,true);

    // The pressure Poisson problem is always singular:
    x_p.plus(-mean_val_pres,0,1,0);

    // The velocity problem is also singular under these cases
    if (theta_alpha == 0.) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            if (geom.isPeriodic(d)) {
                x_u[d].plus(-mean_val_umac[d],0,1,0);
            }
        }
    }
}
