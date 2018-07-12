#include <AMReX_VisMF.H>

#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"

using namespace common;
using namespace gmres;

void GMRES(std::array<MultiFab, AMREX_SPACEDIM>& b_u,
           const MultiFab& b_p,
           std::array<MultiFab, AMREX_SPACEDIM>& x_u,
           MultiFab& x_p,
           const std::array<MultiFab, AMREX_SPACEDIM>& alpha_fc,
           MultiFab& beta,
           std::array<MultiFab, NUM_EDGE>& beta_ed,
           MultiFab& gamma,
           Real theta_alpha,
           const Geometry& geom)
{
    Real temp;
    GMRES(b_u,b_p,x_u,x_p,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom,temp);
}

void GMRES(std::array<MultiFab, AMREX_SPACEDIM>& b_u,
           const MultiFab& b_p,
           std::array<MultiFab, AMREX_SPACEDIM>& x_u,
           MultiFab& x_p,
           const std::array<MultiFab, AMREX_SPACEDIM>& alpha_fc,
           MultiFab& beta,
           std::array<MultiFab, NUM_EDGE>& beta_ed,
           MultiFab& gamma,
           Real theta_alpha,
           const Geometry& geom,
           Real& norm_pre_rhs)
{

    Vector<Real> cs(gmres_max_inner);
    Vector<Real> sn(gmres_max_inner);
    Vector<Real>  y(gmres_max_inner);
    Vector<Real>  s(gmres_max_inner);
    
    Vector<Vector<Real>> H(gmres_max_inner+1,Vector<Real>(gmres_max_inner));

    int iter, total_iter, i_copy; // for looping iteration
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
    AMREX_D_TERM(r_u[0].define(convert(ba,nodal_flag_x), dmap, 1, x_u[0].nGrow());,
                 r_u[1].define(convert(ba,nodal_flag_y), dmap, 1, x_u[0].nGrow());,
                 r_u[2].define(convert(ba,nodal_flag_z), dmap, 1, x_u[0].nGrow()););

    std::array< MultiFab, AMREX_SPACEDIM > w_u;
    AMREX_D_TERM(w_u[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 w_u[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 w_u[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

    std::array< MultiFab, AMREX_SPACEDIM > tmp_u;
    AMREX_D_TERM(tmp_u[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 tmp_u[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 tmp_u[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

    std::array< MultiFab, AMREX_SPACEDIM > V_u;
    AMREX_D_TERM(V_u[0].define(convert(ba,nodal_flag_x), dmap, gmres_max_inner+1, 0);,
                 V_u[1].define(convert(ba,nodal_flag_y), dmap, gmres_max_inner+1, 0);,
                 V_u[2].define(convert(ba,nodal_flag_z), dmap, gmres_max_inner+1, 0););

    // # of ghost cells must match x_p so higher-order stencils can work
    MultiFab r_p  (ba,dmap,                1,x_p.nGrow());
    MultiFab w_p  (ba,dmap,                1,0);
    MultiFab tmp_p(ba,dmap,                1,0);
    MultiFab V_p  (ba,dmap,gmres_max_inner+1,0); // Krylov vectors

    // apply scaling factor
    if (scale_factor != 1.) {
        theta_alpha = theta_alpha*scale_factor;
        // we will solve for scale*x_p so we need to scale the initial guess
        x_p.mult(scale_factor,0,1,x_p.nGrow());
        // scale the rhs:
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            b_u[d].mult(scale_factor,0,1,b_u[d].nGrow());
        }
        // scale the viscosities:
        beta.mult(scale_factor,0,1,beta.nGrow());
        gamma.mult(scale_factor,0,1,gamma.nGrow());
        for (int d=0; d<NUM_EDGE; ++d) {
            beta_ed[d].mult(scale_factor,0,1,beta_ed[d].nGrow());
        }
    }

    // preconditioned norm_b: norm_pre_b
    ApplyPrecon(b_u,b_p,tmp_u,tmp_p,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom);
    StagL2Norm(tmp_u,0,norm_u);
    CCL2Norm(tmp_p,0,norm_p);
    norm_p = p_norm_weight*norm_p;
    norm_pre_b = sqrt(norm_u*norm_u+norm_p*norm_p);

    // calculate the l2 norm of rhs
    StagL2Norm(b_u,0,norm_u);
    CCL2Norm(b_p,0,norm_p);
    norm_p = p_norm_weight*norm_p;
    norm_b = sqrt(norm_u*norm_u+norm_p*norm_p);

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
    iter = 0;


    do {
        // Calculate tmp = Ax
        ApplyMatrix(tmp_u,tmp_p,x_u,x_p,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom);

        // tmp = b - Ax
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Subtract(tmp_u[d],b_u[d],0,0,1,0);
            tmp_u[d].mult(-1.,0,1,0);
        }            
        MultiFab::Subtract(tmp_p,b_p,0,0,1,0);
        tmp_p.mult(-1.,0,1,0);
        
        // un-preconditioned residuals
        StagL2Norm(tmp_u,0,norm_u_noprecon);
        CCL2Norm(tmp_p,0,norm_p_noprecon);
        norm_p_noprecon = p_norm_weight*norm_p_noprecon;
        norm_resid_Stokes=sqrt(norm_u_noprecon*norm_u_noprecon+norm_p_noprecon*norm_p_noprecon);
        if(iter==0) {
            norm_init_Stokes=norm_resid_Stokes;
        }

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

        // solve for r = M^{-1} tmp
        // We should not be counting these toward the number of mg cycles performed
        ApplyPrecon(tmp_u,tmp_p,r_u,r_p,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom);

        // resid = sqrt(dot_product(r, r))
        StagL2Norm(r_u,0,norm_u);
        CCL2Norm(r_p,0,norm_p);
        norm_p = p_norm_weight*norm_p;
        norm_resid = sqrt(norm_u*norm_u+norm_p*norm_p);
        // If first iteration, save the initial preconditioned residual
        if (iter==0) {
            norm_init_resid=norm_resid;
            norm_resid_est=norm_resid;
        }

        if (gmres_verbose >= 3) {
            Print() << "Precond. rel. res. (u,v,p) = " << norm_resid/norm_init_resid << "  " 
                    << norm_u/norm_init_resid << "  " << norm_p/norm_init_resid << std::endl;
        }

        // We need to test the residual now and exit OuterLoop if converged
        if (total_iter >= gmres_max_iter) {
            if (gmres_verbose >= 1) {
                Print() << "GMRES did not converge in max number of total inner iterations: Exiting" 
                        << std::endl;
            }
            break;
        }
        else if (total_iter >= gmres_min_iter) {
            // other options
            if(norm_resid <= gmres_rel_tol*std::min(norm_pre_b, norm_init_resid)) {
                if (gmres_verbose >= 2) {
                    Print() << "GMRES converged: Outer = " << iter << ",  Innter = " << i 
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
            }
            else if (norm_resid <= gmres_abs_tol) {

                if (gmres_verbose >= 2) {
                    Print() << "GMRES converged: Outer = " << iter << ",  Inner = " << i 
                            << " Total=" << total_iter << std::endl;
                }

                break; // exit OuterLoop
            }
        }

        if (iter >= gmres_max_outer) {
            Print() << "GMRES did not converge in max number of outer iterations: Exiting" << std::endl;
            break; // exit OuterLoop
        }
        iter=iter+1;

        // create the first basis in Krylov space
        // V(1) = r / norm(r)
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Copy(V_u[d],r_u[d],0,0,1,0);
            V_u[d].mult(1./norm_resid,0,1,0);
        }
        MultiFab::Copy(V_p,r_p,0,0,1,0);
        V_p.mult(1./norm_resid,0,1,0);

        // s = norm(r) * e_0
        std::fill(s.begin(), s.end(), 0.);
        s[0] = norm_resid;

        ///////////////////////
        // begin inner iteration
        ///////////////////////

        for (i=0; i<gmres_max_inner; ++i) {
            total_iter = total_iter + 1;
            i_copy = i;

            // tmp=A*V(i)
            // we use r_p and r_u as temporaries to hold ith component of V
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                MultiFab::Copy(r_u[d],V_u[d],i,0,1,0);
            }
            MultiFab::Copy(r_p,V_p,i,0,1,0);

            ApplyMatrix(tmp_u,tmp_p,r_u,r_p,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom);

            // w = M^{-1} A*V(i)
            ApplyPrecon(tmp_u,tmp_p,w_u,w_p,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom);

            for (int k=0; k<=i; ++k) {
                // form H(k,i) Hessenberg matrix
                // H(k,i) = dot_product(w, V(k))
                //        = dot_product(w_u, V_u(k))+dot_product(w_p, V_p(k))
                StagInnerProd(w_u,0,V_u,k,inner_prod_vel);
                CCInnerProd(w_p,0,V_p,k,inner_prod_pres);
                H[k][i] = std::accumulate(inner_prod_vel.begin(), inner_prod_vel.end(), 0) +
                          pow(p_norm_weight,2.0)*inner_prod_pres;

                // w = w - H(k,i) * V(k)
                //use tmp_u and tmp_p as temporaries to hold kth component of V(k)
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    MultiFab::Copy(tmp_u[d],V_u[d],k,0,1,0);
                    tmp_u[d].mult(H[k][i],0,1,0);
                    MultiFab::Subtract(w_u[d],tmp_u[d],0,0,1,0);
                }
                MultiFab::Copy(tmp_p,V_p,k,0,1,0);
                tmp_p.mult(H[k][i],0,1,0);
                MultiFab::Subtract(w_p,tmp_p,0,0,1,0);
            }

            // H(i+1,i) = norm(w)
            StagL2Norm(w_u,0,norm_u);
            CCL2Norm(w_p,0,norm_p);
            norm_p = p_norm_weight*norm_p;
            H[i+1][i] = sqrt(norm_u*norm_u+norm_p*norm_p);

            // V(i+1) = w / H(i+1,i)
            if (H[i+1][i] != 0.) {
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    MultiFab::Copy(V_u[d],w_u[d],0,i+1,1,0);
                    V_u[d].mult(1./H[i+1][i],i+1,1,0);
                }
                MultiFab::Copy(V_p,w_p,0,1,1,0);
                V_p.mult(1./H[i+1][i],i+1,1,0);
            }
            else {
                Abort("GMRES.cpp: error in orthogonalization");
            }

            LeastSquares(i,H,cs,sn,s); // solve least square problem
            norm_resid_est = abs(s[i+1]);

            if (gmres_verbose >= 2) {
                Print() << total_iter << ",  est. rel. resid. |Pr|/(Pr0,b)= "
                        << norm_resid_est/norm_init_resid << "  "
                        << norm_resid_est/norm_pre_b << std::endl;
            }

            if (total_iter >= gmres_max_iter) {
                break; // exit InnerLoop
            }
            else if (total_iter >= gmres_min_iter) {
                if ((norm_resid_est <= gmres_rel_tol*std::min(norm_pre_b, norm_init_resid))
                    || (norm_resid_est <= gmres_abs_tol)) {
                    break; // exit InnerLoop
                }
            }

        } // end of inner loop (do i=1,gmres_max_inner)

        // update the solution   ! first, solve for y
        SolveUTriangular(i_copy-1, H, s, y);

        // then, x = x + dot(V(1:i),y(1:i))
        UpdateSol(x_u,x_p,V_u,V_p,y,i_copy);

    } while (true); // end of outer loop (do iter=1,gmres_max_outer)

    // AJN - this is here since I notice epsilon roundoff errors building up
    //       just enough to destroy the asymmetry in time-advancement codes that
    //       ultimately causes lack of convergence in subsequent gmres calls
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        x_u[d].OverrideSync(geom.periodicity());
    }

    // apply scaling factor
    if (scale_factor != 1.) {
        theta_alpha = theta_alpha/scale_factor;
        // the solution we got is scale*x_p
        x_p.mult(1./scale_factor,0,1,x_p.nGrow());
        // unscale the rhs
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            b_u[d].mult(1./scale_factor,0,1,b_u[d].nGrow());
        }
        // unscale the viscosities
        beta.mult(1./scale_factor,0,1,beta.nGrow());
        gamma.mult(1./scale_factor,0,1,gamma.nGrow());
        for (int d=0; d<NUM_EDGE; ++d) {
            beta_ed[d].mult(1./scale_factor,0,1,beta_ed[d].nGrow());
        }
    }

    if (gmres_verbose >= 1) {
        Print() << "Preconditioned GMRES:" << std::endl;
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
        H[k+1][i] = -sn[k]*H[k][i] + cs[k]*H[i+1][i];
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
    else if (abs(b) > abs(a)) {
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
