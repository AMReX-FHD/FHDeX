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
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            b_u[i].mult(scale_factor,0,1,b_u[i].nGrow());
        }
        // scale the viscosities:
        beta.mult(scale_factor,0,1,beta.nGrow());
        gamma.mult(scale_factor,0,1,gamma.nGrow());
        for (int i=0; i<NUM_EDGE; ++i) {
            beta_ed[i].mult(scale_factor,0,1,beta_ed[i].nGrow());
        }
    }

    // preconditioned norm_b: norm_pre_b
    ApplyPrecon(b_u,b_p,tmp_u,tmp_p,alpha_fc,beta,beta_ed,gamma,theta_alpha,geom);
    // call stag_l2_norm(mla,tmp_u,norm_u)
    // call cc_l2_norm(mla,tmp_p,norm_p)
    norm_p = p_norm_weight*norm_p;
    norm_pre_b = sqrt(norm_u*norm_u+norm_p*norm_p);

    //! If norm_b=0 we should return zero as the solution and "return" from this routine
    // It is important to use gmres_abs_tol and not 0 since sometimes due to roundoff we 
    // get a nonzero number that should really be zero
    if (gmres_verbose >= 1) {
        // Useful to print out to give expected scale for gmres_abs_tol
        Print() << "GMRES.cpp: GMRES called with ||rhs||=" << norm_b << std::endl; 
    } 
    if (norm_b <= gmres_abs_tol) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            x_u[i].setVal(0.);
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




    } while (true); // end of outer loop (do iter=1,gmres_max_outer)

}
