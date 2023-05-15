#include "common_functions.H"

#include "gmres_functions.H"
using namespace amrex;

void ApplyMatrix(std::array<MultiFab, AMREX_SPACEDIM> & b_u,
                 MultiFab                             & b_p,
                 std::array<MultiFab, AMREX_SPACEDIM> & x_u,
                 MultiFab                             & x_p,
                 std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                 MultiFab                             & beta,
                 std::array<MultiFab, NUM_EDGE>       & beta_ed,
                 MultiFab                             & gamma,
                 Real                                 & theta_alpha,
                 Geometry                             & geom,
                 int                                  is_inhomogeneous)
{
    Vector<std::array< MultiFab*, AMREX_SPACEDIM >> b_up(1,{AMREX_D_DECL(&b_u[0],&b_u[1],&b_u[2])});
    Vector<std::array< MultiFab*, AMREX_SPACEDIM >> x_up(1,{AMREX_D_DECL(&x_u[0],&x_u[1],&x_u[2])});
    Vector<std::array< MultiFab*, AMREX_SPACEDIM >> alpha_fcp(1,{AMREX_D_DECL(&alpha_fc[0],&alpha_fc[1],&alpha_fc[2])});        
    const Vector<std::array< MultiFab*, NUM_EDGE >> beta_edp(1,{AMREX_D_DECL(&beta_ed[0],&beta_ed[1],&beta_ed[2])});
    
    Vector<MultiFab*> b_pp = {&b_p};
    Vector<MultiFab*> x_pp = {&x_p};
    const Vector<MultiFab*> betap(1,&beta);
    const Vector<MultiFab*> gammap(1,&gamma);
    
    const Vector<Geometry> geomv(1,geom);

    ApplyMatrix(b_up, b_pp, x_up, x_pp, alpha_fcp, betap, beta_edp, gammap, theta_alpha, geomv, is_inhomogeneous,1);
}                 


// This computes A x = b explicitly
void ApplyMatrix(Vector<std::array< MultiFab*, AMREX_SPACEDIM >> & b_u,
                 Vector<MultiFab*>                               & b_p,
                 Vector<std::array< MultiFab*, AMREX_SPACEDIM >> & x_u,
                 Vector<MultiFab*>                               & x_p,
                 Vector<std::array< MultiFab*, AMREX_SPACEDIM >> & alpha_fc,
                 const Vector<MultiFab*>                         & beta,
                 const Vector<std::array< MultiFab*, NUM_EDGE >> & beta_ed,
                 const Vector<MultiFab*>                         & gamma,
                 const Real                                      & theta_alpha,
                 const Vector<Geometry>                          & geom,
                 int                                             is_inhomogeneous,
                 int                                             nlev) {

    BL_PROFILE_VAR("ApplyMatrix()", ApplyMatrix);

    Vector<BoxArray> ba(nlev);
    Vector<DistributionMapping> dmap(nlev);
    Vector<const Real*> dx(nlev);
    for(int lev=0;lev<nlev;lev++)
    {
        const BoxArray & bas = b_p[lev]->boxArray();
        const DistributionMapping & dmaps = b_p[lev]->DistributionMap();
        const Real* dxs = geom[lev].CellSize();
        ba[lev] = bas;
        dmap[lev] = dmaps;
        dx[lev] = dxs;
    }



    // check to make sure x_u and x_p have enough ghost cells
    if (gmres_spatial_order == 2) {
        if (x_u[0][0]->nGrow() < 1) {
            Abort("apply_matrix.f90: x_u needs at least 1 ghost cell");
        }
        if (x_p[0]->nGrow() < 1) {
            Abort("apply_matrix.f90: x_p needs at least 1 ghost cell");
        } else if (gmres_spatial_order == 4) {
            if (x_u[0][0]->nGrow() < 2) {
                Abort("apply_matrix.f90: x_u needs at least 2 ghost cells");
            }
            if (x_p[0]->nGrow() < 2) {
                Abort("apply_matrix.f90: x_p needs at least 2 ghost cells");
            }
        }
    }

//    // fill ghost cells for x_u and x_p
//    // TODO: this should not be done here
    for (int i=0; i<AMREX_SPACEDIM; ++i) {

        MultiFabPhysBCDomainVel(*x_u[0][i], geom[0],i);
        MultiFabPhysBCMacVel(*x_u[0][i], geom[0], i, is_inhomogeneous);
    }
    MultiFabPhysBC(*x_p[0], geom[0], 0, 1, PRES_BC_COMP);

    for(int lev=0;lev<nlev;lev++)
    {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            x_u[lev][i]->FillBoundary(geom[lev].periodicity());
        }
        x_p[lev]->FillBoundary(geom[lev].periodicity());
    }

    // compute b_u = A x_u
    if (gmres_spatial_order == 2) {
        StagApplyOp(geom[0], *beta[0], *gamma[0], *beta_ed[0], *x_u[0], *b_u[0], *alpha_fc[0], dx[0], theta_alpha);
    }
    else if (gmres_spatial_order == 4) {
        Abort("ApplyMatrix.cpp: gmres_spatial_order=4 not supported yet");
    }

//    // compute G x_p and add to b_u
//    if (gmres_spatial_order == 2) {
//        ComputeGrad(x_p, b_u, 0, 0, 1, PRES_BC_COMP, geom, 1);
//    }
//    else if (gmres_spatial_order == 4) {
//        Abort("ApplyMatrix.cpp: gmres_spatial_order=4 not supported yet");
//    }

//    // set b_p = -D x_u
//    ComputeDiv(b_p, x_u, 0, 0, 1, geom, 0);
//    b_p.mult(-1., 0, 1, 0);
}
