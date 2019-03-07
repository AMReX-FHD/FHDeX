#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"


using namespace common;
using namespace gmres;

// This computes A x = b explicitly
void ApplyMatrix(std::array<MultiFab, AMREX_SPACEDIM>& b_u,
                 MultiFab& b_p,
                 std::array<MultiFab, AMREX_SPACEDIM>& x_u,
                 MultiFab& x_p,
                 const std::array<MultiFab, AMREX_SPACEDIM>& alpha_fc,
                 const MultiFab& beta,
                 const std::array<MultiFab, NUM_EDGE>& beta_ed,
                 const MultiFab& gamma,
                 const Real& theta_alpha,
                 const Geometry& geom)
{

    BL_PROFILE_VAR("ApplyMatrix()",ApplyMatrix);

    BoxArray ba = b_p.boxArray();
    DistributionMapping dmap = b_p.DistributionMap();

    const Real* dx = geom.CellSize();

    // check to make sure x_u and x_p have enough ghost cells
    if (gmres_spatial_order == 2) {
        if (x_u[0].nGrow() < 1) {
            Abort("apply_matrix.f90: x_u needs at least 1 ghost cell");
        }
        if (x_p.nGrow() < 1) {
            Abort("apply_matrix.f90: x_p needs at least 1 ghost cell");
        }
        else if (gmres_spatial_order == 4) {
            if (x_u[0].nGrow() < 2) {
                Abort("apply_matrix.f90: x_u needs at least 2 ghost cells");
            }
            if (x_p.nGrow() < 2) {
                Abort("apply_matrix.f90: x_p needs at least 2 ghost cells");
            }
        }
    }

    // fill ghost cells for x_u and x_p
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        x_u[i].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(x_u[i], i, geom);
        MultiFABPhysBCMacVel(x_u[i], i, geom);
    }
    x_p.FillBoundary(geom.periodicity());
    MultiFABPhysBC(x_p, geom);

    std::array< MultiFab, AMREX_SPACEDIM > gx_p;
    AMREX_D_TERM(gx_p[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 gx_p[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 gx_p[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

    // compute b_u = A x_u
    if (gmres_spatial_order == 2) {
        StagApplyOp(beta,gamma,beta_ed,x_u,b_u,alpha_fc,dx,theta_alpha);
    }
    else if (gmres_spatial_order == 4) {
        Abort("ApplyMatrix.cpp: gmres_spatial_order=4 not supported yet");
    }

    // compute G x_p and add to b_u
    if (gmres_spatial_order == 2) {
        ComputeGrad(x_p,gx_p,0,0,1,geom);
    }
    else if (gmres_spatial_order == 4) {
        Abort("ApplyMatrix.cpp: gmres_spatial_order=4 not supported yet");
    }

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        MultiFab::Add(b_u[i],gx_p[i],0,0,1,0);
    }

    // set b_p = -D x_u
    ComputeDiv(b_p,x_u,0,0,1,geom,0);
    b_p.mult(-1.,0,1,0);

}
