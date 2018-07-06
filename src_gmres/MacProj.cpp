#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLLinOp.H>

using namespace common;
using namespace gmres;

void SubtractWeightedGradP(std::array<MultiFab, AMREX_SPACEDIM>& x_u,
                           const std::array<MultiFab, AMREX_SPACEDIM>& alphainv_fc,
                           MultiFab& phi,
                           const Real* dx,
                           const Geometry& geom)
{
    BoxArray ba = phi.boxArray();
    DistributionMapping dmap = phi.DistributionMap();

    std::array< MultiFab, AMREX_SPACEDIM > gradp;
    AMREX_D_TERM(gradp[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 gradp[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 gradp[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

    ComputeGrad(phi,gradp,0,0,1,geom);

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        MultiFab::Multiply(gradp[i],alphainv_fc[i],0,0,1,0);
        MultiFab::Saxpy(x_u[i],-1.,gradp[i],0,0,1,0);        
    }
}

void CCApplyOp(MultiFab& phi,
               MultiFab& Lphi,
               const MultiFab& alpha,
               const std::array<MultiFab, AMREX_SPACEDIM>& beta_fc,
               const Geometry& geom)
{
    int lev=0;

    BoxArray ba = phi.boxArray();
    DistributionMapping dmap = phi.DistributionMap();
    LPInfo info;
    info.setMaxCoarseningLevel(0); // turn off mg coarsening since no actual solve is performed

    MLABecLaplacian mlabec({geom}, {ba}, {dmap}, info);

    // order of stencil
    int stencil_order = 2;
    mlabec.setMaxOrder(stencil_order);

    mlabec.setScalars(0.0, 1.0);

    // build array of boundary conditions needed by MLABecLaplacian
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (geom.isPeriodic(i)) {
            mlmg_lobc[i] = mlmg_hibc[i] = LinOpBCType::Periodic;            
        }
        else {
            Abort("ApplyPrecon only works for periodic");
        }        
    }
    mlabec.setDomainBC(mlmg_lobc,mlmg_hibc);
    mlabec.setLevelBC(lev, &phi);

    // coefficients for solver
    mlabec.setACoeffs(lev,alpha);
    mlabec.setBCoeffs(lev,amrex::GetArrOfConstPtrs(beta_fc));

    MLMG mlmg(mlabec);

    Vector<MultiFab*> Lphi_vec(1);
    Vector<MultiFab*> phi_vec(1);

    Lphi_vec[0] = &Lphi;
    phi_vec[0] = &phi;
   
    mlmg.apply(Lphi_vec,phi_vec);
}
