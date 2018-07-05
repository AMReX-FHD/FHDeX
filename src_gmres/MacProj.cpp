#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"
#include "gmres_namespace.H"

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
