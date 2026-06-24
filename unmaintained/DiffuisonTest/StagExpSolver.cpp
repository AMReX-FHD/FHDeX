#include "gmres_functions.H"

#include "common_functions.H"
#include <AMReX_VisMF.H>

using namespace amrex;

// solve "(theta*alpha*I - L) phi = rhs" using multigrid with Gauss-Seidel relaxation
// if amrex::Math::abs(visc_type) = 1, L = div beta grad
// if amrex::Math::abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
// if amrex::Math::abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
// if visc_type > 1 we assume constant coefficients
// if visc_type < 1 we assume variable coefficients
// beta_cc, and gamma_cc are cell-centered
// alpha_fc, phi_fc, and rhs_fc are face-centered
// beta_ed is nodal (2d) or edge-centered (3d)
// phi_fc must come in initialized to some value, preferably a reasonable guess
void StagExpSolver(const std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
                  const MultiFab& beta_cc,
                  const std::array< MultiFab, NUM_EDGE >& beta_ed,
                  const MultiFab& gamma_cc,
                  std::array< MultiFab, AMREX_SPACEDIM >& phi_fc,
                  const std::array< MultiFab, AMREX_SPACEDIM >& phiorig_fc,
                  const Real& theta,
                  const Geometry& geom)
{

    // get the problem domain and boxarray at level 0
    Box pd_base = geom.Domain();
    BoxArray ba_base = beta_cc.boxArray();

    Box pd = pd_base;
    BoxArray ba(ba_base);
    //////////////////////////////////

    Real weight_lap;
    DistributionMapping dmap = beta_cc.DistributionMap();
    const Real* dx = geom.CellSize();

    std::array< MultiFab, AMREX_SPACEDIM >    phipred_fc;
    std::array< MultiFab, AMREX_SPACEDIM >    Lphipred_fc;
    std::array< MultiFab, AMREX_SPACEDIM >    Lphi_fc;

    AMREX_D_TERM(phipred_fc[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 phipred_fc[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 phipred_fc[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
    AMREX_D_TERM(Lphipred_fc[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 Lphipred_fc[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 Lphipred_fc[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););
    AMREX_D_TERM(Lphi_fc[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);,
                 Lphi_fc[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);,
                 Lphi_fc[2].define(convert(ba,nodal_flag_z), dmap, 1, 1););


    // for (int d=0; d<AMREX_SPACEDIM; ++d) {
    //   MultiFab::Copy(phipred_fc[d],phiorig_fc[d],0,0,1,0);
    //   MultiFab::Copy(phi_fc[d],phiorig_fc[d],0,0,1,0);
    // }

   // phiorig_fc[0].FillBoundary(geom.periodicity());
   // phiorig_fc[1].FillBoundary(geom.periodicity());
   // phiorig_fc[2].FillBoundary(geom.periodicity());

    AMREX_D_TERM(MultiFab::Copy(phipred_fc[0],phiorig_fc[0],0,0,1,0);,
                 MultiFab::Copy(phipred_fc[1],phiorig_fc[1],0,0,1,0);,
                 MultiFab::Copy(phipred_fc[2],phiorig_fc[2],0,0,1,0););
    AMREX_D_TERM(MultiFab::Copy(phi_fc[0],phiorig_fc[0],0,0,1,0);,
                 MultiFab::Copy(phi_fc[1],phiorig_fc[1],0,0,1,0);,
                 MultiFab::Copy(phi_fc[2],phiorig_fc[2],0,0,1,0););

    AMREX_D_TERM(phi_fc[0].FillBoundary(geom.periodicity());,
                 phi_fc[1].FillBoundary(geom.periodicity());,
                 phi_fc[2].FillBoundary(geom.periodicity()););

    AMREX_D_TERM(phipred_fc[0].FillBoundary(geom.periodicity());,
                 phipred_fc[1].FillBoundary(geom.periodicity());,
                 phipred_fc[2].FillBoundary(geom.periodicity()););

    StagApplyOp(geom,beta_cc,gamma_cc,beta_ed,
                phi_fc,Lphipred_fc,alpha_fc,dx,1.);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Add(phipred_fc[d],Lphipred_fc[d],0,0,1,0);
    }

    // VisMF::Write(Lphipred_fc[0],"a_Lphipred0");
    // Abort();

    AMREX_D_TERM(phipred_fc[0].FillBoundary(geom.periodicity());,
                 phipred_fc[1].FillBoundary(geom.periodicity());,
                 phipred_fc[2].FillBoundary(geom.periodicity()););

    StagApplyOp(geom,beta_cc,gamma_cc,beta_ed,
                phipred_fc,Lphi_fc,alpha_fc,dx,1.);
    weight_lap = 0.5;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        Lphipred_fc[d].mult(weight_lap,0,1,0);
        Lphi_fc[d].mult(weight_lap,0,1,0);
        MultiFab::Add(phi_fc[d],Lphipred_fc[d],0,0,1,0);
        MultiFab::Add(phi_fc[d],Lphi_fc[d],0,0,1,0);
    }

    AMREX_D_TERM(phi_fc[0].FillBoundary(geom.periodicity());,
                 phi_fc[1].FillBoundary(geom.periodicity());,
                 phi_fc[2].FillBoundary(geom.periodicity()););

    //////////////////////////////////

    if (stag_mg_verbosity >= 1) {
        Print() << "\nEnd call to stag_exp_solver\n";
    }

}