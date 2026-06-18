#ifndef _StagMGSolver_H_
#define _StagMGSolver_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include "common_functions.H"

using namespace amrex;

class StagMGSolver {

    //////////////////////////////////
    // MultiFabs used in multigrid coarsening

    // face-centered
    // outer vector will be over nlevs_mg; innter array is for face-centered
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > alpha_fc_mg;
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >   rhs_fc_mg;
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >   phi_fc_mg;
    Vector<std::array< MultiFab, AMREX_SPACEDIM > >  Lphi_fc_mg;
    Vector<std::array< MultiFab, AMREX_SPACEDIM > > resid_fc_mg;
    Vector<std::array< MultiFab, NUM_EDGE       > >  beta_ed_mg; // nodal in 2D, edge in 3D

    // cell-centered
    // vector will be over nlevs_mg
    Vector<MultiFab>  beta_cc_mg;
    Vector<MultiFab> gamma_cc_mg;

    // needs to sized to nlevs_mg
    Vector<std::array< Real, AMREX_SPACEDIM > > dx_mg;

    Vector<Geometry> geom_mg;

    int nlevs_mg;

    Box pd_base;
    BoxArray ba_base;
    DistributionMapping dmap;

public:

    StagMGSolver();

    void Define(const BoxArray& ba_in,
                const DistributionMapping& dmap_in,
                const Geometry& geom_in);

    // solve "(theta*alpha*I - L) phi = rhs" using multigrid with Jacobi relaxation
    // if abs(visc_type) = 1, L = div beta grad
    // if abs(visc_type) = 2, L = div [ beta (grad + grad^T) ]
    // if abs(visc_type) = 3, L = div [ beta (grad + grad^T) + I (gamma - (2/3)*beta) div ]
    // if visc_type > 1 we assume constant coefficients
    // if visc_type < 1 we assume variable coefficients
    // beta_cc, and gamma_cc are cell-centered
    // alpha_fc, phi_fc, and rhs_fc are face-centered
    // beta_ed is nodal (2d) or edge-centered (3d)
    // phi_fc must come in initialized to some value, preferably a reasonable guess
    void Solve(const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
               const MultiFab & beta_cc,
               const std::array<MultiFab, NUM_EDGE> & beta_ed,
               const MultiFab & gamma_cc,
               std::array<MultiFab, AMREX_SPACEDIM> & phi_fc,
               const std::array<MultiFab, AMREX_SPACEDIM> & phiorig_fc,
               const Real & theta);


    // compute the number of multigrid levels assuming minwidth is the length of the
    // smallest dimension of the smallest grid at the coarsest multigrid level
    int ComputeNlevsMG(const BoxArray & ba);

    void CCRestriction(MultiFab & phi_c, const MultiFab & phi_f,
                       const Geometry & geom_c);

    void StagRestriction(std::array<MultiFab, AMREX_SPACEDIM> & phi_c,
                         const std::array<MultiFab, AMREX_SPACEDIM > & phi_f,
                         int simple_stencil=0);

    void NodalRestriction(MultiFab & phi_c, const MultiFab & phi_f);

    void EdgeRestriction(std::array<MultiFab, NUM_EDGE> & phi_c,
                         const std::array<MultiFab, NUM_EDGE> & phi_f);

    void StagProlongation(const std::array<MultiFab, AMREX_SPACEDIM> & phi_c_in,
                          std::array<MultiFab, AMREX_SPACEDIM> & phi_f_in);

    void StagMGUpdate(std::array<MultiFab, AMREX_SPACEDIM> & phi_fc,
                      const std::array<MultiFab, AMREX_SPACEDIM> & rhs_fc,
                      const std::array<MultiFab, AMREX_SPACEDIM> & Lphi_fc,
                      const std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                      const MultiFab & beta_cc,
                      const std::array<MultiFab, NUM_EDGE> & beta_ed,
                      const MultiFab & gamma_cc,
                      const Real * dx,
                      const int & color=0);

};

#endif
