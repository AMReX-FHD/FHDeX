#include "electrostatic.H"
#include "common_namespace.H"
#include "common_functions.H"
#include <AMReX_MLMG.H>

using namespace amrex;
using namespace common;

void esSolve(MultiFab& potential, const MultiFab& charge, std::array< MultiFab, AMREX_SPACEDIM >& efield, const std::array< MultiFab, AMREX_SPACEDIM >& external, const Geometry geom)
{

    const BoxArray& ba = potential.boxArray();
    const DistributionMapping& dmap = potential.DistributionMap();


    MLPoisson linop({geom}, {ba}, {dmap});

    linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic,
                                    LinOpBCType::Periodic,
                                    LinOpBCType::Periodic)},
                      {AMREX_D_DECL(LinOpBCType::Periodic,
                                    LinOpBCType::Periodic,
                                    LinOpBCType::Periodic)});

    linop.setLevelBC(0, nullptr);

    MLMG mlmg(linop);
    mlmg.setMaxIter(poisson_max_iter);

    mlmg.setVerbose(poisson_verbose);
    mlmg.setBottomVerbose(poisson_bottom_verbose);

    mlmg.solve({&potential}, {&charge}, poisson_rel_tol, 0.0);

    potential.FillBoundary(geom.periodicity());

    ComputeGrad(potential, efield, 0, 0, 1, geom);

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Add(efield[d], external[d], 0, 0, 1, 0);
        efield[d].FillBoundary(geom.periodicity());
    }

}


