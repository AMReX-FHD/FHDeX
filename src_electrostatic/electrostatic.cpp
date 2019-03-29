#include "electrostatic.H"
#include "common_namespace.H"
#include <AMReX_MLMG.H>

using namespace amrex;
using namespace common;

void poissonSolve(MultiFab& solution, const MultiFab& rhs, const Geometry geom)
{

    const BoxArray& ba = solution.boxArray();
    const DistributionMapping& dmap = solution.DistributionMap();


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

    mlmg.solve({&solution}, {&rhs}, poisson_rel_tol, 0.0);

}

void calculateField(MultiFab& potential, const Geometry geom)
{

    const BoxArray& ba = potential.boxArray();
    const DistributionMapping& dmap = potential.DistributionMap();

}
