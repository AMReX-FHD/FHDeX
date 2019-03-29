#include "electrostatic.H"

using namespace amrex;

void poissonSolve(MultiFab& solution, const MultiFab& rhs )
{


//    MLPoisson linop({geom}, {grids}, {dmap});

//    linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic,
//                                    LinOpBCType::Periodic,
//                                    LinOpBCType::Periodic)},
//                      {AMREX_D_DECL(LinOpBCType::Periodic,
//                                    LinOpBCType::Periodic,
//                                    LinOpBCType::Periodic)});

//    linop.setLevelBC(0, nullptr);

//    linop.setSigma(0, sigma);

//    MLMG mlmg(linop);
//    mlmg.setMaxIter(max_iter);
//    mlmg.setMaxFmgIter(max_fmg_iter);
//    mlmg.setVerbose(verbose);
//    mlmg.setBottomVerbose(bottom_verbose);

//    mlmg.solve({&solution}, {&rhs}, reltol, 0.0);

}
