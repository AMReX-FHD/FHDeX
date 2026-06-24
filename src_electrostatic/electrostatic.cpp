#include "electrostatic.H"
#include "common_functions.H"
#include <AMReX_MLMG.H>

using namespace amrex;

void esSolve(MultiFab& potential, MultiFab& charge,
             std::array< MultiFab, AMREX_SPACEDIM >& efieldCC,
             const std::array< MultiFab, AMREX_SPACEDIM >& external, const Geometry geom)
{
    BL_PROFILE_VAR("esSolve()",esSolve);

    AMREX_D_TERM(efieldCC[0].setVal(0);,
                 efieldCC[1].setVal(0);,
                 efieldCC[2].setVal(0););

    if(es_tog==1 || es_tog==3)
    {

        LinOpBCType lo_linop_bc[3];
        LinOpBCType hi_linop_bc[3];

        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            if (bc_es_lo[i] == -1 && bc_es_hi[i] == -1) {
                lo_linop_bc[i] = LinOpBCType::Periodic;
                hi_linop_bc[i] = LinOpBCType::Periodic;
            }
            if(bc_es_lo[i] == 2)
            {
                lo_linop_bc[i] = LinOpBCType::inhomogNeumann;
//                lo_linop_bc[i] = LinOpBCType::Neumann;
            }
            if(bc_es_hi[i] == 2)
            {
                hi_linop_bc[i] = LinOpBCType::inhomogNeumann;
//                hi_linop_bc[i] = LinOpBCType::Neumann;
            }
            if(bc_es_lo[i] == 1)
            {
                lo_linop_bc[i] = LinOpBCType::Dirichlet;
            }
            if(bc_es_hi[i] == 1)
            {
                hi_linop_bc[i] = LinOpBCType::Dirichlet;
            }
        }

        const BoxArray& ba = charge.boxArray();
        const DistributionMapping& dmap = charge.DistributionMap();

        //create solver opject
        MLPoisson linop({geom}, {ba}, {dmap});

        //set BCs
        linop.setDomainBC({AMREX_D_DECL(lo_linop_bc[0],
                                        lo_linop_bc[1],
                                        lo_linop_bc[2])},
                          {AMREX_D_DECL(hi_linop_bc[0],
                                        hi_linop_bc[1],
                                        hi_linop_bc[2])});

        // fill in ghost cells with Dirichlet/Neumann values
        // the ghost cells will hold the value ON the boundary
        MultiFabPotentialBC_solver(potential,geom);

        // tell MLPoisson about these potentially inhomogeneous BC values
        linop.setLevelBC(0, &potential);

        // this forces the solver to NOT enforce solvability
        // thus if there are Neumann conditions on phi they must
        // be correct or the Poisson solver won't converge
        linop.setEnforceSingularSolvable(false);

        //Multi Level Multi Grid
        MLMG mlmg(linop);

        //Solver parameters
        mlmg.setMaxIter(poisson_max_iter);
        mlmg.setVerbose(poisson_verbose);
        mlmg.setBottomVerbose(poisson_bottom_verbose);

        //Do solve
        mlmg.solve({&potential}, {&charge}, poisson_rel_tol, 0.0);

        potential.FillBoundary(geom.periodicity());
        // set ghost cell values so electric field is calculated properly
        // the ghost cells will hold the values extrapolated to the ghost CC
        MultiFabPotentialBC(potential, geom);

        //Find e field, gradient from cell centers to faces
        ComputeCentredGrad(potential, efieldCC, geom);


    }

    //Add external field on top, then fill boundaries, then setup BCs for peskin interpolation
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        efieldCC[d].mult(-1.0,efieldCC[d].nGrow());
        MultiFab::Add(efieldCC[d], external[d], 0, 0, 1, efieldCC[d].nGrow());
        efieldCC[d].FillBoundary(geom.periodicity());
        MultiFabElectricBC(efieldCC[d], geom);
    }

}


