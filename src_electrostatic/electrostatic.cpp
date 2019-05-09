#include "electrostatic.H"
#include "common_namespace.H"
#include "common_functions.H"
#include <AMReX_MLMG.H>

using namespace amrex;
using namespace common;

void esSolve(MultiFab& potential, MultiFab& charge, std::array< MultiFab, AMREX_SPACEDIM >& efield, const std::array< MultiFab, AMREX_SPACEDIM >& external, const Geometry geom)
{
    AMREX_D_TERM(efield[0].setVal(0);,
                 efield[1].setVal(0);,
                 efield[2].setVal(0););

    if(es_tog==1)
    {

            LinOpBCType lobc[3];
            LinOpBCType hibc[3];

            for (int i=0; i<AMREX_SPACEDIM; ++i) {
                if (bc_lo[i] == -1 && bc_hi[i] == -1) {
                    lobc[i] = LinOpBCType::Periodic;
                    hibc[i] = LinOpBCType::Periodic;
                }
                if(bc_lo[i] == 2)
                {
                    lobc[i] = LinOpBCType::Neumann;
                }
                if(bc_hi[i] == 2)
                {
                    hibc[i] = LinOpBCType::Neumann;
                }
            }


            MultiFABChargeBC(charge, geom); //Adjust spread charge distribtion near boundaries from 

            const BoxArray& ba = potential.boxArray();
            const DistributionMapping& dmap = potential.DistributionMap();

            //create solver opject
            MLPoisson linop({geom}, {ba}, {dmap});

            //set BCs
            linop.setDomainBC({AMREX_D_DECL(lobc[0],
                                            lobc[1],
                                            lobc[2])},
                              {AMREX_D_DECL(hibc[0],
                                            hibc[1],
                                            hibc[2])});

            linop.setLevelBC(0, nullptr);

            //Multi Level Multi Grid
            MLMG mlmg(linop);

            //Solver parameters
            mlmg.setMaxIter(poisson_max_iter);
            mlmg.setVerbose(poisson_verbose);
            mlmg.setBottomVerbose(poisson_bottom_verbose);

            //Do solve
            mlmg.solve({&potential}, {&charge}, poisson_rel_tol, 0.0);

            
            potential.FillBoundary(geom.periodicity());

        //    //Find e field, gradient from cell centers to faces
        //    ComputeGrad(potential, efield, 0, 0, 1, geom);

             //Find e field, gradient from cell centers to faces
             ComputeCentredGrad(potential, efield, geom);
    }

    //Add external field on top, then fill boundaries, then setup BCs for peskin interpolation
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Add(efield[d], external[d], 0, 0, 1, 0);
        efield[d].FillBoundary(geom.periodicity());
        MultiFABElectricBC(efield[d], d, geom);
    }

}


