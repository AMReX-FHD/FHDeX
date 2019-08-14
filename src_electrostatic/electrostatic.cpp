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

    if(es_tog==1 || es_tog==3)
    {

            LinOpBCType lobc[3];
            LinOpBCType hibc[3];

                Print() << "GEOMIN: " << geom << "\n";
                Print() << "BC_ES: " << bc_es_lo[0] << ", " << bc_es_lo[1] << ", " << bc_es_lo[2] << "\n";
                Print() << "BC: " << bc_lo[0] << ", " << bc_lo[1] << ", " << bc_lo[2] << "\n";
// while(true);

            for (int i=0; i<AMREX_SPACEDIM; ++i) {
                if (bc_es_lo[i] == -1 && bc_es_hi[i] == -1) {
                    lobc[i] = LinOpBCType::Periodic;
                    hibc[i] = LinOpBCType::Periodic;
                                Print() << "Here!\n";
                }
                if(bc_es_lo[i] == 2)
                {
                    lobc[i] = LinOpBCType::Neumann;
                }
                if(bc_es_hi[i] == 2)
                {
                    hibc[i] = LinOpBCType::Neumann;
                }
                if(bc_es_lo[i] == 1)
                {                 
                    lobc[i] = LinOpBCType::Dirichlet;
                }
                if(bc_es_hi[i] == 1)
                {
                    hibc[i] = LinOpBCType::Dirichlet;
                }
            }

            //MultiFABChargeBC(charge, geom); //Adjust spread charge distribtion near boundaries from 

            Geometry geomT;

            IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
            IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
            Box domain(dom_lo, dom_hi);

            BoxArray ba;
            ba.define(domain);

            ba.maxSize(IntVect(max_grid_size));

            RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

            Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default

            for (int i=0; i<AMREX_SPACEDIM; ++i) {
                if (bc_es_lo[i] == -1 && bc_es_hi[i] == -1) {
                    is_periodic[i] = 1;
                  
                }
            }

            Print() << "Is periodic: " << is_periodic[0] << ", " << is_periodic[1] <<", " <<  is_periodic[2] << "\n";

            geomT.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
 

            // how boxes are distrubuted among MPI processes
            // AJN needs to be fi
            DistributionMapping dmap(ba);

            MultiFab testCharge(ba, dmap, 1, 2);
            MultiFab testPot(ba, dmap, 1, 2);
            testCharge.setVal(1);
            testPot.setVal(0);


            //const BoxArray& ba = charge.boxArray();
            //const DistributionMapping& dmap = charge.DistributionMap();

            //create solver opject
            MLPoisson linop({geomT}, {ba}, {dmap});

            Print() << "LoBC: " << (int)lobc[1] << "\n";
            Print() << "HiBC: " << (int)hibc[1] << "\n";

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

            Print() << "Pspecs: " << poisson_max_iter << ", " << poisson_verbose << ", " << poisson_bottom_verbose << ", " << poisson_rel_tol << "\n";
            //Do solve
            mlmg.solve({&testPot}, {&testCharge}, poisson_rel_tol, 0.0);

            
            potential.FillBoundary(geom.periodicity());
            MultiFABPotentialBC(potential, geom); //set ghost cell values so electric field is calculated properly

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


