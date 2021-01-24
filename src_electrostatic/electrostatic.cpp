#include "electrostatic.H"
#include "common_functions.H"
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

using namespace amrex;

void esSolve(MultiFab& potential, MultiFab& charge,
             std::array< MultiFab, AMREX_SPACEDIM >& efieldCC,
             const std::array< MultiFab, AMREX_SPACEDIM >& external, const Geometry geom)
{
    AMREX_D_TERM(efieldCC[0].setVal(0);,
                 efieldCC[1].setVal(0);,
                 efieldCC[2].setVal(0););

    if(es_tog==1 || es_tog==3)
    {

        const BoxArray& ba = charge.boxArray();
        const DistributionMapping& dmap = charge.DistributionMap();
        Box dom(geom.Domain());
        
        // Set beta and permittivity_fc equal to permittivity
        std::array< MultiFab, AMREX_SPACEDIM > beta;
        std::array< MultiFab, AMREX_SPACEDIM > permittivity_fc;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            beta[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
            permittivity_fc[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
        }

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            beta[d].setVal(permittivity);
            permittivity_fc[d].setVal(permittivity);
        }

        // Copy charge in the RHS of Poisson Solver
        MultiFab rhs(ba,dmap,1,0);
        MultiFab::Copy(rhs,charge,0,0,1,0);

        // inhomogeneous permittivity on lower y-wall (eps = 0 on dielectric Dirichlet wall) -- set beta = 0
        if (zero_eps_on_wall_type) {
            
            for (MFIter mfi(beta[1]); mfi.isValid(); ++mfi) {

                Box bx = mfi.tilebox();

                const Array4<Real>& data = beta[1].array(mfi);

                if (bx.smallEnd(1) <= dom.smallEnd(1)) {

                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {   

                        if (j == dom.smallEnd(1)) {
                            if ((i <= (int)Math::floor(zero_eps_wall_left_end*(dom.bigEnd(0)+1))) and (k <= (int)Math::floor(zero_eps_wall_left_end*(dom.bigEnd(2)+1)))) { // left side (low x/z)
                                data(i,j,k) = 0.;
                            }
                        }

                        if (j == dom.smallEnd(1)) {
                            if ((i >= (int)Math::floor(zero_eps_wall_right_start*(dom.bigEnd(0)+1))) and (k >= (int)Math::floor(zero_eps_wall_right_start*(dom.bigEnd(2)+1)))) { // right side (high x/z)
                                data(i,j,k) = 0.;
                            }
                        }

                    });

                }

            }

        }

        // compute epsilon*external (where epsilon is stored in permittivity_fc)
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            MultiFab::Multiply(permittivity_fc[i],external[i],0,0,1,0);
        }

        // compute div (epsilon*E_ext) and subtract it from the solver rhs
        // only needed for spatially varying epsilon or external field
        ComputeDiv(rhs,permittivity_fc,0,0,1,geom,-1.);

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

        //create solver opject
        MLABecLaplacian linop({geom}, {ba}, {dmap});
        //MLPoisson linop({geom}, {ba}, {dmap});
 
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

        // set alpha=0, beta=1 (will overwrite beta with epsilon next)
        linop.setScalars(0.0, 1.0);

        // set beta=epsilon
        linop.setBCoeffs(0, amrex::GetArrOfConstPtrs(beta));

        //Multi Level Multi Grid
        MLMG mlmg(linop);

        //Solver parameters
        mlmg.setMaxIter(poisson_max_iter);
        mlmg.setVerbose(poisson_verbose);
        mlmg.setBottomVerbose(poisson_bottom_verbose);
        
        //Do solve -- this give us \phi (as oppowed to -\phi using MLPoisson)
        mlmg.solve({&potential}, {&rhs}, poisson_rel_tol, 0.0);
            
        potential.FillBoundary(geom.periodicity());
        // set ghost cell values so electric field is calculated properly
        // the ghost cells will hold the values extrapolated to the ghost CC
        MultiFabPotentialBC(potential, geom); 

        //Find e field, gradient from cell centers to faces
        //This routine calculates \nabla\phi; we need -\nabla\phi
        //This is done by setting factor to -1
        ComputeCentredGrad(potential, efieldCC, geom, -1);
    }

    //Add external field on top, then fill boundaries, then setup BCs for peskin interpolation
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Add(efieldCC[d], external[d], 0, 0, 1, efieldCC[d].nGrow());
        efieldCC[d].FillBoundary(geom.periodicity());
        MultiFabElectricBC(efieldCC[d], geom);
    }

}


