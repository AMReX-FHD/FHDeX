#include "reactDiff_functions.H"

#include "AMReX_MLMG.H"
#include <AMReX_MLABecLaplacian.H>

// (I - (dt_fac) div D_k grad) n = rhs

void ImplicitDiffusion(MultiFab& n_old,
                       MultiFab& n_new,
                       const MultiFab& rhs,
                       const std::array< MultiFab, AMREX_SPACEDIM >& diff_coef_face,
                       const Geometry& geom,
                       const Real& dt_fac,
                       const Real& time) {

    BoxArray ba = n_old.boxArray();
    DistributionMapping dmap = n_old.DistributionMap();

    // fill n ghost cells
    n_old.FillBoundary(geom.periodicity());
    MultiFabPhysBC(n_old, geom, 0, nspecies, SPEC_BC_COMP, time);

    LPInfo info;

    // operator of the form (ascalar * acoef - bscalar div bcoef grad) phi
    MLABecLaplacian mlabec({geom}, {ba}, {dmap}, info);
    mlabec.setMaxOrder(2);

    // store one component at a time and take L(phi) one component at a time
    MultiFab phi     (ba,dmap,1,1);
    MultiFab rhs_comp(ba,dmap,1,0);

    // build array of boundary conditions needed by MLABecLaplacian
    std::array<LinOpBCType, AMREX_SPACEDIM> lo_mlmg_bc;
    std::array<LinOpBCType, AMREX_SPACEDIM> hi_mlmg_bc;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (bc_mass_lo[idim] == -1 || bc_mass_hi[idim] == -1) {
            if ( !(bc_mass_lo[idim] == -1 && bc_mass_hi[idim] == -1) ) {
                Abort("Both bc_mass_lo and bc_mass_hi must be periodic in a given direction if the other one is");
            }
            lo_mlmg_bc[idim] = LinOpBCType::Periodic;
            hi_mlmg_bc[idim] = LinOpBCType::Periodic;
        }

        if (bc_mass_lo[idim] == 0) {
            lo_mlmg_bc[idim] = LinOpBCType::inhomogNeumann;
        } else if (bc_mass_lo[idim] == 1) {
            lo_mlmg_bc[idim] = LinOpBCType::Dirichlet;
        } else if (bc_mass_lo[idim] != -1) {
            Abort("Invalid bc_mass_lo");
        }

        if (bc_mass_hi[idim] == 0) {
            hi_mlmg_bc[idim] = LinOpBCType::inhomogNeumann;
        } else if (bc_mass_hi[idim] == 1) {
            hi_mlmg_bc[idim] = LinOpBCType::Dirichlet;
        } else if (bc_mass_hi[idim] != -1) {
            Abort("Invalid bc_mass_hi");
        }
    }

    mlabec.setDomainBC(lo_mlmg_bc,hi_mlmg_bc);

    // storage for acoeff and bcoeff in
    // (ascalar * acoeff I - bscalar div bcoeff grad) phi = rhs
    MultiFab acoef(ba,dmap,1,0);
    std::array< MultiFab, AMREX_SPACEDIM > bcoef;
    AMREX_D_TERM(bcoef[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 bcoef[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 bcoef[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

    // set ascalar and bscalar to 1
    mlabec.setScalars(1., 1.);

    // acoeff = 1
    acoef.setVal(1.);
    mlabec.setACoeffs(0, acoef);

    // set bcoeff to dt_fac * D_i
    for (int i=0; i<nspecies; ++i) {

        // load D_fick for species i into bcoef
        // then multiply by dt_fac
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Copy(bcoef[d],diff_coef_face[d],i,0,1,0);
            bcoef[d].mult(dt_fac);
        }
        mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoef));

        // copy in n_old, including ghost cells for boundary conditions, into phi as an initial guess
        // copy in rhs_comp into rhs
        MultiFab::Copy(phi,n_old,i,0,1,1);
        MultiFab::Copy(rhs_comp,rhs,i,0,1,0);

        // tell the operator what the numerical values for physical boundary conditions are
        mlabec.setLevelBC(0, &phi);

        MLMG mlmg(mlabec);

        // solver parameters
        mlmg.setMaxIter(100);
        mlmg.setVerbose(0);
        mlmg.setBottomVerbose(0);

        // do solve
        mlmg.solve({&phi}, {&rhs_comp}, 1.e-10, 0.0);

        MultiFab::Copy(n_new,phi,0,i,1,0);

    }

    n_new.FillBoundary(geom.periodicity());
    MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);

}
