#include "reactDiff_functions.H"

#include "AMReX_MLMG.H"
#include <AMReX_MLABecLaplacian.H>

void DiffusiveNFluxdiv(MultiFab& n_in,
                       MultiFab& diff_fluxdiv,
                       const std::array< MultiFab, AMREX_SPACEDIM >& diff_coef_face,
                       const Geometry& geom,
                       const Real& time) {

    // single cell case set diffusive mass fluxdiv to zero and return
    long cell_count = (AMREX_SPACEDIM==2) ? n_cells[0]*n_cells[1] : n_cells[0]*n_cells[1]*n_cells[2];
    if (cell_count == 1) {
        diff_fluxdiv.setVal(0.);
        return;
    }

    // fill n ghost cells
    n_in.FillBoundary(geom.periodicity());
    MultiFabPhysBC(n_in, geom, 0, nspecies, SPEC_BC_COMP, time);

    BoxArray ba = n_in.boxArray();
    DistributionMapping dmap = n_in.DistributionMap();

    // don't need to set much here for explicit evaluations
    LPInfo info;

    // operator of the form (ascalar * acoef - bscalar div bcoef grad) phi
    MLABecLaplacian mlabec({geom}, {ba}, {dmap}, info);
    mlabec.setMaxOrder(2);

    // store one component at a time and take L(phi) one component at a time
    MultiFab phi (ba,dmap,1,1);
    MultiFab Lphi(ba,dmap,1,0);

    MultiFab acoef(ba,dmap,1,0);
    std::array< MultiFab, AMREX_SPACEDIM > bcoef;
    AMREX_D_TERM(bcoef[0].define(convert(ba,nodal_flag_x), dmap, 1, 0);,
                 bcoef[1].define(convert(ba,nodal_flag_y), dmap, 1, 0);,
                 bcoef[2].define(convert(ba,nodal_flag_z), dmap, 1, 0););

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

    // set acoeff to 0and bcoeff to -1
    mlabec.setScalars(0., -1.);

    acoef.setVal(0.);
    mlabec.setACoeffs(0, acoef);

    for (int i=0; i<nspecies; ++i) {

        // copy ith component of n_in into phi, including ghost cells
        MultiFab::Copy(phi,n_in,i,0,1,1);

        // load D_fick for species i into bcoef
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            MultiFab::Copy(bcoef[d],diff_coef_face[d],i,0,1,0);
        }
        mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoef));

        MLMG mlmg(mlabec);

        mlmg.apply({&Lphi},{&phi});

        MultiFab::Copy(diff_fluxdiv,Lphi,0,i,1,0);

    }

}
