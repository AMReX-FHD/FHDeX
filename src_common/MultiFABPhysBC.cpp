#include "common_functions.H"
#include "MultiFABPhysBC.H"

// Fill 1 ghost cell for pressure based on the velocity boundary conditions.
// We test on bc_vel_lo/hi.  If they are slip or no-slip conditions
// we use homogeneous Neumann conditions
void MultiFABPhysBCPres(MultiFab& data, int sComp, int nComp,
                        const Geometry& geom) {

    if (geom.isAllPeriodic()) {
        return;
    }

    // Physical Domain
    Box dom(geom.Domain());

    // Send BCs to GPU
    GpuArray<int, AMREX_SPACEDIM> bc_lo{AMREX_D_DECL(common::bc_vel_lo[0],
                                                     common::bc_vel_lo[1],
                                                     common::bc_vel_lo[2])};
    GpuArray<int, AMREX_SPACEDIM> bc_hi{AMREX_D_DECL(common::bc_vel_hi[0],
                                                     common::bc_vel_hi[1],
                                                     common::bc_vel_hi[2])};

    for (MFIter mfi(data, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // one ghost cell
        Box bx = mfi.growntilebox(1);

        const Array4<Real>& data_fab = data.array(mfi);

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
        {
            physbc_pres_fab(tbx, dom, data_fab, bc_lo, bc_hi, sComp, nComp);
        });
    }
}

// Set the value of normal velocity on walls to zero
// Set the value of normal ghost cells to the inverse reflection of the interior
// The latter is needed for Peskin kernels and also
// to avoid intermediate NaN propagation in BDS
void MultiFABPhysBCDomainVel(MultiFab& vel, const Geometry& geom, int dim) {

    if (geom.isAllPeriodic()) {
        return;
    }

    // Physical Domain and make sure that the domain index type matches the
    // velocity index type
    Box dom(geom.Domain());
    dom.surroundingNodes(dim);

    int ng = vel.nGrow();
    
    // Send BCs to GPU
    GpuArray<int, AMREX_SPACEDIM> bc_lo{AMREX_D_DECL(common::bc_vel_lo[0],
                                                     common::bc_vel_lo[1],
                                                     common::bc_vel_lo[2])};
    GpuArray<int, AMREX_SPACEDIM> bc_hi{AMREX_D_DECL(common::bc_vel_hi[0],
                                                     common::bc_vel_hi[1],
                                                     common::bc_vel_hi[2])};

    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {

        // Select how much of the ghost region to fill
        Box bx = mfi.growntilebox(ng);

        const Array4<Real>& data_fab = vel.array(mfi);

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
        {
            physbc_domainvel_fab(tbx, dom, data_fab, bc_lo, bc_hi, dim);
        });
    }
}

/* MultiFABPhysBCMacVel

*/

void MultiFABPhysBCMacVel(MultiFab& vel, const Geometry& geom, int dim) {

    if (geom.isAllPeriodic()) {
        return;
    }

    MultiFABPhysBCMacVel(vel, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}

void MultiFABPhysBCMacVel(MultiFab& vel, const IntVect& dim_fill_ghost,
			  const Geometry& geom, int dim) {

    if (geom.isAllPeriodic()) {
        return;
    }

    // Physical Domain and make sure that the domain index type matches the
    // velocity index type
    Box dom(geom.Domain());
    dom.surroundingNodes(dim);

    // Send BCs to GPU
    GpuArray<int, AMREX_SPACEDIM> bc_lo{AMREX_D_DECL(common::bc_vel_lo[0],
                                                     common::bc_vel_lo[1],
                                                     common::bc_vel_lo[2])};
    GpuArray<int, AMREX_SPACEDIM> bc_hi{AMREX_D_DECL(common::bc_vel_hi[0],
                                                     common::bc_vel_hi[1],
                                                     common::bc_vel_hi[2])};

    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {

        // techinically you should only have to grow in the 2 directions that are NOT dim
        // since we are filling transverse velocity ghost cells
        IntVect ngv = vel.nGrowVect();
        Box bx      = mfi.growntilebox(ngv);

        const Array4<Real>& data_fab = vel.array(mfi);
        int n_comp = vel.nComp();

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
        {
            physbc_macvel_fab(tbx, dom, data_fab, bc_lo, bc_hi, 0, n_comp, dim);
        });
    }
}

/* MultiFABElectricBC

*/

void MultiFABElectricBC(MultiFab& data, const Geometry& geom) {
    MultiFABElectricBC(data, IntVect{AMREX_D_DECL(1, 1, 1)}, geom);
}



void MultiFABElectricBC(MultiFab& data, int seq_fill_ghost, const Geometry& geom) {

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABElectricBC(data, fill_ghost, geom);
}



void MultiFABElectricBC(MultiFab& data, const IntVect& dim_fill_ghost,
                        const Geometry& geom) {
    #if (AMREX_SPACEDIM==3 || AMREX_SPACEDIM==2)
    Box dom(geom.Domain());

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();
        fab_electricbc(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_BOX(dom),
                       BL_TO_FORTRAN_FAB(data[mfi]), data.nGrow(),
                       dim_fill_ghost.getVect());
    }
    #endif
}

/* MultiFABPotentialBC

   Note that this currently only operates on 1 layer of ghost cells.
   This routine fills ghost cells with the value extrapolated TO the ghost cell-center
   This is NOT the same as filling the ghost cell with the value on the boundary
   The Poisson solver needs a separate routine to fill ghost cells with the value ON
   the boundary for inhomogeneous Neumann and inhomogeneous Dirichlet; for this we
   use MultiFABPotentialBC_solver

*/
void MultiFABPotentialBC(MultiFab& data, const Geometry& geom) {
    MultiFABPotentialBC(data, IntVect{AMREX_D_DECL(1, 1, 1)}, geom);
}

void MultiFABPotentialBC(MultiFab& data, int seq_fill_ghost, const Geometry& geom) {

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPotentialBC(data, fill_ghost, geom);
}

void MultiFABPotentialBC(MultiFab& data, const IntVect& dim_fill_ghost,
                        const Geometry& geom) {
    #if (AMREX_SPACEDIM==3 || AMREX_SPACEDIM==2)
    Box dom(geom.Domain());

    const Real* dx = geom.CellSize();
    
    for (MFIter mfi(data); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();
        fab_potentialbc(BL_TO_FORTRAN_BOX(bx),
                        BL_TO_FORTRAN_BOX(dom),
                        BL_TO_FORTRAN_FAB(data[mfi]), data.nGrow(),
                        dim_fill_ghost.getVect(),
                        ZFILL(dx));
    }
    #endif
}

/* MultiFABPotentialBC_solver

   Note that this currently only operates on 1 layer of ghost cells.
   This routine fills ghost cells with the value ON the boundary.
   It works for inhomogeneous Neumann and inhomogeneous Dirichlet
   This routine is not to be confused with MultiFABPotentialBC, which fill ghost
   values extrapolated TO the ghost cell-center

*/
void MultiFABPotentialBC_solver(MultiFab& data, const Geometry& geom) {
    MultiFABPotentialBC_solver(data, IntVect{AMREX_D_DECL(1, 1, 1)}, geom);
}

void MultiFABPotentialBC_solver(MultiFab& data, int seq_fill_ghost, const Geometry& geom) {

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPotentialBC_solver(data, fill_ghost, geom);
}

void MultiFABPotentialBC_solver(MultiFab& data, const IntVect& dim_fill_ghost,
                        const Geometry& geom) {
    #if (AMREX_SPACEDIM==3 || AMREX_SPACEDIM==2)
    Box dom(geom.Domain());

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();
        fab_potentialbc_solver(BL_TO_FORTRAN_BOX(bx),
                               BL_TO_FORTRAN_BOX(dom),
                               BL_TO_FORTRAN_FAB(data[mfi]), data.nGrow(),
                               dim_fill_ghost.getVect());
    }
    #endif
}

/* MultiFABPhysBCCharge

*/
void MultiFABPhysBCCharge(MultiFab& data, const Geometry& geom) {
    MultiFABPhysBCCharge(data, IntVect{AMREX_D_DECL(1, 1, 1)}, geom);
}

void MultiFABPhysBCCharge(MultiFab& data, int seq_fill_ghost, const Geometry& geom) {

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBCCharge(data, fill_ghost, geom);
}

void MultiFABPhysBCCharge(MultiFab& data, const IntVect& dim_fill_ghost,
                        const Geometry& geom) {
    #if (AMREX_SPACEDIM==3 || AMREX_SPACEDIM==2)
    Box dom(geom.Domain());

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();
        fab_chargebc(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_BOX(dom),
                       BL_TO_FORTRAN_FAB(data[mfi]), data.nGrow(),
                       dim_fill_ghost.getVect());
    }
    #endif
}

/* MultiFABPhysBCDomainStress

*/
void MultiFABPhysBCDomainStress(MultiFab& stress,
                                const amrex::Geometry& geom, int dim) {
    MultiFABPhysBCDomainStress(stress, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}

void MultiFABPhysBCDomainStress(MultiFab& stress, const IntVect& dim_fill_ghost,
                                const Geometry& geom, int dim) {

    #if (AMREX_SPACEDIM==3 || AMREX_SPACEDIM==2)
    Box dom(geom.Domain());

    for (MFIter mfi(stress); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        fab_physbc_domainstress(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_BOX(dom),
                             BL_TO_FORTRAN_FAB(stress[mfi]), stress.nGrow(),
                             dim_fill_ghost.getVect(), &dim);
    }
    #endif
}

/* MultiFABPhysBCMacStress

*/
void MultiFABPhysBCMacStress(MultiFab& stress, const Geometry& geom, int dim) {
    MultiFABPhysBCMacStress(stress, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}

void MultiFABPhysBCMacStress(MultiFab& stress, const IntVect& dim_fill_ghost,
                             const Geometry& geom, int dim) {

    #if (AMREX_SPACEDIM==3 || AMREX_SPACEDIM==2)
    Box dom(geom.Domain());

    for (MFIter mfi(stress); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();
        fab_physbc_macstress(BL_TO_FORTRAN_BOX(bx),
                          BL_TO_FORTRAN_BOX(dom),
                          BL_TO_FORTRAN_FAB(stress[mfi]), stress.nGrow(),
                          dim_fill_ghost.getVect(), &dim);
    }
    #endif
}
