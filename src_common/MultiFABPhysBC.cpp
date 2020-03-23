#include "common_functions.H"
#include "MultiFABPhysBC.H"

// Fill 1 ghost cell for pressure based on the velocity boundary conditions.
// We test on bc_vel_lo/hi.  If they are slip or no-slip conditions
// we use homogeneous Neumann conditions
void MultiFABPhysBCPres(MultiFab& data, const Geometry& geom) {

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
            physbc_pres_fab(tbx, dom, data_fab, bc_lo, bc_hi);
        });
    }
}

// Set the value of normal velocity on walls to zero
// Set the value of normal ghost cells to the inverse reflection of the interior
// We fill all the ghost cells - they are needed for Perskin kernels and
// to avoid intermediate NaN propagation in BDS
void MultiFABPhysBCDomainVel(MultiFab& vel, const Geometry& geom, int dim) {

    if (geom.isAllPeriodic()) {
        return;
    }

    // grow the domain box by 1 in the dim direction
    Box dom(geom.Domain());

    int ng = vel.nGrow();

    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {

        Box bx = mfi.growntilebox(ng);

        const Array4<Real>& data = vel.array(mfi);

        //___________________________________________________________________________
        // Apply x-physbc to data

        // lo-x faces
        // dim == 0 means we are doing x-velocity on x-faces
        // bc_vel check is to see if we have a wall bc
        // bx/dom comparison is to see if the grid touches a wall
        if ((dim == 0) && (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) && (bx.smallEnd(0) <= dom.smallEnd(0))) {

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) {
                    // set ghost cells to negative of interior value
                    data(i,j,k) = -data(-i,j,k);
                }           
                else if (i == dom.smallEnd(0)){
                    // set normal velocity on boundary to zero
                    data(i,j,k) = 0.;
                }
            });
        }

        // hi-x faces
        if ((dim == 0) && (bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) && (bx.bigEnd(0) >= dom.bigEnd(0)+1)) {

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {        
                if (i > dom.bigEnd(0)+1) {
                    data(i,j,k) = -data(2*dom.bigEnd(0)+2-i,j,k);
                }           
                else if (i == dom.bigEnd(0)+1) {
                    data(i,j,k) = 0.;
                }
            });
        }

#if (AMREX_SPACEDIM >= 2)
        
        //___________________________________________________________________________
        // Apply y-physbc to data

        // lo-y faces
        if ((dim == 1) && (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) && (bx.smallEnd(1) <= dom.smallEnd(1))) {

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {        
                if (j < dom.smallEnd(1)) {
                    data(i,j,k) = -data(i,-j,k);
                }           
                else if (j == dom.smallEnd(1)) {
                    data(i,j,k) = 0.;
                }
            });
        }

        // hi-y faces
        if ((dim == 1) && (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) && (bx.bigEnd(1) >= dom.bigEnd(1)+1)) {

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {        
                if (j > dom.bigEnd(1)+1) {
                    data(i,j,k) = -data(i,2*dom.bigEnd(1)+2-j,k);
                }           
                else if (j == dom.bigEnd(1)+1) {
                    data(i,j,k) = 0.;
                }
            });
        }
#endif
        
#if (AMREX_SPACEDIM >= 3)

        //___________________________________________________________________________
        // Apply z-physbc to data

        // lo-z faces
        if ((dim == 2) && (bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) && (bx.smallEnd(2) <= dom.smallEnd(2))) {

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {        
                if (k < dom.smallEnd(2)) {
                    data(i,j,k) = -data(i,j,-k);
                }           
                else if (k == dom.smallEnd(2)) {
                    data(i,j,k) = 0.;
                }
            });
        }

        // hi-z faces
        if ((dim == 2) && (bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) && (bx.bigEnd(2) >= dom.bigEnd(2)+1)) {

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {        
                if (k > dom.bigEnd(2)+1) {
                    data(i,j,k) = -data(i,j,2*dom.bigEnd(2)+2-k);
                }           
                else if (k == dom.bigEnd(2)+1) {
                    data(i,j,k) = 0.;
                }
            });
        }
#endif
        
    } // end MFIter
}

// Set the value of tranverse ghost cells to +/- the reflection of the interior
// (+ for slip walls, - for no-slip)
// We fill all the ghost cells - they are needed for Perskin kernels
void MultiFABPhysBCMacVel(MultiFab& vel, const Geometry& geom, int dim) {

    if (geom.isAllPeriodic()) {
        return;
    }

    Box dom(geom.Domain());

    int ng = vel.nGrow();

    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {

        Box bx = mfi.growntilebox(ng);

        const Array4<Real>& data = vel.array(mfi);

        //___________________________________________________________________________
        // Apply x-physbc to data
        
        // lo-x faces
        // dim != 0 means we are doing y and z-velocity on x-faces
        // bc_vel check is to see if we have a wall bc
        // bx/dom comparison is to see if the grid touches a wall 
        if ((dim != 0) && (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) && (bx.smallEnd(0) < dom.smallEnd(0))) {
            Real fac = (bc_vel_lo[0] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) {
                    data(i,j,k) = fac*data(-i-1,j,k);
                }
            });
        }

        if ((dim != 0) && (bc_vel_lo[0] == 1 || bc_vel_hi[0] == 2) && (bx.bigEnd(0) > dom.bigEnd(0))) {
            Real fac = (bc_vel_hi[0] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > dom.bigEnd(0)) {
                    data(i,j,k) = fac*data(2*dom.bigEnd(0)-i+1,j,k);
                }
            });
        }
#if (AMREX_SPACEDIM >= 2)

        //___________________________________________________________________________
        // Apply y-physbc to data
        if ((dim != 1) && (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) && (bx.smallEnd(1) < dom.smallEnd(1))) {
            Real fac = (bc_vel_lo[1] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) {
                    data(i,j,k) = fac*data(i,-j-1,k);
                }
            });
        }

        if ((dim != 1) && (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) && (bx.bigEnd(1) > dom.bigEnd(1))) {
            Real fac = (bc_vel_hi[1] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j > dom.bigEnd(1)) {
                    data(i,j,k) = fac*data(i,2*dom.bigEnd(1)-j+1,k);
                }
            });
        }
#endif

        //___________________________________________________________________________
        // Apply z-physbc to data
#if (AMREX_SPACEDIM >= 3)
        if ((dim != 2) && (bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) && (bx.smallEnd(2) < dom.smallEnd(2))) {
            Real fac = (bc_vel_lo[2] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) {
                    data(i,j,k) = fac*data(i,j,-k-1);
                }
            });
        }

        if ((dim != 2) && (bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) && (bx.bigEnd(2) > dom.bigEnd(2))) {
            Real fac = (bc_vel_hi[2] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k > dom.bigEnd(2)) {
                    data(i,j,k) = fac*data(i,j,2*dom.bigEnd(2)-k+1);
                }
            });
        }
#endif
        
    } // end MFIter
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
