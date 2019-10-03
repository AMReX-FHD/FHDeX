#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

#include "MultiFABPhysBC.H"

/* MultiFABPhysBCPres

   Fill ghost cell based on 'pressure' boundary conditions
   We test on bc_vel_lo/hi.  If they are slip or no-slip conditions
   we use homogeneous Neumann conditions

*/

// this version fills ghost cells in all spatial directions
// you can get into trouble accessing uninitilized data using this with mixed bc
// types (wall/wall) at corners when ghost cell data enters uninitialized
void MultiFABPhysBCPres(MultiFab & data, const Geometry & geom) {

    if (geom.isAllPeriodic()) {
        return;
    }

    MultiFABPhysBCPres(data, IntVect{AMREX_D_DECL(1, 1, 1)}, geom);
}

// this version fills ghost cells in spatial directions 0:seq_fill_ghost
void MultiFABPhysBCPres(MultiFab & data, int seq_fill_ghost, const Geometry & geom) {

    if (geom.isAllPeriodic()) {
        return;
    }

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBCPres(data, fill_ghost, geom);
}

// this version fills ghost cells in an arbitrary number of spatial directions
void MultiFABPhysBCPres(MultiFab & data, const IntVect & dim_fill_ghost,
                        const Geometry & geom) {

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

        // Select how much of the ghost region to fill
        IntVect ngv = data.nGrowVect() * dim_fill_ghost;
        Box bx      = mfi.growntilebox(ngv);

        const Array4<Real> & data_fab = data.array(mfi);
        int n_comp = data.nComp();

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
        {
            physbc_pres_fab(tbx, dom, data_fab, bc_lo, bc_hi, 0, n_comp);
        });
    }
}

/* MultiFABPhysBCDomainVel

*/

void MultiFABPhysBCDomainVel(MultiFab & vel, const amrex::Geometry & geom, int dim) {

    if (geom.isAllPeriodic()) {
        return;
    }

    MultiFABPhysBCDomainVel(vel, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}

void MultiFABPhysBCDomainVel(MultiFab & vel, int seq_fill_ghost,
                             const Geometry & geom, int dim) {


//    Abort("MultiFABPhysBC.cpp: Do not call this instance of MultiFABPhysBCDomainVel");
    
    if (geom.isAllPeriodic()) {
        return;
    }

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBCDomainVel(vel, fill_ghost, geom, dim);
}

void MultiFABPhysBCDomainVel(MultiFab & vel, const IntVect & dim_fill_ghost,
                             const Geometry & geom, int dim) {

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

        // Select how much of the ghost region to fill
        IntVect ngv = vel.nGrowVect() * dim_fill_ghost;
        Box bx      = mfi.growntilebox(ngv);

        const Array4<Real> & data_fab = vel.array(mfi);
        int n_comp = vel.nComp();

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
        {
            physbc_domainvel_fab(tbx, dom, data_fab, bc_lo, bc_hi, 0, n_comp, dim);
        });
    }
}

/* MultiFABPhysBCMacVel

*/

void MultiFABPhysBCMacVel(MultiFab & vel, const Geometry & geom, int dim) {

    if (geom.isAllPeriodic()) {
        return;
    }

    MultiFABPhysBCMacVel(vel, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}

void MultiFABPhysBCMacVel(MultiFab & vel, int seq_fill_ghost,
                          const Geometry & geom, int dim) {


//    Abort("MultiFABPhysBC.cpp: Do not call this instance of MultiFABPhysBCMacVel");
    
    if (geom.isAllPeriodic()) {
        return;
    }

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBCMacVel(vel, fill_ghost, geom, dim);
}

void MultiFABPhysBCMacVel(MultiFab & vel, const IntVect & dim_fill_ghost,
			  const Geometry & geom, int dim) {

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

        const Array4<Real> & data_fab = vel.array(mfi);
        int n_comp = vel.nComp();

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
        {
            physbc_macvel_fab(tbx, dom, data_fab, bc_lo, bc_hi, 0, n_comp, dim);
        });
    }
}

/* MultiFABElectricBC

*/

void MultiFABElectricBC(MultiFab & data, const Geometry & geom) {
    MultiFABElectricBC(data, IntVect{AMREX_D_DECL(1, 1, 1)}, geom);
}



void MultiFABElectricBC(MultiFab & data, int seq_fill_ghost, const Geometry & geom) {

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABElectricBC(data, fill_ghost, geom);
}



void MultiFABElectricBC(MultiFab & data, const IntVect & dim_fill_ghost,
                        const Geometry & geom) {
    #if (AMREX_SPACEDIM==3 || AMREX_SPACEDIM==2)
    Box dom(geom.Domain());

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_electricbc(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_BOX(dom),
                       BL_TO_FORTRAN_FAB(data[mfi]), data.nGrow(),
                       dim_fill_ghost.getVect());
    }
    #endif
}

/* MultiFABPotentialBC

   Note that potential currently only operates on 1 layer of ghost cells.

*/
void MultiFABPotentialBC(MultiFab & data, const Geometry & geom) {
    MultiFABPotentialBC(data, IntVect{AMREX_D_DECL(1, 1, 1)}, geom);
}

void MultiFABPotentialBC(MultiFab & data, int seq_fill_ghost, const Geometry & geom) {

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPotentialBC(data, fill_ghost, geom);
}

void MultiFABPotentialBC(MultiFab & data, const IntVect & dim_fill_ghost,
                        const Geometry & geom) {
    #if (AMREX_SPACEDIM==3 || AMREX_SPACEDIM==2)
    Box dom(geom.Domain());

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_potentialbc(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_BOX(dom),
                       BL_TO_FORTRAN_FAB(data[mfi]), data.nGrow(),
                       dim_fill_ghost.getVect());
    }
    #endif
}

/* MultiFABPhysBCCharge

*/
void MultiFABPhysBCCharge(MultiFab & data, const Geometry & geom) {
    MultiFABPhysBCCharge(data, IntVect{AMREX_D_DECL(1, 1, 1)}, geom);
}

void MultiFABPhysBCCharge(MultiFab & data, int seq_fill_ghost, const Geometry & geom) {

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBCCharge(data, fill_ghost, geom);
}

void MultiFABPhysBCCharge(MultiFab & data, const IntVect & dim_fill_ghost,
                        const Geometry & geom) {
    #if (AMREX_SPACEDIM==3 || AMREX_SPACEDIM==2)
    Box dom(geom.Domain());

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_chargebc(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_BOX(dom),
                       BL_TO_FORTRAN_FAB(data[mfi]), data.nGrow(),
                       dim_fill_ghost.getVect());
    }
    #endif
}

/* MultiFABPhysBCDomainStress

*/
void MultiFABPhysBCDomainStress(MultiFab & stress,
                                const amrex::Geometry & geom, int dim) {
    MultiFABPhysBCDomainStress(stress, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}

void MultiFABPhysBCDomainStress(MultiFab & stress, int seq_fill_ghost,
                                const Geometry & geom, int dim) {

    Abort("MultiFABPhysBC.cpp: Do not call this instance of MultiFABPhysBCDomainStress");
        
    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBCDomainStress(stress, fill_ghost, geom, dim);
}

void MultiFABPhysBCDomainStress(MultiFab & stress, const IntVect & dim_fill_ghost,
                                const Geometry & geom, int dim) {

    #if (AMREX_SPACEDIM==3 || AMREX_SPACEDIM==2)
    Box dom(geom.Domain());

        

    for (MFIter mfi(stress); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();

        fab_physbc_domainstress(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_BOX(dom),
                             BL_TO_FORTRAN_FAB(stress[mfi]), stress.nGrow(),
                             dim_fill_ghost.getVect(),&dim);
    }
    #endif
}

/* MultiFABPhysBCMacStress

*/
void MultiFABPhysBCMacStress(MultiFab & stress, const Geometry & geom, int dim) {
    MultiFABPhysBCMacStress(stress, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}

void MultiFABPhysBCMacStress(MultiFab & stress, int seq_fill_ghost,
                             const Geometry & geom, int dim) {

    Abort("MultiFABPhysBC.cpp: Do not call this instance of MultiFABPhysBCMacStress");
    
    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBCMacStress(stress, fill_ghost, geom, dim);
}

void MultiFABPhysBCMacStress(MultiFab & stress, const IntVect & dim_fill_ghost,
                             const Geometry & geom, int dim) {

    #if (AMREX_SPACEDIM==3 || AMREX_SPACEDIM==2)
    Box dom(geom.Domain());

    for (MFIter mfi(stress); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_physbc_macstress(BL_TO_FORTRAN_BOX(bx),
                          BL_TO_FORTRAN_BOX(dom),
                          BL_TO_FORTRAN_FAB(stress[mfi]), stress.nGrow(),
                          dim_fill_ghost.getVect(), &dim);
    }
    #endif
}
