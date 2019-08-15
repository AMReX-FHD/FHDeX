#include "common_functions.H"
#include "common_functions_F.H"



void MultiFABPhysBC(MultiFab & data, const Geometry & geom) {

    if (geom.isAllPeriodic()) {
      return;
    }  
  
    MultiFABPhysBC(data, IntVect{AMREX_D_DECL(1, 1, 1)}, geom);
}



void MultiFABPhysBC(MultiFab & data, int seq_fill_ghost, const Geometry & geom) {

    if (geom.isAllPeriodic()) {
      return;
    }

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBC(data, fill_ghost, geom);
}



void MultiFABPhysBC(MultiFab & data, const IntVect & dim_fill_ghost,
                    const Geometry & geom) {

    if (geom.isAllPeriodic()) {
      return;
    }
    
#if (AMREX_SPACEDIM==2 || AMREX_SPACEDIM==3)
    Box dom(geom.Domain());

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_physbc(BL_TO_FORTRAN_BOX(bx),
                   BL_TO_FORTRAN_BOX(dom),
                   BL_TO_FORTRAN_FAB(data[mfi]), data.nGrow(),
                   dim_fill_ghost.getVect());
    }
#endif
}

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

//Note that potential currently only operates on 1 layer of ghost cells.

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

void MultiFABChargeBC(MultiFab & data, const Geometry & geom) {
    MultiFABChargeBC(data, IntVect{AMREX_D_DECL(1, 1, 1)}, geom);
}



void MultiFABChargeBC(MultiFab & data, int seq_fill_ghost, const Geometry & geom) {

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABChargeBC(data, fill_ghost, geom);
}



void MultiFABChargeBC(MultiFab & data, const IntVect & dim_fill_ghost,
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



void MultiFABPhysBCDomainVel(MultiFab & vel, const amrex::Geometry & geom, int dim) {

    if (geom.isAllPeriodic()) {
      return;
    }
    
    MultiFABPhysBCDomainVel(vel, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}



void MultiFABPhysBCDomainVel(MultiFab & vel, int seq_fill_ghost,
                             const Geometry & geom, int dim) {

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
    
#if (AMREX_SPACEDIM==3 || AMREX_SPACEDIM==2)
    Box dom(geom.Domain());

    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_physbc_domainvel(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_BOX(dom),
                             BL_TO_FORTRAN_FAB(vel[mfi]), vel.nGrow(),
                             dim_fill_ghost.getVect(), &dim);
    }
#endif
}



void MultiFABPhysBCMacVel(MultiFab & vel, const Geometry & geom, int dim) {

    if (geom.isAllPeriodic()) {
      return;
    }
    
    MultiFABPhysBCMacVel(vel, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}



void MultiFABPhysBCMacVel(MultiFab & vel, int seq_fill_ghost,
                          const Geometry & geom, int dim) {

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

    #if (AMREX_SPACEDIM==3 || AMREX_SPACEDIM==2)
    Box dom(geom.Domain());

    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_physbc_macvel(BL_TO_FORTRAN_BOX(bx),
                          BL_TO_FORTRAN_BOX(dom),
                          BL_TO_FORTRAN_FAB(vel[mfi]), vel.nGrow(),
                          dim_fill_ghost.getVect(),&dim);
    }
    #endif
}


void MultiFABPhysBCDomainStress(MultiFab & stress,
                                const amrex::Geometry & geom, int dim) {
    MultiFABPhysBCDomainStress(stress, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}



void MultiFABPhysBCDomainStress(MultiFab & stress, int seq_fill_ghost,
                                const Geometry & geom, int dim) {

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


void MultiFABPhysBCMacStress(MultiFab & stress, const Geometry & geom, int dim) {
    MultiFABPhysBCMacStress(stress, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}



void MultiFABPhysBCMacStress(MultiFab & stress, int seq_fill_ghost,
                             const Geometry & geom, int dim) {

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
