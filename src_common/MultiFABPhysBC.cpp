#include "common_functions.H"
#include "common_functions_F.H"



void MultiFABPhysBC(MultiFab & data, const Geometry & geom) {
    MultiFABPhysBC(data, IntVect{AMREX_D_DECL(1, 1, 1)}, geom);
}



void MultiFABPhysBC(MultiFab & data, int seq_fill_ghost, const Geometry & geom) {

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBC(data, fill_ghost, geom);
}



void MultiFABPhysBC(MultiFab & data, const IntVect & dim_fill_ghost, const Geometry & geom) {

    Box dom(geom.Domain());

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_physbc(BL_TO_FORTRAN_BOX(bx),
                   BL_TO_FORTRAN_BOX(dom),
                   BL_TO_FORTRAN_FAB(data[mfi]), data.nGrow(),
                   dim_fill_ghost.getVect());
    }
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



void MultiFABElectricBC(MultiFab & data, const IntVect & dim_fill_ghost, const Geometry & geom) {

    Box dom(geom.Domain());

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_electricbc(BL_TO_FORTRAN_BOX(bx),
                   BL_TO_FORTRAN_BOX(dom),
                   BL_TO_FORTRAN_FAB(data[mfi]), data.nGrow(),
                   dim_fill_ghost.getVect());
    }
}



void MultiFABPhysBCDomainVel(MultiFab & vel, const amrex::Geometry & geom, int dim) {
    MultiFABPhysBCDomainVel(vel, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}



void MultiFABPhysBCDomainVel(MultiFab & vel, int seq_fill_ghost, const Geometry & geom, int dim) {

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBCDomainVel(vel, fill_ghost, geom, dim);
}



void MultiFABPhysBCDomainVel(MultiFab & vel, const IntVect & dim_fill_ghost, const Geometry & geom, int dim) {

    Box dom(geom.Domain());

    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_physbc_domainvel(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_BOX(dom),
                             BL_TO_FORTRAN_FAB(vel[mfi]), vel.nGrow(),
                             dim_fill_ghost.getVect(), &dim);
    }
}



void MultiFABPhysBCMacVel(MultiFab & vel, const Geometry & geom, int dim) {
    MultiFABPhysBCMacVel(vel, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}



void MultiFABPhysBCMacVel(MultiFab & vel, int seq_fill_ghost, const Geometry & geom, int dim) {

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBCMacVel(vel, fill_ghost, geom, dim);
}



void MultiFABPhysBCMacVel(MultiFab & vel, const IntVect & dim_fill_ghost, const Geometry & geom, int dim) {

    Box dom(geom.Domain());

    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_physbc_macvel(BL_TO_FORTRAN_BOX(bx),
                          BL_TO_FORTRAN_BOX(dom),
                          BL_TO_FORTRAN_FAB(vel[mfi]), vel.nGrow(),
                          dim_fill_ghost.getVect(),&dim);
    }
}


void MultiFABPhysBCDomainStress(MultiFab & stress, const amrex::Geometry & geom, int dim) {
    MultiFABPhysBCDomainStress(stress, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}



void MultiFABPhysBCDomainStress(MultiFab & stress, int seq_fill_ghost, const Geometry & geom, int dim) {

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBCDomainStress(stress, fill_ghost, geom, dim);
}



void MultiFABPhysBCDomainStress(MultiFab & stress, const IntVect & dim_fill_ghost, const Geometry & geom, int dim) {

    Box dom(geom.Domain());

    for (MFIter mfi(stress); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_physbc_domainstress(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_BOX(dom),
                             BL_TO_FORTRAN_FAB(stress[mfi]), stress.nGrow(),
                             dim_fill_ghost.getVect(),&dim);
    }
}


void MultiFABPhysBCMacStress(MultiFab & stress, const Geometry & geom, int dim) {
    MultiFABPhysBCMacStress(stress, IntVect{AMREX_D_DECL(1,1,1)}, geom, dim);
}



void MultiFABPhysBCMacStress(MultiFab & stress, int seq_fill_ghost, const Geometry & geom, int dim) {

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBCMacStress(stress, fill_ghost, geom, dim);
}



void MultiFABPhysBCMacStress(MultiFab & stress, const IntVect & dim_fill_ghost, const Geometry & geom, int dim) {

    Box dom(geom.Domain());

    for (MFIter mfi(stress); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_physbc_macstress(BL_TO_FORTRAN_BOX(bx),
                          BL_TO_FORTRAN_BOX(dom),
                          BL_TO_FORTRAN_FAB(stress[mfi]), stress.nGrow(),
                          dim_fill_ghost.getVect(), &dim);
    }
}
