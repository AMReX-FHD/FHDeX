#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"


void MultiFABPhysBCPres(MultiFab & data, const Geometry & geom) {

    // Call this if the ghost cells are already `setVal` to zero.

    if (geom.isAllPeriodic()) {
        return;
    }

    MultiFABPhysBCPres(data, IntVect{AMREX_D_DECL(1, 1, 1)}, geom);
}



void MultiFABPhysBCPres(MultiFab & data, int seq_fill_ghost, const Geometry & geom) {

    if (geom.isAllPeriodic()) {
        return;
    }

    IntVect fill_ghost{AMREX_D_DECL(0, 0, 0)};
    for(int i=0; i<=seq_fill_ghost; i++)
        fill_ghost[i] = 1;

    MultiFABPhysBCPres(data, fill_ghost, geom);
}



AMREX_GPU_HOST_DEVICE
inline void physbc_pres_fab(const Box & tbx,
                            const Box & dom,
                            const Array4<Real> & data,
                            const GpuArray<int, AMREX_SPACEDIM> & bc_lo,
                            const GpuArray<int, AMREX_SPACEDIM> & bc_hi,
                            int start_cmp, int n_cmp) {

    //___________________________________________________________________________
    // Total work region => the loops below will actually only iterate over
    // cells between tbx and dom
    const Dim3 tlo    = amrex::lbound(tbx);
    const Dim3 thi    = amrex::ubound(tbx);
    const Dim3 dom_lo = amrex::lbound(dom);
    const Dim3 dom_hi = amrex::ubound(dom);


    //___________________________________________________________________________
    // Apply x-physbc to data
    if (((bc_lo[0] == 1) || (bc_lo[0] == 2)) && (tlo.x <= dom_lo.x)) {
        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = tlo.z; k <= thi.z; ++k) {
                for (int j = tlo.y; j <= thi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = tlo.x; i < dom_lo.x; ++i) {
                        int offset = dom_lo.x - i;
                        int i_real = dom_lo.x + offset - 1;
                        data(i, j, k, n) = data(i_real, j, k, n);
                    }
                }
            }
        }
    }

    if (((bc_hi[0] == 1) || (bc_hi[0] == 2)) && (thi.x >= dom_hi.x)){
        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = tlo.z; k <= thi.z; ++k) {
                for (int j = tlo.y; j <= thi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = dom_hi.x + 1; i <= thi.x; ++i) {
                        int offset = i - dom_hi.x;
                        int i_real = dom_hi.x - offset + 1;
                        data(i, j, k, n) = data(i_real, j, k, n);
                    }
                }
            }
        }
    }


    //___________________________________________________________________________
    // Apply y-physbc to data
#if (AMREX_SPACEDIM >= 2)
    if (((bc_lo[1] == 1) || (bc_lo[1] == 2)) && (tlo.y <= dom_lo.y)) {
        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = tlo.z; k <= thi.z; ++k) {
                for (int j = tlo.y; j < dom_lo.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = tlo.x; i <= thi.x; ++i) {
                        int offset = dom_lo.y - j;
                        int j_real = dom_lo.y + offset - 1;
                        data(i, j, k, n) = data(i, j_real, k, n);
                    }
                }
            }
        }
    }

    if (((bc_hi[1] == 1) || (bc_hi[1] == 2)) && (thi.y >= dom_hi.y)) {
        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = tlo.z; k <= thi.z; ++k) {
                for (int j = dom_hi.y + 1; j <= thi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = tlo.x; i <= thi.x; ++i) {
                        int offset = j - dom_hi.y;
                        int j_real = dom_hi.y - offset + 1;
                        data(i, j, k, n) = data(i, j_real, k, n);
                    }
                }
            }
        }
    }
#endif

    //___________________________________________________________________________
    // Apply z-physbc to data
#if (AMREX_SPACEDIM >= 3)
    if (((bc_lo[2] == 1) || (bc_lo[2] == 2)) && (tlo.z <= dom_lo.z)) {
        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = tlo.z; k < dom_lo.z; ++k) {
                for (int j = tlo.y; j <= thi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = tlo.x; i <= thi.x; ++i) {
                        int offset = dom_lo.z - k;
                        int k_real = dom_lo.z + offset - 1;
                        data(i, j, k, n) = data(i, j, k_real, n);
                    }
                }
            }
        }
    }

    if (((bc_hi[2] == 1) || (bc_hi[2] == 2)) && (thi.z >= dom_hi.z)) {
        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = dom_hi.z + 1; k <= thi.z; ++k) {
                for (int j = tlo.y; j <= thi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = tlo.x; i <= thi.x; ++i) {
                        int offset = k - dom_hi.z;
                        int k_real = dom_hi.z - offset + 1;
                        data(i, j, k, n) = data(i, j, k_real, n);
                    }
                }
            }
        }
    }
#endif
}


void MultiFABPhysBCPres(MultiFab & data, const IntVect & dim_fill_ghost,
                        const Geometry & geom) {

    if (geom.isAllPeriodic()) {
        return;
    }

#ifndef GPUBC

#if (AMREX_SPACEDIM==2 || AMREX_SPACEDIM==3)
    // Physical Domain
    Box dom(geom.Domain());

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {

        const Box & bx = mfi.validbox();
        fab_physbc(BL_TO_FORTRAN_BOX(bx),
                   BL_TO_FORTRAN_BOX(dom),
                   BL_TO_FORTRAN_FAB(data[mfi]), data.nGrow(),
                   dim_fill_ghost.getVect());
    }
#endif

#else

    // Physical Domain
    Box dom(geom.Domain());

    // Effective number of ghost cells to iterate over
    int ngc         = data.nGrow();
    IntVect ngc_eff = ngc*dim_fill_ghost;

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



AMREX_GPU_HOST_DEVICE
inline void physbc_domainvel_fab(const Box & tbx,
                                 const Box & dom,
                                 const Array4<Real> & data,
                                 const GpuArray<int, AMREX_SPACEDIM> & bc_lo,
                                 const GpuArray<int, AMREX_SPACEDIM> & bc_hi,
                                 int start_cmp, int n_cmp, int dim) {

    //___________________________________________________________________________
    // Total work region => the loops below will actually only iterate over
    // cells between tbx and dom
    const Dim3 tlo    = amrex::lbound(tbx);
    const Dim3 thi    = amrex::ubound(tbx);
    const Dim3 dom_lo = amrex::lbound(dom);
    const Dim3 dom_hi = amrex::ubound(dom);
    // Compute valid parts only:
    const Dim3 vlo = amrex::elemwiseMax(tlo, dom_lo);
    const Dim3 vhi = amrex::elemwiseMin(thi, dom_hi);



    //___________________________________________________________________________
    // Apply x-physbc to data
    if ((dim == 0) && (bc_lo[0] == 2) && (tlo.x <= dom_lo.x)) {

        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = vlo.z; k <= vhi.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int j = vlo.y; j <= vhi.y; ++j) {
                    data(dom_lo.x, j, k, n) = 0;
                }
            }
        }

        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = tlo.z; k <= thi.z; ++k) {
                for (int j = tlo.y; j <= thi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = tlo.x; i < dom_lo.x; ++i) {
                        int offset = dom_lo.x - i;
                        int i_real = dom_lo.x + offset;
                        data(i, j, k, n) = - data(i_real, j, k, n);
                    }
                }
            }
        }
    }

    if ((dim == 0) && (bc_hi[0] == 2) && (thi.x >= dom_hi.x)) {

        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = vlo.z; k <= vhi.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int j = vlo.y; j <= vhi.y; ++j) {
                    data(dom_hi.x, j, k, n) = 0;
                }
            }
        }

        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = tlo.z; k <= thi.z; ++k) {
                for (int j = tlo.y; j <= thi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = dom_hi.x + 1; i <= thi.x; ++i) {
                        int offset = i - dom_hi.x;
                        int i_real = dom_hi.x - offset;
                        data(i, j, k, n) = - data(i_real, j, k, n);
                    }
                }
            }
        }
    }


    //___________________________________________________________________________
    // Apply y-physbc to data
#if (AMREX_SPACEDIM >= 2)
    if ((dim == 1) && (bc_lo[1] == 2) && (tlo.y <= dom_lo.y)){

        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = vlo.z; k <= vhi.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = vlo.x; i <= vhi.x; ++i) {
                    data(i, dom_lo.y, k, n) = 0;
                }
            }
        }

        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = tlo.z; k <= thi.z; ++k) {
                for (int j = tlo.y; j < dom_lo.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = tlo.x; i <= thi.x; ++i) {
                        int offset = dom_lo.y - j;
                        int j_real = dom_lo.y + offset;
                        data(i, j, k, n) = - data(i, j_real, k, n);
                    }
                }
            }
        }
    }

    if ((dim == 1) && (bc_hi[1] == 2) && (thi.y >= dom_hi.y)) {

        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = vlo.z; k <= vhi.z; ++k) {
                AMREX_PRAGMA_SIMD
                for (int i = vlo.x; i <= vhi.x; ++i) {
                    data(i, dom_hi.y, k, n) = 0;
                }
            }
        }

        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = tlo.z; k <= thi.z; ++k) {
                for (int j = dom_hi.y + 1; j <= thi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = tlo.x; i <= thi.x; ++i) {
                        int offset = j - dom_hi.y;
                        int j_real = dom_hi.y - offset;
                        data(i, j, k, n) = - data(i, j_real, k, n);
                    }
                }
            }
        }
    }
#endif

    //___________________________________________________________________________
    // Apply z-physbc to data
#if (AMREX_SPACEDIM >= 3)
    if ((dim == 2) && (bc_lo[2] == 2) && (tlo.z <= dom_lo.z)) {

        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int j = vlo.y; j <= vhi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = vlo.x; i <= vhi.x; ++i) {
                    data(i, j, dom_lo.z, n) = 0;
                }
            }
        }

        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = tlo.z; k < dom_lo.z; ++k) {
                for (int j = tlo.y; j <= thi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = tlo.x; i <= thi.x; ++i) {
                        int offset = dom_lo.z - k;
                        int k_real = dom_lo.z + offset;
                        data(i, j, k, n) = - data(i, j, k_real, n);
                    }
                }
            }
        }
    }

    if ((dim == 2) && (bc_hi[2] == 2) && (thi.z >= dom_hi.z) ) {

        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int j = vlo.y; j <= vhi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = vlo.x; i <= vhi.x; ++i) {
                    data(i, j, dom_hi.z, n) = 0;
                }
            }
        }

        for (int n = start_cmp; n < start_cmp + n_cmp; ++n) {
            for (int k = dom_hi.z + 1; k <= thi.z; ++k) {
                for (int j = tlo.y; j <= thi.y; ++j) {
                    AMREX_PRAGMA_SIMD
                    for (int i = tlo.x; i <= thi.x; ++i) {
                        int offset = k - dom_hi.z;
                        int k_real = dom_hi.z - offset;
                        data(i, j, k, n) = - data(i, j, k_real, n);
                    }
                }
            }
        }
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

#ifndef GPUBC

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

#else

    // Physical Domain and make sure that the domain index type matches the
    // velocity index type
    Box dom(geom.Domain());
    dom.surroundingNodes(dim);

    // Effective number of ghost cells to iterate over
    int ngc         = vel.nGrow();
    IntVect ngc_eff = ngc*dim_fill_ghost;

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

#endif
}



AMREX_GPU_HOST_DEVICE
inline void physbc_macvel_fab(const Box & tbx,
                              const Box & dom,
                              const Array4<Real> & data,
                              const GpuArray<int, AMREX_SPACEDIM> & bc_lo,
                              const GpuArray<int, AMREX_SPACEDIM> & bc_hi,
                              int dim) {

    //___________________________________________________________________________
    // Total work region => the loops below will actually only iterate over
    // cells between tbx and dom
    const Dim3 tlo    = amrex::lbound(tbx);
    const Dim3 thi    = amrex::ubound(tbx);
    const Dim3 dom_lo = amrex::lbound(dom);
    const Dim3 dom_hi = amrex::ubound(dom);


    //___________________________________________________________________________
    // Apply x-physbc to data
    if ((dim != 0) && (bc_lo[0] == 2) && (tlo.x <= dom_lo.x)) {
        for (int k = tlo.z; k <= thi.z; ++k) {
            for (int j = tlo.y; j <= thi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = tlo.x; i < dom_lo.x; ++i) {
                    int offset = dom_lo.x - i;
                    int i_real = dom_lo.x + offset - 1;
                    data(i, j, k, 0) = - data(i_real, j, k, 0);
                }
            }
        }
    }

    if ((dim != 0) && (bc_hi[0] == 2) && (thi.x >= dom_hi.x)) {
        for (int k = tlo.z; k <= thi.z; ++k) {
            for (int j = tlo.y; j <= thi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = dom_hi.x + 1; i <= thi.x; ++i) {
                    int offset = i - dom_hi.x;
                    int i_real = dom_hi.x - offset + 1;
                    data(i, j, k, 0) = - data(i_real, j, k, 0);
                }
            }
        }
    }


    //___________________________________________________________________________
    // Apply y-physbc to data
#if (AMREX_SPACEDIM >= 2)
    if ((dim != 1) && (bc_lo[1] == 2) && (tlo.y <= dom_lo.y)) {
        for (int k = tlo.z; k <= thi.z; ++k) {
            for (int j = tlo.y; j < dom_lo.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = tlo.x; i <= thi.x; ++i) {
                    int offset = dom_lo.y - j;
                    int j_real = dom_lo.y + offset - 1;
                    data(i, j, k, 0) = - data(i, j_real, k, 0);
                }
            }
        }
    }

    if ((dim != 1) && (bc_hi[1] == 2) && (thi.y >= dom_hi.y)) {
        for (int k = tlo.z; k <= thi.z; ++k) {
            for (int j = dom_hi.y + 1; j <= thi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = tlo.x; i <= thi.x; ++i) {
                    int offset = j - dom_hi.y;
                    int j_real = dom_hi.y - offset + 1;
                    data(i, j, k, 0) = - data(i, j_real, k, 0);
                }
            }
        }
    }
#endif

    //___________________________________________________________________________
    // Apply z-physbc to data
#if (AMREX_SPACEDIM >= 3)
    if ((dim != 2) && (bc_lo[2] == 2) && (tlo.z <= dom_lo.z)) {
        for (int k = tlo.z; k < dom_lo.z; ++k) {
            for (int j = tlo.y; j <= thi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = tlo.x; i <= thi.x; ++i) {
                    int offset = dom_lo.z - k;
                    int k_real = dom_lo.z + offset - 1;
                    data(i, j, k, 0) = - data(i, j, k_real, 0);
                }
            }
        }
    }

    if ((dim != 2) && (bc_hi[2] == 2) && (thi.z >= dom_hi.z)) {
        for (int k = dom_hi.z + 1; k <= thi.z; ++k) {
            for (int j = tlo.y; j <= thi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = tlo.x; i <= thi.x; ++i) {
                    int offset = k - dom_hi.z;
                    int k_real = dom_hi.z - offset + 1;
                    data(i, j, k, 0) = - data(i, j, k_real, 0);
                }
            }
        }
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

#ifndef GPUBC

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

#else

    // Physical Domain and make sure that the domain index type matches the
    // velocity index type
    Box dom(geom.Domain());
    dom.surroundingNodes(dim);

    // Effective number of ghost cells to iterate over
    int ngc         = vel.nGrow();
    IntVect ngc_eff = ngc*dim_fill_ghost;

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

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA(bx, tbx,
        {
            physbc_macvel_fab(tbx, dom, data_fab, bc_lo, bc_hi, dim);
        });
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
