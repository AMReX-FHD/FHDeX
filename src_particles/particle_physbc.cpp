#include "particle_functions.H"

// Folding routine.
// Positive-fold-add or negative-fold-add charge from ghost cells into valid region.
// For Dirichlet potential (bc_es==1), negative-fold-add.
// For Neumann potential (bc_es==2), positive-fold-add.
// Charge is cell-centered with 1 component.
// Note for wall-wall corners, the folding is also done in ghost cells so we get
// corner charges back in.
void MultiFabPhysBCCharge(MultiFab& charge, const Geometry& geom) {

    BL_PROFILE_VAR("MultiFabPhysBCCharge()",MultiFabPhysBCCharge);

#if (AMREX_SPACEDIM >= 2)

    if (geom.isAllPeriodic()) {
        return;
    }

    Box dom(geom.Domain());

    int ng = charge.nGrow();

    for (MFIter mfi(charge); mfi.isValid(); ++mfi) {

        Box bx = mfi.growntilebox(ng);

        const Array4<Real>& data = charge.array(mfi);

        //___________________________________________________________________________
        // Apply x-physbc to data

        // bc_es check is to see if we have a physical boundary condition
        // bx/dom comparison is to see if the grid touches a wall
        if ((bc_es_lo[0] == 1 || bc_es_lo[0] == 2) && (bx.smallEnd(0) <= dom.smallEnd(0))) {
            const Real fac = (bc_es_lo[0] == 1) ? -1. : 1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i >= 0 && i < ng) {
                    data(i,j,k) = data(i,j,k) + fac*data(-i-1,j,k);
                }
            });
        }

        if ((bc_es_lo[0] == 1 || bc_es_hi[0] == 2) && (bx.bigEnd(0) >= dom.bigEnd(0))) {
            const Real fac = (bc_es_hi[0] == 1) ? -1. : 1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i <= dom.bigEnd(0) && i > dom.bigEnd(0)-ng) {
                    data(i,j,k) = data(i,j,k) + fac*data(2*dom.bigEnd(0)-i+1,j,k);
                }
            });
        }

        //___________________________________________________________________________
        // Apply y-physbc to data
        if ((bc_es_lo[1] == 1 || bc_es_lo[1] == 2) && (bx.smallEnd(1) <= dom.smallEnd(1))) {
            const Real fac = (bc_es_lo[1] == 1) ? -1. : 1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j >= 0 && j < ng) {
                    data(i,j,k) = data(i,j,k) + fac*data(i,-j-1,k);
                }
            });
        }

        if ((bc_es_hi[1] == 1 || bc_es_hi[1] == 2) && (bx.bigEnd(1) >= dom.bigEnd(1))) {
            const Real fac = (bc_es_hi[1] == 1) ? -1. : 1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j <= dom.bigEnd(1) && j > dom.bigEnd(1)-ng) {
                    data(i,j,k) = data(i,j,k) + fac*data(i,2*dom.bigEnd(1)-j+1,k);
                }
            });
        }
#endif

        //___________________________________________________________________________
        // Apply z-physbc to data
#if (AMREX_SPACEDIM >= 3)
        if ((bc_es_lo[2] == 1 || bc_es_lo[2] == 2) && (bx.smallEnd(2) <= dom.smallEnd(2))) {
            const Real fac = (bc_es_lo[2] == 1) ? -1. : 1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k >= 0 && k < ng) {
                    data(i,j,k) = data(i,j,k) + fac*data(i,j,-k-1);
                }
            });
        }

        if ((bc_es_hi[2] == 1 || bc_es_hi[2] == 2) && (bx.bigEnd(2) >= dom.bigEnd(2))) {
            const Real fac = (bc_es_hi[2] == 1) ? -1. : 1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k <= dom.bigEnd(2) && k > dom.bigEnd(2)-ng) {
                    data(i,j,k) = data(i,j,k) + fac*data(i,j,2*dom.bigEnd(2)-k+1);
                }
            });
        }
#endif

    } // end MFIter
}

// Folding routine.
// Modifies normal fields on the boundary and the interior.
// Note "stress" is face-centered with direction "dim".
// We test on bc_vel_lo/hi
// 1 =    slip -> leave value on wall alone, add in force from ghost
// 2 = no slip -> set value on wall to zero, add in negative force from ghost cells
void MultiFabPhysBCDomainStress(MultiFab& stress, const Geometry& geom, int dim) {

    BL_PROFILE_VAR("ultiFabPhysBCDomainStress()",ultiFabPhysBCDomainStress);

#if (AMREX_SPACEDIM >= 2)

    if (geom.isAllPeriodic()) {
        return;
    }

    Box dom(geom.Domain());

    int ng = stress.nGrow();

    for (MFIter mfi(stress); mfi.isValid(); ++mfi) {

        Box bx = mfi.growntilebox(ng);

        const Array4<Real>& data = stress.array(mfi);

        //___________________________________________________________________________
        // Apply x-physbc to data

        if (dim == 0) {

            // dim == 0 means we are doing x-stress on x-faces
            // bc_vel check is to see if we have a wall bc
            // bx/dom comparison is to see if the grid touches a wall
            if ((bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) && (bx.smallEnd(0) <= dom.smallEnd(0))) {
                const Real fac = (bc_vel_lo[0] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i == dom.smallEnd(0)) {
                        data(i,j,k) = 0.5*(1.+fac)*data(i,j,k);
                    }
                    else if (i > 0 && i <= ng) {
                        data(i,j,k) += fac*data(-i,j,k);
                    }
                });
            }

            if ((bc_vel_lo[0] == 1 || bc_vel_hi[0] == 2) && (bx.bigEnd(0) >= dom.bigEnd(0)+1)) {
                const Real fac = (bc_vel_hi[0] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i == dom.bigEnd(0)+1) {
                        data(i,j,k) = 0.5*(1.+fac)*data(i,j,k);
                    }
                    else if (i < dom.bigEnd(0)+1 && i >= dom.bigEnd(0)+1-ng) {
                        data(i,j,k) += fac*data(2*dom.bigEnd(0)+2-i,j,k);
                    }
                });
            }

        }

        //___________________________________________________________________________
        // Apply y-physbc to data

        if (dim == 1) {

            if ((bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) && (bx.smallEnd(1) <= dom.smallEnd(1))) {
                const Real fac = (bc_vel_lo[1] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j == dom.smallEnd(1)) {
                        data(i,j,k) = 0.5*(1.+fac)*data(i,j,k);
                    }
                    else if (j > 0 && j <= ng) {
                        data(i,j,k) += fac*data(i,-j,k);
                    }
                });
            }

            if ((bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) && (bx.bigEnd(1) >= dom.bigEnd(1)+1)) {
                const Real fac = (bc_vel_hi[1] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j == dom.bigEnd(1)+1) {
                        data(i,j,k) = 0.5*(1.+fac)*data(i,j,k);
                    }
                    else if (j < dom.bigEnd(1)+1 && j >= dom.bigEnd(1)+1-ng) {
                        data(i,j,k) += fac*data(i,2*dom.bigEnd(1)+2-j,k);
                    }
                });
            }

        }
#endif

        //___________________________________________________________________________
        // Apply z-physbc to data
#if (AMREX_SPACEDIM >= 3)

        if (dim == 2) {

            if ((bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) && (bx.smallEnd(2) <= dom.smallEnd(2))) {
                const Real fac = (bc_vel_lo[2] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k == dom.smallEnd(2)) {
                        data(i,j,k) = 0.5*(1.+fac)*data(i,j,k);
                    }
                    else if (k > 0 && k <= ng) {
                        data(i,j,k) += fac*data(i,j,-k);
                    }
                });
            }

            if ((bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) && (bx.bigEnd(2) >= dom.bigEnd(2)+1)) {
                const Real fac = (bc_vel_hi[2] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k == dom.bigEnd(2)+1) {
                        data(i,j,k) = 0.5*(1.+fac)*data(i,j,k);
                    }
                    else if (k < dom.bigEnd(2)+1 && k >= dom.bigEnd(2)+1-ng) {
                        data(i,j,k) += fac*data(i,j,2*dom.bigEnd(2)+2-k);
                    }
                });
            }

        }
#endif

    } // end MFIter
}

// Folding routine.
// Modifies transverse fields in the interior.
// Note "stress" is face-centered with direction "dim".
// We test on bc_vel_lo/hi
// 1 =    slip -> add in force from ghost
// 2 = no slip -> add in negative force from ghost cells
void MultiFabPhysBCMacStress(MultiFab& stress, const Geometry& geom, int dim) {

    BL_PROFILE_VAR("ultiFabPhysBCMacStress()",ultiFabPhysBCMacStress);

#if (AMREX_SPACEDIM >= 2)

    if (geom.isAllPeriodic()) {
        return;
    }

    Box dom(geom.Domain());

    int ng = stress.nGrow();

    for (MFIter mfi(stress); mfi.isValid(); ++mfi) {

        Box bx = mfi.growntilebox(ng);

        const Array4<Real>& data = stress.array(mfi);

        //___________________________________________________________________________
        // Apply x-physbc to data

        if (dim != 0) {

            // dim != 0 means we are either y- or z-stress on x-faces
            // bc_vel check is to see if we have a wall bc
            // bx/dom comparison is to see if the grid touches a wall
            if ((bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) && (bx.smallEnd(0) <= dom.smallEnd(0))) {
                const Real fac = (bc_vel_lo[0] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i >= 0 && i < ng) {
                        data(i,j,k) += fac*data(-1-i,j,k);
                    }
                });
            }

            if ((bc_vel_lo[0] == 1 || bc_vel_hi[0] == 2) && (bx.bigEnd(0) >= dom.bigEnd(0))) {
                const Real fac = (bc_vel_hi[0] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i <= dom.bigEnd(0) && i > dom.bigEnd(0)-ng) {
                        data(i,j,k) += fac*data(2*dom.bigEnd(0)+1-i,j,k);
                    }
                });
            }

        }

        //___________________________________________________________________________
        // Apply y-physbc to data

        if (dim != 1) {

            if ((bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) && (bx.smallEnd(1) <= dom.smallEnd(1))) {
                const Real fac = (bc_vel_lo[1] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j >= 0 && j < ng) {
                        data(i,j,k) += fac*data(i,-1-j,k);
                    }
                });
            }

            if ((bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) && (bx.bigEnd(1) >= dom.bigEnd(1))) {
                const Real fac = (bc_vel_hi[1] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j <= dom.bigEnd(1) && j > dom.bigEnd(1)-ng) {
                        data(i,j,k) += fac*data(i,2*dom.bigEnd(1)+1-j,k);
                    }
                });
            }

        }
#endif

        //___________________________________________________________________________
        // Apply z-physbc to data
#if (AMREX_SPACEDIM >= 3)

        if (dim != 2) {

            if ((bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) && (bx.smallEnd(2) <= dom.smallEnd(2))) {

                const Real fac = (bc_vel_lo[2] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k >= 0 && k < ng) {
                        data(i,j,k) += fac*data(i,j,-1-k);
                    }
                });
            }

            if ((bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) && (bx.bigEnd(2) >= dom.bigEnd(2))) {
                const Real fac = (bc_vel_hi[2] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k <= dom.bigEnd(2) && k > dom.bigEnd(2)-ng) {
                        data(i,j,k) += fac*data(i,j,2*dom.bigEnd(2)+1-k);
                    }
                });
            }

        }
#endif

    } // end MFIter
}
