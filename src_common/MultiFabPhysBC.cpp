#include "common_functions.H"
#include "InhomogeneousBCVal.H"

// Ghost cell filling routine.
// Fills in ALL ghost cells to the value ON the boundary.
// amrex::BCType::foextrap uses boundary conditions (Neumann) and 1 interior points.
// amrex::BCType::ext_dir copies the supplied Dirichlet condition into the ghost cells.
void MultiFabPhysBC(MultiFab& phi, const Geometry& geom, int scomp, int ncomp, int bccomp, const Real& time) {

    BL_PROFILE_VAR("MultiFabPhysBC()",MultiFabPhysBC);

    // bccomp definitions are in BCPhysToMath.cpp

    if (geom.isAllPeriodic() || phi.nGrow() == 0) {
        return;
    }

    // Physical Domain
    Box dom(geom.Domain());

    GpuArray<Real,3> dx;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        dx[d] = geom.CellSize(d);
    }
    if (AMREX_SPACEDIM == 2) {
        dx[2] = cell_depth;
    }

    int nghost = phi.nGrow();

    // compute mathematical boundary conditions
    Vector<int> bc_lo(AMREX_SPACEDIM);
    Vector<int> bc_hi(AMREX_SPACEDIM);
    BCPhysToMath(bccomp,bc_lo,bc_hi);

    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // one ghost cell
        Box bx = mfi.growntilebox(nghost);

        const Array4<Real>& data = phi.array(mfi);

        //___________________________________________________________________________
        // Apply x-physbc to data

        // lo-x faces
        // bc_vel check is to see if we have a wall bc
        // bx/dom comparison is to see if the grid touches a wall

        int lo = dom.smallEnd(0);
        int hi = dom.bigEnd(0);

        if (bx.smallEnd(0) < lo) {
            Real x = prob_lo[0];
            if (bc_lo[0] == amrex::BCType::foextrap) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (i < lo) {
                        Real y = prob_lo[1] + (j+0.5)*dx[1];
                        Real z = prob_lo[2] + (k+0.5)*dx[2];
                        data(i,j,k,scomp+n) = data(lo,j,k,scomp+n) - 0.5*dx[0]*InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
            else if (bc_lo[0] == amrex::BCType::ext_dir) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (i < lo) {
                        Real y = prob_lo[1] + (j+0.5)*dx[1];
                        Real z = prob_lo[2] + (k+0.5)*dx[2];
                        data(i,j,k,scomp+n) = InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
        }

        if (bx.bigEnd(0) > hi) {
            Real x = prob_hi[0];
            if (bc_hi[0] == amrex::BCType::foextrap) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (i > hi) {
                        Real y = prob_lo[1] + (j+0.5)*dx[1];
                        Real z = prob_lo[2] + (k+0.5)*dx[2];
                        data(i,j,k,scomp+n) = data(hi,j,k,scomp+n) - 0.5*dx[0]*InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
            else if (bc_hi[0] == amrex::BCType::ext_dir) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (i > hi) {
                        Real y = prob_lo[1] + (j+0.5)*dx[1];
                        Real z = prob_lo[2] + (k+0.5)*dx[2];
                        data(i,j,k,scomp+n) = InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
        }

#if (AMREX_SPACEDIM >= 2)
        //___________________________________________________________________________
        // Apply y-physbc to data

        lo = dom.smallEnd(1);
        hi = dom.bigEnd(1);

        if (bx.smallEnd(1) < lo) {
            Real y = prob_lo[1];
            if (bc_lo[1] == amrex::BCType::foextrap) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (j < lo) {
                        Real x = prob_lo[0] + (i+0.5)*dx[0];
                        Real z = prob_lo[2] + (k+0.5)*dx[2];
                        data(i,j,k,scomp+n) = data(i,lo,k,scomp+n) - 0.5*dx[1]*InhomogeneousBCVal(bccomp+n,x,y,z,time);;
                    }
                });
            }
            else if (bc_lo[1] == amrex::BCType::ext_dir) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (j < lo) {
                        Real x = prob_lo[0] + (i+0.5)*dx[0];
                        Real z = prob_lo[2] + (k+0.5)*dx[2];
                        data(i,j,k,scomp+n) = InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
        }

        if (bx.bigEnd(1) > hi) {
            Real y = prob_hi[1];
            if (bc_hi[1] == amrex::BCType::foextrap) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (j > hi) {
                        Real x = prob_lo[0] + (i+0.5)*dx[0];
                        Real z = prob_lo[2] + (k+0.5)*dx[2];
                        data(i,j,k,scomp+n) = data(i,hi,k,scomp+n) - 0.5*dx[1]*InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
            else if (bc_hi[1] == amrex::BCType::ext_dir) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (j > hi) {
                        Real x = prob_lo[0] + (i+0.5)*dx[0];
                        Real z = prob_lo[2] + (k+0.5)*dx[2];
                        data(i,j,k,scomp+n) = InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
        }
#endif

#if (AMREX_SPACEDIM >= 3)
        //___________________________________________________________________________
        // Apply z-physbc to data

        lo = dom.smallEnd(2);
        hi = dom.bigEnd(2);

        if (bx.smallEnd(2) < lo) {
            Real z = prob_lo[2];
            if (bc_lo[2] == amrex::BCType::foextrap) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (k < lo) {
                        Real x = prob_lo[0] + (i+0.5)*dx[0];
                        Real y = prob_lo[1] + (j+0.5)*dx[1];
                        data(i,j,k,scomp+n) = data(i,j,lo,scomp+n) - 0.5*dx[2]*InhomogeneousBCVal(bccomp+n,x,y,z,time);;
                    }
                });
            }
            else if (bc_lo[2] == amrex::BCType::ext_dir) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (k < lo) {
                        Real x = prob_lo[0] + (i+0.5)*dx[0];
                        Real y = prob_lo[1] + (j+0.5)*dx[1];
                        data(i,j,k,scomp+n) = InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
        }

        if (bx.bigEnd(2) > hi) {
            Real z= prob_hi[2];
            if (bc_hi[2] == amrex::BCType::foextrap) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (k > hi) {
                        Real x = prob_lo[0] + (i+0.5)*dx[0];
                        Real y = prob_lo[1] + (j+0.5)*dx[1];
                        data(i,j,k,scomp+n) = data(i,j,hi,scomp+n) - 0.5*dx[2]*InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
            else if (bc_hi[2] == amrex::BCType::ext_dir) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (k > hi) {
                        Real x = prob_lo[0] + (i+0.5)*dx[0];
                        Real y = prob_lo[1] + (j+0.5)*dx[1];
                        data(i,j,k,scomp+n) = InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
        }
#endif

    } // end MFIter
}

// Boundary and ghost cell filling routine for normal velocities.
// Set the value of normal velocity on walls to zero.
// Set the value of normal ghost cells to the inverse reflection of the interior.
// We fill all the ghost cells - they are needed for Perskin kernels and
// to avoid intermediate NaN propagation in BDS.
void MultiFabPhysBCDomainVel(MultiFab& vel, const Geometry& geom, int dim) {

    BL_PROFILE_VAR("MultiFabPhysBCDomainVel()",MultiFabPhysBCDomainVel);

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
                else if (i == dom.smallEnd(0)) {
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

// Ghost cell filling routine for transverse velocities.
// Set the value of tranverse ghost cells to +/- the reflection of the interior
// (+ for slip walls, - for no-slip).
// We fill all the ghost cells - they are needed for Perskin kernels.
void MultiFabPhysBCMacVel(MultiFab& vel, const Geometry& geom, int dim, int is_inhomogeneous) {

    BL_PROFILE_VAR("MultiFabPhysBCMacVel()",MultiFabPhysBCMacVel);

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
            if (bc_vel_lo[0] == 1) { // slip
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        data(i,j,k) = data(-i-1,j,k);
                    }
                });
            } else if (bc_vel_lo[0] == 2) { // no-slip
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        data(i,j,k) = is_inhomogeneous*2.*wallspeed_x_lo[dim] - data(-i-1,j,k);
                    }
                });
            }
        }

        if ((dim != 0) && (bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) && (bx.bigEnd(0) > dom.bigEnd(0))) {
            if (bc_vel_hi[0] == 1) { // slip
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        data(i,j,k) = data(2*dom.bigEnd(0)-i+1,j,k);
                    }
                });
            } else if (bc_vel_hi[0] == 2) { // no-slip
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        data(i,j,k) = is_inhomogeneous*2.*wallspeed_x_hi[dim] - data(2*dom.bigEnd(0)-i+1,j,k);
                    }
                });
            }
        }

#if (AMREX_SPACEDIM >= 2)
        //___________________________________________________________________________
        // Apply y-physbc to data
        if ((dim != 1) && (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) && (bx.smallEnd(1) < dom.smallEnd(1))) {
            if (bc_vel_lo[1] == 1) { // slip
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        data(i,j,k) = data(i,-j-1,k);
                    }
                });
            } else if (bc_vel_lo[1] == 2) { // no-slip
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        data(i,j,k) = is_inhomogeneous*2.*wallspeed_y_lo[dim] - data(i,-j-1,k);
                    }
                });
            }
        }

        if ((dim != 1) && (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) && (bx.bigEnd(1) > dom.bigEnd(1))) {
            if (bc_vel_hi[1] == 1) { // slip
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        data(i,j,k) = data(i,2*dom.bigEnd(1)-j+1,k);
                    }
                });
            } else if (bc_vel_hi[1] == 2) { // no-slip
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        data(i,j,k) = is_inhomogeneous*2.*wallspeed_y_hi[dim] - data(i,2*dom.bigEnd(1)-j+1,k);
                    }
                });
            }
        }
#endif

#if (AMREX_SPACEDIM >= 3)
        //___________________________________________________________________________
        // Apply z-physbc to data
        if ((dim != 2) && (bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) && (bx.smallEnd(2) < dom.smallEnd(2))) {
            if (bc_vel_lo[2] == 1) { // slip
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        data(i,j,k) = data(i,j,-k-1);
                    }
                });
            } else if (bc_vel_lo[2] == 2) { // no-slip
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        data(i,j,k) = is_inhomogeneous*2.*wallspeed_z_lo[dim] - data(i,j,-k-1);
                    }
                });
            }
        }

        if ((dim != 2) && (bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) && (bx.bigEnd(2) > dom.bigEnd(2))) {
            if (bc_vel_hi[2] == 1) { // slip
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        data(i,j,k) = data(i,j,2*dom.bigEnd(2)-k+1);
                    }
                });
            } else if (bc_vel_hi[2] == 2) { // no-slip
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        data(i,j,k) = is_inhomogeneous*2.*wallspeed_z_hi[dim] - data(i,j,2*dom.bigEnd(2)-k+1);
                    }
                });
            }
        }
#endif

    } // end MFIter
}

// Boundary filling routine for normal velocity.
// Set the value of normal velocity on walls to zero.
// Works for slip and no-slip (bc_vel_lo/hi = 1 or 2).
void ZeroEdgevalWalls(std::array<MultiFab, AMREX_SPACEDIM>& edge, const Geometry& geom, int scomp, int ncomp) {

    BL_PROFILE_VAR("ZeroEdgevalWalls()",ZeroEdgevalWalls);

    if (geom.isAllPeriodic()) {
        return;
    }

    Box dom(geom.Domain());

    for (MFIter mfi(edge[0]); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();

        const Array4<Real>& data = edge[0].array(mfi);

        //___________________________________________________________________________
        // Apply x-physbc to data

        // lo-x faces
        // bc_vel check is to see if we have a wall bc
        // bx/dom comparison is to see if the grid touches a wall
        if ((bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) && (bx.smallEnd(0) <= dom.smallEnd(0))) {

            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (i == dom.smallEnd(0)) {
                    // set normal velocity on boundary to zero
                    data(i,j,k,scomp+n) = 0.;
                }
            });
        }

        // hi-x faces
        if ((bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) && (bx.bigEnd(0) >= dom.bigEnd(0)+1)) {

            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (i == dom.bigEnd(0)+1) {
                    data(i,j,k,scomp+n) = 0.;
                }
            });
        }
    }

#if (AMREX_SPACEDIM >= 2)

    for (MFIter mfi(edge[1]); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();

        const Array4<Real>& data = edge[1].array(mfi);

        //___________________________________________________________________________
        // Apply y-physbc to data

        // lo-y faces
        if ((bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) && (bx.smallEnd(1) <= dom.smallEnd(1))) {

            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (j == dom.smallEnd(1)) {
                    data(i,j,k,scomp+n) = 0.;
                }
            });
        }

        // hi-y faces
        if ((bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) && (bx.bigEnd(1) >= dom.bigEnd(1)+1)) {

            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (j == dom.bigEnd(1)+1) {
                    data(i,j,k,scomp+n) = 0.;
                }
            });
        }
    }

#endif

#if (AMREX_SPACEDIM >= 3)

    for (MFIter mfi(edge[2]); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();

        const Array4<Real>& data = edge[2].array(mfi);

        //___________________________________________________________________________
        // Apply z-physbc to data

        // lo-z faces
        if ((bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) && (bx.smallEnd(2) <= dom.smallEnd(2))) {

            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (k == dom.smallEnd(2)) {
                    data(i,j,k,scomp+n) = 0.;
                }
            });
        }

        // hi-z faces
        if ((bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) && (bx.bigEnd(2) >= dom.bigEnd(2)+1)) {

            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (k == dom.bigEnd(2)+1) {
                    data(i,j,k,scomp+n) = 0.;
                }
            });
        }
    }

#endif
}

// Boundary filling routine for normal velocity.
// Set the value of normal velocity on all physical (i.e. non-periodic) walls to zero.
// Checks that bc_vel_lo/hi != -1
void ZeroEdgevalPhysical(std::array<MultiFab, AMREX_SPACEDIM>& edge, const Geometry& geom, int scomp, int ncomp) {

    BL_PROFILE_VAR("ZeroEdgevalPhysical()",ZeroEdgevalPhysical);

    if (geom.isAllPeriodic()) {
        return;
    }

    Box dom(geom.Domain());

    for (MFIter mfi(edge[0]); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();

        const Array4<Real>& data = edge[0].array(mfi);

        //___________________________________________________________________________
        // Apply x-physbc to data

        // lo-x faces
        // bc_vel check is to see if we have a wall bc
        // bx/dom comparison is to see if the grid touches a wall
        if ((bc_vel_lo[0] != -1) && (bx.smallEnd(0) <= dom.smallEnd(0))) {

            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (i == dom.smallEnd(0)) {
                    // set normal velocity on boundary to zero
                    data(i,j,k,scomp+n) = 0.;
                }
            });
        }

        // hi-x faces
        if ((bc_vel_hi[0] != -1) && (bx.bigEnd(0) >= dom.bigEnd(0)+1)) {

            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (i == dom.bigEnd(0)+1) {
                    data(i,j,k,scomp+n) = 0.;
                }
            });
        }
    }

#if (AMREX_SPACEDIM >= 2)

    for (MFIter mfi(edge[1]); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();

        const Array4<Real>& data = edge[1].array(mfi);

        //___________________________________________________________________________
        // Apply y-physbc to data

        // lo-y faces
        if ((bc_vel_lo[1] != -1) && (bx.smallEnd(1) <= dom.smallEnd(1))) {

            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (j == dom.smallEnd(1)) {
                    data(i,j,k,scomp+n) = 0.;
                }
            });
        }

        // hi-y faces
        if ((bc_vel_hi[1] != -1) && (bx.bigEnd(1) >= dom.bigEnd(1)+1)) {

            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (j == dom.bigEnd(1)+1) {
                    data(i,j,k,scomp+n) = 0.;
                }
            });
        }
    }

#endif

#if (AMREX_SPACEDIM >= 3)

    for (MFIter mfi(edge[2]); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();

        const Array4<Real>& data = edge[2].array(mfi);

        //___________________________________________________________________________
        // Apply z-physbc to data

        // lo-z faces
        if ((bc_vel_lo[2] != -1) && (bx.smallEnd(2) <= dom.smallEnd(2))) {

            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (k == dom.smallEnd(2)) {
                    data(i,j,k,scomp+n) = 0.;
                }
            });
        }

        // hi-z faces
        if ((bc_vel_hi[2] != -1) && (bx.bigEnd(2) >= dom.bigEnd(2)+1)) {

            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (k == dom.bigEnd(2)+1) {
                    data(i,j,k,scomp+n) = 0.;
                }
            });
        }
    }

#endif
}

// Ghost cell filling routine.
// Fill all ghost cells for a component of the electric field
// by extrapolating values to the ghost cell-centers (useful for Perskin kernels)
// Note this is not the same as extrapolating values to the boundary.
// (note the input MultiFab is cell-centerd with only 1 component,
//  so this needs to be called in a loop over all directions)
// We test on bc_es_lo/hi, which are bc's for phi (not E)
// 1 = Dirichlet phi -> reflect interior values of E
// 2 = Neumann phi -> reflect and invert interior values of E
void MultiFabElectricBC(MultiFab& efieldCC, const Geometry& geom) {

    BL_PROFILE_VAR("MultiFabElectricBC()",MultiFabElectricBC);

#if (AMREX_SPACEDIM >= 2)

    if (geom.isAllPeriodic()) {
        return;
    }

    Box dom(geom.Domain());

    int ng = efieldCC.nGrow();

    for (MFIter mfi(efieldCC); mfi.isValid(); ++mfi) {

        Box bx = mfi.growntilebox(ng);

        const Array4<Real>& data = efieldCC.array(mfi);

        //___________________________________________________________________________
        // Apply x-physbc to data

        // bc_es check is to see if we have a physical boundary condition
        // bx/dom comparison is to see if the grid touches a wall
        if ((bc_es_lo[0] == 1 || bc_es_lo[0] == 2) && (bx.smallEnd(0) < dom.smallEnd(0))) {
            const Real fac = (bc_es_lo[0] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) {
                    data(i,j,k) = fac*data(-i-1,j,k);
                }
            });
        }

        if ((bc_es_lo[0] == 1 || bc_es_hi[0] == 2) && (bx.bigEnd(0) > dom.bigEnd(0))) {
            const Real fac = (bc_es_hi[0] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > dom.bigEnd(0)) {
                    data(i,j,k) = fac*data(2*dom.bigEnd(0)-i+1,j,k);
                }
            });
        }

        //___________________________________________________________________________
        // Apply y-physbc to data
        if ((bc_es_lo[1] == 1 || bc_es_lo[1] == 2) && (bx.smallEnd(1) < dom.smallEnd(1))) {
            const Real fac = (bc_es_lo[1] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) {
                    data(i,j,k) = fac*data(i,-j-1,k);
                }
            });
        }

        if ((bc_es_hi[1] == 1 || bc_es_hi[1] == 2) && (bx.bigEnd(1) > dom.bigEnd(1))) {
            const Real fac = (bc_es_hi[1] == 1) ? 1. : -1.;
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
        if ((bc_es_lo[2] == 1 || bc_es_lo[2] == 2) && (bx.smallEnd(2) < dom.smallEnd(2))) {
            const Real fac = (bc_es_lo[2] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) {
                    data(i,j,k) = fac*data(i,j,-k-1);
                }
            });
        }

        if ((bc_es_hi[2] == 1 || bc_es_hi[2] == 2) && (bx.bigEnd(2) > dom.bigEnd(2))) {
            const Real fac = (bc_es_hi[2] == 1) ? 1. : -1.;
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

// Ghost cell filling routine.
// This routine fills ghost cells for electric potential with the value
// extrapolated TO the ghost cell-center; we test on bc_es_lo/hi.
// 1 = Dirichlet
// 2 = Neumann
// Uses a 2 point stencil involving boundary value and interior value.
// This is NOT the same as filling the ghost cell with the value on the boundary.
// (not implemented yet)
// This is NOT the same as filling the ghost cell with the Dirichlet or
// Neumann value on the boundary.
// (the Poisson solver needs this; MultifFabPotentialBC_solver())
void MultiFabPotentialBC(MultiFab& phi, const Geometry& geom) {

    BL_PROFILE_VAR("MultiFabPotentialBC()",MultiFabPotentialBC);

#if (AMREX_SPACEDIM >= 2)

    if (geom.isAllPeriodic()) {
        return;
    }

    Box dom(geom.Domain());

    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {

        // fill ONE ghost cell
        Box bx = mfi.growntilebox(1);

        const Array4<Real>& data = phi.array(mfi);

        //___________________________________________________________________________
        // Apply x-physbc to data

        // bx/dom comparison is to see if the grid touches a wall
        // bc_es check is to see if we have a physical boundary condition
        if (bx.smallEnd(0) < dom.smallEnd(0)) {
            if (bc_es_lo[0] == 1) { // Dirichlet
                const Real pot = potential_lo[0];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        data(i,j,k) = -data(i+1,j,k) + 2.*pot;
                    }
                });
            }
            else if (bc_es_lo[0] == 2) { // Neumann
                const Real pot = potential_lo[0];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        data(i,j,k) = data(i+1,j,k) - dx[0]*pot;
                    }
                });
            }
        }

        if (bx.bigEnd(0) > dom.bigEnd(0)) {
            if (bc_es_hi[0] == 1) { // Dirichlet
                const Real pot = potential_hi[0];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        data(i,j,k) = -data(i-1,j,k) + 2.*pot;
                    }
                });
            }
            else if (bc_es_hi[0] == 2) { // Neumann
                const Real pot = potential_hi[0];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        data(i,j,k) = data(i-1,j,k) + dx[0]*pot;
                    }
                });
            }
        }

        //___________________________________________________________________________
        // Apply y-physbc to data

        if (bx.smallEnd(1) < dom.smallEnd(1)) {
            if (bc_es_lo[1] == 1) { // Dirichlet
                const Real pot = potential_lo[1];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        data(i,j,k) = -data(i,j+1,k) + 2.*pot;
                    }
                });
            }
            else if (bc_es_lo[1] == 2) { // Neumann
                const Real pot = potential_lo[1];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        data(i,j,k) = data(i,j+1,k) - dx[1]*pot;
                    }
                });
            }
        }

        if (bx.bigEnd(1) > dom.bigEnd(1)) {
            if (bc_es_hi[1] == 1) { // Dirichlet
                const Real pot = potential_hi[1];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        data(i,j,k) = -data(i,j-1,k) + 2.*pot;
                    }
                });
            }
            else if (bc_es_hi[1] == 2) { // Neumann
                const Real pot = potential_hi[1];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        data(i,j,k) = data(i,j-1,k) + dx[1]*pot;
                    }
                });
            }
        }


#endif
#if (AMREX_SPACEDIM >= 3)

        //___________________________________________________________________________
        // Apply z-physbc to data

        if (bx.smallEnd(2) < dom.smallEnd(2)) {
            if (bc_es_lo[2] == 1) { // Dirichlet
                const Real pot = potential_lo[2];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        data(i,j,k) = -data(i,j,k+1) + 2.*pot;
                    }
                });
            }
            else if (bc_es_lo[2] == 2) { // Neumann
                const Real pot = potential_lo[2];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        data(i,j,k) = data(i,j,k+1) - dx[2]*pot;
                    }
                });
            }
        }

        if (bx.bigEnd(2) > dom.bigEnd(2)) {
            if (bc_es_hi[2] == 1) { // Dirichlet
                const Real pot = potential_hi[2];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        data(i,j,k) = -data(i,j,k-1) + 2.*pot;
                    }
                });
            }
            else if (bc_es_hi[2] == 2) { // Neumann
                const Real pot = potential_hi[2];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        data(i,j,k) = data(i,j,k-1) + dx[2]*pot;
                    }
                });
            }
        }
#endif

    } // end MFIter
}

// Ghost cell filling routine.
// Fill one ghost cell for a component of the electric potential.
// This routine fills the ghost cell with the numerical value of the
// (possibly inhomogeneous) Dirichlet or Neumann value on the boundary.
// This is what the Poisson solver expects.
// This is NOT the same as filling the ghost cell with the value on the boundary.
// (not implemented yet)
// This is NOT the same as filling the ghost cell with the value extrapolated to the ghost cell-center
// (implemented in MultiFabPotentialBC())
void MultiFabPotentialBC_solver(MultiFab& phi, const Geometry& geom) {

    BL_PROFILE_VAR("MultiFabPotentialBC_solver()",MultiFabPotentialBC_solver);

#if (AMREX_SPACEDIM >= 2)

    if (geom.isAllPeriodic()) {
        return;
    }

    Box dom(geom.Domain());

    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {

        // fill ONE ghost cell
        Box bx = mfi.growntilebox(1);

        const Array4<Real>& data = phi.array(mfi);

        //___________________________________________________________________________
        // Apply x-physbc to data

        // bx/dom comparison is to see if the grid touches a wall
        // bc_es check is to see if we have a physical boundary condition
        if (bx.smallEnd(0) < dom.smallEnd(0)) {
            if (bc_es_lo[0] == 1 || bc_es_lo[0] == 2) {
                const Real pot = potential_lo[0];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        data(i,j,k) = pot;
                    }
                });
            }
        }

        if (bx.bigEnd(0) > dom.bigEnd(0)) {
            if (bc_es_hi[0] == 1 || bc_es_hi[0] == 2) {
                const Real pot = potential_hi[0];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        data(i,j,k) = pot;
                    }
                });
            }
        }

        //___________________________________________________________________________
        // Apply y-physbc to data

        if (bx.smallEnd(1) < dom.smallEnd(1)) {
            if (bc_es_lo[1] == 1 || bc_es_lo[1] == 2) {
                const Real pot = potential_lo[1];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        data(i,j,k) = pot;
                    }
                });
            }
        }

        if (bx.bigEnd(1) > dom.bigEnd(1)) {
            if (bc_es_hi[1] == 1 || bc_es_hi[1] == 2) {
                const Real pot = potential_hi[1];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        data(i,j,k) = pot;
                    }
                });
            }
        }


#endif
#if (AMREX_SPACEDIM >= 3)

        //___________________________________________________________________________
        // Apply z-physbc to data

        if (bx.smallEnd(2) < dom.smallEnd(2)) {
            if (bc_es_lo[2] == 1 || bc_es_lo[2] == 2) {
                const Real pot = potential_lo[2];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        data(i,j,k) = pot;
                    }
                });
            }
        }

        if (bx.bigEnd(2) > dom.bigEnd(2)) {
            if (bc_es_hi[2] == 1 || bc_es_hi[2] == 2) {
                const Real pot = potential_hi[2];
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        data(i,j,k) = pot;
                    }
                });
            }
        }
#endif

    } // end MFIter
}

