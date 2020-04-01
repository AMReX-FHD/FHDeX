#include "common_functions.H"


// Fills in all ghost cells to the same value, which is the value AT the boundary.
// FOEXTRAP uses boundary conditions (Neumann) and 1 interior points.
// EXT_DIR copies the supplied Dirichlet condition into the ghost cells.
void MultiFabPhysBC(MultiFab& phi, const Geometry& geom, int varType) {

    if (geom.isAllPeriodic()) {
        return;
    }

    // Physical Domain
    Box dom(geom.Domain());

    int ng = phi.nGrow();
    
    Vector<int> bc_lo(AMREX_SPACEDIM);
    Vector<int> bc_hi(AMREX_SPACEDIM);

    // compute mathematical boundary conditions
    BCPhysToMath(varType,bc_lo,bc_hi);

    for (MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        // one ghost cell
        Box bx = mfi.growntilebox(ng);

        const Array4<Real>& data = phi.array(mfi);

        //___________________________________________________________________________
        // Apply x-physbc to data

        // lo-x faces
        // bc_vel check is to see if we have a wall bc
        // bx/dom comparison is to see if the grid touches a wall
        
        if (bx.smallEnd(0) < dom.smallEnd(0)) {
            if (bc_lo[0] == FOEXTRAP) {
                int lo = dom.smallEnd(0);
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < lo) {
                        data(i,j,k) = data(lo,j,k);
                    }
                });
            }
        }
        
        if (bx.bigEnd(0) > dom.bigEnd(0)) {
            if (bc_hi[0] == FOEXTRAP) {
                int hi = dom.bigEnd(0);
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > hi) {
                        data(i,j,k) = data(hi,j,k);
                    }
                });
            }
        }

#if (AMREX_SPACEDIM >= 2)
        //___________________________________________________________________________
        // Apply y-physbc to data
        if (bx.smallEnd(1) < dom.smallEnd(1)) {
            if (bc_lo[1] == FOEXTRAP) {
                int lo = dom.smallEnd(1);
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < lo) {
                        data(i,j,k) = data(i,lo,k);
                    }
                });
            }
        }

        if (bx.bigEnd(1) > dom.bigEnd(1)) {
            if (bc_hi[1] == FOEXTRAP) {
                int hi = dom.bigEnd(1);
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > hi) {
                        data(i,j,k) = data(i,hi,k);
                    }
                });
            }
        }
#endif

#if (AMREX_SPACEDIM >= 3)
        //___________________________________________________________________________
        // Apply z-physbc to data
        if (bx.smallEnd(2) < dom.smallEnd(2)) {
            if (bc_lo[2] == FOEXTRAP) {
                int lo = dom.smallEnd(2);
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < lo) {
                        data(i,j,k) = data(i,j,lo);
                    }
                });
            }
        }

        if (bx.bigEnd(2) > dom.bigEnd(2)) {
            if (bc_hi[2] == FOEXTRAP) {
                int hi = dom.bigEnd(2);
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > hi) {
                        data(i,j,k) = data(i,j,hi);
                    }
                });
            }
        }
#endif
        
    } // end MFIter
}

// Set the value of normal velocity on walls to zero
// Set the value of normal ghost cells to the inverse reflection of the interior
// We fill all the ghost cells - they are needed for Perskin kernels and
// to avoid intermediate NaN propagation in BDS
void MultiFabPhysBCDomainVel(MultiFab& vel, const Geometry& geom, int dim) {

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

// Set the value of tranverse ghost cells to +/- the reflection of the interior
// (+ for slip walls, - for no-slip)
// We fill all the ghost cells - they are needed for Perskin kernels
void MultiFabPhysBCMacVel(MultiFab& vel, const Geometry& geom, int dim) {

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

// Set the value on walls to zero
void ZeroEdgevalWalls(std::array<MultiFab, AMREX_SPACEDIM>& edge, const Geometry& geom, int scomp, int ncomp) {

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

// Fill all ghost cells for a component of the electric field
// (note efieldCC is cell-centerd with only 1 component,
//  so this needs to be called in a loop over all directions)
// We test on bc_es_lo/hi, which are bc's for phi (not E)
// 1 = Dirichlet phi -> reflect interior values of E
// 2 = Neumann phi -> reflect and invert interior values of E
// this is different from the generic PhysBC routines since it computes
// values extrapolted to the boundary (useful for Perskin kernel interpolation)
// rather than fill ghost cells with values ON the boundary
void MultiFABElectricBC(MultiFab& efieldCC, const Geometry& geom) {
    
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
            Real fac = (bc_es_lo[0] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) {
                    data(i,j,k) = fac*data(-i-1,j,k);
                }
            });
        }

        if ((bc_es_lo[0] == 1 || bc_es_hi[0] == 2) && (bx.bigEnd(0) > dom.bigEnd(0))) {
            Real fac = (bc_es_hi[0] == 1) ? 1. : -1.;
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
            Real fac = (bc_es_lo[1] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) {
                    data(i,j,k) = fac*data(i,-j-1,k);
                }
            });
        }

        if ((bc_es_hi[1] == 1 || bc_es_hi[1] == 2) && (bx.bigEnd(1) > dom.bigEnd(1))) {
            Real fac = (bc_es_hi[1] == 1) ? 1. : -1.;
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
            Real fac = (bc_es_lo[2] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) {
                    data(i,j,k) = fac*data(i,j,-k-1);
                }
            });
        }

        if ((bc_es_hi[2] == 1 || bc_es_hi[2] == 2) && (bx.bigEnd(2) > dom.bigEnd(2))) {
            Real fac = (bc_es_hi[2] == 1) ? 1. : -1.;
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

// Fill one ghost cell for a component of the electric potential
// we test on bc_es_lo/hi
// 2 point stencil involving boundary value and interior value
// 1 = Dirichlet
// 2 = Neumann
// This routine fills ghost cells with the value extrapolated TO the ghost cell-center
// This is NOT the same as filling the ghost cell with the value on the boundary
// The Poisson solver needs a separate routine to fill ghost cells with the value ON
// the boundary for inhomogeneous Neumann and inhomogeneous Dirichlet; for this we
// use MultiFABPotentialBC_solver
void MultiFABPotentialBC(MultiFab& phi, const Geometry& geom) {
#if (AMREX_SPACEDIM >= 2)
    
    if (geom.isAllPeriodic()) {
        return;
    }

    Box dom(geom.Domain());

    const Real* dx_vec  = geom.CellSize();
    GpuArray<Real,AMREX_SPACEDIM> dx{AMREX_D_DECL(dx_vec[0], dx_vec[1], dx_vec[2])};

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
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        data(i,j,k) = -data(i+1,j,k) + 2.*potential_lo[0];
                    }
                });
            }                
            else if (bc_es_lo[0] == 2) { // Neumann
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        data(i,j,k) = data(i+1,j,k) - dx[0]*potential_lo[0];
                    }                    
                });
            }
        }

        if (bx.bigEnd(0) > dom.bigEnd(0)) {
            if (bc_es_hi[0] == 1) { // Dirichlet
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        data(i,j,k) = -data(i-1,j,k) + 2.*potential_hi[0];
                    }
                });
            }                
            else if (bc_es_hi[0] == 2) { // Neumann
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        data(i,j,k) = data(i-1,j,k) + dx[0]*potential_hi[0];
                    }                    
                });
            }
        }

        //___________________________________________________________________________
        // Apply y-physbc to data

        if (bx.smallEnd(1) < dom.smallEnd(1)) {
            if (bc_es_lo[1] == 1) { // Dirichlet
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        data(i,j,k) = -data(i,j+1,k) + 2.*potential_lo[1];
                    }
                });
            }                
            else if (bc_es_lo[1] == 2) { // Neumann
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        data(i,j,k) = data(i,j+1,k) - dx[1]*potential_lo[1];
                    }                    
                });
            }
        }

        if (bx.bigEnd(1) > dom.bigEnd(1)) {
            if (bc_es_hi[1] == 1) { // Dirichlet
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        data(i,j,k) = -data(i,j-1,k) + 2.*potential_hi[1];
                    }
                });
            }                
            else if (bc_es_hi[1] == 2) { // Neumann
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        data(i,j,k) = data(i,j-1,k) + dx[1]*potential_hi[1];
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
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        data(i,j,k) = -data(i,j,k+1) + 2.*potential_lo[2];
                    }
                });
            }                
            else if (bc_es_lo[2] == 2) { // Neumann
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        data(i,j,k) = data(i,j,k+1) - dx[2]*potential_lo[2];
                    }                    
                });
            }
        }

        if (bx.bigEnd(2) > dom.bigEnd(2)) {
            if (bc_es_hi[2] == 1) { // Dirichlet
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        data(i,j,k) = -data(i,j,k-1) + 2.*potential_hi[2];
                    }
                });
            }                
            else if (bc_es_hi[2] == 2) { // Neumann
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        data(i,j,k) = data(i,j,k-1) + dx[2]*potential_hi[2];
                    }                    
                });
            }
        }
#endif
        
    } // end MFIter
}

// Fill one ghost cell for a component of the electric potential
// This routine fills the ghost cell with the value on the boundary, whether it is Dirichlet or Neumann
// It works for inhomogeneous Neumann and inhomogeneous Dirichlet
// This is what the Poisson solver expects
// This routine is not to be confused with MultiFABPotentialBC, which fill ghost
// values extrapolated TO the ghost cell-center
void MultiFABPotentialBC_solver(MultiFab& phi, const Geometry& geom) {
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
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        data(i,j,k) = potential_lo[0];
                    }
                });
            }
        }

        if (bx.bigEnd(0) > dom.bigEnd(0)) {
            if (bc_es_hi[0] == 1 || bc_es_hi[0] == 2) {
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        data(i,j,k) = potential_hi[0];
                    }
                });
            }
        }

        //___________________________________________________________________________
        // Apply y-physbc to data

        if (bx.smallEnd(1) < dom.smallEnd(1)) {
            if (bc_es_lo[1] == 1 || bc_es_lo[1] == 2) {
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        data(i,j,k) = potential_lo[1];
                    }
                });
            }
        }

        if (bx.bigEnd(1) > dom.bigEnd(1)) {
            if (bc_es_hi[1] == 1 || bc_es_hi[1] == 2) {
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        data(i,j,k) = potential_hi[1];
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
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        data(i,j,k) = potential_lo[2];
                    }
                });
            }
        }

        if (bx.bigEnd(2) > dom.bigEnd(2)) {
            if (bc_es_hi[2] == 1 || bc_es_hi[2] == 2) {
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        data(i,j,k) = potential_hi[2];
                    }
                });
            }
        }
#endif
        
    } // end MFIter
}

// positive-fold-add or negative-fold-add charge from ghost cells ghost cells into valid region
// for Dirichlet potential (bc_es==1), negative-fold-add
// for Neumann potential (bc_es==2), positive-fold-add
// charge is cell-centered with 1 component
// note for wall-wall corners, the folding is also done in ghost cells so we get corner charges back in
void MultiFabPhysBCCharge(MultiFab& charge, const Geometry& geom) {
    
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
            Real fac = (bc_es_lo[0] == 1) ? -1. : 1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i >= 0 && i < ng) {
                    data(i,j,k) = data(i,j,k) + fac*data(-i-1,j,k);
                }
            });
        }

        if ((bc_es_lo[0] == 1 || bc_es_hi[0] == 2) && (bx.bigEnd(0) >= dom.bigEnd(0))) {
            Real fac = (bc_es_hi[0] == 1) ? -1. : 1.;
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
            Real fac = (bc_es_lo[1] == 1) ? -1. : 1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j >= 0 && j < ng) {
                    data(i,j,k) = data(i,j,k) + fac*data(i,-j-1,k);
                }
            });
        }

        if ((bc_es_hi[1] == 1 || bc_es_hi[1] == 2) && (bx.bigEnd(1) >= dom.bigEnd(1))) {
            Real fac = (bc_es_hi[1] == 1) ? -1. : 1.;
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
            Real fac = (bc_es_lo[2] == 1) ? -1. : 1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k >= 0 && k < ng) {
                    data(i,j,k) = data(i,j,k) + fac*data(i,j,-k-1);
                }
            });
        }

        if ((bc_es_hi[2] == 1 || bc_es_hi[2] == 2) && (bx.bigEnd(2) >= dom.bigEnd(2))) {
            Real fac = (bc_es_hi[2] == 1) ? -1. : 1.;
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

// Fill all ghost cells for a component of the normal stress
// note "stress" is face-centered with direction "dim"
// We test on bc_vel_lo/hi
// 1 =    slip -> leave value on wall alone, add in force from ghost
// 2 = no slip -> set value on wall to zero, add in negative force from ghost cells
void MultiFabPhysBCDomainStress(MultiFab& stress, const Geometry& geom, int dim) {
    
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
                Real fac = (bc_vel_lo[0] == 1) ? 1. : -1.;
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
                Real fac = (bc_vel_hi[0] == 1) ? 1. : -1.;
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
                Real fac = (bc_vel_lo[1] == 1) ? 1. : -1.;
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
                Real fac = (bc_vel_hi[1] == 1) ? 1. : -1.;
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
                Real fac = (bc_vel_lo[2] == 1) ? 1. : -1.;
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
                Real fac = (bc_vel_hi[2] == 1) ? 1. : -1.;
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

void MultiFabPhysBCMacStress(MultiFab& stress, const Geometry& geom, int dim) {
    
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
                Real fac = (bc_vel_lo[0] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i >= 0 && i < ng) {
                        data(i,j,k) += fac*data(-1-i,j,k);
                    }
                });
            }
            
            if ((bc_vel_lo[0] == 1 || bc_vel_hi[0] == 2) && (bx.bigEnd(0) >= dom.bigEnd(0))) {
                Real fac = (bc_vel_hi[0] == 1) ? 1. : -1.;
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
                Real fac = (bc_vel_lo[1] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j >= 0 && j < ng) {
                        data(i,j,k) += fac*data(i,-1-j,k);
                    }
                });
            }

            if ((bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) && (bx.bigEnd(1) >= dom.bigEnd(1))) {
                Real fac = (bc_vel_hi[1] == 1) ? 1. : -1.;
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

                Real fac = (bc_vel_lo[2] == 1) ? 1. : -1.;
                amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k >= 0 && k < ng) {
                        data(i,j,k) += fac*data(i,j,-1-k);
                    }
                });
            }
            
            if ((bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) && (bx.bigEnd(2) >= dom.bigEnd(2))) {
                Real fac = (bc_vel_hi[2] == 1) ? 1. : -1.;                
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
