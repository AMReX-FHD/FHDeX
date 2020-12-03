#include "compressible_functions.H"
#include "compressible_functions_stag.H"

// Set boundary and ghost cells for staggered compressible code based on BCs
void setBCStag(MultiFab& prim_in, MultiFab& cons_in,
                 std::array< MultiFab, AMREX_SPACEDIM >& cumom_in,
                 std::array< MultiFab, AMREX_SPACEDIM >& vel_in,
                 const amrex::Geometry geom)
{
    BL_PROFILE_VAR("setBCStag()",setBCStag);

    int ng_c = cons_in.nGrow();
    int ng_p = prim_in.nGrow();
    if (ng_c != ng_p) {
        Abort("setBC: prim and cons need the same number of ghost cells");
    }
    int ng_m = cumom_in[0].nGrow();
    int ng_v = vel_in[0].nGrow();
    if (ng_m != ng_v) {
        Abort("setBC: momentum and velocity need the same number of ghost cells");
    }

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        BCMassTempPress(prim_in, geom, i);
        BCMomNormal(cumom_in[i], vel_in[i], geom, i);
        BCMomTrans(cumom_in[i], vel_in[i], geom, i);
    }
    BCRhoRhoE(cons_in, prim_in, cumom_in, geom);
}


// Set mass, pressure and temperature on ghost cells based on BCs
void BCMassTempPress(MultiFab& prim_in,const amrex::Geometry geom,int dim)
{
    BL_PROFILE_VAR("BCMassTempPress()",BCMassTempPress);

    int nspecies_gpu = nspecies;

    Box dom(geom.Domain());
    int ng_p = prim_in.nGrow();

    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng_p);

        const Array4<Real>& prim = prim_in.array(mfi);

        // LO X
        if ((dim == 0) && (bx.smallEnd(0) < dom.smallEnd(0))) {

            int lo = dom.smallEnd(0);

            // mass fractions, wall
            if ( bc_mass_lo[0] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(2*lo-i-1,j,k,n);
                        }
                    }
                });
            }

            // mass fractions, reservoir
            if (bc_mass_lo[0] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        for (int n=0; n<nspecies_gpu; ++n) {
                            prim(i,j,k,6+n)          = 2.*bc_Yk_x_lo[n] - prim(2*lo-i-1,j,k,6+n);
                            prim(i,j,k,6+nspecies+n) = 2.*bc_Xk_x_lo[n] - prim(2*lo-i-1,j,k,6+nspecies+n);
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            if (bc_therm_lo[0] == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(2*lo-i-1,j,k,n);
                        }
                    }
                });
            } else if (bc_therm_lo[0] == 2) { // isothermal
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        prim(i,j,k,4) = -prim(2*lo-i-1,j,k,4) + 2.*t_lo[0];
                        prim(i,j,k,5) = prim(2*lo-i-1,j,k,5);
                    }
                });
            }

        } // end LO X
        
        // HI X
        if ((dim == 0) && (bx.bigEnd(0) > dom.bigEnd(0))) {

            int hi = dom.bigEnd(0);
            
            // mass fractions, wall
            if ( bc_mass_hi[0] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(2*hi-i+1,j,k,n);

                        }
                    }
                });
            }

            // mass fractions, reservoir
            if (bc_mass_hi[0] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        for (int n=0; n<nspecies_gpu; ++n) {
                            prim(i,j,k,6+n)          = 2.*bc_Yk_x_hi[n] - prim(2*hi-i+1,j,k,6+n);
                            prim(i,j,k,6+nspecies+n) = 2.*bc_Xk_x_hi[n] - prim(2*hi-i+1,j,k,6+nspecies+n);
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            if (bc_therm_hi[0] == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(2*hi-i+1,j,k,n);
                        }
                    }
                });
            } else if (bc_therm_hi[0] == 2) { // isothermal
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        prim(i,j,k,4) = -prim(2*hi-i+1,j,k,4) + 2.*t_hi[0];
                        prim(i,j,k,5) = prim(2*hi-i+1,j,k,5);
                    }
                });
            }

        } // end HI X
        
        // LO Y
        if ((dim == 1) && (bx.smallEnd(1) < dom.smallEnd(1))) {

            int lo = dom.smallEnd(1);

            // mass fractions, wall
            if ( bc_mass_lo[1] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(i,2*lo-j-1,k,n);
                        }
                    }
                });
            }

            // mass fractions, reservoir
            if (bc_mass_lo[1] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        for (int n=0; n<nspecies_gpu; ++n) {
                            prim(i,j,k,6+n)          = 2.*bc_Yk_y_lo[n] - prim(i,2*lo-j-1,k,6+n);
                            prim(i,j,k,6+nspecies+n) = 2.*bc_Xk_y_lo[n] - prim(i,2*lo-j-1,k,6+nspecies+n);
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            if (bc_therm_lo[1] == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(i,2*lo-j-1,k,n);
                        }
                    }
                });
            } else if (bc_therm_lo[1] == 2) { // isothermal
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        prim(i,j,k,4) = -prim(i,2*lo-j-1,k,4) + 2.*t_lo[1];
                        prim(i,j,k,5) = prim(i,2*lo-j-1,k,5);
                    }
                });
            }

        } // end LO Y
        

        // HI Y
        if ((dim == 1) && (bx.bigEnd(1) > dom.bigEnd(1))) {

            int hi = dom.bigEnd(1);
            
            // mass fractions, wall
            if ( bc_mass_hi[1] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(i,2*hi-j+1,k,n);

                        }
                    }
                });
            }

            // mass fractions, reservoir
            if (bc_mass_hi[1] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        for (int n=0; n<nspecies_gpu; ++n) {
                            prim(i,j,k,6+n)          = 2.*bc_Yk_y_hi[n] - prim(i,2*hi-j+1,k,6+n);
                            prim(i,j,k,6+nspecies+n) = 2.*bc_Xk_y_hi[n] - prim(i,2*hi-j+1,k,6+nspecies+n);
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            if (bc_therm_hi[1] == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(i,2*hi-j+1,k,n);
                        }
                    }
                });
            } else if (bc_therm_hi[1] == 2) { // isothermal
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        prim(i,j,k,4) = -prim(i,2*hi-j+1,k,4) + 2.*t_hi[1];
                        prim(i,j,k,5) = prim(i,2*hi-j+1,k,5);
                    }
                });
            }

        } // end HI Y
        
        // LO Z
        if ((dim == 2) && (bx.smallEnd(2) < dom.smallEnd(2))) {

            int lo = dom.smallEnd(2);

            // mass fractions, wall
            if ( bc_mass_lo[2] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(i,j,2*lo-k-1,n);
                        }
                    }
                });
            }

            // mass fractions, reservoir
            if (bc_mass_lo[2] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        for (int n=0; n<nspecies_gpu; ++n) {
                            prim(i,j,k,6+n)          = 2.*bc_Yk_z_lo[n] - prim(i,j,2*lo-k-1,6+n);
                            prim(i,j,k,6+nspecies+n) = 2.*bc_Xk_z_lo[n] - prim(i,j,2*lo-k-1,6+nspecies+n);
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            if (bc_therm_lo[2] == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(i,j,2*lo-k-1,n);
                        }
                    }
                });
            } else if (bc_therm_lo[2] == 2) { // isothermal
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        prim(i,j,k,4) = -prim(i,j,2*lo-k-1,4) + 2.*t_lo[2];
                        prim(i,j,k,5) = prim(i,j,2*lo-k-1,5);
                    }
                });
            }

        } // end LO Z
        

        // HI Z
        if ((dim == 2) && (bx.bigEnd(2) > dom.bigEnd(2))) {

            int hi = dom.bigEnd(2);
            
            // mass fractions, wall
            if ( bc_mass_hi[2] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(i,j,2*hi-k+1,n);

                        }
                    }
                });
            }

            // mass fractions, reservoir
            if (bc_mass_hi[2] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        for (int n=0; n<nspecies_gpu; ++n) {
                            prim(i,j,k,6+n)          = 2.*bc_Yk_z_hi[n] - prim(i,j,2*hi-k+1,6+n);
                            prim(i,j,k,6+nspecies+n) = 2.*bc_Xk_z_hi[n] - prim(i,j,2*hi-k+1,6+nspecies+n);
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            if (bc_therm_hi[2] == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(i,j,2*hi-k+1,n);
                        }
                    }
                });
            } else if (bc_therm_hi[2] == 2) { // isothermal
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        prim(i,j,k,4) = -prim(i,j,2*hi-k+1,4) + 2.*t_hi[2];
                        prim(i,j,k,5) = prim(i,j,2*hi-k+1,5);
                    }
                });
            }

        } // end HI Z
    }
}

// Set normal momemntum and velocity on the boundary and ghost cells of the 
// staggered grid based on slip/no-slip BCs
void BCMomNormal(MultiFab& mom_in, MultiFab& vel_in,
                 const amrex::Geometry geom, int dim)
{
    BL_PROFILE_VAR("BCMomNormal()",BCMomNormal);

    Box dom(geom.Domain());
    int ng_v = vel_in.nGrow();

    for ( MFIter mfi(vel_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng_v);
        
        const Array4<Real>& vel = vel_in.array(mfi);
        const Array4<Real>& mom = mom_in.array(mfi);
    
        // LO X
        if ((dim == 0) && (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) && (bx.smallEnd(0) <= dom.smallEnd(0))) {
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) {
                    // set ghost velocity & momentum
                    vel(i,j,k) = -vel(-i,j,k);
                    mom(i,j,k) = -mom(-i,j,k);
                }           
                else if (i == dom.smallEnd(0)) {
                    // set normal velocity & momentum
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }
                
        // HI X
        if ((dim == 0) && (bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) && (bx.bigEnd(0) >= dom.bigEnd(0)+1)) {

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {        
                if (i > dom.bigEnd(0)+1) {
                    // set ghost velocity & momentum
                    vel(i,j,k) = -vel(2*dom.bigEnd(0)+2-i,j,k);
                    mom(i,j,k) = -mom(2*dom.bigEnd(0)+2-i,j,k);
                }           
                else if (i == dom.bigEnd(0)+1) {
                    // set normal velocity & momentum
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }
        
        // LO Y
        if ((dim == 1) && (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) && (bx.smallEnd(1) <= dom.smallEnd(1))) {
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) {
                    // set ghost velocity & momentum
                    vel(i,j,k) = -vel(i,-j,k);
                    mom(i,j,k) = -mom(i,-j,k);
                }           
                else if (j == dom.smallEnd(1)) {
                    // set normal velocity & momentum
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }
                
        // HI Y
        if ((dim == 1) && (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) && (bx.bigEnd(1) >= dom.bigEnd(1)+1)) {

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {        
                if (j > dom.bigEnd(1)+1) {
                    // set ghost velocity & momentum
                    vel(i,j,k) = -vel(i,2*dom.bigEnd(1)+2-j,k);
                    mom(i,j,k) = -mom(i,2*dom.bigEnd(1)+2-j,k);
                }           
                else if (j == dom.bigEnd(1)+1) {
                    // set normal velocity & momentum
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }
        
        // LO Z 
        if ((dim == 2) && (bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) && (bx.smallEnd(2) <= dom.smallEnd(2))) {
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) {
                    // set ghost velocity & momentum
                    vel(i,j,k) = -vel(i,j,-k);
                    mom(i,j,k) = -mom(i,j,-k);
                }           
                else if (k == dom.smallEnd(2)) {
                    // set normal velocity & momentum
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }
                
        // HI Z
        if ((dim == 2) && (bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) && (bx.bigEnd(2) >= dom.bigEnd(2)+1)) {

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {        
                if (k > dom.bigEnd(2)+1) {
                    // set ghost velocity & momentum
                    vel(i,j,k) = -vel(i,j,2*dom.bigEnd(2)+2-k);
                    mom(i,j,k) = -mom(i,j,2*dom.bigEnd(2)+2-k);
                }           
                else if (k == dom.bigEnd(2)+1) {
                    // set normal velocity & momentum
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }
    }
}
        
        
// Set transverse momemntum and velocity on the boundary and ghost cells of the 
// staggered grid based on slip/no-slip BCs
void BCMomTrans(MultiFab& mom_in, MultiFab& vel_in,
                 const amrex::Geometry geom, int dim)
{
    BL_PROFILE_VAR("BCMomTrans()",BCMomTrans);
    
    if (geom.isAllPeriodic()) {
        return;
    }

    Box dom(geom.Domain());
    int ng_v = vel_in.nGrow();

    for ( MFIter mfi(vel_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng_v);
        
        const Array4<Real>& vel = vel_in.array(mfi);
        const Array4<Real>& mom = mom_in.array(mfi);

        // LO X
        if ((dim != 0) && (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) && (bx.smallEnd(0) < dom.smallEnd(0))) {
            const Real fac = (bc_vel_lo[0] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) {
                    vel(i,j,k) = fac*vel(-i-1,j,k);
                    mom(i,j,k) = fac*mom(-i-1,j,k);
                }
            });
        }

        // HI X
        if ((dim != 0) && (bc_vel_lo[0] == 1 || bc_vel_hi[0] == 2) && (bx.bigEnd(0) > dom.bigEnd(0))) {
            const Real fac = (bc_vel_hi[0] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > dom.bigEnd(0)) {
                    vel(i,j,k) = fac*vel(2*dom.bigEnd(0)-i+1,j,k);
                    mom(i,j,k) = fac*mom(2*dom.bigEnd(0)-i+1,j,k);
                }
            });
        }

        // LO Y
        if ((dim != 1) && (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) && (bx.smallEnd(1) < dom.smallEnd(1))) {
            const Real fac = (bc_vel_lo[1] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) {
                    vel(i,j,k) = fac*vel(i,-j-1,k);
                    mom(i,j,k) = fac*mom(i,-j-1,k);
                }
            });
        }

        // HI Y
        if ((dim != 1) && (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) && (bx.bigEnd(1) > dom.bigEnd(1))) {
            const Real fac = (bc_vel_hi[1] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j > dom.bigEnd(1)) {
                    vel(i,j,k) = fac*vel(i,2*dom.bigEnd(1)-j+1,k);
                    mom(i,j,k) = fac*mom(i,2*dom.bigEnd(1)-j+1,k);
                }
            });
        }

        // LO Z
        if ((dim != 2) && (bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) && (bx.smallEnd(2) < dom.smallEnd(2))) {
            const Real fac = (bc_vel_lo[2] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) {
                    vel(i,j,k) = fac*vel(i,j,-k-1);
                    mom(i,j,k) = fac*mom(i,j,-k-1);
                }
            });
        }

        // HI Z
        if ((dim != 2) && (bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) && (bx.bigEnd(2) > dom.bigEnd(2))) {
            const Real fac = (bc_vel_hi[2] == 1) ? 1. : -1.;
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k > dom.bigEnd(2)) {
                    vel(i,j,k) = fac*vel(i,j,2*dom.bigEnd(2)-k+1);
                    mom(i,j,k) = fac*mom(i,j,2*dom.bigEnd(2)-k+1);
                }
            });
        }
    }
}
        

// Set density and energy density on BCs 
void BCRhoRhoE(MultiFab& cons_in, MultiFab& prim_in, 
               std::array< MultiFab, AMREX_SPACEDIM >& cumom_in, 
               const amrex::Geometry geom)
{
    BL_PROFILE_VAR("BCRhoRhoE()",BCRhoRhoE);

    int nprimvars_gpu = nprimvars;
    GpuArray<Real,MAX_SPECIES> molmass_gpu;
    for (int n=0; n<nspecies; ++n) {
        molmass_gpu[n] = molmass[n];
    }
    GpuArray<Real,MAX_SPECIES> hcv_gpu;
    for (int n=0; n<nspecies; ++n) {
        hcv_gpu[n] = hcv[n];
    }
    Real Runiv_gpu = Runiv;
    int nspecies_gpu = nspecies;

    Box dom(geom.Domain());
    int ng_p = prim_in.nGrow();

    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng_p);

        const Array4<Real>& prim = prim_in.array(mfi);
        const Array4<Real>& cons = cons_in.array(mfi);
        AMREX_D_TERM(Array4<Real const> const& momx = cumom_in[0].array(mfi);,
                     Array4<Real const> const& momy = cumom_in[1].array(mfi);,
                     Array4<Real const> const& momz = cumom_in[2].array(mfi););
        
        // LO X
        if ((bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) && (bx.smallEnd(0) < dom.smallEnd(0))) {
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) {
                    
                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies_gpu; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }
                    
                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec,molmass_gpu,Runiv_gpu,nspecies_gpu);
                    GetEnergy(intenergy,fracvec,temp,hcv_gpu,nspecies_gpu);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies_gpu; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }
                    
                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);
                    
                    cons(i,j,k,4) = rho*intenergy + kinenergy; 
                }           
            });
        }
        
        // HI X
        if ((bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) && (bx.bigEnd(0) > dom.bigEnd(0))) {
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > dom.bigEnd(0)) {
                    
                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies_gpu; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }
                    
                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec,molmass_gpu,Runiv_gpu,nspecies_gpu);
                    GetEnergy(intenergy,fracvec,temp,hcv_gpu,nspecies_gpu);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies_gpu; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }
                    
                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);
                    
                    cons(i,j,k,4) = rho*intenergy + kinenergy; 
                }           
            });
        }
        
        // LO Y
        if ((bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) && (bx.smallEnd(1) < dom.smallEnd(1))) {
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) {
                    
                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies_gpu; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }
                    
                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec,molmass_gpu,Runiv_gpu,nspecies_gpu);
                    GetEnergy(intenergy,fracvec,temp,hcv_gpu,nspecies_gpu);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies_gpu; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }
                    
                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);
                    
                    cons(i,j,k,4) = rho*intenergy + kinenergy; 
                }           
            });
        }
        
        // HI Y
        if ((bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) && (bx.bigEnd(1) > dom.bigEnd(1))) {
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j > dom.bigEnd(1)) {
                    
                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies_gpu; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }
                    
                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec,molmass_gpu,Runiv_gpu,nspecies_gpu);
                    GetEnergy(intenergy,fracvec,temp,hcv_gpu,nspecies_gpu);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies_gpu; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }
                    
                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);
                    
                    cons(i,j,k,4) = rho*intenergy + kinenergy; 
                }           
            });
        }
        
        // LO Z 
        if ((bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) && (bx.smallEnd(2) < dom.smallEnd(2))) {
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) {
                    
                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies_gpu; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }
                    
                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec,molmass_gpu,Runiv_gpu,nspecies_gpu);
                    GetEnergy(intenergy,fracvec,temp,hcv_gpu,nspecies_gpu);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies_gpu; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }
                    
                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);
                    
                    cons(i,j,k,4) = rho*intenergy + kinenergy; 
                }           
            });
        }
        
        // HI Z
        if ((bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) && (bx.bigEnd(2) > dom.bigEnd(2))) {
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k > dom.bigEnd(2)) {
                    
                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies_gpu; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }
                    
                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec,molmass_gpu,Runiv_gpu,nspecies_gpu);
                    GetEnergy(intenergy,fracvec,temp,hcv_gpu,nspecies_gpu);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies_gpu; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }
                    
                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);
                    
                    cons(i,j,k,4) = rho*intenergy + kinenergy; 
                }           
            });
        }
    }
}

