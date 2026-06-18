#include "compressible_functions.H"

void SetupBC() {
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_vel_lo[i] == -1) {
            bc_mass_lo[i] = -1;
            bc_mass_hi[i] = -1;
            bc_therm_lo[i] = -1;
            bc_therm_hi[i] = -1;
        }
    }
}

void SetupCWall() {

    Real sumx, sumy;

    // Compute Xk or Yk at the wall, depending on which is defined
    // X walls
    if ((bc_mass_lo[0] == 2) or (bc_mass_lo[0] == 3)) {
       sumx = 0.;
       sumy = 0.;
       for (int ns=0; ns<nspecies; ++ns) {
          sumx = sumx + bc_Xk_x_lo[ns];
          sumy = sumy + bc_Yk_x_lo[ns];
       }
       if (amrex::Math::abs(sumx-1) < 1.e-10) {
           GetMassfrac(bc_Xk_x_lo,bc_Yk_x_lo);;
       } else if (amrex::Math::abs(sumy-1) < 1.e-10) {
          GetMolfrac(bc_Yk_x_lo,bc_Xk_x_lo);
       }
       else {
           Abort("SetupCWall: lo-x; mass or mole fractions do not sum to 1");
       }
    }

    if (bc_mass_lo[0] == 3) {
        // set reservoir pressure equal to inital ambient pressure (for no flow)
        // if t_lo/hi is positive, compute rho_lo/hi (default)
        // if rho_lo/hi is positive, rewrite t_lo/hi (from input script)

        if (p_lo[0] <= 0.0) { // set reservoir pressure to ambient if not specified
            GpuArray<Real,MAX_SPECIES> massvec;
            for (int ns=0;ns<nspecies;++ns) massvec[ns] = rhobar[ns];
            GetPressureGas(p_lo[0],massvec,rho0,T_init[0]);
        }

        if (rho_lo[0] < 0.0) { // specify reservoir density if not specified
            GetDensity(p_lo[0],rho_lo[0],t_lo[0],bc_Yk_x_lo);
        }
        else if (t_lo[0] <= 0.0) { // specify reservoir temperature if not specified
            Real molmix = 0.;
            for (int n=0; n<nspecies; ++n) {
                molmix += bc_Yk_x_lo[n]/molmass[n];
            }
            molmix = 1./molmix;
            t_lo[0] = p_lo[0]*(molmix/Runiv)/rho_lo[0];
        }
    }

    if ((bc_mass_hi[0] == 2) or (bc_mass_hi[0] == 3)) {
       sumx = 0.;
       sumy = 0.;
       for (int ns=0; ns<nspecies; ++ns) {
          sumx = sumx + bc_Xk_x_hi[ns];
          sumy = sumy + bc_Yk_x_hi[ns];
       }
       if (amrex::Math::abs(sumx-1) < 1.e-10) {
          GetMassfrac(bc_Xk_x_hi,bc_Yk_x_hi);
       } else if (amrex::Math::abs(sumy-1) < 1.e-10) {
          GetMolfrac(bc_Yk_x_hi,bc_Xk_x_hi);
       } else {
           Abort("SetupCWall: hi-x; mass or mole fractions do not sum to 1");
       }
    }

    if (bc_mass_hi[0] == 3) {
        // set reservoir pressure equal to inital ambient pressure (for no flow)
        // if t_lo/hi is positive, compute rho_lo/hi (default)
        // if rho_lo/hi is positive, rewrite t_lo/hi (from input script)

        if (p_hi[0] <= 0.0) {  // set reservoir pressure to ambient if not specified
            GpuArray<Real,MAX_SPECIES> massvec;
            for (int ns=0;ns<nspecies;++ns) massvec[ns] = rhobar[ns];
            GetPressureGas(p_hi[0],massvec,rho0,T_init[0]);
        }

        if (rho_hi[0] < 0.0) { // specify reservoir density  if not specified
            GetDensity(p_hi[0],rho_hi[0],t_hi[0],bc_Yk_x_hi);
        }
        else if (t_hi[0] <= 0.0) { // specify reservoir temperature if not specified
            Real molmix = 0.;
            for (int n=0; n<nspecies; ++n) {
                molmix += bc_Yk_x_hi[n]/molmass[n];
            }
            molmix = 1./molmix;
            t_hi[0] = p_hi[0]*(molmix/Runiv)/rho_hi[0];
        }
    }

    // Y walls
    if (bc_mass_lo[1] == 2) {
       sumx = 0.;
       sumy = 0.;
       for (int ns=0; ns<nspecies; ++ns) {
          sumx = sumx + bc_Xk_y_lo[ns];
          sumy = sumy + bc_Yk_y_lo[ns];
       }
       if (amrex::Math::abs(sumx-1) < 1.e-10) {
          GetMassfrac(bc_Xk_y_lo,bc_Yk_y_lo);
       } else if (amrex::Math::abs(sumy-1) < 1.e-10) {
          GetMolfrac(bc_Yk_y_lo,bc_Xk_y_lo);
       } else {
           Abort("SetupCWall: lo-y; mass or mole fractions do not sum to 1");
       }
    }

    if (bc_mass_hi[1] == 2) {
       sumx = 0.;
       sumy = 0.;
       for (int ns=0; ns<nspecies; ++ns) {
          sumx = sumx + bc_Xk_y_hi[ns];
          sumy = sumy + bc_Yk_y_hi[ns];
       }
       if (amrex::Math::abs(sumx-1) < 1.e-10) {
          GetMassfrac(bc_Xk_y_hi,bc_Yk_y_hi);
       } else if (amrex::Math::abs(sumy-1) < 1.e-10) {
          GetMolfrac(bc_Yk_y_hi,bc_Xk_y_hi);
       } else {
           Abort("SetupCWall: hi-y; mass or mole fractions do not sum to 1");
       }
    }

    // Z walls
    if (bc_mass_lo[2] == 2) {
       sumx = 0.;
       sumy = 0.;
       for (int ns=0; ns<nspecies; ++ns) {
          sumx = sumx + bc_Xk_z_lo[ns];
          sumy = sumy + bc_Yk_z_lo[ns];
       }
       if (amrex::Math::abs(sumx-1) < 1.e-10) {
          GetMassfrac(bc_Xk_z_lo,bc_Yk_z_lo);
       } else if (amrex::Math::abs(sumy-1) < 1.e-10) {
          GetMolfrac(bc_Yk_z_lo,bc_Xk_z_lo);
       } else {
           Abort("SetupCWall: lo-z; mass or mole fractions do not sum to 1");
       }
    }

    if (bc_mass_hi[2] == 2) {
       sumx = 0.;
       sumy = 0.;
       for (int ns=0; ns<nspecies; ++ns) {
          sumx = sumx + bc_Xk_z_hi[ns];
          sumy = sumy + bc_Yk_z_hi[ns];
       }
       if (amrex::Math::abs(sumx-1) < 1.e-10) {
          GetMassfrac(bc_Xk_z_hi,bc_Yk_z_hi);
       } else if (amrex::Math::abs(sumy-1) < 1.e-10) {
          GetMolfrac(bc_Yk_z_hi,bc_Xk_z_hi);
       } else {
           Abort("SetupCWall: hi-z; mass or mole fractions do not sum to 1");
       }
    }
}

// in this routine for Dirichlet BC, we set the
// value in the ghost cell to be the value of the
// Dirichlet boundary for temperature, species
// fraction and velocities
void setBC(MultiFab& prim_in, MultiFab& cons_in)
{
    BL_PROFILE_VAR("setBC()",setBC);

    int ng_c = cons_in.nGrow();
    int ng_p = prim_in.nGrow();
    if (ng_c != ng_p) {
        Abort("setBC: prim and cons need the same number of ghost cells");
    }

    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& gbx = mfi.growntilebox(ng_c);

        const Array4<Real> cons = cons_in.array(mfi);
        const Array4<Real> prim = prim_in.array(mfi);

        // LO X
        if (gbx.smallEnd(0) < 0) {

            int lo = 0;

            // mass fractions, wall
            if ( bc_mass_lo[0] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < 0) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(2*lo-i-1,j,k,n);
                        }
                    }
                });
            }

            // mass fractions, concentration
            if (bc_mass_lo[0] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < 0) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_x_lo[n];
                            prim(i,j,k,6+nspecies+n) = bc_Xk_x_lo[n];
                        }
                    }
                });
            }

            // mass fractions, density, temperature and pressure in the reservoir
            if (bc_mass_lo[0] == 3 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < 0) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_x_lo[n]; // set ghost cell equal to reservoir mass fraction
                            prim(i,j,k,6+nspecies+n) = bc_Xk_x_lo[n]; // set ghost cell equal to reservoir mole fraction

                            prim(i,j,k,0) = rho_lo[0]; // set ghost cell equal to reservoir density
                            cons(i,j,k,0) = rho_lo[0]; // set ghost cell equal to reservoir density

                            prim(i,j,k,4) = t_lo[0]; // set ghost cell equal to reservoir temperature
                            prim(i,j,k,5) = p_lo[0]; // set ghost cell equal to reservoir pressure
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            else if (bc_therm_lo[0] == 1) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < 0) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(2*lo-i-1,j,k,n);
                        }
                    }
                });
            } else if (bc_therm_lo[0] == 2) { // isothermal
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < 0) {
                        prim(i,j,k,4) = t_lo[0];
                        prim(i,j,k,5) = prim(2*lo-i-1,j,k,5);
                    }
                });
            }

            // momentum, velocity, rho, rhoY, rhoE
            if (bc_mass_lo[0] == 3) { // reservoir
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < 0) {

                        cons(i,j,k,1) = 0.0;
                        cons(i,j,k,2) = 0.0;
                        cons(i,j,k,3) = 0.0;

                        prim(i,j,k,1) = 0.0;
                        prim(i,j,k,2) = 0.0;
                        prim(i,j,k,3) = 0.0;

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        // Real pt = prim(i,j,k,5);
                        Real rho = prim(i,j,k,0);
                        Real intenergy;

                        GetEnergy(intenergy,fracvec,temp);

                        // total density depends on temperature
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            }
            else if (bc_vel_lo[0] == 1) { // slip
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < 0) {

                        cons(i,j,k,1) = 0.0;
                        cons(i,j,k,2) = cons(2*lo-i-1,j,k,2);
                        cons(i,j,k,3) = cons(2*lo-i-1,j,k,3);

                        prim(i,j,k,1) = 0.0;
                        prim(i,j,k,2) = prim(2*lo-i-1,j,k,2);
                        prim(i,j,k,3) = prim(2*lo-i-1,j,k,3);

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        Real pt = prim(i,j,k,5);
                        Real rho;
                        Real intenergy;

                        GetDensity(pt,rho,temp,fracvec);
                        GetEnergy(intenergy,fracvec,temp);

                        // total density depends on temperature
                        prim(i,j,k,0) = rho;
                        cons(i,j,k,0) = rho;
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            } else if (bc_vel_lo[0] == 2) { // no slip
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < 0) {

                        prim(i,j,k,1) = 0.0;
                        prim(i,j,k,2) = 0.0;
                        prim(i,j,k,3) = 0.0;

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        Real pt = prim(i,j,k,5);
                        Real rho;
                        Real intenergy;

                        GetDensity(pt,rho,temp,fracvec);
                        GetEnergy(intenergy,fracvec,temp);

                        prim(i,j,k,0) = rho;
                        cons(i,j,k,0) = rho;
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        cons(i,j,k,1) = rho*prim(i,j,k,1);
                        cons(i,j,k,2) = rho*prim(i,j,k,2);
                        cons(i,j,k,3) = rho*prim(i,j,k,3);

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            }
        } // end LO X


        // HI X
        if (gbx.bigEnd(0) > n_cells[0]-1) {

            int hi = n_cells[0]-1;

            // mass fractions, wall
            if ( bc_mass_hi[0] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i >= n_cells[0]) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(2*hi-i+1,j,k,n);

                        }
                    }
                });
            }

            // mass fractions, concentration
            if (bc_mass_hi[0] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > n_cells[0]-1) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_x_hi[n];
                            prim(i,j,k,6+nspecies+n) = bc_Xk_x_hi[n];
                        }
                    }
                });
            }

            // mass fractions, density, temperature and pressure in the reservoir
            if (bc_mass_hi[0] == 3 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > n_cells[0]-1) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_x_hi[n]; // set ghost cell equal to reservoir mass fraction
                            prim(i,j,k,6+nspecies+n) = bc_Xk_x_hi[n]; // set ghost cell equal to reservoir mole fraction

                            prim(i,j,k,0) = rho_hi[0]; // set ghost cell equal to reservoir density
                            cons(i,j,k,0) = rho_hi[0]; // set ghost cell equal to reservoir density

                            prim(i,j,k,4) = t_hi[0]; // set ghost cell equal to reservoir temperature
                            prim(i,j,k,5) = p_hi[0]; // set ghost cell equal to reservoir pressure
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            else if (bc_therm_hi[0] == 1) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > n_cells[0]-1) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(2*hi-i+1,j,k,n);
                        }
                    }
                });
            } else if (bc_therm_hi[0] == 2) { // isothermal
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > n_cells[0]-1) {
                        prim(i,j,k,4) = t_hi[0];
                        prim(i,j,k,5) = prim(2*hi-i+1,j,k,5);
                    }
                });
            }

            // momentum, velocity, rho, rhoY, rhoE
            if (bc_mass_hi[0] == 3) { // reservoir
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > n_cells[0]-1) {


                        cons(i,j,k,1) = 0.0;
                        cons(i,j,k,2) = 0.0;
                        cons(i,j,k,3) = 0.0;

                        prim(i,j,k,1) = 0.0;
                        prim(i,j,k,2) = 0.0;
                        prim(i,j,k,3) = 0.0;

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        // Real pt = prim(i,j,k,5);
                        Real rho = prim(i,j,k,0);
                        Real intenergy;

                        GetEnergy(intenergy,fracvec,temp);

                        // total density depends on temperature
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            } else if (bc_vel_hi[0] == 1) { // slip
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > n_cells[0]-1) {

                        cons(i,j,k,1) = 0.0;
                        cons(i,j,k,2) = cons(2*hi-i+1,j,k,2);
                        cons(i,j,k,3) = cons(2*hi-i+1,j,k,3);

                        prim(i,j,k,1) = 0.0;
                        prim(i,j,k,2) = prim(2*hi-i+1,j,k,2);
                        prim(i,j,k,3) = prim(2*hi-i+1,j,k,3);

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        Real pt = prim(i,j,k,5);
                        Real rho;
                        Real intenergy;

                        GetDensity(pt,rho,temp,fracvec);
                        GetEnergy(intenergy,fracvec,temp);

                        // total density depends on temperature
                        prim(i,j,k,0) = rho;
                        cons(i,j,k,0) = rho;
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            } else if (bc_vel_hi[0] == 2) { // no slip
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > n_cells[0]-1) {

                        prim(i,j,k,1) = 0.0;
                        prim(i,j,k,2) = 0.0;
                        prim(i,j,k,3) = 0.0;

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        Real pt = prim(i,j,k,5);
                        Real rho;
                        Real intenergy;

                        GetDensity(pt,rho,temp,fracvec);
                        GetEnergy(intenergy,fracvec,temp);

                        prim(i,j,k,0) = rho;
                        cons(i,j,k,0) = rho;
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        cons(i,j,k,1) = rho*prim(i,j,k,1);
                        cons(i,j,k,2) = rho*prim(i,j,k,2);
                        cons(i,j,k,3) = rho*prim(i,j,k,3);

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            }
        } // end HI X

        // LO Y
        if (gbx.smallEnd(1) < 0) {

            int lo = 0;

            // mass fractions, wall
            if ( bc_mass_lo[1] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < 0) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(i,2*lo-j-1,k,n);
                        }
                    }
                });
            }

            // mass fractions, concentration
            if (bc_mass_lo[1] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < 0) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_y_lo[n];
                            prim(i,j,k,6+nspecies+n) = bc_Xk_y_lo[n];
                        }
                    }
                });
            }

            // mass fractions, density, temperature and pressure in the reservoir
            if (bc_mass_lo[1] == 3 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < 0) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_y_lo[n]; // set ghost cell equal to reservoir mass fraction
                            prim(i,j,k,6+nspecies+n) = bc_Xk_y_lo[n]; // set ghost cell equal to reservoir mole fraction

                            prim(i,j,k,0) = rho_lo[1]; // set ghost cell equal to reservoir density
                            cons(i,j,k,0) = rho_lo[1]; // set ghost cell equal to reservoir density

                            prim(i,j,k,4) = t_lo[1]; // set ghost cell equal to reservoir temperature
                            prim(i,j,k,5) = p_lo[1]; // set ghost cell equal to reservoir pressure
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            if (bc_therm_lo[1] == 1) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < 0) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(i,2*lo-j-1,k,n);
                        }
                    }
                });
            } else if (bc_therm_lo[1] == 2) { // isothermal
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < 0) {
                        prim(i,j,k,4) = t_lo[1];
                        prim(i,j,k,5) = prim(i,2*lo-j-1,k,5);
                    }
                });
            }

            // momentum, velocity, rho, rhoY, rhoE
            if (bc_mass_lo[1] == 3) { // reservoir

            } else if (bc_vel_lo[1] == 1) { // slip
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < 0) {

                        cons(i,j,k,1) = cons(i,2*lo-j-1,k,1);
                        cons(i,j,k,2) = 0.0;
                        cons(i,j,k,3) = cons(i,2*lo-j-1,k,3);

                        prim(i,j,k,1) = prim(i,2*lo-j-1,k,1);
                        prim(i,j,k,2) = 0.0;
                        prim(i,j,k,3) = prim(i,2*lo-j-1,k,3);

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        Real pt = prim(i,j,k,5);
                        Real rho;
                        Real intenergy;

                        GetDensity(pt,rho,temp,fracvec);
                        GetEnergy(intenergy,fracvec,temp);

                        // total density depends on temperature
                        prim(i,j,k,0) = rho;
                        cons(i,j,k,0) = rho;
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            } else if (bc_vel_lo[1] == 2) { // no slip
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < 0) {

                        prim(i,j,k,1) = 0.0;
                        prim(i,j,k,2) = 0.0;
                        prim(i,j,k,3) = 0.0;

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        Real pt = prim(i,j,k,5);
                        Real rho;
                        Real intenergy;

                        GetDensity(pt,rho,temp,fracvec);
                        GetEnergy(intenergy,fracvec,temp);

                        prim(i,j,k,0) = rho;
                        cons(i,j,k,0) = rho;
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        cons(i,j,k,1) = rho*prim(i,j,k,1);
                        cons(i,j,k,2) = rho*prim(i,j,k,2);
                        cons(i,j,k,3) = rho*prim(i,j,k,3);

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            }
        } // end LO Y


        // HI Y
        if (gbx.bigEnd(1) > n_cells[1]-1) {

            int hi = n_cells[1]-1;

            // mass fractions, wall
            if ( bc_mass_hi[1] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j >= n_cells[1]) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(i,2*hi-j+1,k,n);

                        }
                    }
                });
            }

            // mass fractions, concentration
            if (bc_mass_hi[1] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > n_cells[1]-1) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_y_hi[n];
                            prim(i,j,k,6+nspecies+n) = bc_Xk_y_hi[n];
                        }
                    }
                });
            }

            // mass fractions, density, temperature and pressure in the reservoir
            if (bc_mass_hi[1] == 3 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > n_cells[1]-1) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_y_hi[n]; // set ghost cell equal to reservoir mass fraction
                            prim(i,j,k,6+nspecies+n) = bc_Xk_y_hi[n]; // set ghost cell equal to reservoir mole fraction

                            prim(i,j,k,0) = rho_hi[1]; // set ghost cell equal to reservoir density
                            cons(i,j,k,0) = rho_hi[1]; // set ghost cell equal to reservoir density

                            prim(i,j,k,4) = t_hi[1]; // set ghost cell equal to reservoir temperature
                            prim(i,j,k,5) = p_hi[1]; // set ghost cell equal to reservoir pressure
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            if (bc_therm_hi[1] == 1) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > n_cells[1]-1) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(i,2*hi-j+1,k,n);
                        }
                    }
                });
            } else if (bc_therm_hi[1] == 2) { // isothermal
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > n_cells[1]-1) {
                        prim(i,j,k,4) = t_hi[1];
                        prim(i,j,k,5) = prim(i,2*hi-j+1,k,5);
                    }
                });
            }

            // momentum, velocity, rho, rhoY, rhoE
            if (bc_mass_hi[1] == 3) { // reservoir

            } else if (bc_vel_hi[1] == 1) { // slip
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > n_cells[1]-1) {

                        cons(i,j,k,1) = cons(i,2*hi-j+1,k,1);
                        cons(i,j,k,2) = 0.0;
                        cons(i,j,k,3) = cons(i,2*hi-j+1,k,3);

                        prim(i,j,k,1) = prim(i,2*hi-j+1,k,1);
                        prim(i,j,k,2) = 0.0;
                        prim(i,j,k,3) = prim(i,2*hi-j+1,k,3);

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        Real pt = prim(i,j,k,5);
                        Real rho;
                        Real intenergy;

                        GetDensity(pt,rho,temp,fracvec);
                        GetEnergy(intenergy,fracvec,temp);

                        // total density depends on temperature
                        prim(i,j,k,0) = rho;
                        cons(i,j,k,0) = rho;
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            } else if (bc_vel_hi[1] == 2) { // no slip

                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > n_cells[1]-1) {

                        prim(i,j,k,1) = 0.0;
                        prim(i,j,k,2) = 0.0;
                        prim(i,j,k,3) = 0.0;

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        Real pt = prim(i,j,k,5);
                        Real rho;
                        Real intenergy;

                        GetDensity(pt,rho,temp,fracvec);
                        GetEnergy(intenergy,fracvec,temp);

                        prim(i,j,k,0) = rho;
                        cons(i,j,k,0) = rho;
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        cons(i,j,k,1) = rho*prim(i,j,k,1);
                        cons(i,j,k,2) = rho*prim(i,j,k,2);
                        cons(i,j,k,3) = rho*prim(i,j,k,3);

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            }
        } // end HI Y

        // LO Z
        if (gbx.smallEnd(2) < 0) {

            int lo = 0;

            // mass fractions, wall
            if ( bc_mass_lo[2] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < 0) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(i,j,2*lo-k-1,n);
                        }
                    }
                });
            }

            // mass fractions, concentration
            if (bc_mass_lo[2] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < 0) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_z_lo[n];
                            prim(i,j,k,6+nspecies+n) = bc_Xk_z_lo[n];
                        }
                    }
                });
            }

            // mass fractions, density, temperature and pressure in the reservoir
            if (bc_mass_lo[2] == 3 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < 0) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_z_lo[n]; // set ghost cell equal to reservoir mass fraction
                            prim(i,j,k,6+nspecies+n) = bc_Xk_z_lo[n]; // set ghost cell equal to reservoir mole fraction

                            prim(i,j,k,0) = rho_lo[2]; // set ghost cell equal to reservoir density
                            cons(i,j,k,0) = rho_lo[2]; // set ghost cell equal to reservoir density

                            prim(i,j,k,4) = t_lo[2]; // set ghost cell equal to reservoir temperature
                            prim(i,j,k,5) = p_lo[2]; // set ghost cell equal to reservoir pressure
                        }
                    }
                });
            }


            // temperature and pressure, adiabatic
            if (bc_therm_lo[2] == 1) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < 0) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(i,j,2*lo-k-1,n);
                        }
                    }
                });
            } else if (bc_therm_lo[2] == 2) { // isothermal
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < 0) {
                        prim(i,j,k,4) = t_lo[2];
                        prim(i,j,k,5) = prim(i,j,2*lo-k-1,5);
                    }
                });
            }

            // momentum, velocity, rho, rhoY, rhoE
            if (bc_mass_lo[2] == 3) { // reservoir

            } else if (bc_vel_lo[2] == 1) { // slip
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < 0) {

                        cons(i,j,k,1) = cons(i,j,2*lo-k-1,1);
                        cons(i,j,k,2) = cons(i,j,2*lo-k-1,2);
                        cons(i,j,k,3) = 0.0;

                        prim(i,j,k,1) = prim(i,j,2*lo-k-1,1);
                        prim(i,j,k,2) = prim(i,j,2*lo-k-1,2);
                        prim(i,j,k,3) = 0.0;;

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        Real pt = prim(i,j,k,5);
                        Real rho;
                        Real intenergy;

                        GetDensity(pt,rho,temp,fracvec);
                        GetEnergy(intenergy,fracvec,temp);

                        // total density depends on temperature
                        prim(i,j,k,0) = rho;
                        cons(i,j,k,0) = rho;
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            } else if (bc_vel_lo[2] == 2) { // no slip
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < 0) {

                        prim(i,j,k,1) = 0.0;
                        prim(i,j,k,2) = 0.0;
                        prim(i,j,k,3) = 0.0;

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        Real pt = prim(i,j,k,5);
                        Real rho;
                        Real intenergy;

                        GetDensity(pt,rho,temp,fracvec);
                        GetEnergy(intenergy,fracvec,temp);

                        prim(i,j,k,0) = rho;
                        cons(i,j,k,0) = rho;
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        cons(i,j,k,1) = rho*prim(i,j,k,1);
                        cons(i,j,k,2) = rho*prim(i,j,k,2);
                        cons(i,j,k,3) = rho*prim(i,j,k,3);

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            }
        } // end LO Z


        // HI Z
        if (gbx.bigEnd(2) > n_cells[2]-1) {

            int hi = n_cells[2]-1;

            // mass fractions, wall
            if ( bc_mass_hi[2] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k >= n_cells[2]) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(i,j,2*hi-k+1,n);

                        }
                    }
                });
            }

            // mass fractions, concentration
            if (bc_mass_hi[2] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > n_cells[2]-1) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_z_hi[n];
                            prim(i,j,k,6+nspecies+n) = bc_Xk_z_hi[n];
                        }
                    }
                });
            }

            // mass fractions, density, temperature and pressure in the reservoir
            if (bc_mass_hi[2] == 3 && algorithm_type == 2) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > n_cells[2]-1) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_z_hi[n]; // set ghost cell equal to reservoir mass fraction
                            prim(i,j,k,6+nspecies+n) = bc_Xk_z_hi[n]; // set ghost cell equal to reservoir mole fraction

                            prim(i,j,k,0) = rho_hi[2]; // set ghost cell equal to reservoir density
                            cons(i,j,k,0) = rho_hi[2]; // set ghost cell equal to reservoir density

                            prim(i,j,k,4) = t_hi[1]; // set ghost cell equal to reservoir temperature
                            prim(i,j,k,5) = p_hi[1]; // set ghost cell equal to reservoir pressure
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            if (bc_therm_hi[2] == 1) {
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > n_cells[2]-1) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(i,j,2*hi-k+1,n);
                        }
                    }
                });
            } else if (bc_therm_hi[2] == 2) { // isothermal
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > n_cells[2]-1) {
                        prim(i,j,k,4) = t_hi[2];
                        prim(i,j,k,5) = prim(i,j,2*hi-k+1,5);
                    }
                });
            }

            // momentum, velocity, rho, rhoY, rhoE
            if (bc_mass_hi[2] == 3) { // reservoir

            } else if (bc_vel_hi[2] == 1) { // slip
                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > n_cells[2]-1) {

                        cons(i,j,k,1) = cons(i,j,2*hi-k+1,1);
                        cons(i,j,k,2) = cons(i,j,2*hi-k+1,2);
                        cons(i,j,k,3) = 0.0;

                        prim(i,j,k,1) = prim(i,j,2*hi-k+1,1);
                        prim(i,j,k,2) = prim(i,j,2*hi-k+1,2);
                        prim(i,j,k,3) = 0.0;

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        Real pt = prim(i,j,k,5);
                        Real rho;
                        Real intenergy;

                        GetDensity(pt,rho,temp,fracvec);
                        GetEnergy(intenergy,fracvec,temp);

                        // total density depends on temperature
                        prim(i,j,k,0) = rho;
                        cons(i,j,k,0) = rho;
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            } else if (bc_vel_hi[2] == 2) { // no slip

                amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > n_cells[2]-1) {

                        prim(i,j,k,1) = 0.0;
                        prim(i,j,k,2) = 0.0;
                        prim(i,j,k,3) = 0.0;

                        // thermal & species (+pressure) BCs must be enforced first
                        GpuArray<Real,MAX_SPECIES> fracvec;
                        for (int n=0; n<nspecies; ++n) {
                            fracvec[n] = prim(i,j,k,6+n);
                        }
                        Real temp = prim(i,j,k,4);
                        Real pt = prim(i,j,k,5);
                        Real rho;
                        Real intenergy;

                        GetDensity(pt,rho,temp,fracvec);
                        GetEnergy(intenergy,fracvec,temp);

                        prim(i,j,k,0) = rho;
                        cons(i,j,k,0) = rho;
                        if (algorithm_type == 2) {
                            for (int n=0; n<nspecies; ++n) {
                                cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                            }
                        }

                        cons(i,j,k,1) = rho*prim(i,j,k,1);
                        cons(i,j,k,2) = rho*prim(i,j,k,2);
                        cons(i,j,k,3) = rho*prim(i,j,k,3);

                        // must be last BC enforced: depends on rho, vel, & temp
                        cons(i,j,k,4) = rho*intenergy + 0.5*rho*(prim(i,j,k,1)*prim(i,j,k,1) +
                                                                 prim(i,j,k,2)*prim(i,j,k,2) +
                                                                 prim(i,j,k,3)*prim(i,j,k,3));
                    }
                });
            }
        } // end HI Z
    }
}

// set species and total density flux to zero for wall boundary conditions
void BCWallSpeciesFlux(std::array< MultiFab, AMREX_SPACEDIM >& faceflux, const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("BCWallSpeciesFlux()",BCWallSpeciesFlux);

    // LO X
    if (bc_mass_lo[0] == 1) {

        // domain grown nodally based on faceflux[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xlo = amrex::bdryNode(dom_x, Orientation(0, Orientation::low));

        for (MFIter mfi(faceflux[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xlo;
            Array4<Real> const& flux = (faceflux[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) = 0.;
                    }
                    // density
                    flux(i,j,k,0) = 0.;
                    // Dufour
                    flux(i,j,k,nvars+3) = 0.;
                });
            }
        }
    }
    // HI X
    if (bc_mass_hi[0] == 1) {

        // domain grown nodally based on faceflux[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux[0].ixType());

        // this is the x-hi domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xhi = amrex::bdryNode(dom_x, Orientation(0, Orientation::high));

        for (MFIter mfi(faceflux[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xhi;
            Array4<Real> const& flux = (faceflux[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) = 0.;
                    }
                    // density
                    flux(i,j,k,0) = 0.;
                    // Dufour
                    flux(i,j,k,nvars+3) = 0.;
                });
            }
        }
    }
    // LO Y
    if (bc_mass_lo[1] == 1) {

        // domain grown nodally based on faceflux[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux[1].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_ylo = amrex::bdryNode(dom_y, Orientation(1, Orientation::low));

        for (MFIter mfi(faceflux[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_ylo;
            Array4<Real> const& flux = (faceflux[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) = 0.;
                    }
                    // density
                    flux(i,j,k,0) = 0.;
                    // Dufour
                    flux(i,j,k,nvars+3) = 0.;
                });
            }
        }
    }
    // HI Y
    if (bc_mass_hi[1] == 1) {

        // domain grown nodally based on faceflux[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux[1].ixType());

        // this is the y-hi domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yhi = amrex::bdryNode(dom_y, Orientation(1, Orientation::high));

        for (MFIter mfi(faceflux[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yhi;
            Array4<Real> const& flux = (faceflux[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) = 0.;
                    }
                    // density
                    flux(i,j,k,0) = 0.;
                    // Dufour
                    flux(i,j,k,nvars+3) = 0.;
                });
            }
        }
    }
    // LO Z
    if (bc_mass_lo[2] == 1) {

        // domain grown nodally based on faceflux[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux[2].ixType());

        // this is the z-lo domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zlo = amrex::bdryNode(dom_z, Orientation(2, Orientation::low));

        for (MFIter mfi(faceflux[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zlo;
            Array4<Real> const& flux = (faceflux[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) = 0.;
                    }
                    // density
                    flux(i,j,k,0) = 0.;
                    // Dufour
                    flux(i,j,k,nvars+3) = 0.;
                });
            }
        }
    }
    // HI Z
    if (bc_mass_hi[2] == 1) {

        // domain grown nodally based on faceflux[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux[2].ixType());

        // this is the z-hi domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zhi = amrex::bdryNode(dom_z, Orientation(2, Orientation::high));

        for (MFIter mfi(faceflux[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zhi;
            Array4<Real> const& flux = (faceflux[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) = 0.;
                    }
                    // density
                    flux(i,j,k,0) = 0.;
                    // Dufour
                    flux(i,j,k,nvars+3) = 0.;
                });
            }
        }
    }
}

void StochFlux(std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in,
               const amrex::Geometry& geom) {


    BL_PROFILE_VAR("StochFlux()",StochFlux);

    // First we do mass boundary conditions (species fluxes reside on faces)
    // LO X
    if (bc_mass_lo[0] == 1 || bc_mass_lo[0] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_lo[0] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_in[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xlo = amrex::bdryNode(dom_x, Orientation(0, Orientation::low));

        for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xlo;
            Array4<Real> const& flux = (faceflux_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) *= factor;
                    }
                    // Set Dufour as well
                    flux(i,j,k,nvars+3) *= factor;
                });
            }
        }
    }
    // HI X
    if (bc_mass_hi[0] == 1 || bc_mass_hi[0] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_hi[0] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_in[0].ixType());

        // this is the x-hi domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xhi = amrex::bdryNode(dom_x, Orientation(0, Orientation::high));

        for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xhi;
            Array4<Real> const& flux = (faceflux_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) *= factor;
                    }
                    // Set Dufour as well
                    flux(i,j,k,nvars+3) *= factor;
                });
            }
        }
    }
    // LO Y
    if (bc_mass_lo[1] == 1 || bc_mass_lo[1] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_lo[1] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_in[1].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_ylo = amrex::bdryNode(dom_y, Orientation(1, Orientation::low));

        for (MFIter mfi(faceflux_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_ylo;
            Array4<Real> const& flux = (faceflux_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) *= factor;
                    }
                    // Set Dufour as well
                    flux(i,j,k,nvars+3) *= factor;
                });
            }
        }
    }
    // HI Y
    if (bc_mass_hi[1] == 1 || bc_mass_hi[1] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_hi[1] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_in[1].ixType());

        // this is the y-hi domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yhi = amrex::bdryNode(dom_y, Orientation(1, Orientation::high));

        for (MFIter mfi(faceflux_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yhi;
            Array4<Real> const& flux = (faceflux_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) *= factor;
                    }
                    // Set Dufour as well
                    flux(i,j,k,nvars+3) *= factor;
                });
            }
        }
    }
    // LO Z
    if (bc_mass_lo[2] == 1 || bc_mass_lo[2] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_lo[2] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_in[2].ixType());

        // this is the z-lo domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zlo = amrex::bdryNode(dom_z, Orientation(2, Orientation::low));

        for (MFIter mfi(faceflux_in[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zlo;
            Array4<Real> const& flux = (faceflux_in[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) *= factor;
                    }
                    // Set Dufour as well
                    flux(i,j,k,nvars+3) *= factor;
                });
            }
        }
    }
    // HI Z
    if (bc_mass_hi[2] == 1 || bc_mass_hi[2] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = reservoir   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_hi[2] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_in[2].ixType());

        // this is the z-hi domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zhi = amrex::bdryNode(dom_z, Orientation(2, Orientation::high));

        for (MFIter mfi(faceflux_in[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zhi;
            Array4<Real> const& flux = (faceflux_in[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) *= factor;
                    }
                    // Set Dufour as well
                    flux(i,j,k,nvars+3) *= factor;
                });
            }
        }
    }

    // Next we do thermal boundary conditions (energy fluxes reside on faces)
    // LO X
    if (bc_therm_lo[0] == 1 || bc_therm_lo[0] == 2) {

        // 1 = adiabatic        : multiply fluxes on wall by 0
        // 2 = isothermal   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_therm_lo[0] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_lo[0] == 3) factor = 1.0;

        // domain grown nodally based on faceflux_in[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_in[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xlo = amrex::bdryNode(dom_x, Orientation(0, Orientation::low));

        for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xlo;
            Array4<Real> const& flux = (faceflux_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // heat flux
                    flux(i,j,k,nvars) *= factor;
                });
            }
        }
    }
    // HI X
    if (bc_therm_hi[0] == 1 || bc_therm_hi[0] == 2) {

        // 1 = adiabatic        : multiply fluxes on wall by 0
        // 2 = isothermal   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_therm_hi[0] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_hi[0] == 3) factor = 1.0;

        // domain grown nodally based on faceflux_in[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_in[0].ixType());

        // this is the x-hi domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xhi = amrex::bdryNode(dom_x, Orientation(0, Orientation::high));

        for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xhi;
            Array4<Real> const& flux = (faceflux_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // heat flux
                    flux(i,j,k,nvars) *= factor;
                });
            }
        }
    }
    // LO Y
    if (bc_therm_lo[1] == 1 || bc_therm_lo[1] == 2) {

        // 1 = adiabatic        : multiply fluxes on wall by 0
        // 2 = isothermal   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_therm_lo[1] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_in[1].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_ylo = amrex::bdryNode(dom_y, Orientation(1, Orientation::low));

        for (MFIter mfi(faceflux_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_ylo;
            Array4<Real> const& flux = (faceflux_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // heat flux
                    flux(i,j,k,nvars) *= factor;
                });
            }
        }
    }
    // HI Y
    if (bc_therm_hi[1] == 1 || bc_therm_hi[1] == 2) {

        // 1 = adiabatic        : multiply fluxes on wall by 0
        // 2 = isothermal   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_therm_hi[1] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_in[1].ixType());

        // this is the y-hi domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yhi = amrex::bdryNode(dom_y, Orientation(1, Orientation::high));

        for (MFIter mfi(faceflux_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yhi;
            Array4<Real> const& flux = (faceflux_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // heat flux
                    flux(i,j,k,nvars) *= factor;
                });
            }
        }
    }
    // LO Z
    if (bc_therm_lo[2] == 1 || bc_therm_lo[2] == 2) {

        // 1 = adiabatic        : multiply fluxes on wall by 0
        // 2 = isothermal   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_therm_lo[2] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_in[2].ixType());

        // this is the z-lo domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zlo = amrex::bdryNode(dom_z, Orientation(2, Orientation::low));

        for (MFIter mfi(faceflux_in[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zlo;
            Array4<Real> const& flux = (faceflux_in[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // heat flux
                    flux(i,j,k,nvars) *= factor;
                });
            }
        }
    }
    // HI Z
    if (bc_therm_hi[2] == 1 || bc_therm_hi[2] == 2) {

        // 1 = adiabatic        : multiply fluxes on wall by 0
        // 2 = isothermal   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_therm_hi[2] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_in[2].ixType());

        // this is the z-hi domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zhi = amrex::bdryNode(dom_z, Orientation(2, Orientation::high));

        for (MFIter mfi(faceflux_in[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zhi;
            Array4<Real> const& flux = (faceflux_in[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    flux(i,j,k,nvars) *= factor;
                });
            }
        }
    }


    Real sqrtTwo = sqrt(2.);

    // Last we do velocity boundary conditions
    // LO X
    if (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) {

        // always multiply normal velocity by sqrt(2)
        // 1 = slip wall:    multiply transverse velocity by 0
        // 2 = no-slip wall: multiply transverse velocity by sqrt(2)
        Real factor = (bc_vel_lo[0] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_in[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xlo = amrex::bdryNode(dom_x, Orientation(0, Orientation::low));

        for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xlo;
            Array4<Real> const& flux = (faceflux_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, AMREX_SPACEDIM, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (n == 0) {
                        // normal velocity
                        flux(i,j,k,1+n) *= sqrtTwo;
                        // viscous heating (diagonal & shear)
                        flux(i,j,k,nvars+1) *= sqrtTwo;
                        flux(i,j,k,nvars+2) *= factor;
                    } else {
                        // transverse velocity
                        flux(i,j,k,1+n) *= factor;
                    }
                });
            }
        }
    }
    // HI X
    if (bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) {

        // always multiply normal velocity by sqrt(2)
        // 1 = slip wall:    multiply transverse velocity by 0
        // 2 = no-slip wall: multiply transverse velocity by sqrt(2)
        Real factor = (bc_vel_hi[0] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_in[0].ixType());

        // this is the x-hi domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xhi = amrex::bdryNode(dom_x, Orientation(0, Orientation::high));

        for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xhi;
            Array4<Real> const& flux = (faceflux_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, AMREX_SPACEDIM, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (n == 0) {
                        // normal velocity
                        flux(i,j,k,1+n) *= sqrtTwo;
                        // viscous heating (diagonal & shear)
                        flux(i,j,k,nvars+1) *= sqrtTwo;
                        flux(i,j,k,nvars+2) *= factor;
                    } else {
                        // transverse velocity
                        flux(i,j,k,1+n) *= factor;
                    }
                });
            }
        }
    }
    // LO Y
    if (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) {

        // always multiply normal velocity by sqrt(2)
        // 1 = slip wall:    multiply transverse velocity by 0
        // 2 = no-slip wall: multiply transverse velocity by sqrt(2)
        Real factor = (bc_vel_lo[1] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_in[1].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_ylo = amrex::bdryNode(dom_y, Orientation(1, Orientation::low));

        for (MFIter mfi(faceflux_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_ylo;
            Array4<Real> const& flux = (faceflux_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, AMREX_SPACEDIM, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (n == 1) {
                        // normal velocity
                        flux(i,j,k,1+n) *= sqrtTwo;
                        // viscous heating (diagonal & shear)
                        flux(i,j,k,nvars+1) *= sqrtTwo;
                        flux(i,j,k,nvars+2) *= factor;
                    } else {
                        // transverse velocity
                        flux(i,j,k,1+n) *= factor;
                    }
                });
            }
        }
    }
    // HI Y
    if (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) {

        // always multiply normal velocity by sqrt(2)
        // 1 = slip wall:    multiply transverse velocity by 0
        // 2 = no-slip wall: multiply transverse velocity by sqrt(2)
        Real factor = (bc_vel_hi[1] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_in[1].ixType());

        // this is the y-hi domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yhi = amrex::bdryNode(dom_y, Orientation(1, Orientation::high));

        for (MFIter mfi(faceflux_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yhi;
            Array4<Real> const& flux = (faceflux_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, AMREX_SPACEDIM, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (n == 1) {
                        // normal velocity
                        flux(i,j,k,1+n) *= sqrtTwo;
                        // viscous heating (diagonal & shear)
                        flux(i,j,k,nvars+1) *= sqrtTwo;
                        flux(i,j,k,nvars+2) *= factor;
                    } else {
                        // transverse velocity
                        flux(i,j,k,1+n) *= factor;
                    }
                });
            }
        }
    }
    // LO Z
    if (bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) {

        // always multiply normal velocity by sqrt(2)
        // 1 = slip wall:    multiply transverse velocity by 0
        // 2 = no-slip wall: multiply transverse velocity by sqrt(2)
        Real factor = (bc_vel_lo[2] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_in[2].ixType());

        // this is the z-lo domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zlo = amrex::bdryNode(dom_z, Orientation(2, Orientation::low));

        for (MFIter mfi(faceflux_in[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zlo;
            Array4<Real> const& flux = (faceflux_in[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, AMREX_SPACEDIM, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (n == 2) {
                        // normal velocity
                        flux(i,j,k,1+n) *= sqrtTwo;
                        // viscous heating (diagonal & shear)
                        flux(i,j,k,nvars+1) *= sqrtTwo;
                        flux(i,j,k,nvars+2) *= factor;
                    } else {
                        // transverse velocity
                        flux(i,j,k,1+n) *= factor;
                    }
                });
            }
        }
    }
    // HI Z
    if (bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) {

        Real factor = (bc_vel_hi[2] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_in[2].ixType());

        // this is the z-hi domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zhi = amrex::bdryNode(dom_z, Orientation(2, Orientation::high));

        for (MFIter mfi(faceflux_in[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zhi;
            Array4<Real> const& flux = (faceflux_in[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, AMREX_SPACEDIM, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (n == 2) {
                        // normal velocity
                        flux(i,j,k,1+n) *= sqrtTwo;
                        // viscous heating (diagonal & shear)
                        flux(i,j,k,nvars+1) *= sqrtTwo;
                        flux(i,j,k,nvars+2) *= factor;
                    } else {
                        // transverse velocity
                        flux(i,j,k,1+n) *= factor;
                    }
                });
            }
        }
    }
}

void MembraneFlux(std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in,
                  const amrex::Geometry& /*geom*/) {

    BL_PROFILE_VAR("MembraneFlux()",MembraneFlux);

    // Loop over boxes
    for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {

        const Array4<Real> & xflux = (faceflux_in[0]).array(mfi);

        // since the MFIter is built on a nodal MultiFab we need to build the
        // nodal tileboxes for each direction in this way
        Box bx_x = mfi.tilebox(nodal_flag_x);

        if (bx_x.smallEnd(0) == membrane_cell || bx_x.bigEnd(0) == membrane_cell) {

            amrex::ParallelFor(bx_x, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                if (i == membrane_cell) {
                    xflux(i,j,k,1) = 0.;
                    xflux(i,j,k,2) = 0.;
                    xflux(i,j,k,3) = 0.;
                    xflux(i,j,k,4) = 0.;
                }

             });

        }
    }
}
