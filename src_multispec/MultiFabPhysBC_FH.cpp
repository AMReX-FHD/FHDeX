#include "multispec_functions.H"
#include "InhomogeneousBCVal.H"

// Ghost cell filling routine.
// Fills in ALL ghost cells to the value ON the boundary.
// FOEXTRAP uses boundary conditions (Neumann) and 1 interior points.
// EXT_DIR copies the supplied Dirichlet condition into the ghost cells.
void MultiFabPhysBCFH(MultiFab& phi, const Geometry& geom, int scomp, int ncomp, const Real& scale, const Real& time) {

    BL_PROFILE_VAR("MultiFabPhysBCFH()",MultiFabPhysBCFH);
    
    // bccomp definitions are in BCPhysToMath.cpp

//    amrex::Print() << " entering special FH bc routine " << std::endl;
    
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
    int bccomp = SPEC_BC_COMP;

    amrex::Real kappa_coeff = fh_kappa(0,1);
    amrex::Real ce = fh_ce;
    amrex::Real omce = 1.-ce;
    amrex::Real hack = 1./(1.-2.*ce);

//    amrex::Print() << "scale and coeff " << scale << " " << kappa_coeff << " " << scale*kappa_coeff << " " << fh_tension <<std::endl;

    //Real coeff = 6.*fh_tension/(kappa_coeff*scale);
    Real coeff = 3.*fh_tension*hack/(kappa_coeff*scale);

//    amrex::Print() << " coeff " << coeff << std::endl;

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
            if (bc_mass_lo[0] == 1) {
//                amrex::Print() << " entering special FH bc routine lo " <<  bc_mass_lo[0] << std::endl;
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (i < lo) {
                        Real y = prob_lo[1] + (j+0.5)*dx[1];
                        Real z = prob_lo[2] + (k+0.5)*dx[2];
                        data(i,j,k,scomp+n) = data(lo,j,k,scomp+n) - dx[0]*InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
            else if (bc_mass_lo[0] == 4) {
//                amrex::Print() << " entering special FH bc routine lo " <<  bc_mass_lo[0] << std::endl;
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < lo) {
                        amrex::Real scratch = (data(lo,j,k,scomp+0)-ce)/(omce-ce);
                        data(i,j,k,scomp+0) = data(lo,j,k,scomp+0) + dx[0]*std::cos(contact_angle_lo[0])*scratch*(1.-scratch)*coeff;
                        //data(i,j,k,scomp+0) = data(lo,j,k,scomp+0) + dx[0]*std::cos(contact_angle_lo[0])*data(lo,j,k,scomp+0)*data(lo,j,k,scomp+1)*coeff;
                        data(i,j,k,scomp+0) = amrex::min(1.,amrex::max(0.,data(i,j,k,scomp+0)));
                        data(i,j,k,scomp+1) = 1.-data(i,j,k,scomp+0);
//                        amrex::Print() << " Fh data left " << j << " " << data(i,j,k,scomp) <<" " << data(i,j,k,scomp+1) <<  std::endl;

//                        if(j == 18){
//                            Real foo = std::cos(contact_angle_lo[0])* coeff;
//                            amrex::Print() << "coeff and cos " << coeff << " " << contact_angle_lo[0] << " " << foo << std::endl;
//                        }

                        //data(i,j,k,scomp+1) = data(lo,j,k,scomp+1) + 0.5*dx[0]*data(lo,j,k,scomp+0)*data(lo,j,k,scomp+1);
                    }
                });

            }
            else if ( bc_mass_lo[0] != -1) {
               amrex::Print{} << "Should not be here lo x " << bc_mass_lo[0]  << std::endl;
               Abort("Wrong parameter in FH boundary condition");
            }
        }
        
        if (bx.bigEnd(0) > hi) {
            Real x = prob_hi[0];
            if (bc_mass_hi[0] == 1) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (i > hi) {
                        Real y = prob_lo[1] + (j+0.5)*dx[1];
                        Real z = prob_lo[2] + (k+0.5)*dx[2];
                        data(i,j,k,scomp+n) = data(hi,j,k,scomp+n) - dx[0]*InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
            else if (bc_mass_hi[0] == 4) {
//                amrex::Print() << " entering special FH bc routine hi " <<  bc_mass_hi[0] << std::endl;
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > hi) {
                        amrex::Real scratch = (data(hi,j,k,scomp+0)-ce)/(omce-ce);
                        data(i,j,k,scomp+0) = data(hi,j,k,scomp+0) + dx[0]*std::cos(contact_angle_hi[0])*scratch*(1.-scratch)*coeff;
                        //data(i,j,k,scomp+0) = data(hi,j,k,scomp+0) + dx[0]*std::cos(contact_angle_hi[0])*data(hi,j,k,scomp+0)*data(hi,j,k,scomp+1)*coeff;
                        data(i,j,k,scomp+0) = amrex::min(1.,amrex::max(0.,data(i,j,k,scomp+0)));
                        data(i,j,k,scomp+1) = 1.-data(i,j,k,scomp+0);
//                        amrex::Print() << " Fh data right " << j << " " << data(i,j,k,scomp) << " " << data(i,j,k,scomp+1) << std::endl;

                        //data(i,j,k,scomp+1) = data(lo,j,k,scomp+1) + 0.5*dx[0]*data(lo,j,k,scomp+0)*data(lo,j,k,scomp+1);
                    }
                });

            }
            else if ( bc_mass_hi[0] != -1) {
               amrex::Print{} << "Should not be here hi x " << bc_mass_hi[0]  << std::endl;
               Abort("Wrong parameter in FH boundary condition");
            }
        }

#if (AMREX_SPACEDIM >= 2)
        //___________________________________________________________________________
        // Apply y-physbc to data

        lo = dom.smallEnd(1);
        hi = dom.bigEnd(1);
        
        if (bx.smallEnd(1) < lo) {
            Real y = prob_lo[1];
            if (bc_mass_lo[1] == 1) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (j < lo) {
                        Real x = prob_lo[0] + (i+0.5)*dx[0];
                        Real z = prob_lo[2] + (k+0.5)*dx[2];
                        data(i,j,k,scomp+n) = data(i,lo,k,scomp+n) - dx[1]*InhomogeneousBCVal(bccomp+n,x,y,z,time);;
                    }
                });
            }
            else if (bc_mass_lo[1] == 4) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < lo) {
                        //amrex::Real scratch = (data(i,lo,k,scomp+0)-.0348115)/(.9651885-.0348115);
                        //data(i,j,k,scomp+0) = data(i,lo,k,scomp+0) + dx[1]*std::cos(contact_angle_lo[1])*data(i,lo,k,scomp+0)*data(i,lo,k,scomp+1)*coeff;
                        amrex::Real scratch = (data(i,lo,k,scomp+0)-ce)/(omce-ce);
                        data(i,j,k,scomp+0) = data(i,lo,k,scomp+0) + dx[1]*std::cos(contact_angle_lo[1])*scratch*(1.-scratch)*coeff;
                        data(i,j,k,scomp+0) = amrex::min(1.,amrex::max(0.,data(i,j,k,scomp+0)));
                        data(i,j,k,scomp+1) = 1.-data(i,j,k,scomp+0);

//                        if(j == lo-1 && i > 60 && i < 70){
//                           amrex::Print() << " normal derivative " << cos(contact_angle_lo[1])*data(i,lo,k,scomp+0)*data(i,lo,k,scomp+1)*coeff << std::endl;
//                        }
                        //data(i,j,k,scomp+1) = data(lo,j,k,scomp+1) + 0.5*dx[0]*data(lo,j,k,scomp+0)*data(lo,j,k,scomp+1);
                    }
                });

            }
            else if ( bc_mass_lo[1] != -1) {
               amrex::Print{} << "Should not be here lo y " << bc_mass_lo[1]  << std::endl;
               Abort("Wrong parameter in FH boundary condition");
            }
        }

        if (bx.bigEnd(1) > hi) {
            Real y = prob_hi[1];
            if (bc_mass_hi[1] == 1) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (j > hi) {
                        Real x = prob_lo[0] + (i+0.5)*dx[0];
                        Real z = prob_lo[2] + (k+0.5)*dx[2];
                        data(i,j,k,scomp+n) = data(i,hi,k,scomp+n) - dx[1]*InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
            else if (bc_mass_hi[1] == 4) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > hi) {
                        amrex::Real scratch = (data(i,hi,k,scomp+0)-ce)/(omce-ce);
                        data(i,j,k,scomp+0) = data(i,hi,k,scomp+0) + dx[1]*std::cos(contact_angle_hi[1])*scratch*(1.-scratch)*coeff;
                        //data(i,j,k,scomp+0) = data(i,hi,k,scomp+0) + dx[1]*std::cos(contact_angle_hi[1])*data(i,hi,k,scomp+0)*data(i,hi,k,scomp+1)*coeff;
                        data(i,j,k,scomp+0) = amrex::min(1.,amrex::max(0.,data(i,j,k,scomp+0)));
                        data(i,j,k,scomp+1) = 1.-data(i,j,k,scomp+0);

                        //data(i,j,k,scomp+1) = data(lo,j,k,scomp+1) + 0.5*dx[0]*data(lo,j,k,scomp+0)*data(lo,j,k,scomp+1);
                    }
                });

            }
            else if ( bc_mass_hi[1] != -1) {
               amrex::Print{} << "Should not be here hi y " << bc_mass_hi[1]  << std::endl;
               Abort("Wrong parameter in FH boundary condition");
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
            if (bc_mass_lo[2] == 1) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (k < lo) {
                        Real x = prob_lo[0] + (i+0.5)*dx[0];
                        Real y = prob_lo[1] + (j+0.5)*dx[1];
                        data(i,j,k,scomp+n) = data(i,j,lo,scomp+n) - dx[2]*InhomogeneousBCVal(bccomp+n,x,y,z,time);;
                    }
                });
            }
            else if (bc_mass_lo[2] == 4) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < lo) {
                        amrex::Real scratch = (data(i,j,lo,scomp+0)-ce)/(omce-ce);
                        data(i,j,k,scomp+0) = data(i,j,lo,scomp+0) + dx[2]*std::cos(contact_angle_lo[2])*scratch*(1.-scratch)*coeff;
                        //data(i,j,k,scomp+0) = data(i,j,lo,scomp+0) + dx[2]*std::cos(contact_angle_lo[2])*data(i,j,lo,scomp+0)*data(i,j,lo,scomp+1)*coeff;
                        data(i,j,k,scomp+0) = amrex::min(1.,amrex::max(0.,data(i,j,k,scomp+0)));
                        data(i,j,k,scomp+1) = 1.-data(i,j,k,scomp+0);

                        //data(i,j,k,scomp+1) = data(lo,j,k,scomp+1) + 0.5*dx[0]*data(lo,j,k,scomp+0)*data(lo,j,k,scomp+1);
                    }
                });

            }
            else if ( bc_mass_lo[2] != -1) {
               amrex::Print{} << "Should not be here " << std::endl;
               Abort("Wrong parameter in FH boundary condition");
            }
        }

        if (bx.bigEnd(2) > hi) {
            Real z= prob_hi[2];
            if (bc_mass_hi[2] == 1) {
                amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (k > hi) {
                        Real x = prob_lo[0] + (i+0.5)*dx[0];
                        Real y = prob_lo[1] + (j+0.5)*dx[1];
                        data(i,j,k,scomp+n) = data(i,j,hi,scomp+n) - dx[2]*InhomogeneousBCVal(bccomp+n,x,y,z,time);
                    }
                });
            }
            else if (bc_mass_hi[2] == 4) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < hi) {
                        amrex::Real scratch = (data(i,j,hi,scomp+0)-ce)/(omce-ce);
                        data(i,j,k,scomp+0) = data(i,j,hi,scomp+0) + dx[2]*std::cos(contact_angle_hi[2])*scratch*(1.-scratch)*coeff;
                        //data(i,j,k,scomp+0) = data(i,j,hi,scomp+0) + dx[2]*std::cos(contact_angle_hi[2])*data(i,j,hi,scomp+0)*data(i,j,hi,scomp+1)*coeff;
                        data(i,j,k,scomp+0) = amrex::min(1.,amrex::max(0.,data(i,j,k,scomp+0)));
                        data(i,j,k,scomp+1) = 1.-data(i,j,k,scomp+0);

                        //data(i,j,k,scomp+1) = data(lo,j,k,scomp+1) + 0.5*dx[0]*data(lo,j,k,scomp+0)*data(lo,j,k,scomp+1);
                    }
                });

            }
            else if ( bc_mass_hi[2] != -1) {
               amrex::Print{} << "Should not be here " << std::endl;
               Abort("Wrong parameter in FH boundary condition");
            }
        }
#endif
        
    } // end MFIter
}

