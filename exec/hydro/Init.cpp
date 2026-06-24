#include "hydro_test_functions.H"

using namespace amrex;

void InitVel(std::array< MultiFab, AMREX_SPACEDIM >& umac,
             const Geometry& geom) {

    // restart from a possibly-coarser plotfile
    if (prob_type == -1) {

        MultiFab plot_init;

        plot_init_file += "/Level_0/Cell";

        // read in the (coarse) plotfile
        VisMF::Read(plot_init, plot_init_file);

        // read in BoxArray from (coarse) plotfile
        BoxArray ba_init_file = plot_init.boxArray();

        // create a single Box that spans the domain of the (coarse) plotfile
        Box bx_init_file = ba_init_file.minimalBox();

        // create BoxArray with this one single box and a corresponding DistrubtionMap
        BoxArray ba_onegrid(bx_init_file);
        DistributionMapping dm_onegrid(ba_onegrid);

        // create MultiFabs with one (coarse) box on one MPI rank
        // cell-centered
        MultiFab mf_onegrid(ba_onegrid, dm_onegrid, 1, 0);
        // face-centered
        std::array< MultiFab, AMREX_SPACEDIM > umac_onegrid;
        AMREX_D_TERM(umac_onegrid[0].define(convert(ba_onegrid,nodal_flag_x), dm_onegrid, 1, 0);,
                     umac_onegrid[1].define(convert(ba_onegrid,nodal_flag_y), dm_onegrid, 1, 0);,
                     umac_onegrid[2].define(convert(ba_onegrid,nodal_flag_z), dm_onegrid, 1, 0););

        // check to make sure refinement ratio is valid in the x-direction
        int hi_x = bx_init_file.bigEnd(0) + 1;
        int rr = n_cells[0] / hi_x;
        if (n_cells[0] % hi_x != 0) {
            Abort("Init.cpp prob_type 1; plotfile is not evenly divisible by n_cells in the x direction");
        }

        // check remaining dimensions
        for (int i=1; i<AMREX_SPACEDIM; ++i) {
            int hi = bx_init_file.bigEnd(i);
            if (n_cells[i] / hi != rr) {
                Abort("Init.cpp prob_type 1; refinement ratio not the same in all directions");
            }
            else if (n_cells[i] & hi != 0) {
                Abort("Init.cpp prob_type 1; plotfile is not evenly divisible by n_cells in y or z direction");
            }
        }

        Print() << "Initializing from a plotfile coarser by a factor of " << rr << std::endl;

        for (int i=0; i<AMREX_SPACEDIM; ++i) {

            // parallel copy (coarse) plotfile data to (coarse) MF with one grid
            mf_onegrid.ParallelCopy(plot_init,i+3,0,1);

            // shift velocities onto faces
            ShiftCCToFace_onegrid(umac_onegrid[i],0,mf_onegrid,0,1);

            // obtain BoxArray and DistributionMap from the final output MF
            BoxArray ba_final = umac[i].boxArray();
            DistributionMapping dm_final = umac[i].DistributionMap();

            // make a coarsened version of the BoxArray from the final output MF
            BoxArray ba_final_coarse = ba_final.coarsen(rr);

            // error checking to make sure coarsened domains match
            Box bx_final_coarse1 = ba_final_coarse.minimalBox();
            Box bx_final_coarse2 = umac_onegrid[i].boxArray()[0];
            if (bx_final_coarse1 != bx_final_coarse2) {
                Abort("Init: prob_type=-1; coarsened box does not equal read-in box.  Check refinement ratio?");
            }

            // allocate a face-centered MultiFab with the
            // coarsened final BoxArray and DistributionMap
            MultiFab umac_coarse(ba_final_coarse,dm_final,1,0);

            // parallel copy data from one (coarse) grid into MFs with same distribution
            // map as the final output MF but on the coarsened BoxArray
            umac_coarse.ParallelCopy(umac_onegrid[i],0,0,1);

            // inject the (coarse) distributed data into the final output MF
            for (MFIter mfi(umac[i],TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                Array4<Real      > const& fine = umac[i].array(mfi);
                Array4<Real const> const& crse = umac_coarse.array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    fine(i,j,k) = crse(i/rr,j/rr,k/rr);
                });
            }
        }

        return;
    }


    Real zshft = (AMREX_SPACEDIM == 2) ? 0. : 0.5;

    GpuArray<Real,AMREX_SPACEDIM> reallo = geom.ProbLoArray();
    GpuArray<Real,AMREX_SPACEDIM> realhi = geom.ProbHiArray();
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    GpuArray<Real,AMREX_SPACEDIM> center;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        center[d] = 0.5*(realhi[d]-reallo[d]);
    }

    // IC parameters
    Real L_hlf = center[0];
    // k1 & k2 determine steepness of velocity profile:
    Real k1 = 0.1*L_hlf;
    // k1 = 1d-6*L_hlf
    Real k2 = k1;
    Real k1_inv = 1./k1;
    Real k2_inv = 1./k2;

    // Vortex:
    // [r_a r_b] defines radial bounds of velocity bump:
    Real r_a = 0.35*L_hlf;
    Real r_b = L_hlf - r_a;

    // Kelvin-Helmholtz:
    Real pi = 3.1415926535897932;
    Real freq = 02.*pi/L_hlf;
    Real amp = 2.0e-3*L_hlf;
    // amp = 2.0d-1*L_hlf
    Real width1 = L_hlf/2.;

    if (prob_type == 0) {
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            umac[d].setVal(0.);
        }
        return;
    }

    for ( MFIter mfi(umac[0]); mfi.isValid(); ++mfi ) {

        AMREX_D_TERM(const Array4<Real> & u = (umac[0]).array(mfi);,
                     const Array4<Real> & v = (umac[1]).array(mfi);,
                     const Array4<Real> & w = (umac[2]).array(mfi););

        // since the MFIter is built on a nodal MultiFab we need to build the
        // nodal tileboxes for each direction in this way
        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x);,
                     Box bx_y = mfi.tilebox(nodal_flag_y);,
                     Box bx_z = mfi.tilebox(nodal_flag_z););


        amrex::ParallelFor(bx_x,
                           bx_y,
#if (AMREX_SPACEDIM == 3)
                           bx_z,
#endif
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            GpuArray<Real,AMREX_SPACEDIM> itVec;

            AMREX_D_TERM(itVec[0] = i*dx[0];,
                         itVec[1] = (j+0.5)*dx[1];,
                         itVec[2] = (k+zshft)*dx[2];);

            GpuArray<Real,AMREX_SPACEDIM> relpos;
            Real rad2 = 0.;
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                relpos[d] = reallo[d] + itVec[d] - center[d];
            }

            // only sum first two velocities
            for (int d=0; d<2; ++d) {
                rad2 += relpos[d]*relpos[d];
            }

            Real rad = std::sqrt(rad2);

            // note prob_type == 0 handled above
            if (prob_type == 1) {

                // Multiply velocity magnitude by sin(theta)
                u(i,j,k) = 0.25*(1.+std::tanh(k1_inv*(rad-r_a)))*(1.+std::tanh(k2_inv*(r_b-rad)))
                    *(relpos[1]/rad);

            } else if (prob_type == 2) {

                Real perturb = amp*sin(freq*relpos[0]);
                u(i,j,k) = 0.25*(1.+std::tanh(k1_inv*(relpos[1] - (-width1/2.+perturb))))
                    *(1.+std::tanh(k2_inv*((width1/2.+perturb) - relpos[1])));

            } else if (prob_type == 3) {

                u(i,j,k) = 0.25*(1.+std::tanh(k1_inv*(relpos[1] - (-width1/2.))))
                    *(1.+std::tanh(k2_inv*((width1/2.) - relpos[1])));

            }
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real,AMREX_SPACEDIM> itVec;

            AMREX_D_TERM(itVec[0] = (i+0.5)*dx[0];,
                         itVec[1] = j*dx[1];,
                         itVec[2] = (k+zshft)*dx[2];);

            GpuArray<Real,AMREX_SPACEDIM> relpos;
            Real rad2 = 0.;
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                relpos[d] = reallo[d] + itVec[d] - center[d];
            }

            // only sum first two velocities
            for (int d=0; d<2; ++d) {
                rad2 += relpos[d]*relpos[d];
            }

            Real rad = std::sqrt(rad2);

            // note prob_type == 0 handled above
            if (prob_type == 1) {

                // Multiply velocity magnitude by -cos(theta)
                v(i,j,k) = 0.25*(1.+std::tanh(k1_inv*(rad-r_a)))*(1.+std::tanh(k2_inv*(r_b-rad)))
                    *(-relpos[0]/rad);

            } else if (prob_type == 2) {

                Real perturb = amp*sin(freq*relpos[0]);
                Real slope = amp*freq*cos(freq*relpos[0]);
                Real fun_ptrb = 0.25*(1.*tanh(k1_inv*(relpos[1] - (-width1/2.+perturb))))
                    *(1.+tanh(k2_inv*((width1/2.+perturb) - relpos[1])));
                v(i,j,k) = slope*fun_ptrb;
            } else if (prob_type == 3) {
                v(i,j,k) = 0.;
            }
        }
#if (AMREX_SPACEDIM == 3)
                         , [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            w(i,j,k) = 0.;
        }
#endif
        );
    }
}
