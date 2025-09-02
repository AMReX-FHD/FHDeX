#include <main_driver.H>
#include <main_driver_F.H>

#include <hydro_functions.H>
#include <StochMomFlux.H>
#include <common_functions.H>

#include <gmres_functions.H>

#include <ib_functions.H>

#include <immbdy_namespace.H>
#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>

#include <IBMarkerContainer.H>
//#include <IBMarkerMD.H>

#include "chrono"

#include <bits/stdc++.h>
using namespace std;

using namespace std::chrono;
using namespace amrex;
using namespace immbdy;
//using namespace immbdy_md;
using namespace ib_flagellum;


// argv contains the name of the inputs file entered at the command line
void main_driver(const char * argv) {

    BL_PROFILE_VAR("main_driver()", main_driver);

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();


    //___________________________________________________________________________
    // Load parameters from inputs file, and initialize global parameters
    std::string inputs_file = argv;

    read_immbdy_namelist(inputs_file.c_str(), inputs_file.size() + 1);

    // copy contents of F90 modules to C++ namespaces NOTE: any changes to
    // global settings in fortran/c++ after this point need to be synchronized
    InitializeCommonNamespace();
    InitializeImmbdyNamespace();
    InitializeIBFlagellumNamespace();


    //___________________________________________________________________________
    // Set boundary conditions

    // is the problem periodic?
    Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
    for (int i=0; i<AMREX_SPACEDIM; ++i)
        if (bc_vel_lo[i] <= -1 && bc_vel_hi[i] <= -1)
            is_periodic[i] = 1;


    //___________________________________________________________________________
    // Make BoxArray, DistributionMapping, and Geometry

    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(             0,              0,              0));
        IntVect dom_hi(AMREX_D_DECL(n_cells[0] - 1, n_cells[1] - 1, n_cells[2] - 1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);

        // Break up boxarray "ba" into chunks no larger than "max_grid_size"
        // along a direction note we are converting "Vector<int> max_grid_size"
        // to an IntVect
        ba.maxSize(IntVect(max_grid_size));

        // This defines the physical box, [-1, 1] in each direction
        RealBox real_box({AMREX_D_DECL(prob_lo[0], prob_lo[1], prob_lo[2])},
                         {AMREX_D_DECL(prob_hi[0], prob_hi[1], prob_hi[2])});

        // This defines a Geometry object
        geom.define(domain, & real_box, CoordSys::cartesian, is_periodic.data());
    }

    // how boxes are distrubuted among MPI processes
    DistributionMapping dmap(ba);


    //___________________________________________________________________________
    // Cell size, and time step
    Real dt         = fixed_dt;
    Real dtinv      = 1.0 / dt;
    const Real * dx = geom.CellSize();


    //___________________________________________________________________________
    // Initialize immersed boundaries
    // Make sure that the nghost (last argument) is big enough!

    BL_PROFILE_VAR("main_create_markers", create_markers);

    // Find the optimal number of ghost cells for the IBMarkerContainer
    Real min_dx = dx[0];
    for (int d=1; d<AMREX_SPACEDIM; ++d)
        min_dx = amrex::min(min_dx, dx[d]);

    int ib_nghost = 1;
    Print() << "Initializing IBMarkerContainer with "
            << ib_nghost << " ghost cells" << std::endl;
    // Initialize immersed boundary container
    IBMarkerContainer ib_mc(geom, dmap, ba, ib_nghost);

    for (int i_ib=0; i_ib < n_immbdy; ++i_ib) {

        if (n_marker[i_ib] <= 0) continue;

        int N  = n_marker[i_ib];
        Real L = ib_flagellum::length[i_ib];

        Real l_link = L/(N-1);

        const RealVect & x_0 = offset_0[i_ib];

        Print() << "Initializing flagellum:" << std::endl;
        Print() << "N=      " << N           << std::endl;
        Print() << "L=      " << L           << std::endl;
        Print() << "l_link= " << l_link      << std::endl;
        Print() << "x_0=    " << x_0         << std::endl;

        // using fourier modes => first two nodes reserved as "anchor"
        int N_markers = immbdy::contains_fourier ? N+1 : N;

        Vector<RealVect> marker_positions(N_markers);
        for (int i=0; i<marker_positions.size(); ++i) {
            Real x = x_0[0] + i*l_link;
            // Compute periodic offset. Will work as long as winding number = 1
            Real x_period = x < geom.ProbHi(0) ? x : x - geom.ProbLength(0);

            marker_positions[i] = RealVect{x_period, x_0[1], x_0[2]};
        }

        Vector<Real> marker_radii(N_markers);
        for (int i=0; i<marker_radii.size(); ++i) marker_radii[i] = 4*l_link;

        ib_mc.InitList(0, marker_radii, marker_positions, i_ib);
    }

    ib_mc.UpdatePIDMap();
    ib_mc.fillNeighbors();
    ib_mc.PrintMarkerData(0);
    BL_PROFILE_VAR_STOP(create_markers);

    int size = ParallelDescriptor::NProcs();
    auto & num_ids = ib_mc.getNumIDs();
    auto & cpu_offset = ib_mc.getCPUOffset();

    Print() << "num_ids = ";
    for (auto & i:num_ids) Print() << i << " ";
    Print() << std::endl;

    Print() << "cpu_offset = ";
    for (auto & i:cpu_offset) Print() << i << " ";
    Print() << std::endl;

    Vector<Real> pos_x(ib_mc.getTotalNumIDs());
    ib_mc.PullDown(0, pos_x, -1);
    Vector<Real> pos_y(ib_mc.getTotalNumIDs());
    ib_mc.PullDown(0, pos_y, -2);
    Vector<Real> pos_z(ib_mc.getTotalNumIDs());
    ib_mc.PullDown(0, pos_z, -3);

    Print() << "pos_x = ";
    for (auto & i:pos_x) Print() << i << " ";
    Print() << std::endl;

    Print() << "pos_y = ";
    for (auto & i:pos_y) Print() << i << " ";
    Print() << std::endl;

    Print() << "pos_z = ";
    for (auto & i:pos_z) Print() << i << " ";
    Print() << std::endl;

    Vector<int> ids(ib_mc.getTotalNumIDs());
    ib_mc.PullDownInt(0, ids, -1);
    Vector<int> cpus(ib_mc.getTotalNumIDs());
    ib_mc.PullDownInt(0, cpus, -2);

    Print() << "ids = ";
    for (auto & i:ids) Print() << i << " ";
    Print() << std::endl;

    Print() << "cpus = ";
    for (auto & i:cpus) Print() << i << " ";
    Print() << std::endl;


    //id_1 records actual ids on a given ib flagellum
    Vector<int> id1s(ib_mc.getTotalNumIDs());
    ib_mc.PullDownInt(0, id1s, IBMInt::id_1);
    //cpu_1 records the order of ib flagellum when generated
    Vector<int> ibs(ib_mc.getTotalNumIDs());
    ib_mc.PullDownInt(0, ibs, IBMInt::cpu_1);

    Print() << "real_id = ";
    for (auto & i:id1s) Print() << i << " ";
    Print() << std::endl;

    Print() << "number of flagellum = ";
    for (auto & i:ibs) Print() << i << " ";
    Print() << std::endl;

     // Get sorted ibs list
    vector<pair<int, int>> sorted_ibs = ib_mc.get_sorted_map();

    //       vector<pair<int, int>> sorted_ibs;

    //        for (int i = 0; i < ibs.size(); ++i) {
    //            sorted_ibs.push_back(make_pair(ibs[i], i));
    //        }

    //        sort(sorted_ibs.begin(), sorted_ibs.end());

        Print() << "Flagellum number\t" << "index in PullDown Vector" << std::endl;
    for (int i = 0; i < ibs.size(); i++) {
            Print() << sorted_ibs[i].first << "\t" << sorted_ibs[i].second << std::endl;
    }


    //Vectors for storing forces in each direction
    Vector<Real> fx(ib_mc.getTotalNumIDs());
    Vector<Real> fy(ib_mc.getTotalNumIDs());
    Vector<Real> fz(ib_mc.getTotalNumIDs());

    int index_start = 0;

    for (int i_ib=0; i_ib < n_immbdy; ++i_ib) {

        if (n_marker[i_ib] <= 0) continue;

        int N       = ib_flagellum::n_marker[i_ib];
        Real L      = ib_flagellum::length[i_ib];
        Real l_link = L/(N-1);

        Real k_spr  = ib_flagellum::k_spring[i_ib];
        Real k_driv = ib_flagellum::k_driving[i_ib];

        for (int ind=index_start; ind < index_start+N; ++ind ) {    //going through the sorted ibs index

            //Getting index for the current marker in the PullDown Vectors
            int i_c = sorted_ibs[ind].second;

            if(sorted_ibs[ind].first != i_ib) Abort("Mismatched flagella detected in predictor!!! flee for your lunch!");

            if(ind>index_start and ind<index_start+N-1) {   //exclude the first and last marker on the flagellum that doesn't have either next or previous marker
                //Getting indexes for the previous/minus and next/plus markers in the PullDown Vectors

                int i_p = sorted_ibs[ind+1].second;
                int i_m = sorted_ibs[ind-1].second;

                RealVect      pos = {pos_x[i_c],   pos_y[i_c],   pos_z[i_c]};
                RealVect next_pos = {pos_x[i_p],   pos_y[i_p],   pos_z[i_p]};
                RealVect prev_pos = {pos_x[i_m],   pos_y[i_m],   pos_z[i_m]};

                RealVect r_p = next_pos - pos, r_m = pos - prev_pos;

                Real lp_m = r_m.vectorLength(),         lp_p = r_p.vectorLength();
                Real fm_0 = k_spr * (lp_m-l_link)/lp_m, fp_0 = k_spr * (lp_p-l_link)/lp_p;

                //update spring forces between previous/minus and current markers
                fx[i_m] += fm_0 * r_m[1];   fy[i_m] += fm_0 * r_m[2];   fz[i_m] += fm_0 * r_m[3];
                fx[i_c] -= fm_0 * r_m[1];   fy[i_c] -= fm_0 * r_m[2];   fz[i_c] -= fm_0 * r_m[3];

                //update spring forces between next/plus and current markers
                fx[i_c] += fp_0 * r_p[1];   fy[i_c] += fp_0 * r_p[2];   fz[i_c] += fp_0 * r_p[3];
                fx[i_p] -= fp_0 * r_p[1];   fy[i_p] -= fp_0 * r_p[2];   fz[i_p] -= fp_0 * r_p[3];
            }

            // update bending forces for curent, minus/prev, and next/plus
            if(immbdy::contains_fourier) {  //for periodic waveform of flagellum in Fourier series
                Vector<RealVect> marker_positions(N); // =  equil_pos(i_ib, 0, geom);
                //int marker_seq_id                 = mark.idata(IBMInt::id_1);
                RealVect target_pos               = marker_positions[ids[i_c]];

                fx[i_c] += k_driv*(target_pos[1] - pos_x[i_c]);
                fy[i_c] += k_driv*(target_pos[2] - pos_y[i_c]);
                fz[i_c] += k_driv*(target_pos[3] - pos_z[i_c]);

            } else {                        // for harmonic waveform of Taylor lines

                if(ind>index_start and ind<index_start+N-1) {   //exclude the first and last markers on the flagellum that doesn't have either next or previous marker

                    int i_p = sorted_ibs[ind+1].second;
                    int i_m = sorted_ibs[ind-1].second;

                    RealVect      pos = {pos_x[i_c], pos_y[i_c], pos_z[i_c]};
                    RealVect next_pos = {pos_x[i_p], pos_y[i_p], pos_z[i_p]};
                    RealVect prev_pos = {pos_x[i_m], pos_y[i_m], pos_z[i_m]};

                    const RealVect & r = pos, & r_m = prev_pos, & r_p = next_pos;

                    // Set bending forces to zero before passing to 'driving_f'
                    RealVect f_p = RealVect{AMREX_D_DECL(0., 0., 0.)};
                    RealVect f   = RealVect{AMREX_D_DECL(0., 0., 0.)};
                    RealVect f_m = RealVect{AMREX_D_DECL(0., 0., 0.)};

                    // calling the active bending force calculation

//                    Real th = theta(
//                                    driv_amp, time, i_ib, ids[i_c] - 1
//                                   );
//                    driving_f(f, f_p, f_m, r, r_p, r_m, driv_u, th, k_driv);

                    // updating the force on the minus, current, and plus particles.
                    fx[i_m] += f_m[1];   fy[i_m] += f_m[2];   fz[i_m] += f_m[3];
                    fx[i_c] += f[1];     fy[i_c] += f[2];     fz[i_c] += f[3];
                    fx[i_p] += f_p[1];   fy[i_p] += f_p[2];   fz[i_p] += f_p[3];
                }
            }
        }
        index_start += N;
    }

    Print() << "force in x = ";
    for (auto & i:fx) Print() << i << " ";
    Print() << std::endl;

    Print() << "force in y = ";
    for (auto & i:fy) Print() << i << " ";
    Print() << std::endl;

    Print() << "force in z = ";
    for (auto & i:fz) Print() << i << " ";
    Print() << std::endl;


//   for (int i=0; i<ib_mc.getTotalNumIDs(); i++) {
//        std::vector<int>::iterator gx = std::find(sorted_ibs.begin(), sorted_ibs.end(), i_ib);

   //actual index in PullDown vectors
//    int i_c = sorted.ibs[std::distance(sorted_ibs.begin(), gx) + id].second;

//    Print() << "Adding forces to particles..." << std::endl;




    // Just for fun, print out the max runtime
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    Print() << "Run time = " << stop_time << std::endl;
}