#include <main_driver.H>

#include <AMReX_ParallelDescriptor.H>

#include <immbdy_namespace.H>

#include <IBMarkerMD.H>

// #include <bits/stdc++.h>

using namespace amrex;

using namespace immbdy;
using namespace immbdy_md;
using namespace ib_flagellum;

void constrain_ibm_marker(IBMarkerContainer & ib_mc, int ib_lev, int component) {

    BL_PROFILE_VAR("constrain_ibm_marker", TIMER_CONSTRAIN_IBM_MARKER);

    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        TileIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        long np = ib_mc.GetParticles(ib_lev).at(index).numParticles();

        for (int i = 0; i < np; ++i) {

            ParticleType & mark = markers[i];
            // Zero component only
            mark.rdata(component) = 0.;
        }
    }

    BL_PROFILE_VAR_STOP(TIMER_CONSTRAIN_IBM_MARKER);
}



void anchor_first_marker(IBMarkerContainer & ib_mc, int ib_lev, int component) {

    BL_PROFILE_VAR("anchor_first_marker", TIMER_ANCHOR_FIRST_MARKER);

    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        TileIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();

        long np = ib_mc.GetParticles(ib_lev).at(index).numParticles();

        for (int i = 0; i < np; ++i) {

            ParticleType & mark = markers[i];

            if((mark.idata(IBMInt::id_1) == 0)||(mark.idata(IBMInt::id_1) == 1))
                for (int d=0; d<AMREX_SPACEDIM; ++d) mark.rdata(component + d) = 0.;
        }
    }

    BL_PROFILE_VAR_STOP(TIMER_ANCHOR_FIRST_MARKER);
}


Real theta(Real amp_ramp, Real time, int i_ib, int index_marker) {

    if(immbdy::contains_fourier) {

        // First two nodes reserved as "anchor"
        // index_marker = amrex::max(0, index_marker-1);

        int N                 = chlamy_flagellum::N[i_ib][index_marker];
        int coef_len          = ib_flagellum::fourier_coef_len;
        Vector<Real> & a_coef = chlamy_flagellum::a_coef[i_ib][index_marker];
        Vector<Real> & b_coef = chlamy_flagellum::b_coef[i_ib][index_marker];

        Real frequency = ib_flagellum::frequency[i_ib];

        Real dt = 1./N;
        Real th = 0;

        Real k_fact = 2*M_PI/N;
        Real t_unit = frequency*time/dt;
        for (int i=0; i < coef_len; ++i) {
            Real k = k_fact * i;
            th += a_coef[i]*cos(k*t_unit);
            th -= b_coef[i]*sin(k*t_unit);
        }

        return amp_ramp*th/N;

    } else {

        int  N          = ib_flagellum::n_marker[i_ib];
        Real L          = ib_flagellum::length[i_ib];
        Real wavelength = ib_flagellum::wavelength[i_ib];
        Real frequency  = ib_flagellum::frequency[i_ib];
        Real amplitude  = ib_flagellum::amplitude[i_ib];
        Real l_link     = L/(N-1);

        Real theta = l_link*amp_ramp*amplitude*sin(
                2*M_PI*(frequency*time + index_marker*l_link/wavelength)
                + M_PI/2 + M_PI*(l_link/wavelength)
            );

        return theta;
    }
}



void update_ibm_marker(const RealVect & driv_u, Real driv_amp, Real time,
                       IBMarkerContainer & ib_mc, int ib_lev,
                       int component, bool pred_pos,
                       const Geometry & geom) {

    BL_PROFILE_VAR("update_ibm_marker", UpdateForces);

    // PullDown all particle positions onto one processor
    int size = ParallelDescriptor::NProcs();
    auto & num_ids = ib_mc.getNumIDs();
    auto & cpu_offset = ib_mc.getCPUOffset();

    Vector<Real> pos_x(ib_mc.getTotalNumIDs());
    ib_mc.PullDown(0, pos_x, -1);
    Vector<Real> pos_y(ib_mc.getTotalNumIDs());
    ib_mc.PullDown(0, pos_y, -2);
    Vector<Real> pos_z(ib_mc.getTotalNumIDs());
    ib_mc.PullDown(0, pos_z, -3);

    //id_1 records marker ids on a given ib
    Vector<int> ids(ib_mc.getTotalNumIDs());
    ib_mc.PullDownInt(0, ids, IBMInt::id_1);
    //cpu_1 records the id of each of ib
    Vector<int> ibs(ib_mc.getTotalNumIDs());
    ib_mc.PullDownInt(0, ibs, IBMInt::cpu_1);

    Print() << "maker_id = ";
    for (auto & i:ids) Print() << i << " ";
    Print() << std::endl;

    Print() << "ib_id = ";
    for (auto & i:ibs) Print() << i << " ";
    Print() << std::endl;

    //Vectors for storing forces in each direction
    Vector<Real> fx(ib_mc.getTotalNumIDs());
    for (auto & x:fx) x = 0.;
    Vector<Real> fy(ib_mc.getTotalNumIDs());
    for (auto & x:fy) x = 0.;
    Vector<Real> fz(ib_mc.getTotalNumIDs());
    for (auto & x:fz) x = 0.;

    // Get sorted ibs list
    std::vector<std::tuple<int, int, int>> sorted_ibs = ib_mc.get_sorted_map();
    std::vector<int> reduced_ibs = ib_mc.get_reduced_map();

    int index_start = 0;

    for (int i_ib=0; i_ib < n_immbdy; ++i_ib) {

        if (n_marker[i_ib] <= 0) continue;

        int N       = ib_flagellum::n_marker[i_ib];
        Real L      = ib_flagellum::length[i_ib];
        Real l_link = L/(N-1);

        Real k_spr  = ib_flagellum::k_spring[i_ib];
        Real k_driv = ib_flagellum::k_driving[i_ib];

        for (int ind = index_start; ind <(index_start + N - 1); ++ind){
            int i_0 = IBMarkerContainer::storage_idx(sorted_ibs[ind]);
            int i_1 = IBMarkerContainer::storage_idx(sorted_ibs[ind+1]);

            if ((pos_x[i_1] - pos_x[i_0]) < -geom.ProbLength(0)/2) {
                pos_x[i_1] = pos_x[i_1] + geom.ProbLength(0);
            }

            if ((pos_x[i_1] - pos_x[i_0]) > geom.ProbLength(0)/2) {
                pos_x[i_1] = pos_x[i_1] - geom.ProbLength(0);
            }

            if ((pos_y[i_1] - pos_y[i_0]) < -geom.ProbLength(1)/2) {
                pos_y[i_1] = pos_y[i_1] + geom.ProbLength(1);
            }

            if ((pos_y[i_1] - pos_y[i_0]) > geom.ProbLength(1)/2) {
                pos_y[i_1] = pos_y[i_1] - geom.ProbLength(1);
            }

            if ((pos_z[i_1] - pos_z[i_0]) < -geom.ProbLength(2)/2) {
                pos_z[i_1] = pos_z[i_1] + geom.ProbLength(2);
            }

            if ((pos_z[i_1] - pos_z[i_0]) > geom.ProbLength(2)/2) {
                pos_z[i_1] = pos_z[i_1] - geom.ProbLength(2);
            }
        }


        for (int ind=index_start; ind < index_start+N; ++ind ) {    //going through the sorted ibs index

                //Getting index for the current marker in the PullDown Vectors
                int i_c = IBMarkerContainer::storage_idx(sorted_ibs[ind]);

                if(IBMarkerContainer::immbdy_idx(sorted_ibs[ind]) != i_ib) {
                    Print() << IBMarkerContainer::immbdy_idx(sorted_ibs[ind])
                            << " i_ib = " << i_ib << std::endl;
                    Abort("Mismatched flagella detected in predictor!!!");
                }

                if(ind>index_start and ind<index_start+N-1) {
                    //exclude the first and last marker on the flagellum that
                    //doesn't have either next or previous marker Getting
                    //indexes for the previous/minus and next/plus markers in
                    //the PullDown Vectors

                    int i_p = IBMarkerContainer::storage_idx(sorted_ibs[ind+1]);
                    int i_m = IBMarkerContainer::storage_idx(sorted_ibs[ind-1]);

                    RealVect      pos = {pos_x[i_c],   pos_y[i_c],   pos_z[i_c]};
                    RealVect next_pos = {pos_x[i_p],   pos_y[i_p],   pos_z[i_p]};
                    RealVect prev_pos = {pos_x[i_m],   pos_y[i_m],   pos_z[i_m]};

                    RealVect r_p = next_pos - pos, r_m = pos - prev_pos;

                    Real lp_m = r_m.vectorLength(),         lp_p = r_p.vectorLength();
                    Real fm_0 = k_spr * (lp_m-l_link)/lp_m, fp_0 = k_spr * (lp_p-l_link)/lp_p;

                    Print() << "Updating spring forces..." << std::endl;

                    //update spring forces between previous/minus and current markers
                    fx[i_m] += fm_0 * r_m[0]; fy[i_m] += fm_0 * r_m[1]; fz[i_m] += fm_0 * r_m[2];
                    fx[i_c] -= fm_0 * r_m[0]; fy[i_c] -= fm_0 * r_m[1]; fz[i_c] -= fm_0 * r_m[2];

                    //update spring forces between next/plus and current markers
                    fx[i_c] += fp_0 * r_p[0]; fy[i_c] += fp_0 * r_p[1]; fz[i_c] += fp_0 * r_p[2];
                    fx[i_p] -= fp_0 * r_p[0]; fy[i_p] -= fp_0 * r_p[1]; fz[i_p] -= fp_0 * r_p[2];
                }

                // update bending forces for curent, minus/prev, and next/plus
                if(immbdy::contains_fourier) {
                    //for periodic waveform of flagellum in Fourier series
                    Vector<RealVect> marker_positions = equil_pos(i_ib, time, geom);
                    RealVect target_pos = marker_positions[ids[i_c]];

                    fx[i_c] += k_driv*(target_pos[0] - pos_x[i_c]);
                    fy[i_c] += k_driv*(target_pos[1] - pos_y[i_c]);
                    fz[i_c] += k_driv*(target_pos[2] - pos_z[i_c]);
                } else {
                    // for harmonic waveform of Taylor lines
                    if(ind > index_start && ind < index_start+N-1) {
                        //exclude the first and last markers on the flagellum
                        //that doesn't have either next or previous marker

                        Print() << "Updating bending forces..." << std::endl;

                        int i_p = IBMarkerContainer::storage_idx(sorted_ibs[ind+1]);
                        int i_m = IBMarkerContainer::storage_idx(sorted_ibs[ind-1]);

                        RealVect      pos = {pos_x[i_c], pos_y[i_c], pos_z[i_c]};
                        RealVect next_pos = {pos_x[i_p], pos_y[i_p], pos_z[i_p]};
                        RealVect prev_pos = {pos_x[i_m], pos_y[i_m], pos_z[i_m]};

                        const RealVect & r = pos, & r_m = prev_pos, & r_p = next_pos;

                        // Set bending forces to zero before passing to 'driving_f'
                        RealVect f_p = RealVect{AMREX_D_DECL(0., 0., 0.)};
                        RealVect f   = RealVect{AMREX_D_DECL(0., 0., 0.)};
                        RealVect f_m = RealVect{AMREX_D_DECL(0., 0., 0.)};

                        // calling the active bending force calculation

                        Real th = theta(driv_amp, time, i_ib, ids[i_c] - 1);
                        driving_f(f, f_p, f_m, r, r_p, r_m, driv_u, th, k_driv);

                        // updating the force on the minus, current, and plus particles.
                        fx[i_m] += f_m[0];   fy[i_m] += f_m[1];   fz[i_m] += f_m[2];
                        fx[i_c] += f[0];     fy[i_c] += f[1];     fz[i_c] += f[2];
                        fx[i_p] += f_p[0];   fy[i_p] += f_p[1];   fz[i_p] += f_p[2];
                    }
                }
        }
        index_start += N;
    }

    // Fianlly, Iterating through all markers and add forces
    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();
        long np = ib_mc.GetParticles(ib_lev).at(index).numParticles();

        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            ParticleType & mark = markers[m_index.first];

            int id      = mark.idata(IBMInt::id_1);   //id # on the flagellum
            int i_ib    = mark.idata(IBMInt::cpu_1);  // flagellum #
            int N       = ib_flagellum::n_marker[i_ib];

            //int i_c = id + std::distance(sorted_ibs.begin(), std::find_if(sorted_ibs.begin(), sorted_ibs.end(), [&](const auto& pair) { return pair.first == i_ib; }));
            int i_c = IBMarkerContainer::storage_idx(sorted_ibs[id + reduced_ibs[i_ib]]);
            Print() << "Adding forces to particles..." << std::endl;
            mark.rdata(IBMReal::forcex) += fx[i_c];
            mark.rdata(IBMReal::forcey) += fy[i_c];
            mark.rdata(IBMReal::forcez) += fz[i_c];
        }
    }
    BL_PROFILE_VAR_STOP(UpdateForces);
};


void yeax_ibm_marker(Real mot, IBMarkerContainer & ib_mc, int ib_lev,
                     int component_src, int component_dest) {

    BL_PROFILE_VAR("move_ibm_marker", move_ibm_marker);

    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        TileIndex index(pti.index(), pti.LocalTileIndex());
        AoS & markers = ib_mc.GetParticles(ib_lev).at(index).GetArrayOfStructs();
        long np = ib_mc.GetParticles(ib_lev).at(index).numParticles();

        for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

            ParticleType & mark = markers[m_index.first];

            for (int d=0; d<AMREX_SPACEDIM; ++d)
                 mark.rdata(component_dest + d) += mot * mark.rdata(component_src + d);

        }
    }

    BL_PROFILE_VAR_STOP(move_ibm_marker);
};



Vector<RealVect> equil_pos(
        int i_ib, Real time, const Geometry & geom
    ) {
    // TODO: make this function work on general planes and orientation -- using
    // the 3D rotation matrix defined in IBMarkerMD.cpp

    int N  = ib_flagellum::n_marker[i_ib];
    Real L = ib_flagellum::length[i_ib];

    Real l_link = L/(N-1);

    const RealVect & x_0 = offset_0[i_ib];

    // using fourier modes => first two nodes reserved as "anchor"
    int N_markers = immbdy::contains_fourier ? N+1 : N;

    Vector<Real> thetas(N_markers - 2);
    for (int i = 0; i < N_markers - 2; ++i) {
        Real th = theta(1, time, i_ib, i);
        thetas[i] = th;
    }

    Real x = x_0[0];
    Real y = x_0[1];
    Real z = x_0[2];

    // TODO: implement general orientation vector
    Real tx = 1.;
    Real ty = 0.;

    Vector<RealVect> marker_positions(N_markers);
    marker_positions[0] = RealVect{x, y, z};
    for (int i=1; i<marker_positions.size()-1; ++i) {

        // TODO: generalize to 3D
        // Real x = x_0[0] + i*l_link;
        Real nx, ny;
        next_node_z(nx, ny, x, y, tx, ty, l_link);
        x = nx;
        y = ny;

        // Compute periodic offset. Will work as long as winding number = 1
        Real x_period = x;
        Real y_period = y;
        if (geom.isPeriodic(0))
            x_period = x < geom.ProbHi(0) ? x : x - geom.ProbLength(0);
        if (geom.isPeriodic(1))
            y_period = y < geom.ProbHi(1) ? y : y - geom.ProbLength(1);

        // marker_positions[i] = RealVect{x_period, x_0[1], x_0[2]};
        marker_positions[i] = RealVect{x_period, y_period, z};

        // TODO: Generalize to 3D
        Real rx, ry;
        rotate_z(rx, ry, tx, ty, thetas[i-1]);
        tx = rx;
        ty = ry;
    }

    Real nx, ny;
    next_node_z(nx, ny, x, y, tx, ty, l_link);
    // Compute periodic offset. Will work as long as winding number = 1
    Real x_period = nx;
    Real y_period = ny;
    if (geom.isPeriodic(0))
        x_period = x < geom.ProbHi(0) ? nx : nx - geom.ProbLength(0);
    if (geom.isPeriodic(1))
        y_period = y < geom.ProbHi(1) ? ny : ny - geom.ProbLength(1);
    // marker_positions[i] = RealVect{x_period, x_0[1], x_0[2]};
    marker_positions[N_markers-1] = RealVect{x_period, y_period, z};

    return marker_positions;
}