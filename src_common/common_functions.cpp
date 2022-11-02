#include "common_functions.H"

#include <AMReX_ParmParse.H>


void InitializeCommonNamespace() {

    BL_PROFILE_VAR("InitializeCommonNamespace()",InitializeCommonNameSpace);
    
    nodal_flag_dir.resize(AMREX_SPACEDIM);
    nodal_flag_edge.resize(AMREX_SPACEDIM);

    for (int i=0; i<AMREX_SPACEDIM; ++i) {

        //_______________________________________________________________________
        // Designates data on nodes
        nodal_flag[i] = 1;

        //_______________________________________________________________________
        // Designates data on faces
        nodal_flag_x[i] = int(i==0);
        nodal_flag_y[i] = int(i==1);
        nodal_flag_z[i] = int(i==2);

        // Enable indexing flags above in loops
        AMREX_D_TERM(nodal_flag_dir[0][i] = nodal_flag_x[i];,
                     nodal_flag_dir[1][i] = nodal_flag_y[i];,
                     nodal_flag_dir[2][i] = nodal_flag_z[i];);

        //_______________________________________________________________________
        // Designates data on edges
        nodal_flag_xy[i] = int(i==0 || i==1);
        nodal_flag_xz[i] = int(i==0 || i==2);
        nodal_flag_yz[i] = int(i==1 || i==2);

        // Enable indexing flags above in loops
        AMREX_D_TERM(nodal_flag_edge[0][i] = nodal_flag_xy[i];,
                     nodal_flag_edge[1][i] = nodal_flag_xz[i];,
                     nodal_flag_edge[2][i] = nodal_flag_yz[i];);
    }

    max_grid_size.resize(AMREX_SPACEDIM);
    max_grid_size_structfact.resize(AMREX_SPACEDIM);
    max_particle_tile_size.resize(AMREX_SPACEDIM);
    dof.resize(MAX_SPECIES);
    u_init.resize(2);

    // boundary condition flags
    //bc_vel_lo.resize(AMREX_SPACEDIM);
    //bc_vel_hi.resize(AMREX_SPACEDIM);
    //bc_es_lo.resize(AMREX_SPACEDIM);
    //bc_es_hi.resize(AMREX_SPACEDIM);
    //bc_mass_lo.resize(AMREX_SPACEDIM);
    //bc_mass_hi.resize(AMREX_SPACEDIM);
    //bc_therm_lo.resize(AMREX_SPACEDIM);
    //bc_therm_hi.resize(AMREX_SPACEDIM);

    wallspeed_lo.resize((AMREX_SPACEDIM-1)*AMREX_SPACEDIM);
    wallspeed_hi.resize((AMREX_SPACEDIM-1)*AMREX_SPACEDIM);

    max_grid_projection.resize(AMREX_SPACEDIM-1);

    density_weights.resize(MAX_SPECIES);
    shift_cc_to_boundary.resize(AMREX_SPACEDIM*LOHI);

    mass.resize(MAX_SPECIES);
    nfrac.resize(MAX_SPECIES);
    particle_count.resize(MAX_SPECIES);
    p_move_tog.resize(MAX_SPECIES);
    p_force_tog.resize(MAX_SPECIES);
    particle_n0.resize(MAX_SPECIES);

    qval.resize(MAX_SPECIES);

    diff.resize(MAX_SPECIES);

    //eamp.resize(3);    
    //efreq.resize(3);
    //ephase.resize(3);

    char temp_plot_base_name[128];
    char temp_chk_base_name[128];

    initialize_common_namespace(prob_lo.begin(), prob_hi.begin(), n_cells.data(),
                                max_grid_size.dataPtr(), max_grid_size_structfact.dataPtr(),
                                max_particle_tile_size.dataPtr(), &cell_depth, ngc.getVect(),
                                &nvars, &nprimvars,
                                &membrane_cell, &cross_cell, &transmission,
                                qval.dataPtr(), pkernel_fluid.begin(), pkernel_es.begin(),
                                eskernel_fluid.begin(), eskernel_beta.begin(),
                                &fixed_dt, &cfl, &rfd_delta, &max_step,
                                &plot_int, &plot_stag, temp_plot_base_name, 128,
                                &chk_int, temp_chk_base_name, 128,
                                &prob_type, &restart, &particle_restart, &print_int, &project_eos_int,
                                grav.data(), &nspecies, molmass.data(), diameter.data(),
                                dof.dataPtr(), hcv.data(), hcp.data(),
                                rhobar.data(),
                                &rho0, &variance_coef_mom, &variance_coef_mass, &k_B, &Runiv,
                                T_init.begin(),
                                &algorithm_type,  &advection_type,
                                &barodiffusion_type, &use_bl_rng, &seed,
                                &seed_momentum, &seed_diffusion, &seed_reaction,
                                &seed_init_mass,
                                &seed_init_momentum, &visc_coef, &visc_type,
                                &filtering_width, &stoch_stress_form, u_init.dataPtr(),
                                &perturb_width, &smoothing_width, &initial_variance_mom,
                                &initial_variance_mass, &domega,
                                bc_vel_lo.data(), bc_vel_hi.data(),
                                bc_es_lo.data(), bc_es_hi.data(),
                                bc_mass_lo.data(), bc_mass_hi.data(),
                                bc_therm_lo.data(), bc_therm_hi.data(),
                                p_lo.data(), p_hi.data(),
                                t_lo.data(), t_hi.data(),
                                bc_Yk_x_lo.data(), bc_Yk_x_hi.data(),
                                bc_Yk_y_lo.data(), bc_Yk_y_hi.data(),
                                bc_Yk_z_lo.data(), bc_Yk_z_hi.data(),
                                bc_Xk_x_lo.data(), bc_Xk_x_hi.data(),
                                bc_Xk_y_lo.data(), bc_Xk_y_hi.data(),
                                bc_Xk_z_lo.data(), bc_Xk_z_hi.data(),
                                wallspeed_lo.dataPtr(), wallspeed_hi.dataPtr(),
                                potential_lo.data(), potential_hi.data(),
                                &struct_fact_int, &radialdist_int, &cartdist_int,
                                &n_steps_skip, &binSize, &searchDist,
				&project_dir, &slicepoint, max_grid_projection.dataPtr(),
                                &histogram_unit,
                                density_weights.dataPtr(), shift_cc_to_boundary.dataPtr(),
                                &particle_placement, particle_count.dataPtr(),
                                p_move_tog.dataPtr(), p_force_tog.dataPtr(),
                                p_int_tog.begin(),p_int_tog_wall.begin(), &particle_neff,
                                particle_n0.dataPtr(), mass.dataPtr(), nfrac.dataPtr(),
                                &permittivity,
                                &induced_charge_eo, &zero_eps_on_wall_type, &zero_eps_wall_left_end, &zero_eps_wall_right_start,
                                &wall_mob,rmin.begin(),rmax.begin(), eepsilon.begin(), sigma.begin(),rmin_wall.begin(),rmax_wall.begin(), eepsilon_wall.begin(), sigma_wall.begin(),
                                &poisson_verbose, &poisson_bottom_verbose, &poisson_max_iter,
                                &poisson_rel_tol, &particle_grid_refine, &es_grid_refine,
                                diff.dataPtr(), &all_dry, &fluid_tog, &es_tog, &drag_tog, &move_tog, &rfd_tog,
                                &dry_move_tog, &sr_tog, &graphene_tog, &crange, &thermostat_tog, &zero_net_force, &body_force,
                                &images, eamp.data(), efreq.data(), ephase.data(), F_body.data(),
                                &plot_ascii, &plot_means, &plot_vars,
                                &solve_chem, &diffcoeff, &scaling_factor,
                                &source_strength, &regrid_int, &do_reflux, &particle_motion,
                                &turb_a, &turb_b, &turbForcing);

    plot_base_name = temp_plot_base_name;
    chk_base_name = temp_chk_base_name;

    ParmParse pp;

    // read in from command line
    pp.query("max_step",max_step);

    // copy value into fortran namelist
    set_max_step(&max_step);

    // read in from command line
    pp.query("domega",domega);

    // copy value into fortran namelist
    set_domega(&domega);
    
}
