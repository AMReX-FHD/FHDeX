#ifndef _PARTICLE_DATA_H_
#define _PARTICLE_DATA_H_

#include <string>
#include <iostream>

// #include <IndexDefines.H>
#include <AMReX_Particles.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Gpu.H>

#include <StochasticPC.H>

/**
 * Container holding many of the particle-related data and options
 */

struct ParticleData {
  public:
    void init_particle_params(int max_level)
    {
        amrex::ParmParse pp(pp_prefix);

        use_particles = 0;
        pp.query("use_particles", use_particles);

        if (use_particles) {
            if (max_level != 1) {
                amrex::Abort("When using particles you must have max_level == 1");
            }
        }

    }

    void init_particles(amrex::ParGDBBase* a_gdb, const amrex::BoxArray& grown_fba,
                        amrex::MultiFab& phi_crse, amrex::MultiFab& phi_fine)
    {
        // Initialize stochastic particles and phi at level 1
        if (use_particles) {
            const amrex::Geometry& geom_fine = a_gdb->Geom(1);
            const auto dx_fine = geom_fine.CellSize();
#if (AMREX_SPACEDIM == 2)
            amrex::Real cell_vol_fine = dx_fine[0]*dx_fine[1];
#else
            amrex::Real cell_vol_fine = dx_fine[0]*dx_fine[1]*dx_fine[2];
#endif

            // **********************************************************************
            // Define the particle class with the default gdb, then change the particle
            //        boxArray (and DM) to be the grown ba and associated dm
            // **********************************************************************
            stochastic_particles = std::make_unique<StochasticPC>(a_gdb);
            stochastic_particles->SetParticleBoxArray(1,grown_fba);
            amrex::DistributionMapping grown_dm { grown_fba, amrex::ParallelDescriptor::NProcs() };
            stochastic_particles->SetParticleDistributionMap(1,grown_dm);

            // **********************************************************************
            // Fill a temporary phi_grown with phi from "coarse" level
            // **********************************************************************
            amrex::MultiFab phi_grown(grown_fba, grown_dm, 1, 0);
            phi_crse.FillBoundary(geom_fine.periodicity());
            phi_grown.ParallelCopy(phi_crse,0,0,1,1,0,geom_fine.periodicity());

            // **********************************************************************
            // Use this grown boxArray to hold the particles
            // **********************************************************************
            stochastic_particles->InitParticles(phi_grown);

            amrex::Print() << "Before Removal we have  " << stochastic_particles->TotalNumberOfParticles() <<
                              " stochastic particles." << std::endl;

            // **********************************************************************
            // Remove particles not in the real fine boxArray
            // **********************************************************************
            stochastic_particles->RemoveParticlesNotInBA(phi_fine.boxArray());

            // **********************************************************************
            // Redistribute so every particle is now in the "correct" place
            // **********************************************************************
            stochastic_particles->Redistribute();
            amrex::Print() << "Initialized " << stochastic_particles->TotalNumberOfParticles() <<
                              " stochastic particles." << std::endl;

            // **********************************************************************
            // Define phi at level 1 from the particle count
            // **********************************************************************
            int lev = 1;
            phi_fine.setVal(0.0);
            stochastic_particles->Increment(phi_fine, lev);
            phi_fine.mult(1./cell_vol_fine);
            amrex::Print() << "Used the particle locations to define phi at level 1 " << std::endl;
        }
    }

    void regrid_particles(const amrex::BoxArray& new_fine_ba, const amrex::BoxArray& new_ungrown_ba,
                          const amrex::BoxArray& old_fine_ba, amrex::MultiFab& phi_fine)
    {
        amrex::DistributionMapping new_ungrown_dm { new_ungrown_ba, amrex::ParallelDescriptor::NProcs() };
        stochastic_particles->SetParticleBoxArray(1, new_ungrown_ba);
        stochastic_particles->SetParticleDistributionMap(1, new_ungrown_dm);

        stochastic_particles->Redistribute();
        stochastic_particles->RemoveParticlesAtLevel(0);
        amrex::Real factor = .9;
        stochastic_particles->AddParticles(phi_fine, old_fine_ba, factor);
        phi_fine.setVal(0.0);
        const int lev = 1;
        const auto dx = stochastic_particles->Geom(lev).CellSizeArray();

#if (AMREX_SPACEDIM == 2)
    const amrex::Real cell_vol = dx[0]*dx[1];
#else
    const amrex::Real cell_vol = dx[0]*dx[1]*dx[2];
#endif
        stochastic_particles->Increment(phi_fine, lev);
        phi_fine.mult(1./cell_vol);

        amrex::DistributionMapping new_dm { new_fine_ba, amrex::ParallelDescriptor::NProcs() };
        stochastic_particles->SetParticleBoxArray(1, new_fine_ba);
        stochastic_particles->SetParticleDistributionMap(1, new_dm);
    }

    void writePlotFile (const std::string& a_fname, amrex::MultiFab const& phi_fine) const
    {
        if (use_particles) {
            bool is_checkpoint = false;
            amrex::MultiFab phi_grown(stochastic_particles->ParticleBoxArray(1),
                                      stochastic_particles->ParticleDistributionMap(1),1,0);
            const amrex::Geometry& geom_fine = stochastic_particles->GetParGDB()->Geom(1);
            phi_grown.ParallelCopy(phi_fine,0,0,1,1,0,geom_fine.periodicity());
            stochastic_particles->ColorParticlesWithPhi(phi_grown);
            stochastic_particles->Checkpoint( a_fname, "particles", is_checkpoint, particle_varnames );
        }
    }

    void Checkpoint (const std::string& a_fname) const
    {
        if (use_particles) {
            bool is_checkpoint = true;
            stochastic_particles->Checkpoint(a_fname, "stochastic_particles", is_checkpoint, particle_varnames);
        }
    }

    void Restart(amrex::ParGDBBase* a_gdb, std::string& restart_file)
    {
        if (use_particles) {
            stochastic_particles = std::make_unique<StochasticPC>(a_gdb);
            std::string particle_file("stochastic_particles");
            stochastic_particles->Restart(restart_file, particle_file);
        }
    }

    void advance_particles(int lev, amrex::Real dt_lev, amrex::Real cell_vol,
                           amrex::MultiFab& phi_crse_old,
                           amrex::MultiFab& phi_crse_new,
                           amrex::MultiFab& phi_fine)
    {
        // Update stochastic particles
        if (use_particles) {
            AMREX_ALWAYS_ASSERT(lev == 1);
            const amrex::Geometry& geom_fine = stochastic_particles->GetParGDB()->Geom(lev);

            // **********************************************************************
            // Fill a temporary phi_grown with phi from "coarse" level
            // **********************************************************************
            // amrex::Print() << "PARTICLE BOXARRAY " << stochastic_particles->ParticleBoxArray(0) << std::endl;
            // amrex::Print() << "PARTICLE BOXARRAY " << stochastic_particles->ParticleBoxArray(1) << std::endl;
            amrex::MultiFab phi_grown(stochastic_particles->ParticleBoxArray(1),
                                      stochastic_particles->ParticleDistributionMap(1),1,0);
            phi_crse_old.FillBoundary(geom_fine.periodicity());
            phi_grown.ParallelCopy(phi_crse_old,0,0,1,1,0,geom_fine.periodicity());

            // **********************************************************************
            // Create new particles in first layer of ghost cells around the level 1 grids
            // **********************************************************************
            amrex::BoxArray ba_to_exclude(phi_fine.boxArray());

            amrex::Print() << "Adding particles in advance " << std::endl;

            amrex::Real factor = .9;
            stochastic_particles->AddParticles(phi_grown, ba_to_exclude,factor);

            // **********************************************************************
            // Advect the particles in the valid region and in the ghost region
            // **********************************************************************
            stochastic_particles->AdvectWithRandomWalk(lev, dt_lev);

            // **********************************************************************
            // Account for particles that went from crse --> fine BEFORE Redistribute
            // **********************************************************************
            amrex::MultiFab phi_grown2(phi_grown.boxArray(),phi_grown.DistributionMap(),1,1);
            phi_grown2.setVal(0.);

            stochastic_particles->RefluxCrseToFine(phi_fine.boxArray(), phi_grown2);

            phi_grown2.mult(1./cell_vol,0,1,1);

            phi_crse_new.ParallelAdd(phi_grown2, 0, 0, 1, amrex::IntVect::TheUnitVector(), amrex::IntVect::TheUnitVector());

            // **********************************************************************
            // Put particles on the correct levels and in the right boxes
            // **********************************************************************
            stochastic_particles->Redistribute();

            // **********************************************************************
            // Remove particles on the coarsest level
            // **********************************************************************

            //amrex::Print() << "Before Removal we have  " << stochastic_particles->TotalNumberOfParticles() <<
            //                  " stochastic particles." << std::endl;

            int crse_lev = 0;
            stochastic_particles->RemoveParticlesAtLevel(crse_lev);

            //amrex::Print() << " After Removal we have  " << stochastic_particles->TotalNumberOfParticles() <<
            //                  " stochastic particles." << std::endl;

            // **********************************************************************
            // Account for particles that went from fine --> crse
            // **********************************************************************

            phi_grown.setVal(0.);

            stochastic_particles->RefluxFineToCrse(phi_fine.boxArray(), phi_grown);

            phi_grown.mult(1./cell_vol);
#if (AMREX_SPACEDIM == 2)
            print_state(phi_grown, amrex::IntVect(31,16));
#else
            print_state(phi_grown, amrex::IntVect(31,16,0));
#endif

            phi_crse_new.ParallelAdd(phi_grown,0,0,1);

            // **********************************************************************
            // Remove particles that are not in the valid boxArray used by phi_fine
            // **********************************************************************
            stochastic_particles->RemoveParticlesNotInBA(phi_fine.boxArray());

            // **********************************************************************
            // Define phi at level 1 from the particle count
            // **********************************************************************
            phi_grown.setVal(0.);
            stochastic_particles->Increment(phi_grown, lev);

            // *****************************************************************************
            // Overwrite phi at level 0 from phi at level 1 only in the valid level 1 region
            // *****************************************************************************
            phi_fine.ParallelCopy(phi_grown,0,0,1);

            phi_fine.mult(1./cell_vol);
            amrex::Print() << "Used the particle locations to define phi at level 1 " << std::endl;
        }
    }

    void Redistribute()
    {
        if (use_particles) {
            stochastic_particles->Redistribute();
        }
    }

    std::string pp_prefix {"amr"};

    int use_particles;
    std::unique_ptr<StochasticPC> stochastic_particles;
    amrex::Vector<std::string> particle_varnames = {"xold", "yold", "zold"};
};
#endif
