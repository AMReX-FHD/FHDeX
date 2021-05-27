#include "INS_functions.H"

#include "common_namespace_declarations.H"

#include "species.H"
#include "paramPlane.H"

#include "StructFact.H"

#include "particle_functions.H"

#include "chrono"

using namespace std::chrono;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{
	// timer for total simulation time
	Real strt_time = ParallelDescriptor::second();

	std::string inputs_file = argv;

	// read in parameters from inputs file into F90 modules
	// we use "+1" because of amrex_string_c_to_f expects a null char termination
	read_common_namelist(inputs_file.c_str(),inputs_file.size()+1);


	// copy contents of F90 modules to C++ namespaces
	InitializeCommonNamespace();
	InitializeGmresNamespace();

	int step = 1;
	Real time = 0.;
	int statsCount = 1;

	/*
		Terms prepended with a 'C' are related to the particle grid; only used for finding neighbor lists
      Those with 'P' are for the electostatic grid.
      Those without are for the fluid grid.
      The particle grid and es grids are created as a corsening or refinement of the fluid grid.
	*/

	// BoxArray for the particles (collision grid)
	BoxArray ba;
	// N_cells -< collision cell count on each dim
    
	// Box for the fluid
	IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
	IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
	Box domain(dom_lo, dom_hi);
    
	// how boxes are distrubuted among MPI processes
	DistributionMapping dmap;

	// MFs for storing particle statistics
	// A lot of these relate to gas kinetics, but many are still useful so leave in for now.
	MultiFab particleMeans;
	MultiFab particleVars;
	MultiFab particleInstant;
	
	// MFs for DSMC
	MultiFab selectionsCell; // define(..,..,max_species*max_species)
   MultiFab vrmax
    
	if (restart < 0) {
		if (seed > 0) {
			// initializes the seed for C++ random number calls
			InitRandom(seed+ParallelDescriptor::MyProc());
		} else if (seed == 0) {
			// initializes the seed for C++ random number calls based on the clock
			auto now = time_point_cast<nanoseconds>(system_clock::now());
			int randSeed = now.time_since_epoch().count();
			// broadcast the same root seed to all processors
			ParallelDescriptor::Bcast(&randSeed,1,ParallelDescriptor::IOProcessorNumber());
			InitRandom(randSeed+ParallelDescriptor::MyProc());
		} else {
			Abort("Must supply non-negative seed");
		}

		// Initialize the boxarray "ba" from the single box "bx"
		ba.define(domain);

		// Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
		// note we are converting "Vector<int> max_grid_size" to an IntVect
		ba.maxSize(IntVect(max_grid_size));

		// how boxes are distrubuted among MPI processes
		dmap.define(ba);
  
		// Statistics
		particleMeans.define(ba, dmap, 8+nspecies, 0);
		particleMeans.setVal(0.);

		particleVars.define(ba, dmap, 8+nspecies, 0);
		particleVars.setVal(0.);

		particleInstant.define(ba, dmap, 8+nspecies, 0);
		particleInstant.setVal(0.);
		
		// Collision Cell Vars
		selectionsCell.define(ba, dmap, nspecies*nspecies, 0);
		selectionsCell.setVal(0.);
		
		vrmax.define(ba, dmap, nspecies, 0);
		vrmax.setVal(0.); // we will reset this later in initParticles
	}
	else {
		// restart from checkpoint
	}

	Vector<int> is_periodic (AMREX_SPACEDIM,0);
	for (int i=0; i<AMREX_SPACEDIM; ++i) {
		if (bc_vel_lo[i] == -1 && bc_vel_hi[i] == -1) {
			is_periodic [i] = -1;
		}
	}

	// This defines a Geometry object
	RealBox realDomain({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
		{AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

	Geometry geom (domain ,&realDomain,CoordSys::cartesian,is_periodic.  data());
	const Real* dx = geom.CellSize();

	Real dt = fixed_dt;
	// we will want to define a variable timestep based on the max granular temperature
 
    std::ifstream planeFile("paramplanes.dat");
    int fileCount;
    planeFile >> fileCount;
    planeFile.close();

    int paramPlaneCount = fileCount+6;
    paramPlane paramPlaneList[paramPlaneCount];
    BuildParamplanes(paramPlaneList,paramPlaneCount,realDomain.lo(),realDomain.hi());

   // IBMarkerContainerBase default behaviour is to do tiling. Turn off here:


	// Particle tile size
	Vector<int> ts(BL_SPACEDIM);
    
	for (int d=0; d<AMREX_SPACEDIM; ++d) {        
		if (max_particle_tile_size[d] > 0) {
			ts[d] = max_particle_tile_size[d];
		}
		else {
			ts[d] = max_grid_size[d];
		}
	}

	ParmParse pp ("particles");
	pp.addarr("tile_size", ts);

	//int num_neighbor_cells = 4; replaced by input var
	//Particles! Build on geom & box array for collision cells/ poisson grid?

	int cRange = 0;

	FhdParticleContainer particles(geom, dmap, ba, cRange);

	if (restart < 0 && particle_restart < 0) {
		// Pass Multifab to set vrmax
		particles.InitParticles(Real T_init, MutliFab& vrmax);
		
	}
	else {
        //load from checkpoint
	}

	// cell centered real coordinates - es grid
	MultiFab RealCenteredCoords;
	RealCenteredCoords.define(ba, dmap, AMREX_SPACEDIM, 0);

	//FindCenterCoords(RealCenteredCoords, geom);
	Real init_time = ParallelDescriptor::second() - strt_time;
	ParallelDescriptor::ReduceRealMax(init_time);
	amrex::Print() << "Initialization time = " << init_time << " seconds " << std::endl;

	for (int istep=step; istep<=max_step; ++istep) {
		// timer for time step
		Real time1 = ParallelDescriptor::second();
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Collide Particles
		if (collide_tog != 0) {
			particles.CollideParticles(MultiFab& selectionsCell, dt);
		}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Move Particles
		// total particle move (1=single step, 2=midpoint)
		if (move_tog != 0) {
			particles.MoveParticlesCPP(dt, paramPlaneList, paramPlaneCount);

            // reset statistics after step n_steps_skip
            // if n_steps_skip is negative, we use it as an interval
            if ((n_steps_skip > 0 && istep == n_steps_skip) ||
                (n_steps_skip < 0 && istep%n_steps_skip == 0) ) {

                
            }
            else {
                
            }

            Print() << "Finish move.\n";
        }

        particles.EvaluateStats(particleInstant, particleMeans, particleVars, dt,statsCount);
        statsCount++;

        if (istep%plot_int == 0) {
            
            WritePlotFile(istep, time, geom, particleInstant, particleMeans, particleVars, particles);

        }
 
        if ((n_steps_skip > 0 && istep == n_steps_skip) ||
            (n_steps_skip < 0 && istep%n_steps_skip == 0) ) {
            
            particleMeans.setVal(0.0);
            particleVars.setVal(0.0);

            Print() << "Resetting stat collection.\n";

            statsCount = 1;
        }
 
        // timer for time step
        Real time2 = ParallelDescriptor::second() - time1;
        ParallelDescriptor::ReduceRealMax(time2);
        amrex::Print() << "Advanced step " << istep << " in " << time2 << " seconds\n";
        
        time = time + dt;
        // MultiFab memory usage
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        amrex::Long min_fab_megabytes  = amrex::TotalBytesAllocatedInFabsHWM()/1048576;
        amrex::Long max_fab_megabytes  = min_fab_megabytes;

        ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

        amrex::Print() << "High-water FAB megabyte spread across MPI nodes: ["
                       << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";

        min_fab_megabytes  = amrex::TotalBytesAllocatedInFabs()/1048576;
        max_fab_megabytes  = min_fab_megabytes;

        ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
        ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

        amrex::Print() << "Curent     FAB megabyte spread across MPI nodes: ["
                       << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";
        
    }
    ///////////////////////////////////////////
        //test change
    // timer for total simulation time
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << " seconds" << std::endl;

}
