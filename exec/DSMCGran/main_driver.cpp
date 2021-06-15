#include "INS_functions.H"

#include "common_namespace_declarations.H"

#include "gmres_namespace_declarations.H"

#include "species.H"
#include "paramPlane.H"

#include "StructFact.H"

#include "StochMomFlux.H"

#include "hydro_functions.H"

#include "electrostatic.H"

#include "particle_functions.H"

#include "chrono"

#include "iostream"
#include "fstream"

using namespace std::chrono;
using namespace std;
// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{
	// timer for total simulation time
	Real strt_time = ParallelDescriptor::second();

	std::string inputs_file = argv;

	// read in parameters from inputs file into F90 modules
	// we use "+1" because of amrex_string_c_to_f expects a null char termination
	read_common_namelist(inputs_file.c_str(),inputs_file.size()+1);

	// amrex::Print() << "HERE \n";
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
	
	// Box for collision cells (if different)
	// not yet implemented
	// normally would use tiles
	// IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
	// IntVect dom_hi(AMREX_D_DECL(n_cells[0]*-1, n_cells[1]-1, n_cells[2]-1));
	// Box domain(dom_lo, dom_hi);
    
	// how boxes are distrubuted among MPI processes
	DistributionMapping dmap;
	
	MultiFab cuInst;
	MultiFab cuMean;
	MultiFab cuDel;
	MultiFab covar;
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
		/*
		0 	- n - number density
		1	- rho - mass density
		2 	- Jx - x-mom density
		3 	- Jy - y-mom density
		4 	- Jz - z-mom density
		5 	- tau_xx - xx shear stress
		6 	- tau_xy - xy shear stress
		7 	- tau_xz - xz shear stress
		8 	- tau_yy - yy shear stress
		9  - tau_yz - yz shear stress
		10 - tau_zz - zz shear stress
		11 - T_gran - granular temperature
		12 - X_i - mole fraction for spec. i
		13 - Y_i - mass fraction for spec. i
		14 - u_i - x-vel for spec. i
		15 - v_i - y-vel for spec. i
		16 - w_i - z-vel for spec. i
		17 - uu_i 
		18 - uv_i
		19 - uw_i
		20 - vv_i
		21 - vw_i
		22 - ww_i
		23 - T_gran_i
		... (repeat for each add. species)
		*/
		int statsz = (nspecies+1)*12;
		cuInst.define(ba, dmap, statsz, 0); cuInst.setVal(0.);
		cuMean.define(ba, dmap, statsz, 0); cuMean.setVal(0.);
		cuDel.define(ba, dmap, statsz, 0);  cuDel.setVal(0.);
		
		// Variances
		// Each has (nspecies)(nspecies+1)*0.5 data points
		/*
		0  - drho.drho
		1  - drho.dJx
		2  - drho.dJy
		3  - drho.dJz
		4  - drho.dE
		5  - dJx.dJx
		6  - dJx.dJy
		7  - dJx.dJz
		8  - dJy.dJy
		9  - dJy.dJz
		10 - dJz.dJz
		11 - dJx.dE
		12 - dJy.dE
		13 - dJz.dE
		14 - dE.dE
		*/

		/*
		0 	- drho.drho
		1	- drho.du
		2 	- drho.dv
		3 	- drho.dw
		4 	- drho.dTg
		5 	- du.du
		6 	- du.dv
		7 	- du.dw
		8 	- dv.dv
		9  - dv.dw
		10 - dw.dw
		11 - du.dTg
		12 - dv.dTg
		13 - dw.dTg
		14 - dTg.dTg
		*/
		// Add multifabs for variances
		int nnspec = std::ceil((double)nspecies*(nspecies-1)*0.5);
		//MultiFab covar(ba, dmap, (nnspec+1)*15, 0);
		covar.define(ba, dmap, 15, 0);
		
		// just track ones you want
	} else {
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

	Geometry geom (domain ,&realDomain,CoordSys::cartesian,is_periodic.data());
	const Real* dx = geom.CellSize();

	Real dt = fixed_dt;
	// we will want to define a variable timestep based on the max granular temperature
 
	int paramPlaneCount = 6;
	paramPlane paramPlaneList[paramPlaneCount];
	BuildParamplanes(paramPlaneList,paramPlaneCount,realDomain.lo(),realDomain.hi());

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


	int cRange = 0;
	FhdParticleContainer particles(geom, dmap, ba, cRange);
	if (restart < 0 && particle_restart < 0) {
		// Collision Cell Vars
		particles.mfselect.define(ba, dmap, nspecies*nspecies, 0);
		particles.mfselect.setVal(0.);
		
		particles.mfphi.define(ba, dmap, nspecies, 0);
		particles.mfphi.setVal(0.);
		
		particles.mfvrmax.define(ba, dmap, nspecies*nspecies, 0);
		particles.mfvrmax.setVal(0.);
		
		//particles.mfgrantemp.define(ba,dmap,nspecies,0);
		//particles.mfgrantemp.setVal(0.);
		
		particles.InitParticles();

		particles.InitCollisionCells();
		
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
		amrex::Print() << "Time step: " << istep << "\n";
		particles.CalcSelections(dt);
		particles.CollideParticles(dt);
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
        }

		  cuInst.setVal(0.);
        particles.EvaluateStats(cuInst, cuMean, cuDel, covar, statsCount);
        statsCount++;
 
        if ((n_steps_skip > 0 && istep == n_steps_skip) ||
            (n_steps_skip < 0 && istep%n_steps_skip == 0) ) {
            
            //particleMeans.setVal(0.0);
            //particleVars.setVal(0.0);

            // Print() << "Resetting stat collection.\n";

            statsCount = 1;
        }
 
        // timer for time step
        /*Real time2 = ParallelDescriptor::second() - time1;
        ParallelDescriptor::ReduceRealMax(time2);
        amrex::Print() << "Advanced step " << istep << " in " << time2 << " seconds\n";*/
    }


    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << " seconds" << std::endl;


	 if(particle_input<0) {particles.OutputParticles();} // initial condition
	 //particles.OutputParticles();
}
