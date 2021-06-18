#include "INS_functions.H"
#include "common_namespace_declarations.H"
#include "gmres_namespace_declarations.H"
#include "species.H"
#include "paramPlane.H"
#include "StructFact.H"
#include "particle_functions.H"
#include "chrono"
#include "iostream"
#include "fstream"
#include "DsmcParticleContainer.H"
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

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

	InitializeCommonNamespace();
	InitializeGmresNamespace();

	BoxArray ba;
	IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
	IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
	Box domain(dom_lo, dom_hi);	
	DistributionMapping dmap;

	// Will likely need to redo this
	if (restart < 0) {
		if (seed > 0) {
			InitRandom(seed+ParallelDescriptor::MyProc());
		} else if (seed == 0) {
			auto now = time_point_cast<nanoseconds>(system_clock::now());
			int randSeed = now.time_since_epoch().count();
			// broadcast the same root seed to all processors
			ParallelDescriptor::Bcast(&randSeed,1,ParallelDescriptor::IOProcessorNumber());
			InitRandom(randSeed+ParallelDescriptor::MyProc());
		} else {
			Abort("Must supply non-negative seed");
		}

		ba.define(domain);
		ba.maxSize(IntVect(max_grid_size));
		dmap.define(ba);
		
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

	// Currently overwritten later
	Real dt = fixed_dt;

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
	
	/*
		Conserved Vars:
		0  - rho = (1/V) += m
		1  - Jx  = (1/V) += mu
		2  - Jy  = (1/V) += mv
		3  - Jz  = (1/V) += mw
		4  - K   = (1/V) += m|v|^2
		... (repeat for each species)
	*/
	
	/*
	   Primitive Vars:
		0	- n   (X_ns)
		1  - rho (Y_ns)
		2  - u   (u_ns)
		3  - v   (v_ns)
		4  - w   (w_ns)
		5  - uu  (uu_ns)
		6  - uv  (uv_ns)
		7  - uw  (uw_ns)
		8  - vv  (vv_ns)
		9  - vw  (vw_ns)
		10 - ww  (ww_ns)
		11 - T   (T_ns)
		12 - P   (P_ns)
		13 - E   (E_ns)
		... (repeat for each species)
	*/   	
   MultiFab cuConInst,  cuConMeans;
	MultiFab cuPrimInst, cuPrimMeans;
	MultiFab cuVars;
	
	int ncon  = (nspecies+1)*5;
	cuConInst.define(ba, dmap, ncon, 0);    cuConInst.setVal(0.);
	cuConMeans.define(ba, dmap, ncon, 0);   cuConMeans.setVal(0.);
	
	int nprim = (nspecies+1)*14;
	cuPrimInst.define(ba, dmap, nprim, 0);   cuPrimInst.setVal(0.);
	cuPrimMeans.define(ba, dmap, nprim, 0);  cuPrimMeans.setVal(0.);
	cuVars.define(ba, dmap, ncon+nprim, 0);       cuVars.setVal(0.);
		
	// Covariances
	/*
		0  - drho.dJx
		1  - drho.dJy
		2  - drho.dJz
		3  - drho.dT
		4  - drho.d(rho*E)
		5  - dJx.dJy
		6  - dJx.dJz
		7  - dJy.dJz
		8  - dJx.d(rho*E)
		9  - dJy.d(rho*E)
		10 - dJz.d(rho*E)
		11 - drho.du
		12 - drho.dv
		13 - drho.dw
		14 - du.dv
		15 - du.dw
		16 - dv.dw
		17 - drho.dT
		18 - du.dT
		19 - dv.dT
		20 - dw.dT
	*/
	
	// Add multifabs for variances
	int ncovar = 21;
	MultiFab coVars(ba, dmap, ncovar, 0);   coVars.setVal(0.);
   
   //////////////////////////////////////
   // Structure Factor Setup
   //////////////////////////////////////
      
   // Output all primitives for structure factor
   int nvarstruct = nprim;
	const Real* dx = geom.CellSize();
	int nstruct = std::ceil((double)nvarstruct*(nvarstruct+1)/2);
	// scale SF results by inverse cell volume
   Vector<Real> var_scaling(nstruct);
   for (int d=0; d<var_scaling.size(); ++d) {var_scaling[d] = 1./(dx[0]*dx[1]*dx[2]);}

   // Structure Factor labels
	Vector< std::string > cu_struct_names(nvarstruct);
	int cnt = 0;
	std::string varname;
	cu_struct_names[cnt++] = "n";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("X",ispec,2);
   }
   cu_struct_names[cnt++] = "rho";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("rho",ispec,2);
   }
   cu_struct_names[cnt++] = "u";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("u",ispec,2);
   }
   cu_struct_names[cnt++] = "v";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("v",ispec,2);
   }
   cu_struct_names[cnt++] = "w";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("w",ispec,2);
   }
   cu_struct_names[cnt++] = "uu";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("uu",ispec,2);
   }
   cu_struct_names[cnt++] = "uv";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("uv",ispec,2);
   }
   cu_struct_names[cnt++] = "uw";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("uw",ispec,2);
   }
	cu_struct_names[cnt++] = "vv";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("vv",ispec,2);
   }
   cu_struct_names[cnt++] = "vw";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("vw",ispec,2);
   }
   cu_struct_names[cnt++] = "ww";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("ww",ispec,2);
   }
   cu_struct_names[cnt++] = "T";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("T",ispec,2);
   }
   cu_struct_names[cnt++] = "P";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("P",ispec,2);
   }
   cu_struct_names[cnt++] = "E";
	for (int ispec=0; ispec<nspecies; ispec++) {
     	cu_struct_names[cnt++] = amrex::Concatenate("E",ispec,2);
   }
   
   // Structure Factor
	StructFact structFactPrim  (ba, dmap, cu_struct_names, var_scaling);
	MultiFab   structFactPrimMF(ba, dmap,      nvarstruct,           0);
	
	if (restart < 0 && particle_restart < 0) {
		// Collision Cell Vars
		particles.mfselect.define(ba, dmap, nspecies*nspecies, 0);
		particles.mfselect.setVal(0.);
		
		particles.mfphi.define(ba, dmap, nspecies, 0);
		particles.mfphi.setVal(0.);
		
		particles.mfvrmax.define(ba, dmap, nspecies*nspecies, 0);
		particles.mfvrmax.setVal(0.);
		// overwrite dt
		particles.InitParticles(dt);

		particles.InitCollisionCells();
		amrex::Print() << "Overwritten dt so Courant number <1: " << dt << "\n";
	}
	
	else {

	}

	Real init_time = ParallelDescriptor::second() - strt_time;
	ParallelDescriptor::ReduceRealMax(init_time);
	amrex::Print() << "Initialization time = " << init_time << " seconds " << std::endl;

	// Frequency of plot and checkpoints
	int plot_step = 10;
	//int check_step = 10;
	Real time = 0.;
	int statsCount = 1;
	for (int istep=0; istep<=max_step; ++istep) {
		Real tbegin = ParallelDescriptor::second();
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Checkpoints
		//if (plot_int > 0 && istep > 0 && istep%plot_step == 0) {
           
		//}

		//if (istep==1 && restart > 0) {
      //     ReadCheckPoint(istep, statsCount, time, geom,
      //     						cu, cuMeans, cuVars, coVars);//,
           //						particles.mfselect, particles.mfphi, particles.mfvrmax);
		//}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Collide Particles
		//amrex::Print() << "Collisions per particles per species: " << 
		//	particles.CountedCollision[
		particles.CalcSelections(dt);
		particles.CollideParticles(dt);
		
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Move Particles
		particles.MoveParticlesCPP(dt, paramPlaneList, paramPlaneCount);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Stats
		cuConInst.setVal(0.); cuPrimInst.setVal(0.);
		particles.EvaluateStats(cuConInst,
										cuConMeans,
		                        cuPrimInst,
		                        cuPrimMeans,
		                        cuVars,
		                        coVars,
		                        statsCount++);
		                        
		if(istep%n_steps_skip == 0 && istep > 0) {
			MultiFab::Copy(structFactPrimMF,cuPrimInst,0,0,nvarstruct,0);
			structFactPrim.FortStructure(structFactPrimMF,geom,fft_type);
		}
		
		if(istep%plot_int == 0) {
			particles.writePlotFile(cuConInst,
											  cuConMeans,
		                        	  cuPrimInst,
		                        	  cuPrimMeans,
		                        	  cuVars,
		                        	  coVars,
		                        	  geom,
		                        	  time,
		                        	  istep);
		}
		
		//if(istep%struct_fact_int ==0 && istep>0) {structFactPrim.WritePlotFile(istep,time,geom,"plt_SF_prim");}
 
		Real tend = ParallelDescriptor::second() - tbegin;
		ParallelDescriptor::ReduceRealMax(tend);
		amrex::Print() << "Advanced step " << istep << " in " << tend << " seconds\n";
	}

	Real stop_time = ParallelDescriptor::second() - strt_time;
	ParallelDescriptor::ReduceRealMax(stop_time);
	amrex::Print() << "Run time = " << stop_time << " seconds" << std::endl;


	//if(particle_input<0 && ParallelDescriptor::MyProc() == 0) {particles.OutputParticles();} // initial condition
}
