#include "INS_functions.H"
#include "common_functions.H"
#include "species.H"
#include "DsmcParticleContainer.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;
using namespace std;

void updateTimeStep(const MultiFab& mfcuInst,
						const MultiFab& mfcuMeans,
						const MultiFab& mfcuVars,
						const MultiFab& mfprimInst,
						const MultiFab& mfprimMeans,
						const MultiFab& mfprimVars,
						const MultiFab& mfcoVars,
						const Geometry& geom,
						Real time,
						int step) {
	BL_PROFILE_VAR("updateTimeStep()",writePlotFile);

	const Real* dx = Geom(lev).CellSize();

	Real umax = 0., vmax = 0., wmax = 0.;
	Real dtmax = 0.0;
	for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {

		const int grid_id = pti.index();
		const int tile_id = pti.LocalTileIndex();
		const Box& tile_box  = pti.tilebox();

		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();
		const long np = particles.numParticles();

		for (int i = 0; i < np; ++ i) {
			ParticleType & part = particles[i];
			Real u = part.rData(FHD_realData::velx);
			Real v = part.rData(FHD_realData::vely);
			Real w = part.rData(FHD_realData::velz);
			umax = std::max(umax,u);
			vmax = std::max(vmax,v);
			wmax = std::max(wmax,w);
		}
	}
	
	dtmax  = umax/dx[0];
	dtmax += vmax/dx[1];
	dtmax += wmax/dx[2];
	dtmax  = 0.2/dtmax; // Courant number of 0.2
	// Want to compare dt across each processor and choose largest one (gather?)
	if(ParallelDescriptor::MyProc() == 0) {
		//ParallelDescriptor::Gather()
		if(dtmax>dt) { 
			dt = dtmax;
			ParallelDescriptor::Bcast(&dt,1,ParallelDescriptor::IOProcessorNumber());
		}
	}
}
