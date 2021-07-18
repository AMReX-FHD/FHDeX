#include "INS_functions.H"
#include "common_functions.H"
#include "DsmcParticleContainer.H"
#include <sstream>
#include <string>
#include <fstream>

using namespace std;

void FhdParticleContainer::InitParticles(Real & dt) {
	const int lev = 0;
	const Geometry& geom = Geom(lev);

	int pcount = 0;

	bool proc0_enter = true;
    
	// Search for max relative speed
	// ... estimate as double the mag of a single particle speed
	std::array<Real, 3> vpart = {0., 0., 0.};
	Real spd, csx;
	Real csxmax = 0., spdmax = 0.;
	Real umax = 0, vmax = 0, wmax = 0;
	maxDiam = 0.;
	
	//tTg = 0;
	for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi) {
		// take tile/box
		const Box& tile_box  = mfi.tilebox();
		const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
		const int grid_id = mfi.index();
		const int tile_id = mfi.LocalTileIndex();
		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();
		
		//Assuming tile=box for now, i.e. no tiling.
		IntVect smallEnd = tile_box.smallEnd();
		IntVect bigEnd = tile_box.bigEnd();
		
		if(ParallelDescriptor::MyProc() == 0 && mfi.LocalTileIndex() == 0 && proc0_enter) {
			proc0_enter = false;
			// Loop over species (mean diameters)
			for(int ispec=0; ispec < nspecies; ispec++) {
				Real dmean = diameter[ispec];
				Real dstd  = qval[ispec];
				Real rhomean = mass[ispec]/(pi_usr*pow(dmean,3.0)/6.0);
				for (int ipart=0; ipart<properties[ispec].total;ipart++) {
					ParticleType p;
					p.id()  = ParticleType::NextID();
					p.cpu() = ParallelDescriptor::MyProc();
					p.idata(FHD_intData::sorted) = -1;
					p.idata(FHD_intData::i) = -100;
					p.idata(FHD_intData::j) = -100;
					p.idata(FHD_intData::k) = -100;
					p.rdata(FHD_realData::timeFrac) = 1;
					p.idata(FHD_intData::species) = ispec;
					p.idata(FHD_intData::species_change) = ispec;
					
					// Normal Distribution
					Real diam = -1;
					while(diam<0) {
						diam = amrex::RandomNormal(dmean,dstd);
					}
					Real volp = pi_usr*pow(diam,3.0)/6.0;
					
					p.rdata(FHD_realData::radius) = diam*0.5;
					p.rdata(FHD_realData::mass) = rhomean*volp;
					p.rdata(FHD_realData::R) = k_B/p.rdata(FHD_realData::mass);

          p.pos(0) = prob_lo[0] + amrex::Random()*(prob_hi[0]-prob_lo[0]);
          p.pos(1) = prob_lo[1] + amrex::Random()*(prob_hi[1]-prob_lo[1]);
          p.pos(2) = prob_lo[2] + amrex::Random()*(prob_hi[2]-prob_lo[2]);
          
          Real stdev = sqrt(T_init[0]*p.rdata(FHD_realData::R));
          vpart[0] = stdev*amrex::RandomNormal(0.,1.);
					vpart[1] = stdev*amrex::RandomNormal(0.,1.);
					vpart[2] = stdev*amrex::RandomNormal(0.,1.);

					p.rdata(FHD_realData::velx) = vpart[0];
					p.rdata(FHD_realData::vely) = vpart[1];
					p.rdata(FHD_realData::velz) = vpart[2];
					spd = sqrt(pow(vpart[0],2)+pow(vpart[1],2)+pow(vpart[2],2));
					csx = pi_usr*pow(diam,2.0);
				  
					if(csx>csxmax){ csxmax = csx; }
					if(spd>spdmax){ spdmax = spd; }
					// For calculating timstep from Courant number
					if(vpart[0]>umax) { umax=vpart[0]; }
					if(vpart[1]>vmax) { vmax=vpart[1]; }
					if(vpart[2]>wmax) { wmax=vpart[2]; }

					particle_tile.push_back(p);
					
				}
			}
		}
	}
	// Set guess of max relative velocity
	Real cvrmax = csxmax*spdmax;
	ParallelDescriptor::Bcast(&cvrmax,1,ParallelDescriptor::IOProcessorNumber());
	mfvrmax.setVal(cvrmax);

	if(ParallelDescriptor::MyProc() == 0) {
		dt  = umax*n_cells[0]/(prob_hi[0]-prob_lo[0]);
		dt += vmax*n_cells[1]/(prob_hi[1]-prob_lo[1]);
		dt += wmax*n_cells[2]/(prob_hi[2]-prob_lo[2]);
		dt  = 0.2/dt; // Courant number of 0.2
	}
	ParallelDescriptor::Bcast(&dt,1,ParallelDescriptor::IOProcessorNumber());
	amrex::Print() << "My dt " << dt << "\n";

	Redistribute();
	SortParticles();
}

void FhdParticleContainer::ReInitParticles() {
	const int lev = 0;
	const Geometry& geom = Geom(lev);

	int pcount = 0;
	bool proc0_enter = true;

	std::array<Real, 3> vpart = {0.,0.,0.};
	Real u,v,w;
	Real spd, spdmax = 0.;
	Real csx, csxmax = 0.;

	for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi) {
		const Box& tile_box = mfi.tilebox();
		const RealBox tile_realbox{tile_box, geom.CellSize(), geom.ProbLo()};
		const int grid_id = mfi.index();
		const int tile_id = mfi.LocalTileIndex();
		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();
		const long np = particles.numParticles();

		for (int i=0; i<np; ++i) {
			ParticleType & part = particles[i];
			part.id() = ParticleType::NextID();
			part.cpu() = ParallelDescriptor::MyProc();
			part.idata(FHD_intData::sorted) = -1;
			part.idata(FHD_intData::i) = -100;
			part.idata(FHD_intData::j) = -100;
			part.idata(FHD_intData::k) = -100;

			part.rdata(FHD_realData::timeFrac) = 1;
			Real u = part.rdata(FHD_realData::velx);
			Real v = part.rdata(FHD_realData::vely);
			Real w = part.rdata(FHD_realData::velz);
			spd = sqrt(pow(u,2.0)+pow(v,2.0)+pow(w,2.0));
			Real diam = 2.0*part.rdata(FHD_realData::radius);
			csx = pi_usr*pow(diam,2.0);
			if(csx>csxmax) { csxmax=csx; }
			if(spd>spdmax) { spdmax=spd; }
		}
	}
	Real cvrmax;
	cvrmax = spdmax*csxmax;
	ParallelDescriptor::Bcast(&cvrmax,1,ParallelDescriptor::IOProcessorNumber());
	mfvrmax.setVal(cvrmax);

	Redistribute();
	SortParticles();
}
