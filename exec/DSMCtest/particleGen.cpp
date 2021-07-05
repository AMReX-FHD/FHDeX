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
	Real spdmax = 0;
	Real umax = 0, vmax = 0, wmax = 0;
	Real spd;
	Real u[nspecies],v[nspecies],w[nspecies];
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
			// If full particle input provided
			if(particle_input>0) {
				std::ifstream particleFile("particles.dat");
				while(true) {
					ParticleType p;
					p.id()  = ParticleType::NextID();
					p.cpu() = ParallelDescriptor::MyProc();
					p.idata(FHD_intData::sorted) = -1;
					particleFile >> p.pos(0);
					particleFile >> p.pos(1);
          particleFile >> p.pos(2);
          particleFile >> vpart[0];
          particleFile >> vpart[1];
          particleFile >> vpart[2];
          p.rdata(FHD_realData::velx) = vpart[0];
          p.rdata(FHD_realData::vely) = vpart[1];
          p.rdata(FHD_realData::velz) = vpart[2];
          
          spd = sqrt(pow(vpart[0],2)+pow(vpart[1],2)+pow(vpart[2],2));
          if(spd>spdmax){ spdmax=spd; }
					// For calculating timstep from Courant number
					if(vpart[0]>umax) { umax=vpart[0]; }
					if(vpart[1]>vmax) { vmax=vpart[1]; }
					if(vpart[2]>wmax) { wmax=vpart[2]; }
          
          particleFile >> p.idata(FHD_intData::species);
          int ispec = p.idata(FHD_intData::species);
          p.rdata(FHD_realData::R) = k_B/properties[ispec].mass;
          p.rdata(FHD_realData::boostx) = 0;
					p.rdata(FHD_realData::boosty) = 0;
					p.rdata(FHD_realData::boostz) = 0;
          if( particleFile.eof() ) break;
				}
				particleFile.close();
			} else {
				// Initialize to bulk velocities
				for(int i_spec=0; i_spec < nspecies; i_spec++) {
					u[i_spec] = 0.0; v[i_spec] = 0.0; w[i_spec] = 0.0;
				}
				for(int i_spec=0; i_spec < nspecies; i_spec++) {
					// Standard deviation of velocity at temperature T_init
					Real R     = k_B/properties[i_spec].mass;
					Real stdev = sqrt(T_init[i_spec]*R);
					for (int i_part=0; i_part<properties[i_spec].total;i_part++) {
						ParticleType p;
						p.id()  = ParticleType::NextID();
						p.cpu() = ParallelDescriptor::MyProc();
						p.idata(FHD_intData::sorted) = -1;
						p.idata(FHD_intData::species) = i_spec;
						p.rdata(FHD_realData::R) = R;	
            //p.pos(0) = prob_lo[0] + amrex::Random()*(prob_hi[0]-prob_lo[0]);
            p.pos(0) = (prob_lo[0] + properties[i_spec].radius) + amrex::Random()*(prob_hi[0]-(prob_lo[0] + properties[i_spec].radius));
            // p.pos(1) = prob_lo[1] + amrex::Random()*(prob_hi[1]-prob_lo[1]);
            p.pos(1) = (prob_lo[1] + properties[i_spec].radius) + amrex::Random()*(prob_hi[1]-(prob_lo[1] + properties[i_spec].radius));
            // p.pos(2) = prob_lo[2] + amrex::Random()*(prob_hi[2]-prob_lo[2]);
            p.pos(2) = (prob_lo[2] + properties[i_spec].radius) + amrex::Random()*(prob_hi[2]-(prob_lo[2] + properties[i_spec].radius));
						vpart[0] = stdev*amrex::RandomNormal(0.,1.);
						vpart[1] = stdev*amrex::RandomNormal(0.,1.);
						vpart[2] = stdev*amrex::RandomNormal(0.,1.);

						p.rdata(FHD_realData::velx) = vpart[0]; u[p.idata(FHD_intData::species)] += vpart[0];
						p.rdata(FHD_realData::vely) = vpart[1]; v[p.idata(FHD_intData::species)] += vpart[1];
						p.rdata(FHD_realData::velz) = vpart[2]; w[p.idata(FHD_intData::species)] += vpart[2];
						spd = sqrt(pow(vpart[0],2)+pow(vpart[1],2)+pow(vpart[2],2));
					
						if(spd>spdmax){ spdmax=spd; }
						// For calculating timstep from Courant number
						if(vpart[0]>umax) { umax=vpart[0]; }
						if(vpart[1]>vmax) { vmax=vpart[1]; }
						if(vpart[2]>wmax) { wmax=vpart[2]; }

						p.rdata(FHD_realData::boostx) = 0;
						p.rdata(FHD_realData::boosty) = 0;
						p.rdata(FHD_realData::boostz) = 0;
						particle_tile.push_back(p);
						pcount++;
					}
				}
			}
			
			int nstart = 0;
			for(int i_spec=0; i_spec < nspecies; i_spec++) {			
				// Zero out the bulk velocities
				u[i_spec] = u[i_spec]/properties[i_spec].total;
				v[i_spec] = v[i_spec]/properties[i_spec].total;
				w[i_spec] = w[i_spec]/properties[i_spec].total;
				for (int i_part=nstart; i_part<nstart+properties[i_spec].total;i_part++) {
					ParticleType & p = particles[i_part];
					p.rdata(FHD_realData::velx) = p.rdata(FHD_realData::velx) - u[i_spec];
					p.rdata(FHD_realData::vely) = p.rdata(FHD_realData::vely) - v[i_spec];
					p.rdata(FHD_realData::velz) = p.rdata(FHD_realData::velz) - w[i_spec];
				}
				nstart += properties[i_spec].total;
			}
		}
	}
	// Set guess of max relative velocity
	ParallelDescriptor::Bcast(&spdmax,1,ParallelDescriptor::IOProcessorNumber());
	mfvrmax.setVal(spdmax);

	// Calculate global timestep
	// Assume IO processor is 0
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
