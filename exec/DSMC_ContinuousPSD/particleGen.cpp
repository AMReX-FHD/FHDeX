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
	Real cvrmax = 0, spd, csx;
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
				for(int i_spec=0; i_spec < nspecies; i_spec++) {
					// Standard deviation of velocity at temperature T_init
					Real R     = k_B/properties[i_spec].mass;
					Real stdev = sqrt(T_init[i_spec]*R);
					for (int i_part=0; i_part<properties[i_spec].total;i_part++) {
						ParticleType p;
						p.id()  = ParticleType::NextID();
						p.cpu() = ParallelDescriptor::MyProc();
						p.idata(FHD_intData::sorted) = -1;
						p.idata(FHD_intData::species) = 0;
						p.rdata(FHD_realData::mass) = properties[i_spec].mass;
						p.rdata(FHD_realData::R) = k_B/p.rdata(FHD_realData::mass);
						p.rdata(FHD_realData::radius) = properties[i_spec].diameter*0.5;
						
            p.pos(0) = prob_lo[0] + amrex::Random()*(prob_hi[0]-prob_lo[0]);
            p.pos(1) = prob_lo[1] + amrex::Random()*(prob_hi[1]-prob_lo[1]);
            p.pos(2) = prob_lo[2] + amrex::Random()*(prob_hi[2]-prob_lo[2]);
            vpart[0] = stdev*amrex::RandomNormal(0.,1.);
						vpart[1] = stdev*amrex::RandomNormal(0.,1.);
						vpart[2] = stdev*amrex::RandomNormal(0.,1.);

						p.rdata(FHD_realData::velx) = vpart[0];
						p.rdata(FHD_realData::vely) = vpart[1];
						p.rdata(FHD_realData::velz) = vpart[2];
						spd = sqrt(pow(vpart[0],2)+pow(vpart[1],2)+pow(vpart[2],2));
						csx = pi_usr*pow(2.0*p.rdata(FHD_realData::radius),2.0);
					  
						if((csx*spd)>cvrmax){ cvrmax=csx*spd; }
						// For calculating timstep from Courant number
						if(vpart[0]>umax) { umax=vpart[0]; }
						if(vpart[1]>vmax) { vmax=vpart[1]; }
						if(vpart[2]>wmax) { wmax=vpart[2]; }

						particle_tile.push_back(p);
					}
				}
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
