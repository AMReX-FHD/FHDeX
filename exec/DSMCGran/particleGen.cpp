#include "INS_functions.H"
#include "common_functions.H"
#include "DsmcParticleContainer.H"
#include <sstream>
#include <string>
#include <fstream>

using namespace std;

void FhdParticleContainer::InitParticles() {
	const int lev = 0;
	const Geometry& geom = Geom(lev);

	int pcount = 0;

	bool proc0_enter = true;
    
	// Search for max relative speed
	// ... estimate as double the mag of a single particle speed
	std::array<Real, 3> vpart = {0., 0., 0.};
	Real spdmax = 0.;
	Real spd;
	Real stdev;
	Real u,v,w;
	
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
			std::ifstream particleFile("particles.dat");
			Real vmean;
			for(int i_spec=0; i_spec < nspecies; i_spec++) {
				u = 0; v = 0; w = 0;
				stdev = sqrt(T_init[i_spec]);
				for (int i_part=0; i_part<properties[i_spec].total;i_part++) {
					ParticleType p;
					p.id()  = ParticleType::NextID();
					p.cpu() = ParallelDescriptor::MyProc();
					p.idata(FHD_intData::sorted) = -1;
					p.idata(FHD_intData::species) = i_spec;

					vpart[0] = stdev*amrex::RandomNormal(0.,1.);
					vpart[1] = stdev*amrex::RandomNormal(0.,1.);
					vpart[2] = stdev*amrex::RandomNormal(0.,1.);
					
					if(particle_input > 0) {
               	particleFile >> p.pos(0);
               	particleFile >> p.pos(1);
               	particleFile >> p.pos(2);
               	// Overwrite velocities if provided
               	particleFile >> vpart[0];
               	particleFile >> vpart[1];
               	particleFile >> vpart[2];
               	particleFile >> p.idata(FHD_intData::species);
					} else if(particle_placement > 0) {
               	particleFile >> p.pos(0);                       
               	particleFile >> p.pos(1);
               	particleFile >> p.pos(2);
               } else {
               	p.pos(0) = prob_lo[0] + amrex::Random()*(prob_hi[0]-prob_lo[0]);
               	p.pos(1) = prob_lo[1] + amrex::Random()*(prob_hi[1]-prob_lo[1]);
               	p.pos(2) = prob_lo[2] + amrex::Random()*(prob_hi[2]-prob_lo[2]);
					}

					spd = sqrt(vpart[0]*vpart[0]+vpart[1]*vpart[1]+vpart[2]*vpart[2]);
					if(spd>spdmax){ spdmax = spd; }
					
					p.rdata(FHD_realData::velx) = vpart[0]; u += vpart[0];
					p.rdata(FHD_realData::vely) = vpart[1]; v += vpart[1];
					p.rdata(FHD_realData::velz) = vpart[2]; w += vpart[2];

					p.rdata(FHD_realData::boostx) = 0;
					p.rdata(FHD_realData::boosty) = 0;
					p.rdata(FHD_realData::boostz) = 0;
                    
					particle_tile.push_back(p);

					pcount++;
				}
				// Zero out the bulk velocities
				u = u/properties[i_spec].total;
				v = v/properties[i_spec].total;
				w = w/properties[i_spec].total;
				for (int i_part=0; i_part<properties[i_spec].total;i_part++) {
					ParticleType p = particles[i_part];
					p.rdata(FHD_realData::velx) = p.rdata(FHD_realData::velx) - u;
					p.rdata(FHD_realData::vely) = p.rdata(FHD_realData::vely) - v;
					p.rdata(FHD_realData::velz) = p.rdata(FHD_realData::velz) - w;
				}
			}
			particleFile.close();
		}
	}
	mfvrmax.setVal(spdmax);

	Redistribute();
	SortParticles();
}
