#include "DsmcParticleContainer.H"
using namespace std;
void FhdParticleContainer::EvaluateStatsParticles(MultiFab& mfPartInst,
						MultiFab& mfPartMeans) {
	BL_PROFILE_VAR("EvaluateStatsParticles()",EvaluateStats);
	const Real osteps = 1.0/steps;
	const int stepsMinusOne = steps-1;

	// TODO: Add Heat Fluxes
	const int lev = 0;
	mfPartInst.setVal(0.);
	for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {
		const int grid_id = pti.index();
		const int tile_id = pti.LocalTileIndex();
		const Box& tile_box  = pti.tilebox();
		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();
		IntVect smallEnd = tile_box.smallEnd();
		IntVect bigEnd = tile_box.bigEnd();

		// Rader, Gallis (2006) Conduction DSMC
		/*
			Particle Vars:
			0  - S   = Running total number of particles (samples)
			1  - rho
			2  - Jx
			3  - Jy
			4  - Jz
			5  - K
			6  - c
			7  - (c^2)*u
			8  - (c^2)*v
			9  - (c^2)*w
			10 - T
			11 - qx
			12 - qy
			13 - qz
		*/

		int npart  = 9*(nspecies+1);
		Array4<Real> partInst    = mfPartInst[pti].array();
		Array4<Real> partMeans   = mfPartMeans[pti].array();

		//////////////////////////////////////
		// Primitve and Conserved Instantaneous Values
		//////////////////////////////////////

		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);
			int ipart = 9;
			Real cv  = 0.;

			for (int l=0; l<nspecies; l++) {
				const long np_spec = m_cell_vectors[l][grid_id][imap].size();
				Real mass = properties[l].mass*properties[l].Neff;
				Real moV  = properties[l].mass*ocollisionCellVol;
				partInst(i,j,k,ipart+0)   = np_spec;
				partInst(i,j,k,0)        += np_spec;
				partMeans(i,j,k,ipart+0) += np_spec;
				partMeans(i,j,k,0)       += np_spec;

				partInst(i,j,k,ipart+1)   = np_spec*moV;
				partInst(i,j,k,1)        += np_spec*moV;
				partMeans(i,j,k,ipart+1) += np_spec*moV;
				partMeans(i,j,k,1)       += np_spec*moV;

				// Read particle data
				for (int m=0; m<np_spec; m++) {
					int pind = m_cell_vectors[l][grid_id][imap][m];
					ParticleType ptemp = particles[pind];
					ParticleType & p = ptemp;
					// ParticleType & p = particles[pind];
					Real u = p.rdata(FHD_realData::velx);
					Real v = p.rdata(FHD_realData::vely);
					Real w = p.rdata(FHD_realData::velz);

					// Jx, Jy, Jz
					partInst(i,j,k,ipart+2) += (mass*u); partInst(i,j,k,2) += (mass*u);
					partInst(i,j,k,ipart+3) += (mass*v); partInst(i,j,k,3) += (mass*v);
					partInst(i,j,k,ipart+4) += (mass*w); partInst(i,j,k,4) += (mass*w);

					partMeans(i,j,k,ipart+2) += (mass*u); partMeans(i,j,k,2) += (mass*u);
					partMeans(i,j,k,ipart+3) += (mass*v); partMeans(i,j,k,3) += (mass*v);
					partMeans(i,j,k,ipart+4) += (mass*w); partMeans(i,j,k,4) += (mass*w);
					
					Real spdsq = pow(u,2)+pow(v,2)+pow(w,2);
					Real spd = pow(spdsq,0.5);
					
					// K
					partInst(i,j,k,ipart+5)  += (mass*spd); partInst(i,j,k,5) += (mass*spd);
					partMeans(i,j,k,ipart+5) += (mass*spd); partMeans(i,j,k,5) += (mass*spd);
					
					// qx, qy, qz
					partInst(i,j,k,ipart+7)  += (spdsq*u); partInst(i,j,k,7) += partInst(i,j,k,ipart+7);
					partInst(i,j,k,ipart+8)  += (spdsq*v); partInst(i,j,k,8) += partInst(i,j,k,ipart+8);
					partInst(i,j,k,ipart+9)  += (spdsq*w); partInst(i,j,k,9) += partInst(i,j,k,ipart+9);
					
					partMeans(i,j,k,ipart+7) += partInst(i,j,k,ipart+7); partMeans(i,j,k,7) += partInst(i,j,k,ipart+7);
					partMeans(i,j,k,ipart+8) += partInst(i,j,k,ipart+8); partMeans(i,j,k,8) += partInst(i,j,k,ipart+8);
					partMeans(i,j,k,ipart+9) += partInst(i,j,k,ipart+9); partMeans(i,j,k,9) += partInst(i,j,k,ipart+9);
				}

				// T
				partInst(i,j,k,ipart+6)  += ((partInst(i,j,k,ipart+5)
					- pow(partInst(i,j,k,ipart+2),2)/mass - pow(partInst(i,j,k,ipart+3),2)/mass
					- pow(partInst(i,j,k,ipart+4),2)/mass)/(3.0*k_B));
				partInst(i,j,k,6)        += partInst(i,j,k,ipart+6);
				partMeans(i,j,k,ipart+4) += partInst(i,j,k,ipart+6);
				partMeans(i,j,k,4)       += partInst(i,j,k,ipart+6);

				// qx, qy, qz
				partInst(i,j,k,ip
				
				ipart += 9;
			}
		});
	}
}

