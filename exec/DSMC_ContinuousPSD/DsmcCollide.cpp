#include "DsmcParticleContainer.H"
#include <math.h>

void FhdParticleContainer::InitCollisionCells() {
	BL_PROFILE_VAR("InitCollisionCells()",InitCollisionCells);
	interproperties[0].alpha = 1.0; //alpha_pp[cnt];

	const int lev = 0;
	for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {
		const int grid_id = pti.index();
		const int tile_id = pti.LocalTileIndex();
		const Box& tile_box  = pti.tilebox();
		
		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();
		
		// Convert MultiFabs -> arrays
		const Array4<Real> & arrselect = mfselect.array(pti);
		
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			arrselect(i,j,k,0) = 0.0;
			
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);
		});
	}
}

// Compute selections here
void FhdParticleContainer::CalcSelections(Real dt) {
	BL_PROFILE_VAR("CalcSelections()",CalcSelections);
	int lev = 0;
	mfselect.setVal(0.0);
	for(MFIter mfi(mfvrmax); mfi.isValid(); ++mfi) {
		const Box& tile_box  = mfi.tilebox();
		const int grid_id = mfi.index();
		const int tile_id = mfi.LocalTileIndex();
		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();

		// Convert MultiFabs -> arrays
		const Array4<Real> & arrvrmax = mfvrmax.array(mfi);
		const Array4<Real> & arrselect = mfselect.array(mfi);
		// anything defined outside of parallelfor is read-only
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {			
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);		

			const long np = m_cell_vectors[0][grid_id][imap].size();
			const Real vrmax = arrvrmax(i,j,k,0); // includes cross section
			Real NSel = 2.0*particle_neff*np*(np-1)*vrmax*ocollisionCellVol*dt;
			arrselect(i,j,k,0) = std::floor(NSel + amrex::Random());
		});
	}
}

void FhdParticleContainer::CollideParticles(Real dt) {
	BL_PROFILE_VAR("CollideParticles()",CollideParticles);
	int lev = 0;
	for(MFIter mfi(mfvrmax); mfi.isValid(); ++mfi) {
		const Box& tile_box  = mfi.tilebox();
		const int grid_id = mfi.index();
		const int tile_id = mfi.LocalTileIndex();
		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();

		// Convert MultiFabs -> arrays
		const Array4<Real> & arrvrmax = mfvrmax.array(mfi);
		const Array4<Real> & arrselect = mfselect.array(mfi);

		const long np = particles.numParticles();
		
		IntVect smallEnd = tile_box.smallEnd();
		IntVect bigEnd = tile_box.bigEnd();
		
		for (int i = smallEnd[0]; i <= bigEnd[0]; i++) {
		for (int j = smallEnd[1]; j <= bigEnd[1]; j++) {
		for (int k = smallEnd[2]; k <= bigEnd[2]; k++) {
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);

			Real NSel[nspecies*nspecies], totalSel, selrun;
			int pindxi, pindxj; // index of randomly sampled particles
			
			RealVect eij, vreij;
			Real eijmag, vreijmag;
			RealVect vi, vj, vij;
			Real massi, massj, massij;
			Real radj, radi, dij;
			Real csx, vrmag;
			Real cvrmag, cvrmax, cvreijmag;
			
			totalSel = (int)arrselect(i,j,k,0);;
			const long np = m_cell_vectors[0][grid_id][imap].size();
			cvrmax = arrvrmax(i,j,k,0); // includes cross section
			
			for (int isel = 0; isel<totalSel; isel++) {
				pindxi = floor(amrex::Random()*np);
				pindxj = floor(amrex::Random()*np);
				pindxi = m_cell_vectors[0][grid_id][imap][pindxi];
				pindxj = m_cell_vectors[0][grid_id][imap][pindxj];
				ParticleType & parti = particles[pindxi];
				ParticleType & partj = particles[pindxj];
				massi = parti.rdata(FHD_realData::mass);
				massj = partj.rdata(FHD_realData::mass);
				radi  = parti.rdata(FHD_realData::radius);
				radj  = partj.rdata(FHD_realData::radius);

				massij = massi+massj; dij = radi+radj;

				vi[0] = parti.rdata(FHD_realData::velx);
				vi[1] = parti.rdata(FHD_realData::vely);
				vi[2] = parti.rdata(FHD_realData::velz);

				vj[0] = partj.rdata(FHD_realData::velx);
				vj[1] = partj.rdata(FHD_realData::vely);
				vj[2] = partj.rdata(FHD_realData::velz);

				vrmag = sqrt(pow(vij[0],2)+pow(vij[1],2)+pow(vij[2],2));
				csx = pi_usr*pow(dij,2.0);
				cvrmag = csx*vrmag;
				if(cvrmag>cvrmax) {
					cvrmax = cvrmag;
					arrvrmax(i,j,k,0) = cvrmax;
				}

				vij[0] = vi[0]-vj[0]; vij[1] = vi[1]-vj[1]; vij[2] = vi[2]-vj[2];
				Real eijmag = 0;
				while(eijmag < 1.0e-12) {
					eij[0] = amrex::RandomNormal(0.,1.);
					eij[1] = amrex::RandomNormal(0.,1.);
					eij[2] = amrex::RandomNormal(0.,1.);
					eijmag = pow(eij[0],2)+pow(eij[1],2)+pow(eij[2],2);
					eijmag = pow(eijmag,0.5);
				}
				for(int idim=0; idim<3; idim++) { eij[idim] /= eijmag; }
				vreijmag = vij[0]*eij[0]+vij[1]*eij[1]+vij[2]*eij[2];
				cvreijmag = vreijmag*csx;

				if(amrex::Math::abs(cvreijmag)>cvrmax*amrex::Random()) {

					vreijmag = vreijmag*(1.0+interproperties[0].alpha)/massij;
					vreij[0] = vreijmag*eij[0];
					vreij[1] = vreijmag*eij[1];
					vreij[2] = vreijmag*eij[2];
							
					// Update velocities
					parti.rdata(FHD_realData::velx) = vi[0] - vreij[0]*massj;
					parti.rdata(FHD_realData::vely) = vi[1] - vreij[1]*massj;
					parti.rdata(FHD_realData::velz) = vi[2] - vreij[2]*massj;
					partj.rdata(FHD_realData::velx) = vj[0] + vreij[0]*massi;
					partj.rdata(FHD_realData::vely) = vj[1] + vreij[1]*massi;
					partj.rdata(FHD_realData::velz) = vj[2] + vreij[2]*massi;
				}		
			}
		}
		}
		}
	}
}
