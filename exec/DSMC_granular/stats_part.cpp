#include "DsmcParticleContainer.H"
using namespace std;
void FhdParticleContainer::EvaluateStatsPart(MultiFab& mfvmom)
{
	BL_PROFILE_VAR("EvaluateStatsPart()",EvaluateStats);
	const int lev = 0; 
  for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {
      const int grid_id = pti.index();
      const int tile_id = pti.LocalTileIndex();
      const Box& tile_box  = pti.tilebox();
      auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
      auto& particles = particle_tile.GetArrayOfStructs();
      IntVect smallEnd = tile_box.smallEnd();
      IntVect bigEnd = tile_box.bigEnd();
    /*
      Velocity Moments:
      0  - Nsample
      1  - Pk_(1,1)
      2  - Pk_(1,2)
      3  - Pk_(1,3)
      4  - Pk_(2,2)
      5  - Pk_(2,3)
      6  - Pk_(3,3)
      7  - qx
      8  - qy
      9  - qz
    */

    Array4<Real> vmom = mfvmom[pti].array();

    amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept 
    {
        const IntVect& iv = {i,j,k};
        long imap = tile_box.index(iv);

				// Only works with 1 specie right now
        //for (int l=0; l<nspecies; l++) {
        const long np_spec = m_cell_vectors[0][grid_id][imap].size();
				
				RealVect vbulk = {0.,0.,0.};
        for (int m=0; m<np_spec; m++) {
          int pind = m_cell_vectors[0][grid_id][imap][m];
          ParticleType ptemp = particles[pind];
          ParticleType & p = ptemp;
          Real u = p.rdata(FHD_realData::velx);
          Real v = p.rdata(FHD_realData::vely);
          Real w = p.rdata(FHD_realData::velz);

          vbulk[0] += u;
          vbulk[1] += v;
          vbulk[2] += w;
        }
        vbulk[0] /= np_spec;
        vbulk[1] /= np_spec;
        vbulk[2] /= np_spec;

				RealVect c = {0.,0.,0.};
				long Nsample = vmom(i,j,k,0);
        for(int m=1; m<=9; m++)
        {
					vmom(i,j,k,m) *= Nsample;
				}					

        for (int m=0; m<np_spec; m++) {
          int pind = m_cell_vectors[0][grid_id][imap][m];
          ParticleType ptemp = particles[pind];
          ParticleType & p = ptemp;
          Real u = p.rdata(FHD_realData::velx);
          Real v = p.rdata(FHD_realData::vely);
          Real w = p.rdata(FHD_realData::velz);

          c[0] = u - vbulk[0];
          c[1] = v - vbulk[1];
          c[2] = w - vbulk[2];

          vmom(i,j,k,1) = vmom(i,j,k,1)+(c[0]*c[0]);
          vmom(i,j,k,2) = vmom(i,j,k,2)+(c[0]*c[1]);
          vmom(i,j,k,3) = vmom(i,j,k,3)+(c[0]*c[2]);
          vmom(i,j,k,4) = vmom(i,j,k,4)+(c[1]*c[1]);
          vmom(i,j,k,5) = vmom(i,j,k,5)+(c[1]*c[2]);
          vmom(i,j,k,6) = vmom(i,j,k,6)+(c[2]*c[2]);
          
          Real spdsq = c[0]*c[0]+c[1]*c[1]+c[2]*c[2];
					spdsq = spdsq*0.5;
          vmom(i,j,k,7) = vmom(i,j,k,7)+spdsq*c[0];
          vmom(i,j,k,8) = vmom(i,j,k,8)+spdsq*c[1];
          vmom(i,j,k,9) = vmom(i,j,k,9)+spdsq*c[2];
        }

        vmom(i,j,k,0) += np_spec;
        Nsample = vmom(i,j,k,0);
        for(int m=1; m<=9; m++)
        {
					vmom(i,j,k,m) /= Nsample;
				}
				//Print() << "T: " << mass[0]*(vmom(i,j,k,1) + vmom(i,j,k,4) + vmom(i,j,k,6))/(3.0*k_B) << "\n"; 
        //}
    });
	}
}
