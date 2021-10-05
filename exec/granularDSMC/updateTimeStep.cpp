#include "INS_functions.H"
#include "common_functions.H"
#include "species.H"
#include "DsmcParticleContainer.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;
using namespace std;

void FhdParticleContainer::updateTimeStep(const Geometry& geom, Real& dt) {
	BL_PROFILE_VAR("updateTimeStep()",writePlotFile);

	const int lev = 0;
	const Real* dx = Geom(lev).CellSize();

	Real dtmin;
	for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {
		const int grid_id = pti.index();
		const int tile_id = pti.LocalTileIndex();
		const Box& tile_box  = pti.tilebox();

    auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
    auto& particles = particle_tile.GetArrayOfStructs();

		IntVect smallEnd = tile_box.smallEnd();
		IntVect bigEnd = tile_box.bigEnd();

		dtmin = -1;
		for (int i = smallEnd[0]; i <= bigEnd[0]; i++) {
		for (int j = smallEnd[1]; j <= bigEnd[1]; j++) {
		for (int k = smallEnd[2]; k <= bigEnd[2]; k++) {
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);
			for (int ispec=0; ispec<nspecies; ispec++) {
				const long np = m_cell_vectors[ispec][grid_id][imap].size();
		    RealVect vbulk = {0.,0.,0.};
		    
		    for (int m=0; m<np; m++) {
		      int pind = m_cell_vectors[ispec][grid_id][imap][m];
		      ParticleType ptemp = particles[pind];
		      ParticleType & p = ptemp;

		      vbulk[0] += p.rdata(FHD_realData::velx);
		      vbulk[1] += p.rdata(FHD_realData::vely);
		      vbulk[2] += p.rdata(FHD_realData::velz);
		    }
		    vbulk[0] /= np;
		    vbulk[1] /= np;
		    vbulk[2] /= np;

				RealVect vc = {0.,0.,0.};
				Real T = 0.;
		    for (int m=0; m<np; m++) {
		      int pind = m_cell_vectors[ispec][grid_id][imap][m];
		      ParticleType ptemp = particles[pind];
		      ParticleType & p = ptemp;
		      Real u = p.rdata(FHD_realData::velx);
		      Real v = p.rdata(FHD_realData::vely);
		      Real w = p.rdata(FHD_realData::velz);

		      vc[0] = u - vbulk[0];
		      vc[1] = v - vbulk[1];
		      vc[2] = w - vbulk[2];

		      Real spdsq = vc[0]*vc[0]+vc[1]*vc[1]+vc[2]*vc[2];
					T += spdsq;
		    }
		    Real mass = properties[ispec].mass;
		    T = T*mass/np;
				T *= (2.0/3.0);

				Real vmean = sqrt(T);
				Real olam = sqrt(2.0)*np*pi_usr*
					pow(diameter[0],2.0)*particle_neff*ocollisionCellVol;
				Real lam = 1.0/olam;
				Real dt_cell = lam/vmean;
				if(dtmin<=0)
				{
					dtmin = dt_cell;
				}
				else
				{
					dtmin = std::min(dtmin,dt_cell);
				}
			}
		}
		}
		}
	}

	dtmin = dtmin/5.0;
	dtmin = std::min(dt,dtmin);
	ParallelDescriptor::ReduceRealMin(dtmin);
	if(dtmin<dt)
	{
		dt = dtmin;
		// for geometry
		fixed_dt = dt;
		ParallelDescriptor::Bcast(&fixed_dt,1,ParallelDescriptor::IOProcessorNumber());
		amrex::Print() << "Time step updated to " << dt << "\n";
	}
}
