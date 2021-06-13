#include "DsmcParticleContainer.H"

// #include "particle_functions_K.H"
#include "paramplane_functions_K.H"
#include <math.h>
using namespace std;
// Same as FhdParticleContainer.H?
FhdParticleContainer::FhdParticleContainer(const Geometry & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba,
                              int ncells)
    : NeighborParticleContainer<FHD_realData::count, FHD_intData::count> (geom, dmap, ba, ncells)
{
	BL_PROFILE_VAR("FhdParticleContainer()",FhdParticleContainer);

	realParticles = 0;
	simParticles = 0;
	 
	totalCollisionCells = n_cells[0]*n_cells[1]*n_cells[2];
	domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*(prob_hi[2] - prob_lo[2]);

	collisionCellVol = domainVol/totalCollisionCells;
	ocollisionCellVol = 1/collisionCellVol;
	for(int i=0;i<nspecies;i++) {
		properties[i].mass = mass[i];
		properties[i].radius = diameter[i]/2.0;
		properties[i].partVol = pow(diameter[i],3)*pi_usr/6;
      properties[i].part2cellVol = properties[i].partVol*ocollisionCellVol;
		properties[i].Neff = particle_neff; // assume F_n is same for each
      
      // Overwrite particle_count
      if( particle_count[i] >= 0 ) {
      	
      	properties[i].total = particle_count[i];
      	properties[i].n0 = particle_neff*properties[i].total/domainVol;
      	
      	amrex::Print() <<  "Species "<< i << " count " << properties[i].total << "\n";
      	amrex::Print() <<  "Species "<< i << " n0 " << properties[i].total << "\n";
      	
      } else if( phi_domain[i] >= 0 ) {
      
			properties[i].total = (int)amrex::Math::ceil(
      		(phi_domain[i]*domainVol)/(properties[i].partVol*properties[i].Neff) );
      	properties[i].n0 = properties[i].total/domainVol;

      	amrex::Print() <<  "Species "<< i << " count " << properties[i].total << "\n";
      	amrex::Print() <<  "Species "<< i << " n0 " << properties[i].total << "\n";

      } else {
      
      	properties[i].total = (int)amrex::Math::ceil(particle_n0[i]*domainVol/particle_neff);
      	properties[i].n0 = properties[i].total/domainVol;
      
      	amrex::Print() <<  "Species "<< i << " count " << properties[i].total << "\n";
      	amrex::Print() <<  "Species "<< i << " n0 " << properties[i].total << "\n";
      }
      
      realParticles = realParticles + properties[i].total*particle_neff;
      simParticles = simParticles + properties[i].total;
      amrex::Print() << "Particles per cell for species: " << i << " is " << properties[i].total/totalCollisionCells << "\n";
	}
   
   int indx, ij_spec;
   int cnt = 0;
	for(int i_spec=0;i_spec<nspecies;i_spec++) {
		for(int j_spec=0;j_spec<nspecies;j_spec++) {
			ij_spec = getSpeciesIndex(i_spec,j_spec);
    		interproperties[ij_spec].alpha = alpha_pp[cnt]; // need to find better way to input this
    		//amrex::Print() << "i: " << i_spec << " j: " << j_spec << " ij: " << ij_spec 
    		//	<< " alpha: " << interproperties[ij_spec].alpha << "\n";
    		interproperties[ij_spec].csx = pow(properties[i_spec].radius+properties[j_spec].radius,2)*pi_usr;
    		countedCollisions[ij_spec] = 0;
    		expectedCollisions[ij_spec] = 0;
    		NSel_spec[ij_spec] = 0;
    		cnt++;
		}
	}

	for (int d=0; d<AMREX_SPACEDIM; ++d) {
		domSize[d] = prob_hi[d] - prob_lo[d];
	}

	const int lev=0;
	for(MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi) {
		const Box& box = mfi.validbox();
		const int grid_id = mfi.index();
		for(int i=0;i<nspecies;i++){
			m_cell_vectors[i][grid_id].resize(box.numPts());
		}
	}

	//initTg = -1;
	//tTg = 0;
}

void FhdParticleContainer::MoveParticlesCPP(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount)
{
    BL_PROFILE_VAR("MoveParticlesCPP()", MoveParticlesCPP);

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();
    const Real* phi = Geom(lev).ProbHi();

    int        np_tile = 0 ,       np_proc = 0 ; // particle count
    Real    moves_tile = 0.,    moves_proc = 0.; // total moves
    Real  maxspeed_proc = 0.; // max speed
    Real  maxdist_proc = 0.; // max displacement (fraction of radius)

    Real adj = 0.99999;
    Real adjalt = 2.0*(1.0-0.99999);
    Real runtime, inttime;
    int intsurf, intside, push;

    Real maxspeed = 0;
    Real maxdist = 0;

    long moves = 0;
    int reDist = 0;
    

    for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {

        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const long np = particles.numParticles();
        
        // Update solid fraction here ( and at sort_particles() )
        const Array4<Real> & arrphi = mfphi.array(pti);

        Box bx  = pti.tilebox();
        IntVect myLo = bx.smallEnd();
        IntVect myHi = bx.bigEnd();

        np_proc += np;

        for (int i = 0; i < np; ++ i) {
              ParticleType & part = particles[i];

                if(part.id() == 332)
                {
                       //Print() << "stated moving particle " << i << "\n";
                }


              Real speed = 0;

              for (int d=0; d<AMREX_SPACEDIM; ++d)
              {                   
                  speed += part.rdata(FHD_realData::velx + d)*part.rdata(FHD_realData::velx + d);
              }

              if(speed > maxspeed)
              {
                  maxspeed = speed;
              }

              moves++;

              runtime = dt*part.rdata(FHD_realData::timeFrac);

              while(runtime > 0)
              {

                  find_inter_gpu(part, runtime, paramPlaneList, paramPlaneCount, &intsurf, &inttime, &intside, ZFILL(plo), ZFILL(phi));

                  for (int d=0; d<AMREX_SPACEDIM; ++d)
                  {
                      part.pos(d) += inttime * part.rdata(FHD_realData::velx + d)*adj;
                  }

                  runtime = runtime - inttime;

                  if(intsurf > 0)
                  {
                      const paramPlane& surf = paramPlaneList[intsurf-1];//find_inter indexes from 1 to maintain compatablity with fortran version
      
                      Real posAlt[3];

                      for (int d=0; d<AMREX_SPACEDIM; ++d)
                      {
                          posAlt[d] = inttime * part.rdata(FHD_realData::velx + d)*adjalt;
                      }

                      Real dummy = 1;
                       //Print() << "surf: " << intsurf-1 << "\n";
                      app_bc_gpu(&surf, part, intside, domSize, &push, &runtime, dummy);
                       //Print() << "rt: " << runtime << "\n";

                      if(push == 1)
                      {
                          for (int d=0; d<AMREX_SPACEDIM; ++d)
                          {
                              part.pos(d) += part.pos(d) + posAlt[d];
                          }
                      }
                  }

              }
                if(part.id() == 332)
                {
              //Print() << "finished move particle " << i << "\n";
}
              part.rdata(FHD_realData::timeFrac) = 1;


            int cell[3];
            cell[0] = (int)floor((part.pos(0)-plo[0])/dx[0]);
            cell[1] = (int)floor((part.pos(1)-plo[1])/dx[1]);
            cell[2] = (int)floor((part.pos(2)-plo[2])/dx[2]);

            //n << "current pos: " << part.pos(0) << ", " << part.pos(1) << ", " << part.pos(2) << "\n";

//            if((cell[0] < myLo[0]) || (cell[1] < myLo[1]) || (cell[2] < myLo[2]) || (cell[0] > myHi[0]) || (cell[1] > myHi[1]) || (cell[2] > myHi[2]))
//            {
//                reDist++;
//            }

            if((part.idata(FHD_intData::i) != cell[0]) || (part.idata(FHD_intData::j) != cell[1]) || (part.idata(FHD_intData::k) != cell[2]) || part.id() < 0)
            {
                //remove particle from old cell


                IntVect iv(part.idata(FHD_intData::i), part.idata(FHD_intData::j), part.idata(FHD_intData::k));
                long imap = tile_box.index(iv);
                int lastIndex = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].size() - 1;
                int lastPart = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap][lastIndex];
                int newIndex = part.idata(FHD_intData::sorted);

                m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap][newIndex] = lastPart;
                m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].pop_back();

                particles[lastPart].idata(FHD_intData::sorted) = newIndex;

                part.idata(FHD_intData::sorted) = -1;

                // reduce solid fraction
                int ispec = part.idata(FHD_intData::species);
                arrphi(iv[0],iv[1],iv[2],ispec) =
                	arrphi(iv[0],iv[1],iv[2],ispec) - properties[ispec].part2cellVol*properties[ispec].Neff;

                //Print() << "Removed\n";
            }

              //Print() << "finished sorting particle " << i << "\n";
        }

        maxspeed_proc = amrex::max(maxspeed_proc, maxspeed);
        maxdist_proc  = amrex::max(maxdist_proc, maxdist);
    }

    // gather statistics
    ParallelDescriptor::ReduceIntSum(np_proc);
    ParallelDescriptor::ReduceRealSum(moves_proc);
    ParallelDescriptor::ReduceRealMax(maxspeed_proc);
    ParallelDescriptor::ReduceRealMax(maxdist_proc);
    ParallelDescriptor::ReduceIntSum(reDist);

    // write out global diagnostics
    // if (ParallelDescriptor::IOProcessor()) {
    //    Print() << "I see " << np_proc << " particles\n";
    //    Print() << reDist << " particles to be redistributed.\n";
    //    Print() <<"Maximum observed speed: " << sqrt(maxspeed_proc) << "\n";
    //    Print() <<"Maximum observed displacement (fraction of radius): " << maxdist_proc << "\n";
    //}

    Redistribute();
    SortParticles();
}

void FhdParticleContainer::InitCollisionCells() {	
	const int lev = 0;
	for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {
		const int grid_id = pti.index();
		const int tile_id = pti.LocalTileIndex();
		const Box& tile_box  = pti.tilebox();
		
		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();
		
		// Convert MultiFabs -> arrays
		const Array4<Real> & arrvrmax = mfvrmax.array(pti);
		const Array4<Real> & arrphi = mfphi.array(pti);
		const Array4<Real> & arrselect = mfselect.array(pti);
		
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

			int ij_spec;
			for (int i_spec=0; i_spec<nspecies; i_spec++) {
				for (int j_spec = i_spec; j_spec < nspecies; j_spec++) {
					ij_spec = getSpeciesIndex(i_spec,j_spec);
					arrselect(i,j,k,ij_spec) = 0.0;
				}
			}
					
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);

			for (int i_spec=0; i_spec<nspecies; i_spec++) {
				arrphi(i,j,k,i_spec) = m_cell_vectors[i_spec][grid_id][imap].size()
					*properties[i_spec].part2cellVol*properties[i_spec].Neff;
			}
		});
	}
}


// Compute selections here
void FhdParticleContainer::CalcSelections(Real dt) {
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
		const Array4<Real> & arrphi = mfphi.array(mfi);
		const Array4<Real> & arrselect = mfselect.array(mfi);		
		// anything defined outside of parallelfor is read-only
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {       
			int ij_spec;
			long np_i, np_j;
			
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);		

			Real vrmax;
			Real NSel;
			Real phi1, phi2, chi0;
			Real crossSection;

			for (int i_spec = 0; i_spec < nspecies; i_spec++) {
				for (int j_spec = i_spec; j_spec < nspecies; j_spec++) {
					ij_spec = getSpeciesIndex(i_spec,j_spec);
					np_i = m_cell_vectors[i_spec][grid_id][imap].size();
					np_j = m_cell_vectors[j_spec][grid_id][imap].size();
					phi1 = arrphi(i,j,k,i_spec);
					phi2 = arrphi(i,j,k,j_spec);
					// comment out if expecting dilute
					chi0 = g0_Ma_Ahmadi(i_spec,j_spec, phi1, phi2);
					vrmax = arrvrmax(i,j,k,i_spec);
					crossSection = interproperties[ij_spec].csx;
					
					NSel = 2.0*particle_neff*np_i*np_j*crossSection*vrmax*ocollisionCellVol*chi0*dt;
					if(i_spec==j_spec) {NSel = NSel*0.5;}
					arrselect(i,j,k,ij_spec) = std::floor(NSel + amrex::Random());
					NSel_spec[ij_spec] += arrselect(i,j,k,ij_spec);
				}
			}
		});

	}

}

void FhdParticleContainer::CollideParticles(Real dt) {
	int lev = 0;
	//tTg = 0.0;
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
		// may be better if written with AMREX_FOR_1D
		//amrex::ParallelForRNG(tile_box,
		//	[=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept {
		
		IntVect smallEnd = tile_box.smallEnd();
		IntVect bigEnd = tile_box.bigEnd();
		
		for (int i = smallEnd[0]; i <= bigEnd[0]; i++) {
		for (int j = smallEnd[1]; j <= bigEnd[1]; j++) {
		for (int k = smallEnd[2]; k <= bigEnd[2]; k++) {
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);

			Real NSel[nspecies*nspecies], totalSel, selrun;
			int np[nspecies];
			int pindxi, pindxj; // index of randomly sampled particles
			int ij_spec;
			
			RealVect eij, vreij;
			Real phi, theta, eijmag;
			RealVect vi, vj, vij;
			RealVect vijpost, boost;
			Real oboostmag;
			Real massi, massj, massij;
			Real vrmag, vrmax, vreijmag;
			
			totalSel = 0;
			for (int i_spec = 0; i_spec<nspecies; i_spec++) {
				np[i_spec] = m_cell_vectors[i_spec][grid_id][imap].size();
				for (int j_spec = i_spec; j_spec < nspecies; j_spec++) {
					ij_spec = getSpeciesIndex(i_spec,j_spec);
					NSel[ij_spec] = (int)arrselect(i,j,k,ij_spec);
					totalSel += NSel[ij_spec];
				}
			}
			int speci, specj, specij;
			while (totalSel>0) {
				Real RR = amrex::Random();
				bool spec_select = false;
				selrun = 0;
				speci = -1; specj = -1; specij = -1;
				for(int i_spec=0;i_spec<nspecies;i_spec++) { 
					for(int j_spec=i_spec;j_spec<nspecies;j_spec++) { 
						int ij_spec = getSpeciesIndex(i_spec,j_spec);
						selrun += NSel[ij_spec];
						if(selrun/totalSel>RR && !spec_select) {
							spec_select = true;
							specij = ij_spec;
							speci = i_spec; specj = j_spec;
							NSel[ij_spec] -= 1;
						}
					}
				}
				totalSel--;
				massi = properties[speci].mass;
				massj = properties[specj].mass;
				massij = properties[speci].mass + properties[specj].mass;
				vrmax = arrvrmax(i,j,k,specij);
				pindxi = floor(amrex::Random()*np[speci]);
				pindxj = floor(amrex::Random()*np[specj]);
				pindxi = m_cell_vectors[speci][grid_id][imap][pindxi];
				pindxj = m_cell_vectors[specj][grid_id][imap][pindxj];
				ParticleType &	parti = particles[pindxi];
				ParticleType & partj = particles[pindxj];
						
				vi[0] = parti.rdata(FHD_realData::velx);
				vi[1] = parti.rdata(FHD_realData::vely);
				vi[2] = parti.rdata(FHD_realData::velz);

				vj[0] = partj.rdata(FHD_realData::velx);
				vj[1] = partj.rdata(FHD_realData::vely);
				vj[2] = partj.rdata(FHD_realData::velz);
						
				vij[0] = vi[0]-vj[0]; vij[1] = vi[1]-vj[1]; vij[2] = vi[2]-vj[2];
				vrmag = sqrt(pow(vij[0],2)+pow(vij[1],2)+pow(vij[2],2));
				if(vrmag>vrmax) {vrmax = vrmag*1.05; arrvrmax(i,j,k,ij_spec) = vrmax;}

				theta = 2.0*pi_usr*amrex::Random();
				phi = std::acos(2.0*amrex::Random()-1.0);
				eij[0] = std::sin(theta)*std::cos(phi);
				eij[1] = std::sin(theta)*std::sin(phi);
				eij[2] = std::cos(theta);
				vreijmag = vij[0]*eij[0]+vij[1]*eij[1]+vij[2]*eij[2];					
				if(amrex::Math::abs(vreijmag)>vrmax*amrex::Random()) {
					countedCollisions[specij] += 2;

					//virial.push_back(vreijmag);
					vreijmag = vreijmag*(1.0+interproperties[specij].alpha)/massij;
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
							
					// Boosted Velocities
					vijpost[0] = vi[0]-vj[0];
					vijpost[1] = vi[1]-vj[1];
					vijpost[2] = vi[2]-vj[2];
							
					boost[0] = vijpost[0]-vij[0];
					boost[1] = vijpost[1]-vij[1];
					boost[2] = vijpost[2]-vij[2];
					oboostmag = (properties[speci].radius+properties[specj].radius)/
						(sqrt(pow(boost[0],2)+pow(boost[1],2)+pow(boost[2],2))*dt);
					boost[0] = boost[0]*oboostmag;
					boost[1] = boost[1]*oboostmag;
					boost[2] = boost[2]*oboostmag;
							
					parti.rdata(FHD_realData::boostx) = boost[0];
					parti.rdata(FHD_realData::boosty) = boost[1];
					parti.rdata(FHD_realData::boostz) = boost[2];
					partj.rdata(FHD_realData::boostx) = -boost[0];
					partj.rdata(FHD_realData::boosty) = -boost[1];
					partj.rdata(FHD_realData::boostz) = -boost[2];
				}		
			}
		}
		}
		}
	}
}

void FhdParticleContainer::OutputParticles() {
	string tTgFile = "particles.dat";
	ofstream myfile;
	myfile.open(tTgFile);
	int lev = 0;
	for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {
		const int grid_id = pti.index();
		const int tile_id = pti.LocalTileIndex();
		const Box& tile_box  = pti.tilebox();

		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();
		const long np = particles.numParticles();

		for (int i = 0; i < np; ++i) {
			ParticleType & part = particles[i];
			myfile << fixed << setprecision(5) << part.pos(0) << "  " << part.pos(1) << "  " << part.pos(2) << " "
				<< part.rdata(FHD_realData::velx) << " " << part.rdata(FHD_realData::vely) << " "
				<< part.rdata(FHD_realData::velz) << " " << part.idata(FHD_intData::species) << "\n"; 
		}
	}
	myfile.close();	
}

void FhdParticleContainer::EvaluateStats(MultiFab& particleInstant, MultiFab& particleMeans,
                                         MultiFab& particleVars, const Real delt, int steps) {
	BL_PROFILE_VAR("EvaluateStats()",EvaluateStats);
    
	const int lev = 0;
    
	BoxArray ba = particleMeans.boxArray();
	long cellcount = ba.numPts();

	const Real* dx = Geom(lev).CellSize();
	const Real dxInv = 1.0/dx[0];
	const Real cellVolInv = 1.0/(dx[0]*dx[0]*dx[0]); // this is recorded in ocollisionCellVol

	const Real stepsInv = 1.0/steps;
	const int stepsMinusOne = steps-1;

	// zero instantaneous values
	particleInstant.setVal(0.);
   
	for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) {
		PairIndex index(pti.index(), pti.LocalTileIndex());
		const int np = this->GetParticles(lev)[index].numRealParticles();
		auto& plev = this->GetParticles(lev);
		auto& ptile = plev[index];
		auto& aos   = ptile.GetArrayOfStructs();
		const Box& tile_box  = pti.tilebox();
		ParticleType* particles = aos().dataPtr();

		GpuArray<int, 3> bx_lo = {tile_box.loVect()[0], tile_box.loVect()[1], tile_box.loVect()[2]};
		GpuArray<int, 3> bx_hi = {tile_box.hiVect()[0], tile_box.hiVect()[1], tile_box.hiVect()[2]};

		Array4<Real> part_inst = particleInstant[pti].array();
		Array4<Real> part_mean = particleMeans[pti].array();
		Array4<Real> part_var = particleVars[pti].array();

		AMREX_FOR_1D( np, ni,
		{
			ParticleType & part = particles[ni];

			int i = floor(part.pos(0)*dxInv);
			int j = floor(part.pos(1)*dxInv);
			int k = floor(part.pos(2)*dxInv);
            
			amrex::Gpu::Atomic::Add(&part_inst(i,j,k,0), 1.0);

			for(int l=0;l<nspecies;l++) {

			}

		});

		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			part_inst(i,j,k,1) = part_inst(i,j,k,1)*cellVolInv;         
		});

		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

			part_mean(i,j,k,1)  = (part_mean(i,j,k,1)*stepsMinusOne + part_inst(i,j,k,1))*stepsInv;

			for(int l=0;l<nspecies;l++)
			{

			}
             
		});

		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

            Real del1 = part_inst(i,j,k,1) - part_mean(i,j,k,1);

            part_var(i,j,k,1)  = (part_var(i,j,k,1)*stepsMinusOne + del1*del1)*stepsInv;            
		});
		
		// Granular Temperature
		/*
		const long grid_id = pti.index();
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
					
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);
			Real spdsq;
			ParticleType ptemp;
			long pindx;

			for (int i_spec=0; i_spec<nspecies; i_spec++) {
				// phi = np * particle volume * collision cell volume * Neff
				long npcell = m_cell_vectors[i_spec][grid_id][imap].size();
				for (int i_part=0; i_part<npcell; i_part++) {
					pindx = m_cell_vectors[i_spec][grid_id][imap][i_part];
					ptemp =  particles[pindx];
					ParticleType & p = ptemp;
					spdsq = sqrt(pow(p.rdata(FHD_realData::velx),2)+
						pow(p.rdata(FHD_realData::vely),2)+pow(p.rdata(FHD_realData::velz),2));
					tTg = tTg + properties[i_spec].mass*spdsq;
				}
			}
		
		});*/
		
	}
}



// Likely excessive
// Relates lower diagonal matrix indices to those of 1D array
// Problematic when nspecies = max_species (last index will cause seg fault)
int FhdParticleContainer::getSpeciesIndex(int species1, int species2) {
	if(species1<species2){
		return species2+(nspecies-1)*species1;
	} else {
		return species1+(nspecies-1)*species2;
	}
}

void FhdParticleContainer::Source(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount) {
    int lev = 0;
    bool proc0_enter = true;
    //Do this all on rank 0 for now
    for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi)
    {
        if(ParallelDescriptor::MyProc() == 0 && proc0_enter)
        {
            proc0_enter = false;

            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

            for(int i = 0; i< paramPlaneCount; i++)
            {
                if(paramPlaneList[i].sourceLeft == 1)
                {
                    for(int j = 0; j< nspecies; j++)
                    {
                        Real density = paramPlaneList[i].densityLeft[j];
                        Real temp = paramPlaneList[i].temperatureLeft;
                        Real area = paramPlaneList[i].area;

                        Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
                        Real fluxVar = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;

                        Real totalFlux = dt*fluxMean + sqrt(dt*fluxVar)*amrex::RandomNormal(0.,1.);

                        //Print() << "Flux mean " << dt*fluxMean << ", flux sd " << sqrt(dt*fluxVar) << "\n";

                        int totalFluxInt =  (int)floor(totalFlux);
                        Real totalFluxLeftOver = totalFlux - totalFluxInt;

                        if(amrex::Random() < totalFluxLeftOver)
                        {
                            totalFluxInt++;
                        }

                        Print() << "Surface " << i << " generating " << totalFluxInt << " of species " << j << "\n";

                        for(int k=0;k<totalFluxInt;k++)
                        {
                            Real uCoord = amrex::Random()*paramPlaneList[i].uTop;
                            Real vCoord = amrex::Random()*paramPlaneList[i].vTop;

                            ParticleType p;
                            p.id() = ParticleType::NextID();

                            p.cpu() = ParallelDescriptor::MyProc();
                            p.idata(FHD_intData::sorted) = -1;

                            p.idata(FHD_intData::species) = j;

                            p.pos(0) = paramPlaneList[i].x0 + paramPlaneList[i].ux*uCoord + paramPlaneList[i].vx*vCoord;
                            p.pos(1) = paramPlaneList[i].y0 + paramPlaneList[i].uy*uCoord + paramPlaneList[i].vy*vCoord;
                            p.pos(2) = paramPlaneList[i].z0 + paramPlaneList[i].uz*uCoord + paramPlaneList[i].vz*vCoord;

                            //Print() << "origin: " << paramPlaneList[i].x0 << ", " << paramPlaneList[i].y0 << ", " << paramPlaneList[i].z0 << ", uCoord: " << uCoord << "\n";

                            //move the particle slightly off the surface so it doesn't intersect it when it moves
                            p.pos(0) = p.pos(0) + uCoord*0.00000001*paramPlaneList[i].lnx;
                            p.pos(1) = p.pos(1) + uCoord*0.00000001*paramPlaneList[i].lny;
                            p.pos(2) = p.pos(2) + uCoord*0.00000001*paramPlaneList[i].lnz;

                            p.rdata(FHD_realData::boostx) = 0;
                            p.rdata(FHD_realData::boosty) = 0;
                            p.rdata(FHD_realData::boostz) = 0;

                            p.idata(FHD_intData::i) = -100;
                            p.idata(FHD_intData::j) = -100;
                            p.idata(FHD_intData::k) = -100;

                            p.rdata(FHD_realData::R) = properties[j].R;
                            p.rdata(FHD_realData::timeFrac) = amrex::Random();

                            Real srt = sqrt(p.rdata(FHD_realData::R)*temp);

                            p.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.);
                            p.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.);
                            p.rdata(FHD_realData::velz) = sqrt(2)*srt*sqrt(-log(amrex::Random()));

                            const paramPlane surf = paramPlaneList[i];

                            rotation(surf.cosThetaLeft, surf.sinThetaLeft, surf.cosPhiLeft, surf.sinPhiLeft, &p.rdata(FHD_realData::velx), &p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz));

                            //Print() << "Pushing back " << p.id() << ", pos: " << p.pos(0) << ", " << p.pos(1) << ", " << p.pos(2) << "\n";
                            //Print() << "Pushing back " << p.id() << ", vel: " << p.rdata(FHD_realData::velx) << ", " << p.rdata(FHD_realData::vely) << ", " << p.rdata(FHD_realData::velz) << "\n";

                            particle_tile.push_back(p);

                            if(p.id() == 332)
                            {
                                Print() << "Generated particle!\n";
                            }

                        }


                    }

                }

            }

        }
    }

    Redistribute();
    SortParticles();

}

void FhdParticleContainer::SortParticles() {
   // Print() << "SORTING\n";
    int lev = 0;
    for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {

        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        const long np = particles.numParticles();

        const Array4<Real> & arrphi = mfphi.array(pti);

        for (int i = 0; i < np; ++ i)
        {
            ParticleType & part = particles[i];

            //Print() << "Checking " << i << ", " << part.id() << ", " << part.idata(FHD_intData::sorted) << "\n";
            //Print() << "vel is " << part.rdata(FHD_realData::velx) << "\n";
            if(part.idata(FHD_intData::sorted) == -1)
            {
                const IntVect& iv = this->Index(part, lev);

               // cout << "part " << i << " is in cell " << iv[0] << ", " << iv[1] << ", " << iv[2] << "\n";

                part.idata(FHD_intData::i) = iv[0];
                part.idata(FHD_intData::j) = iv[1];
                part.idata(FHD_intData::k) = iv[2];

                //Print() << "cell recorded\n";

                long imap = tile_box.index(iv);

                //Print() << "map built\n";
               // Print() << "pti is: " << pti.index() << "\n";
               // Print() << "spec is: " << part.idata(FHD_intData::species) << "\n";

                part.idata(FHD_intData::sorted) = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].size();

               // Print() << "size is: " << part.idata(FHD_intData::sorted) << "\n";

                m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].push_back(i);
                
					 int ispec = part.idata(FHD_intData::species);
					 arrphi(iv[0],iv[1],iv[2],ispec) =
					    arrphi(iv[0],iv[1],iv[2],ispec) + properties[ispec].part2cellVol*properties[ispec].Neff;

                if(part.id() == 332)
                {
                    //Print() << "Adding to " << iv[0] << ", " << iv[1] << ", " << iv[2] << "\n";
                }

            }

        }
    }
}

// Use to approximate pair correlation func with radial distr. func
Real FhdParticleContainer::g0_Ma_Ahmadi(int ispec, int jspec, Real iphi, Real jphi){
	const Real C1 = 2.5, C2 = 4.5904;
	const Real C3 = 4.515439, CPOW = 0.67802;
	Real chi;
	if(ispec==jspec) {
		// Ma Ahmadi 1990s ish
		Real numer = 1 + C1*iphi + C2*pow(iphi,2) + C3*pow(iphi,3);
		numer = numer * 4.0*iphi;
		Real phiRatio = iphi/phi_max;
		Real denom = pow(1.0 - pow(phiRatio,3),CPOW);
		chi = 1+numer/denom;
		return std::min(chi,chi_max);
	} else {
		// Mansoori et al (1971) Equ. Thermo. Prop. of the Mix. of HS
		// Agreement looks good up to 0.50 solid fraction from paper
		// Will likely underestimate near packing
		// Also allows packing fractions past the max packing fraction
		const Real phiTotal = iphi + jphi;
		const Real irad = properties[ispec].radius;
		const Real jrad = properties[jspec].radius;
		Real a1 = 1.0/(1.0-phiTotal);
		Real a2 = 3.0*irad*jrad/(irad+jrad);
		Real eta2 = (iphi/irad) + (jphi/jrad);
		eta2 = eta2/pow(1.0-phiTotal,2);
		a2 = a2 * eta2;
		Real a3 = 2.0 * pow((irad*jrad)/(irad+jrad),2);
		a3 = a3 * pow(eta2,2) / pow(1.0-phiTotal,3);
		chi = a1 + a2 + a3;
		return std::min(chi,chi_max);
	}
}
