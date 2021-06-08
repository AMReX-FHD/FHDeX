#include "DsmcParticleContainer.H"

// #include "particle_functions_K.H"
#include "paramplane_functions_K.H"
#include <math.h>

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
	// amrex::Print() << "Max Grid Size: " << max_grid_size[0]*max_grid_size[1]*max_grid_size[2] << "\n";
	domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*(prob_hi[2] - prob_lo[2]);

	collisionCellVol = domainVol/totalCollisionCells;
	ocollisionCellVol = 1/collisionCellVol;
	for(int i=0;i<nspecies;i++) {
		properties[i].mass = mass[i];
		properties[i].radius = diameter[i]/2.0;
		properties[i].partVol = pow(diameter[i],3)*pi_usr/6;
		//amrex::Print() << "Particle Volume:" << properties[i].partVol << "\n";
      properties[i].part2cellVol = properties[i].partVol*ocollisionCellVol;
		properties[i].Neff = particle_neff; // assume F_n is same for each
      
      // Overwrite particle_count
      properties[i].total = (phi_domain[i]*domainVol)/properties[i].partVol;
      realParticles = realParticles + properties[i].total*properties[i].Neff;
      //amrex::Print() << "Measured Phi: " << properties[i].total*properties[i].partVol/domainVol << "\n";
      properties[i].total = properties[i].total/properties[i].Neff;
      //amrex::Print() << "Measured Phi: " << properties[i].total*properties[i].partVol/domainVol << "\n";
      simParticles = simParticles + properties[i].total;
      //amrex::Print() << "Npi: " << properties[i].total << "\n";
	}
   
   int indx;
   int ij_spec;
   int cnt = 0;
	for(int i_spec=0;i_spec<nspecies;i_spec++) {
		for(int j_spec=i_spec;j_spec<nspecies;j_spec++) {
			ij_spec = getSpeciesIndex(i_spec,j_spec);
    		interproperties[ij_spec].inel = (1+alpha_pp[cnt])/6.0; // need to find better way to input this
    		interproperties[ij_spec].csx = pow(properties[i_spec].radius+properties[j_spec].radius,2)*pi_usr;
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

              runtime = dt;

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

                      app_bc_gpu(&surf, part, intside, domSize, &push, dummy, dummy);

                      if(push == 1)
                      {
                          for (int d=0; d<AMREX_SPACEDIM; ++d)
                          {
                              part.pos(d) += part.pos(d) + posAlt[d];
                          }
                      }
                  }

              }


            int cell[3];
            cell[0] = (int)floor((part.pos(0)-plo[0])/dx[0]);
            cell[1] = (int)floor((part.pos(1)-plo[1])/dx[1]);
            cell[2] = (int)floor((part.pos(2)-plo[2])/dx[2]);

            //cout << "current cell: " << cell[0] << ", " << cell[1] << ", " << cell[2] << "\cout";
            //n << "current pos: " << part.pos(0) << ", " << part.pos(1) << ", " << part.pos(2) << "\n";

            if((cell[0] < myLo[0]) || (cell[1] < myLo[1]) || (cell[2] < myLo[2]) || (cell[0] > myHi[0]) || (cell[1] > myHi[1]) || (cell[2] > myHi[2]))
            {
                reDist++;
            }

            if((part.idata(FHD_intData::i) != cell[0]) || (part.idata(FHD_intData::j) != cell[1]) || (part.idata(FHD_intData::k) != cell[2]))
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
                //int ispec = part.idata(FHD_intData::species);
                //arrphi(iv[0],iv[1],iv[2],ispec) =
                //	arrphi(iv[0],iv[1],iv[2],ispec) - properties[ispec].partVol*properties[ispec].Neff*ocollisionCellVol;
                
            }
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

void FhdParticleContainer::SortParticles() {
	int lev = 0;
	for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {
		const int grid_id = pti.index();
		const int tile_id = pti.LocalTileIndex();
		const Box& tile_box  = pti.tilebox();

		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();
		const long np = particles.numParticles();
		
      // Update solid fraction here ( and at MoveParticlesCPP() )
      const Array4<Real> & arrphi = mfphi.array(pti);

		for (int i = 0; i < np; ++ i) {
			ParticleType & part = particles[i];
			if(part.idata(FHD_intData::sorted) == -1) {
				const IntVect& iv = this->Index(part, lev);

				part.idata(FHD_intData::i) = iv[0];
				part.idata(FHD_intData::j) = iv[1];
				part.idata(FHD_intData::k) = iv[2];

				long imap = tile_box.index(iv);
				//cout << "part " << i << " is in cell " << iv[0] << ", " << iv[1] << ", " << iv[2] << ", adding to element " << m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].size() << "\n";

				part.idata(FHD_intData::sorted) = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].size();

				m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].push_back(i);
				
				// increase solid fraction in cell
				//int ispec = part.idata(FHD_intData::species);
				//arrphi(iv[0],iv[1],iv[2],ispec) =
				//	arrphi(iv[0],iv[1],iv[2],ispec) + properties[ispec].partVol*properties[ispec].Neff*ocollisionCellVol;
			}
		}
	}
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
		
		// if we have scalar in the particle container, can update and keep track
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

			
			// Initialize mfs
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
				// phi = np * particle volume * collision cell volume * Neff
				arrphi(i,j,k,i_spec) = m_cell_vectors[i_spec][grid_id][imap].size()
					*properties[i_spec].partVol*ocollisionCellVol*properties[i_spec].Neff;
			}
		});
	}
}


// Compute selections here
void FhdParticleContainer::CalcSelections(Real dt) {
	int lev = 0;
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
			
			// Update number of selections
			// Note that the inner loop goes from 0 to nspecies but the selections calculated
			// ... for same species pair are stored in the same index.
			// i.e. 1-2 and 2-1 both calculated and stored at same index ij_spec
			for (int i_spec = 0; i_spec < nspecies; i_spec++) {
				for (int j_spec = 0; j_spec < nspecies; j_spec++) {
					ij_spec = getSpeciesIndex(i_spec,j_spec);
					np_i = m_cell_vectors[i_spec][grid_id][imap].size();
					np_j = m_cell_vectors[j_spec][grid_id][imap].size();
					phi1 = arrphi(i,j,k,i_spec);
					phi2 = arrphi(i,j,k,j_spec);
					// comment out if expecting dilute
					chi0 = g0_Ma_Ahmadi(i_spec,j_spec, phi1, phi2);
					vrmax = arrvrmax(i,j,k,i_spec);
					crossSection = interproperties[ij_spec].csx;
					
					//amrex::Print() << "Old Selections: " << arrselect(i,j,k,ij_spec) << "\n";
					NSel = 0.5*np_i*np_j*crossSection*vrmax*ocollisionCellVol*chi0*dt;
					// Currently running total
					arrselect(i,j,k,ij_spec) = arrselect(i,j,k,ij_spec) + NSel;
					//amrex::Print() << "New Selections: " << arrselect(i,j,k,ij_spec) << "\n";
					// Needs to be tested!
					// NSel = std::floor(Nsel + amrex::random());
					// arrselect(i,j,k,spec_indx) = NSel;
				}
			}
			// check that solid fractions are being updated over time
			//for (int i_spec = 0; i_spec < nspecies; i_spec++) {
			//	Real phimf = arrphi(i,j,k,i_spec);
			//	Real phinp = m_cell_vectors[i_spec][grid_id][imap].size()
			//		*properties[i_spec].partVol*properties[i_spec].Neff*ocollisionCellVol;
			//	amrex::Print() << "Phi MF: " << phimf << " Phi MCell: " << phinp << "\n";
			//}
		});
	}
}

void FhdParticleContainer::CollideParticles() {
	/*
	int lev = 0;
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
		const Array4<Real> & arrnspec = mfnspec.array(mfi);
		//Print() << arrvrmax << "\n";
		
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
		
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);
			bool pairFound;

			// Garcia: Random selections - no species order
			/*
			int spec_indx;
			int totalSel = 0;
			int Nsel_spec[nspecies][3];
			for (int i_spec=0; i < nspecies; i_spec++) {
				for(int j_spec=i_spec; j < nspecies; j_spec++) {
					spec_indx = getSpeciesIndex(i_spec,j_spec);
					Nsel_spec[spec_indx][0] = i_spec;
					Nsel_spec[spec_indx][1] = j_spec;
					Nsel_spec[spec_indx][2] = arrselect(i,j,k,spec_indx);
					totalSel = totalSel + arrselect(i,j,k,spec_indx);
				}
			}

			Real RR;
			int selSum;
			int i_spec, j_spec;
			for (int isel=0; isel<totalSel; isel++) {
				RR = amrex::Random();
				selSum = 0.;
				spec_indx = 0;
				while((float)selSum/totalSel<RR) {
					selSum = selSum + Nsel_spec[spec_indx][2];
					spec_indx++;
				}
				spec_indx--;
				i_spec = Nsel_spec[spec_indx][0];
				j_spec = Nsel_spec[spec_indx][1];
				Nsel_spec[spec_indx][2] = Nsel_spec[spec_indx][2] - 1;
				
				
			}
			//
			
			
			
			
			int NSel; // Number of selections this time step
			int pindxi, pindxj; // index of randomly sampled particles
			ParticleType pi, pj;
			int np_i, np_j; // number of particles in collision cell for species i and j
			
			RealVect eij, vreij;
			RealVect vi, vj, vij;
			Real vrmag, vrmax, vreijmag;
			// Loops through species pairs
			for (int i_spec = 0; i_spec < nspecies; i_spec++) {
				for (int j_spec = i_spec; j_spec < nspecies; j_spec++) {
				spec_indx = getSpeciesIndex(i_spec,j_spec);
				NSel = floor(arrselect(i,j,k,spec_indx));
				arrselect(i,j,k,spec_indx) = arrselect(i,j,k,spec_indx) - NSel;
				np_i = m_cell_vectors[i_spec][grid_id][imap].size();
				np_j = m_cell_vectors[j_spec][grid_id][imap].size();
				np_total = np_i + np_j; // initially sample from all particles
				vrmax = arrvrmax(i,j,k,spec_indx);
				// Loop through selections
				for (int isel = 0; isel < NSel; isel++) {
					pindxi = floor(amrex::Random()*np_i);
					pindxj = floor(amrex::Random()*np_j);
					pindxi = m_cell_vectors[i_spec][grid_id][imap][pindxi];
					pindxj = m_cell_vectors[j_spec][grid_id][imap][pindxj];
					pi = particles[pindx1];
					pj = particles[pindx2];

					vi[0] = pi.rdata(FHD_realData::velx);
					vi[1] = pi.rdata(FHD_realData::vely);
					vi[2] = pi.rdata(FHD_realData::velz);
					vj[0] = pj.rdata(FHD_realData::velx);
					vj[1] = pj.rdata(FHD_realData::vely);
					vj[2] = pj.rdata(FHD_realData::velz);
					vij[0] = vj[0]-vi[0];
					vij[1] = vj[1]-vi[1];
					vij[2] = vj[2]-vi[2];
					// vrmag = sqrt(vr.dotProduct(vr));
					vrmag = sqrt(vr[0]*vr[0]+vr[1]*vr[1]+vr[2]*vr[2]);
					// If relative speed greater than max relative speed, replace
					if(vrmag>vrmax) {
						vrmax = vrmag*1.2;
						arrvrmax(i,j,k,spec_indx) = vrmax;
					}
					// later want to reject non-approaching
					if(vrmag>vrmax*amrex::Random()) {
						// sample random unit vector at impact
						for (int dir = 0; dir < 3; dir++) {
							eij[dir] = amrex::Random();
						}
						vreijmag = vr[0]*eij[0]+vr[1]*eij[1]+vr[2]*eij[2]; // dot_product
						vreijmag = vreijmag*interproperties[spec_indx].inel;
						//vreij = vr.dotProduct(eij) * eij * interproperties[spec_indx].inel;;
						vreij = vreijmag*eij;
						vel2 = vel2 + vreij;
						vel1 = vel1 - vreij;
						// add boosted velocity calculations here
					}
				}
			}
		}
	}
	*/
}

void FhdParticleContainer::EvaluateStats(MultiFab& particleInstant, MultiFab& particleMeans,
                                         MultiFab& particleVars, const Real delt, int steps) {
	BL_PROFILE_VAR("EvaluateStats()",EvaluateStats);
    
	const int lev = 0;
    
	BoxArray ba = particleMeans.boxArray();
	long cellcount = ba.numPts();

	const Real* dx = Geom(lev).CellSize();
	const Real dxInv = 1.0/dx[0];
	const Real cellVolInv = 1.0/(dx[0]*dx[0]*dx[0]);

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
	}
}
