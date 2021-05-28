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
		amrex::Print() << "Particle Volume:" << properties[i].partVol << "\n";
      properties[i].part2cellVol = properties[i].partVol*ocollisionCellVol;
		properties[i].Neff = particle_neff; // assume F_n is same for each
      
      // Overwrite particle_count
      properties[i].total = (phi_domain[i]*domainVol)/properties[i].partVol;
		//amrex::Print() << "Phi: " << phi_domain[i] << "\n";
		//amrex::Print() << "Part Vol: " << properties[i].partVol << "\n";
		//amrex::Print() << "Domain Vol: " << domainVol << "\n";
		//amrex::Print() << "Np: " << properties[i].total << "\n";
		//amrex::Print() << "Phi: " << phi_domain[i] << "\n";
      //amrex::Print() << "\n";
      realParticles = realParticles + properties[i].total*properties[i].Neff;
      //amrex::Print() << "Measured Phi: " << properties[i].total*properties[i].partVol/domainVol << "\n";
      properties[i].total = properties[i].total/properties[i].Neff;
      //amrex::Print() << "Measured Phi: " << properties[i].total*properties[i].partVol/domainVol << "\n";
      simParticles = simParticles + properties[i].total;
      amrex::Print() << "Npi: " << properties[i].total << "\n";
	}
   
   int indx;
   int cnt = 0;
	for(int i_spec=0;i_spec<nspecies;i_spec++) {
		for(int j_spec=i_spec;j_spec<nspecies;j_spec++) {
			indx = getSpeciesIndex(i_spec,j_spec);
    		interproperties[indx].inel = (1+alpha_pp[cnt])/6.0;
    		// amrex::Print() << "Alpha: " << alpha_pp[cnt] << "\n";
    		interproperties[indx].csx = (diameter[i_spec] + diameter[j_spec])*0.5*pi_usr;
    		// amrex::Print() << "i: " << i_spec << " j: " << j_spec << "\n";
    		// amrex::Print() << "InElastic: " << interproperties[indx].inel << "\n";
    		// amrex::Print() << "CSX: " << interproperties[indx].csx << "\n";
    		// amrex::Print() << "\n";
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
		m_cell_vectors[grid_id].resize(box.numPts());
	}
}

// Likely excessive
// Relates lower diagonal matrix indices to those of 1D array
int FhdParticleContainer::getSpeciesIndex(int species1, int species2) {
	if(species1<species2){
		return species1+(nspecies-1)*species2;
	} else {
		return species2+(nspecies-1)*species1;
	}
}

// Use to approximate pair correlation fuunc with radial distr. func
Real FhdParticleContainer::g0_Ma_Ahmadi(int species1, int species2, Real phi1, Real phi2){
	const double C1 = 2.5, C2 = 4.5904;
	const double C3 = 4.515439, CPOW = 0.67802;

	if(species1==species2) {
		Real phiTotal = phi1 + phi2;
		Real numer = 1 + C1*phiTotal + C2*pow(phiTotal,2) + C3*pow(phiTotal,3);
		numer = numer * 4.0*phiTotal;
		Real phiRatio = phiTotal/phi_max;
		Real denom = pow(1.0 - pow(phiRatio,3),CPOW);
		return std::min(1+numer/denom,chi_max);
	} else {
	 // add Lebowitz RDF here
	}
	return 0;
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

		Box bx  = pti.tilebox();
		IntVect myLo = bx.smallEnd();
		IntVect myHi = bx.bigEnd();

		np_proc += np;

		for (int i = 0; i < np; ++ i) {
			ParticleType & part = particles[i];
			Real speed = 0;

			for (int d=0; d<AMREX_SPACEDIM; ++d) {                   
				speed += part.rdata(FHD_realData::velx + d)*part.rdata(FHD_realData::velx + d);
			}

			if(speed > maxspeed) {
				maxspeed = speed;
			}

			moves++;
			runtime = dt;

			while(runtime > 0) {
				find_inter_gpu(part, runtime, paramPlaneList, paramPlaneCount, &intsurf, &inttime, &intside, ZFILL(plo), ZFILL(phi));

				for (int d=0; d<AMREX_SPACEDIM; ++d) {
					part.pos(d) += inttime * part.rdata(FHD_realData::velx + d)*adj;
				}
				runtime = runtime - inttime;

				if(intsurf > 0) {
					const paramPlane& surf = paramPlaneList[intsurf-1];
					Real posAlt[3];
					for (int d=0; d<AMREX_SPACEDIM; ++d) {
						posAlt[d] = inttime * part.rdata(FHD_realData::velx + d)*adjalt;
					}
					Real dummy = 1;
					app_bc_gpu(&surf, part, intside, domSize, &push, dummy, dummy);
					if(push == 1) {
						for (int d=0; d<AMREX_SPACEDIM; ++d) {
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

			if((cell[0] < myLo[0]) || (cell[1] < myLo[1]) || (cell[2] < myLo[2])) {
				reDist++;
			} 
			else if((cell[0] > myHi[0]) || (cell[1] > myHi[1]) || (cell[2] > myHi[2])) {
				reDist++;
			}

			if((part.idata(FHD_intData::i) != cell[0]) || (part.idata(FHD_intData::j) != cell[1]) || (part.idata(FHD_intData::k) != cell[2])) {
				//remove particle from old cell
				IntVect iv(part.idata(FHD_intData::i), part.idata(FHD_intData::j), part.idata(FHD_intData::k));
				long imap = tile_box.index(iv);
				int lastIndex = m_cell_vectors[pti.index()][imap].size() - 1;
				int lastPart = m_cell_vectors[pti.index()][imap][lastIndex];
				int newIndex = part.idata(FHD_intData::sorted);

				m_cell_vectors[pti.index()][imap][newIndex] = m_cell_vectors[pti.index()][imap][lastIndex];
				m_cell_vectors[pti.index()][imap].pop_back();
            part.idata(FHD_intData::sorted) = -1; // indicates need to be sorted
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

		for (int i = 0; i < np; ++ i) {
			ParticleType & part = particles[i];
			if(part.idata(FHD_intData::sorted) == -1) {
				const IntVect& iv = this->Index(part, lev);

				part.idata(FHD_intData::i) = iv[0];
				part.idata(FHD_intData::j) = iv[1];
				part.idata(FHD_intData::k) = iv[2];

				long imap = tile_box.index(iv);

			   // size gives num of particles in the coll cell
				part.idata(FHD_intData::sorted) = m_cell_vectors[pti.index()][imap].size();
				m_cell_vectors[pti.index()][imap].push_back(i);
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
		// const long np = particles.numParticles();
		
		// Convert MultiFabs -> arrays
		const Array4<Real> & arrvrmax = mfvrmax.array(pti);
		const Array4<Real> & arrphi = mfphi.array(pti);
		const Array4<Real> & arrselect = mfselect.array(pti);
		const Array4<Real> & arrnspec = mfnspec.array(pti);
		
		// if we have scalar in the particle container, can update and keep track
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

			// Initialize mfs
			int ij_spec;
			for (int i_spec=0; i_spec<nspecies; i_spec++) {
				arrnspec(i,j,k,i_spec) = 0.0;	
				for (int j_spec = i_spec; j_spec < nspecies; j_spec++) {
					ij_spec = getSpeciesIndex(i_spec,j_spec);
					// Print() << ij_spec << "\n";
					arrselect(i,j,k,ij_spec) = 0.0;
				}
			}
			
			ParticleType p;
			int p_spec, pindx;
			
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);
			const long np = m_cell_vectors[grid_id][imap].size();
			for (int l=0; l<np; l++) {
				pindx = m_cell_vectors[grid_id][imap][l];
				ParticleType p = particles[pindx];
				p_spec = p.idata(FHD_intData::species);
				arrnspec(i,j,k,p_spec) = arrnspec(i,j,k,p_spec) + 1; // merge daniel's
			}
			// amrex::Print() << "Np: " << arrnspec(i,j,k,0) << "\n";
			for (int i_spec=0; i_spec<nspecies; i_spec++) {
				// phi = np * particle volume * collision cell volume * Neff
				arrphi(i,j,k,i_spec) = arrnspec(i,j,k,i_spec)
					*properties[i_spec].partVol*ocollisionCellVol*properties[i_spec].Neff;
				// amrex::Print() << "Phi: " << arrphi(i,j,k,i_spec) << "\n";
			}
		});
		
	}
}


// Compute selections here
void FhdParticleContainer::CalcCollisionCells(Real dt) {
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
		
		int spec_indx;
		Real NSel; // Number of selections this time step
		Real phi1, phi2, chi0; // radial distribution function (enhance collision frequency)
		Real crossSection;
		int np_i, np_j; // number of particles in collision cell for species i and j
		
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {       
		
			Real vrmax;
			// Update number of selections
			// Note that the inner loop goes from 0 to nspecies but the selections calculated
			// ... are stored in the same index.
			for (int i_spec = 0; i_spec < nspecies; i_spec++) {
				for (int j_spec = 0; j_spec < nspecies; j_spec++) {
					spec_indx = getSpeciesIndex(i_spec,j_spec);
					np_i = arrnspec(i,j,k,spec_indx);
					np_j = arrnspec(i,j,k,spec_indx);
					phi1 = arrphi(i,j,k,i_spec);
					phi2 = arrphi(i,j,k,j_spec);
					// comment out if expecting dilute
					chi0 = g0_Ma_Ahmadi(i_spec,j_spec, phi1, phi2);
					vrmax = arrvrmax(i,j,k,spec_indx);
					crossSection = interproperties[spec_indx].csx;
					NSel = 0.5*np_i*np_j*crossSection*vrmax*ocollisionCellVol*chi0*dt;
					arrselect(i,j,k,spec_indx) = arrselect(i,j,k,spec_indx) + NSel;
					// Needs to be tested!
					// NSel = std::floor(Nsel + amrex::random());
					// arrselect(i,j,k,spec_indx) = NSel;
				}
			}
		});
	}
}

void FhdParticleContainer::CollideParticles() {
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
		
			RealVect eij, vreij;
			RealVect vel1, vel2, vr;
			Real vrmag, vrmax, vreijmag;

			// Garic: Random selections - no species order
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
			
			
			
			
			int NSel; // Number of selections this time step
			int pindx1, pindx2; // index of randomly sampled particles
			int p1_spec, p2_spec; // species of random particles
			ParticleType p1, p2;
			int np_i, np_j, np_total; // number of particles in collision cell for species i and j
			// Loops through species pairs
			for (int i_spec = 0; i_spec < nspecies; i_spec++) {
				for (int j_spec = i_spec; j_spec < nspecies; j_spec++) {
				spec_indx = getSpeciesIndex(i_spec,j_spec);
				NSel = floor(arrselect(i,j,k,spec_indx));
				arrselect(i,j,k,spec_indx) = arrselect(i,j,k,spec_indx) - NSel;
				np_i = arrnspec(i,j,k,spec_indx);
				np_j = arrnspec(i,j,k,spec_indx);
				np_total = np_i + np_j; // initially sample from all particles
				vrmax = arrvrmax(i,j,k,spec_indx);
				// Loop through selections
				for (int isel = 0; isel < NSel; isel++) {
					pairFound = false;
					while(!pairFound) {
						// p1 = floor(amrex::Random()*np_i);
						// p2 = floor(amrex::Random()*np_j);
						pindx1 = floor(amrex::Random()*np_total);
						pindx2 = floor(amrex::Random()*np_total);
						// p1 = m_cell_vectors[i_spec][grid_id][imap][p1];
						// p2 = m_cell_vectors[j_spec][grid_id][imap][p2];
						pindx1 = m_cell_vectors[grid_id][imap][pindx1];
						pindx2 = m_cell_vectors[grid_id][imap][pindx2];
						p1 = particles[pindx1];
						p2 = particles[pindx2];
						p1_spec = p1.idata(FHD_intData::species);
						p2_spec = p2.idata(FHD_intData::species);
						if(p1_spec == i_spec && p2_spec == j_spec) {
							pairFound = true;
						}
					}
					vel1[0] = p1.rdata(FHD_realData::velx);
					vel1[1] = p1.rdata(FHD_realData::vely);
					vel1[2] = p1.rdata(FHD_realData::velz);
					vel2[0] = p2.rdata(FHD_realData::velx);
					vel2[1] = p2.rdata(FHD_realData::vely);
					vel2[2] = p2.rdata(FHD_realData::velz);
					vr[0] = vel2[0]-vel1[0];
					vr[1] = vel2[1]-vel1[1];
					vr[2] = vel2[2]-vel1[2];
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
