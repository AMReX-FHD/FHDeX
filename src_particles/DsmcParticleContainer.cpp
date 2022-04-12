#include "DsmcParticleContainer.H"

// #include "particle_functions_K.H"
#include "paramplane_functions_K.H"
#include <math.h>
using namespace std;
FhdParticleContainer::FhdParticleContainer(const Geometry & geom, const DistributionMapping & dmap,
	const BoxArray & ba, int ncells)
	: NeighborParticleContainer<FHD_realData::count, FHD_intData::count> (geom, dmap, ba, ncells)
{
	BL_PROFILE_VAR("FhdParticleContainer()",FhdParticleContainer);

	realParticles = 0;
	simParticles = 0;
	 
	totalCollisionCells = n_cells[0]*n_cells[1]*n_cells[2];
	domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*(prob_hi[2] - prob_lo[2]);

	collisionCellVol = domainVol/totalCollisionCells;
	ocollisionCellVol = 1.0/collisionCellVol;
	
	bool rho_defined = true;
	if(rho0<0)
	{
		rho0=0.;
		rho_defined = false;
	}
	
	// total = simulated total (total*particle_neff = real number of particles)
	// n0 = real number density
	for(int i=0;i<nspecies;i++)
	{
		properties[i].mass = mass[i];
		properties[i].radius = diameter[i]/2.0;
		properties[i].partVol = pow(diameter[i],3.0)*pi_usr/6.0;
		properties[i].part2cellVol = properties[i].partVol*ocollisionCellVol;
		properties[i].Neff = particle_neff; //   real / simulated
		properties[i].R = k_B/properties[i].mass;

		if(particle_count[i]>=0)
		{
			//Print() << "pcnt\n";
			properties[i].total = particle_count[i];
			properties[i].n0 = particle_neff*properties[i].total/domainVol;
			particle_n0[i] = properties[i].n0;

			amrex::Print() <<  "Species "<< i << " count " << properties[i].total << "\n";
			amrex::Print() <<  "Species "<< i << " n0 " << properties[i].n0 << "\n";
			amrex::Print() <<  "Species " << i << " rho0 " << properties[i].n0*properties[i].mass << "\n";

			rho0 += properties[i].n0*properties[i].mass;
			Yk0[i] = properties[i].n0*properties[i].mass;
		}
		else if(phi_domain[i]>=0)
		{
			//Print() << "phi\n";
			properties[i].total = (int)amrex::Math::ceil(
				(phi_domain[i]*domainVol)/(properties[i].partVol*particle_neff) );
			properties[i].n0 = particle_neff*properties[i].total/domainVol;
			particle_n0[i] = properties[i].n0;

			amrex::Print() <<  "Species "<< i << " count " << properties[i].total << "\n";
			amrex::Print() <<  "Species "<< i << " n0 " << properties[i].n0 << "\n";
			amrex::Print() <<  "Species " << i << " rho0 " << properties[i].n0*properties[i].mass << "\n";

			rho0 += properties[i].n0*properties[i].mass;
			Yk0[i] = properties[i].n0*properties[i].mass;
		}
		else if(rho_defined > 0)
		{
			//Print() << "rho\n";
			Real Yktot = 0.0;
			for(int j=0; j<nspecies; j++) { Yktot += Yk0[j]; }
			for(int j=0; j<nspecies; j++) {	Yk0[j] /= Yktot; }

			Real rhop = properties[i].mass/domainVol;
			properties[i].total = std::ceil((rho0*Yk0[i])/(rhop*particle_neff));
			properties[i].n0 = particle_neff*properties[i].total/domainVol;
			particle_n0[i] = properties[i].n0;

//            Real specRho = rho0*Yk0[i];
//            properties[i].total = std::ceil((specRho/properties[i].mass)*domainVol/particle_neff);
//            properties[i].n0 = particle_neff*properties[i].total/properties[i].mass;
//            particle_n0[i] = properties[i].n0;
            
			amrex::Print() <<  "Species "<< i << " count " << properties[i].total << "\n";
			amrex::Print() <<  "Species "<< i << " n0 " << properties[i].n0 << "\n";
			amrex::Print() <<  "Species " << i << " rho0 " << properties[i].n0*properties[i].mass << ", from " << rho0*Yk0[i] << "\n";
			for(int j=0; j<nspecies; j++)
			{
				Yk0[j] *= rho0;
			}
		}
		else
		{
			//Print() << "n0\n";
			amrex::Print() << "n0: " << particle_n0[i] << "\n";
			properties[i].total = (int)amrex::Math::ceil(particle_n0[i]
				*domainVol/particle_neff);
			properties[i].n0 = particle_neff*properties[i].total/domainVol;
			particle_n0[i] = properties[i].n0;

			amrex::Print() <<  "Species " << i << " count " << properties[i].total << "\n";
			amrex::Print() <<  "Species " << i << " n0 " << properties[i].n0 << "\n";
			amrex::Print() <<  "Species " << i << " rho0 " << properties[i].n0*properties[i].mass << "\n";

			rho0 += properties[i].n0*properties[i].mass;
			Yk0[i] = properties[i].n0*properties[i].mass;
		}
		realParticles = realParticles + properties[i].total*particle_neff;
		simParticles = simParticles + properties[i].total;
		amrex::Print() << "Particles per cell for species " << i 
			<< " is " << properties[i].total/totalCollisionCells << "\n";
	}
	amrex::Print() <<  "Rho0: " << rho0 << "\n";
	amrex::Print() <<  "Total n0: " << realParticles/domainVol << "\n";
	
	for(int i=0; i<nspecies ; i++)
	{
		Yk0[i] = Yk0[i]/rho0;
	}
	
	for (int d=0; d<AMREX_SPACEDIM; ++d)
	{
		domSize[d] = prob_hi[d] - prob_lo[d];
	}

	const int lev=0;
	for(MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi)
	{
		const Box& box = mfi.validbox();
		const int grid_id = mfi.index();

		for(int i=0;i<nspecies;i++)
		{
			m_cell_vectors[i][grid_id].resize(box.numPts());
		}
	}
}

void FhdParticleContainer::MoveParticlesCPP(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount)
{
	BL_PROFILE_VAR("MoveParticlesCPP()", MoveParticlesCPP);

	const int lev = 0;
	const GpuArray<Real, 3> dx = Geom(lev).CellSizeArray();
	const GpuArray<Real, 3> plo = Geom(lev).ProbLoArray();
	const GpuArray<Real, 3> phi = Geom(lev).ProbHiArray();

	int np_tile = 0, np_proc = 0;

	Real adj = 0.99999;
	Real adjalt = 2.0*(1.0-0.99999);
	Real runtime, inttime;
	int intsurf, intside, push;

    int delCount[MAX_SPECIES];
    for(int i = 0; i<MAX_SPECIES;i++)
    {
        delCount[i]=0;
    }
	
	int totalParts = 0;
    amrex::RandomEngine engine;

	for (FhdParIter pti(* this, lev); pti.isValid(); ++pti)
	{
		const int grid_id = pti.index();
		const int tile_id = pti.LocalTileIndex();
		const Box& tile_box  = pti.tilebox();

		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();
		const long np = particles.numParticles();

		Box bx  = pti.tilebox();
		IntVect myLo = bx.smallEnd();
		IntVect myHi = bx.bigEnd();
		
		//cout << "Rank " << ParallelDescriptor::MyProc() << " sees " << np << " particles\n";
		
		totalParts += np;

		for (int i = 0; i < np; i++)
		{
			ParticleType & part = particles[i];
			runtime = dt*part.rdata(FHD_realData::timeFrac);
			while(runtime > 0)
			{
				find_inter_gpu(part, runtime, paramPlaneList, paramPlaneCount,
					&intsurf, &inttime, &intside, ZFILL(plo), ZFILL(phi));

				for (int d=0; d<AMREX_SPACEDIM; ++d)
				{
					part.pos(d) += inttime * part.rdata(FHD_realData::velx + d)*adj;
				}

				runtime = runtime - inttime;

				if(intsurf > 0)
				{
					//find_inter indexes from 1 to maintain compatablity with fortran version
					const paramPlane& surf = paramPlaneList[intsurf-1];
      
					Real posAlt[3];

					for (int d=0; d<AMREX_SPACEDIM; ++d)
					{
						posAlt[d] = inttime * part.rdata(FHD_realData::velx + d)*adjalt;
					}

					Real dummy = 1;
                    //cout << "particle " << part.id() << " intersected " << intsurf << " with vel " << part.rdata(FHD_realData::velx) << ", " << part.rdata(FHD_realData::vely) << ", " << part.rdata(FHD_realData::velz) << endl;
					app_bc_gpu(&surf, part, intside, domSize, &push, &runtime, dummy, engine);
                    if(part.id() == -1)
                    {
                        delCount[part.idata(FHD_intData::species)]++;
                    }
					if(push == 1)
					{
						for (int d=0; d<AMREX_SPACEDIM; ++d)
						{
							part.pos(d) += part.pos(d) + posAlt[d];
						}
					}
				}
			}

			part.rdata(FHD_realData::timeFrac) = 1;

			int cell[3];
			cell[0] = (int)floor((part.pos(0)-plo[0])/dx[0]);
			cell[1] = (int)floor((part.pos(1)-plo[1])/dx[1]);
			cell[2] = (int)floor((part.pos(2)-plo[2])/dx[2]);

			if((part.idata(FHD_intData::i) != cell[0]) || (part.idata(FHD_intData::j) != cell[1]) ||
				(part.idata(FHD_intData::k) != cell[2]) || part.id() < 0 || part.idata(FHD_intData::newSpecies) != -1)
			{
				IntVect iv(part.idata(FHD_intData::i), part.idata(FHD_intData::j), part.idata(FHD_intData::k));
				long imap = tile_box.index(iv);

				int lastIndex = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].size() - 1;
				int lastPart = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap][lastIndex];
				int newIndex = part.idata(FHD_intData::sorted);

				m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap][newIndex] = lastPart;
				m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].pop_back();

				particles[lastPart].idata(FHD_intData::sorted) = newIndex;

				part.idata(FHD_intData::sorted) = -1;
				if(part.idata(FHD_intData::newSpecies) != -1)
				{
				    part.idata(FHD_intData::species) = part.idata(FHD_intData::newSpecies);
				    part.idata(FHD_intData::newSpecies) = -1;
				}
			}
		}
	}
	
//    for(int i = 0; i<MAX_SPECIES;i++)
//    {
//        int temp = delCount[i];
//        ParallelDescriptor::ReduceIntSum(temp);
//        if(temp != 0)
//        {
//            Print() << "Deleted " << temp << " of species " << i << "\n";
//        }
//    }

    ParallelDescriptor::ReduceIntSum(totalParts);
	//Print() << "Total particles: " << totalParts << "\n";
	
	Redistribute();
	SortParticles();
}

void FhdParticleContainer::MovePhononsCPP(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount, const int step)
{
	BL_PROFILE_VAR("MoveParticlesCPP()", MoveParticlesCPP);

	const int lev = 0;
	const GpuArray<Real, 3> dx = Geom(lev).CellSizeArray();
	const GpuArray<Real, 3> plo = Geom(lev).ProbLoArray();
	const GpuArray<Real, 3> phi = Geom(lev).ProbHiArray();

	int np_tile = 0, np_proc = 0, scatterCount = 0, count = 0, specCount = 0;

	Real adj = 0.99999;
	Real adjalt = 2.0*(1.0-0.99999);
	Real runtime, inttime;
	int intsurf, intside, push;

    int delCount[MAX_SPECIES];
    for(int i = 0; i<MAX_SPECIES;i++)
    {
        delCount[i]=0;
    }
	
	int totalParts = 0;
    amrex::RandomEngine engine;

	for (FhdParIter pti(* this, lev); pti.isValid(); ++pti)
	{
		const int grid_id = pti.index();
		const int tile_id = pti.LocalTileIndex();
		const Box& tile_box  = pti.tilebox();

		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();
		const long np = particles.numParticles();

		Box bx  = pti.tilebox();
		IntVect myLo = bx.smallEnd();
		IntVect myHi = bx.bigEnd();
		
		//cout << "Rank " << ParallelDescriptor::MyProc() << " sees " << np << " particles\n";
		
		totalParts += np;

		for (int i = 0; i < np; i++)
		{
			ParticleType & part = particles[i];
			runtime = dt*part.rdata(FHD_realData::timeFrac);
			while(runtime > 0)
			{
				find_inter_gpu(part, runtime, paramPlaneList, paramPlaneCount,
					&intsurf, &inttime, &intside, ZFILL(plo), ZFILL(phi));
				
				Real tauImpurityInv = pow(part.rdata(FHD_realData::omega),4)/tau_i;
				Real tauTAInv = part.rdata(FHD_realData::omega)*pow(T_init[0],4)/tau_ta;
				Real tauLAInv = pow(part.rdata(FHD_realData::omega),2)*pow(T_init[0],3)/tau_la;
				Real tauNormalInv = (2.0*tauTAInv+tauLAInv)/3.0;
				Real tauInv = tauImpurityInv + tauNormalInv;
				
				Real scatterTime = -log(amrex::Random(engine))/tauInv;
				
                if(scatterTime > inttime)
                {
                    runtime = runtime - inttime;
                    
                    for (int d=0; d<AMREX_SPACEDIM; ++d)
				    {
					    part.pos(d) += inttime * part.rdata(FHD_realData::velx + d)*adj;
				    }
                    if(intsurf > 0)
				    {
					    //find_inter indexes from 1 to maintain compatablity with fortran version
					    const paramPlane& surf = paramPlaneList[intsurf-1];
          
					    Real posAlt[3];

					    for (int d=0; d<AMREX_SPACEDIM; ++d)
					    {
						    posAlt[d] = inttime * part.rdata(FHD_realData::velx + d)*adjalt;
					    }
					    
					    app_bc_phonon_gpu(&surf, part, intside, domSize, &push, &runtime, step, &count, &specCount, engine);
                        if(part.id() == -1)
                        {
                            delCount[part.idata(FHD_intData::species)]++;
                        }
					    if(push == 1)
					    {
						    for (int d=0; d<AMREX_SPACEDIM; ++d)
						    {
							    part.pos(d) += part.pos(d) + posAlt[d];
						    }
					    }
				    }				    
				    
				    
                }else
                {
                    runtime = runtime - scatterTime;
                    scatterCount++;
                    for (int d=0; d<AMREX_SPACEDIM; ++d)
				    {
					    part.pos(d) += scatterTime * part.rdata(FHD_realData::velx + d)*adj;
				    }
				    
				    randomSphere(&part.rdata(FHD_realData::velx),&part.rdata(FHD_realData::vely), &part.rdata(FHD_realData::velz), engine);				    
                }
                
			}
			
			part.rdata(FHD_realData::timeFrac) = 1.0;

			int cell[3];
			cell[0] = (int)floor((part.pos(0)-plo[0])/dx[0]);
			cell[1] = (int)floor((part.pos(1)-plo[1])/dx[1]);
			cell[2] = (int)floor((part.pos(2)-plo[2])/dx[2]);

			if((part.idata(FHD_intData::i) != cell[0]) || (part.idata(FHD_intData::j) != cell[1]) ||
				(part.idata(FHD_intData::k) != cell[2]) || part.id() < 0)
			{
				IntVect iv(part.idata(FHD_intData::i), part.idata(FHD_intData::j), part.idata(FHD_intData::k));
				long imap = tile_box.index(iv);

				int lastIndex = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].size() - 1;
				int lastPart = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap][lastIndex];
				int newIndex = part.idata(FHD_intData::sorted);

				m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap][newIndex] = lastPart;
				m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].pop_back();

				particles[lastPart].idata(FHD_intData::sorted) = newIndex;

				part.idata(FHD_intData::sorted) = -1;
			}
		}
	}
	
//    for(int i = 0; i<MAX_SPECIES;i++)
//    {
//        int temp = delCount[i];
//        ParallelDescriptor::ReduceIntSum(temp);
//        if(temp != 0)
//        {
//            Print() << "Deleted " << temp << " of species " << i << "\n";
//        }
//    }

    ParallelDescriptor::ReduceIntSum(scatterCount);
    ParallelDescriptor::ReduceIntSum(totalParts);
    ParallelDescriptor::ReduceIntSum(count);
    ParallelDescriptor::ReduceIntSum(specCount);
	Print() << "Total particles: " << totalParts << "\n";
    Print() << "Internal scattering events: " << scatterCount << "\n";
    if(count != 0)
    {
        Print() << "Fraction of boundary interactions which were specular: " << (double)specCount/((double)count) << "\n";
    }
	
	Redistribute();
	SortParticles();
}

void FhdParticleContainer::SortParticles()
{
	int lev = 0;

	const Real* dx = Geom(lev).CellSize();
	const Real* plo = Geom(lev).ProbLo();
	const Real* phi = Geom(lev).ProbHi();

	for (FhdParIter pti(* this, lev); pti.isValid(); ++pti)
	{
		const int grid_id = pti.index();
		const int tile_id = pti.LocalTileIndex();
		const Box& tile_box  = pti.tilebox();

		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();
		const long np = particles.numParticles();

		for (int i = 0; i < np; ++ i)
		{
			ParticleType & part = particles[i];
			IntVect iv ={0,0,0};

			if(part.idata(FHD_intData::sorted) == -1)
			{
				iv[0] = (int)floor((part.pos(0)-plo[0])/dx[0]);
				iv[1] = (int)floor((part.pos(1)-plo[1])/dx[1]);
				iv[2] = (int)floor((part.pos(2)-plo[2])/dx[2]);
				part.idata(FHD_intData::i) = iv[0];
				part.idata(FHD_intData::j) = iv[1];
				part.idata(FHD_intData::k) = iv[2];
				long imap = tile_box.index(iv);
				part.idata(FHD_intData::sorted) = m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].size();
				m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap].push_back(i);
			}
			else
			{
				iv[0] = part.idata(FHD_intData::i);
				iv[1] = part.idata(FHD_intData::j);
				iv[2] = part.idata(FHD_intData::k);
				long imap = tile_box.index(iv);
				m_cell_vectors[part.idata(FHD_intData::species)][pti.index()][imap][part.idata(FHD_intData::sorted)] = i;
			}
		}
	}
}


//void FhdParticleContainer::SpecChange(FhdParticleContainer::ParticleType& part) {
//	int lev = 0;
//	bool proc_enter = true;
//	
//	const Real* dx = Geom(lev).CellSize();
//	Real smallNumber = dx[0];
//	if(dx[1] < smallNumber){smallNumber = dx[1];}
//	if(dx[2] < smallNumber){smallNumber = dx[2];}
//	smallNumber = smallNumber*0.000001;
//	
//	for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi)
//	{
//		if(proc_enter && ParallelDescriptor::MyProc()==0)
//		{
//			proc_enter = false;//Make sure this runs only once incase of tiling

//			const int grid_id = mfi.index();
//			const int tile_id = mfi.LocalTileIndex();
//			auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

//			for(int i = 0; i< paramPlaneCount; i++)
//			{
//				
//				if(paramPlaneList[i].sourceRight == 1)
//				{
//					for(int j=nspecies-1; j>=0; j--)
//					{
//						Real density = paramPlaneList[i].densityRight[j];						
//						Real temp = paramPlaneList[i].temperatureRight;
////						Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
//						Real area = paramPlaneList[i].area;
//						Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;

//						//cout << "R: " << properties[j].R << endl;
//						//cout << "fluxmean right species " << j << " surface " << i << ": " << fluxMean << "\n";
//												
////						Real elapsedTime = -log(amrex::Random())/fluxMean;
////                        int totalFluxInt = 0;
////						while(elapsedTime < dt)
////						{
////						    totalFluxInt++;
////						    elapsedTime += -log(amrex::Random())/fluxMean;
////						}

//                        int totalFluxInt = amrex::RandomPoisson(dt*fluxMean);

////						Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
////						Real fluxVar = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
////						
////						Real totalFlux = dt*fluxMean + sqrt(dt*fluxVar)*amrex::RandomNormal(0.,1.);
////						totalFlux = std::max(totalFlux,0.);

////						int totalFluxInt =  (int)floor(totalFlux);
////						Real totalFluxLeftOver = totalFlux - totalFluxInt;

////						if(amrex::Random() < totalFluxLeftOver)
////						{
////							totalFluxInt++;
////						}
//						//Print() << "Surface " << i << " generating " << totalFluxInt << " of species " << j << " on the right.\n";

//						for(int k=0;k<totalFluxInt;k++)
//						{
//							Real uCoord = amrex::Random()*paramPlaneList[i].uTop;
//							Real vCoord = amrex::Random()*paramPlaneList[i].vTop;

//							ParticleType p;
//							p.id() = ParticleType::NextID();

//							p.cpu() = ParallelDescriptor::MyProc();
//							p.idata(FHD_intData::sorted) = -1;

//							p.idata(FHD_intData::species) = j;

//							p.pos(0) = paramPlaneList[i].x0 + paramPlaneList[i].ux*uCoord + paramPlaneList[i].vx*vCoord;
//							p.pos(1) = paramPlaneList[i].y0 + paramPlaneList[i].uy*uCoord + paramPlaneList[i].vy*vCoord;
//							p.pos(2) = paramPlaneList[i].z0 + paramPlaneList[i].uz*uCoord + paramPlaneList[i].vz*vCoord;

//							p.pos(0) = p.pos(0) + smallNumber*paramPlaneList[i].rnx;
//							p.pos(1) = p.pos(1) + smallNumber*paramPlaneList[i].rny;
//							p.pos(2) = p.pos(2) + smallNumber*paramPlaneList[i].rnz;
//							
//							p.rdata(FHD_realData::boostx) = 0;
//							p.rdata(FHD_realData::boosty) = 0;
//							p.rdata(FHD_realData::boostz) = 0;

//							p.idata(FHD_intData::i) = -100;
//							p.idata(FHD_intData::j) = -100;
//							p.idata(FHD_intData::k) = -100;

//							p.rdata(FHD_realData::R) = properties[j].R;
//							p.rdata(FHD_realData::timeFrac) = amrex::Random();

//							Real srt = sqrt(p.rdata(FHD_realData::R)*temp);
//							p.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.);
//							p.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.);
//							p.rdata(FHD_realData::velz) = sqrt(2)*srt*sqrt(-log(amrex::Random()));

//							const paramPlane surf = paramPlaneList[i];

//							rotation(surf.cosThetaRight, surf.sinThetaRight, surf.cosPhiRight, surf.sinPhiRight,
//								&p.rdata(FHD_realData::velx), &p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz));
//							
////							cout << p.rdata(FHD_realData::velx) << "," << p.rdata(FHD_realData::vely) << ", " << p.rdata(FHD_realData::velz) << "\n";
////    						cout << p.pos(0) << "," << p.pos(1) << ", " << p.pos(2) << "\n";

//							particle_tile.push_back(p);
//						}
////                        ParallelDescriptor::ReduceIntSum(totalFluxInt);
////                        if(totalFluxInt != 0)
////                        {
////                            Print() << "Surface " << i << " generated " << totalFluxInt << " of species " << j << "\n";
////                        }
//					}
//				}else if(paramPlaneList[i].sourceLeft == 1)
//				{
//					for(int j=nspecies-1; j>=0; j--)
//					{
//						Real density = paramPlaneList[i].densityLeft[j];
//						Real temp = paramPlaneList[i].temperatureLeft;
////						Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
//						Real area = paramPlaneList[i].area;
//						Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;

//						//cout << "R: " << properties[j].R << endl;
//						//cout << "fluxmean left species " << j << " surface " << i << ": " << fluxMean << "\n";
//												
////						Real elapsedTime = -log(amrex::Random())/fluxMean;
////                        int totalFluxInt = 0;
////						while(elapsedTime < dt)
////						{
////						    totalFluxInt++;
////						    elapsedTime += -log(amrex::Random())/fluxMean;
////						}

//                        int totalFluxInt = amrex::RandomPoisson(dt*fluxMean);


////						Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
////						Real fluxVar = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;

////						Real totalFlux = dt*fluxMean + sqrt(dt*fluxVar)*amrex::RandomNormal(0.,1.);
////						totalFlux = std::max(totalFlux,0.);

////						int totalFluxInt =  (int)floor(totalFlux);
////						Real totalFluxLeftOver = totalFlux - totalFluxInt;

////						if(amrex::Random() < totalFluxLeftOver)
////						{
////							totalFluxInt++;
////						}
//						
//						//

//						for(int k=0;k<totalFluxInt;k++)
//						{
//							Real uCoord = amrex::Random()*paramPlaneList[i].uTop;
//							Real vCoord = amrex::Random()*paramPlaneList[i].vTop;

//							ParticleType p;
//							p.id() = ParticleType::NextID();

//							p.cpu() = ParallelDescriptor::MyProc();
//							p.idata(FHD_intData::sorted) = -1;

//							p.idata(FHD_intData::species) = j;

//							p.pos(0) = paramPlaneList[i].x0 + paramPlaneList[i].ux*uCoord + paramPlaneList[i].vx*vCoord;
//							p.pos(1) = paramPlaneList[i].y0 + paramPlaneList[i].uy*uCoord + paramPlaneList[i].vy*vCoord;
//							p.pos(2) = paramPlaneList[i].z0 + paramPlaneList[i].uz*uCoord + paramPlaneList[i].vz*vCoord;

//							//move the particle slightly off the surface so it doesn't intersect it when it moves
//							p.pos(0) = p.pos(0) + smallNumber*paramPlaneList[i].lnx;
//							p.pos(1) = p.pos(1) + smallNumber*paramPlaneList[i].lny;
//							p.pos(2) = p.pos(2) + smallNumber*paramPlaneList[i].lnz;

//							p.rdata(FHD_realData::boostx) = 0;
//							p.rdata(FHD_realData::boosty) = 0;
//							p.rdata(FHD_realData::boostz) = 0;

//							p.idata(FHD_intData::i) = -100;
//							p.idata(FHD_intData::j) = -100;
//							p.idata(FHD_intData::k) = -100;

//							p.rdata(FHD_realData::R) = properties[j].R;
//							p.rdata(FHD_realData::timeFrac) = amrex::Random();

//							Real srt = sqrt(p.rdata(FHD_realData::R)*temp);
//							p.rdata(FHD_realData::velx) = srt*amrex::RandomNormal(0.,1.);
//							p.rdata(FHD_realData::vely) = srt*amrex::RandomNormal(0.,1.);
//							p.rdata(FHD_realData::velz) = sqrt(2)*srt*sqrt(-log(amrex::Random()));

//							const paramPlane surf = paramPlaneList[i];
//							
//	
//							rotation(surf.cosThetaLeft, surf.sinThetaLeft, surf.cosPhiLeft, surf.sinPhiLeft,
//                                                &p.rdata(FHD_realData::velx),&p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz));


//							particle_tile.push_back(p);
//						}
////                        ParallelDescriptor::ReduceIntSum(totalFluxInt);
////                        if(totalFluxInt != 0)
////                        {
////                            Print() << "Surface " << i << " generated " << totalFluxInt << " of species " << j << "\n";
////                        }
//                        
//					}
//				}
//			}
//		}
//	}
//	Redistribute();
//	SortParticles();
//}


void FhdParticleContainer::Source(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount) {
	int lev = 0;
	bool proc_enter = true;
	
	const Real* dx = Geom(lev).CellSize();
	Real smallNumber = dx[0];
	if(dx[1] < smallNumber){smallNumber = dx[1];}
	if(dx[2] < smallNumber){smallNumber = dx[2];}
	smallNumber = smallNumber*0.000001;
	
	for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi)
	{
		if(proc_enter)
		{
			proc_enter = false;//Make sure this runs only once incase of tiling

			const int grid_id = mfi.index();
			const int tile_id = mfi.LocalTileIndex();
			auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

			for(int i = 0; i< paramPlaneCount; i++)
			{
				
				if(paramPlaneList[i].sourceRight == 1)
				{
					for(int j=nspecies-1; j>=0; j--)
					{
						Real density = paramPlaneList[i].densityRight[j];						
						Real temp = paramPlaneList[i].temperatureRight;
						Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
						Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;

						//cout << "R: " << properties[j].R << endl;
						//cout << "fluxmean right species " << j << " surface " << i << ": " << fluxMean << "\n";
												
//						Real elapsedTime = -log(amrex::Random())/fluxMean;
//                        int totalFluxInt = 0;
//						while(elapsedTime < dt)
//						{
//						    totalFluxInt++;
//						    elapsedTime += -log(amrex::Random())/fluxMean;
//						}

                        int totalFluxInt = amrex::RandomPoisson(dt*fluxMean);

//						Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
//						Real fluxVar = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
//						
//						Real totalFlux = dt*fluxMean + sqrt(dt*fluxVar)*amrex::RandomNormal(0.,1.);
//						totalFlux = std::max(totalFlux,0.);

//						int totalFluxInt =  (int)floor(totalFlux);
//						Real totalFluxLeftOver = totalFlux - totalFluxInt;

//						if(amrex::Random() < totalFluxLeftOver)
//						{
//							totalFluxInt++;
//						}
						//Print() << "Surface " << i << " generating " << totalFluxInt << " of species " << j << " on the right.\n";

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

							p.pos(0) = p.pos(0) + smallNumber*paramPlaneList[i].rnx;
							p.pos(1) = p.pos(1) + smallNumber*paramPlaneList[i].rny;
							p.pos(2) = p.pos(2) + smallNumber*paramPlaneList[i].rnz;
							
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

							rotation(surf.cosThetaRight, surf.sinThetaRight, surf.cosPhiRight, surf.sinPhiRight,
								&p.rdata(FHD_realData::velx), &p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz));
							
//							cout << p.rdata(FHD_realData::velx) << "," << p.rdata(FHD_realData::vely) << ", " << p.rdata(FHD_realData::velz) << "\n";
//    						cout << p.pos(0) << "," << p.pos(1) << ", " << p.pos(2) << "\n";

							particle_tile.push_back(p);
						}
//                        ParallelDescriptor::ReduceIntSum(totalFluxInt);
//                        if(totalFluxInt != 0)
//                        {
//                            Print() << "Surface " << i << " generated " << totalFluxInt << " of species " << j << "\n";
//                        }
					}
				}else if(paramPlaneList[i].sourceLeft == 1)
				{
					for(int j=nspecies-1; j>=0; j--)
					{
						Real density = paramPlaneList[i].densityLeft[j];
						Real temp = paramPlaneList[i].temperatureLeft;
						Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
						Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;

						//cout << "R: " << properties[j].R << endl;
						//cout << "fluxmean left species " << j << " surface " << i << ": " << fluxMean << "\n";
												
//						Real elapsedTime = -log(amrex::Random())/fluxMean;
//                        int totalFluxInt = 0;
//						while(elapsedTime < dt)
//						{
//						    totalFluxInt++;
//						    elapsedTime += -log(amrex::Random())/fluxMean;
//						}

                        int totalFluxInt = amrex::RandomPoisson(dt*fluxMean);


//						Real fluxMean = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;
//						Real fluxVar = density*area*sqrt(properties[j].R*temp/(2.0*M_PI))/particle_neff;

//						Real totalFlux = dt*fluxMean + sqrt(dt*fluxVar)*amrex::RandomNormal(0.,1.);
//						totalFlux = std::max(totalFlux,0.);

//						int totalFluxInt =  (int)floor(totalFlux);
//						Real totalFluxLeftOver = totalFlux - totalFluxInt;

//						if(amrex::Random() < totalFluxLeftOver)
//						{
//							totalFluxInt++;
//						}
						
						//

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

							//move the particle slightly off the surface so it doesn't intersect it when it moves
							p.pos(0) = p.pos(0) + smallNumber*paramPlaneList[i].lnx;
							p.pos(1) = p.pos(1) + smallNumber*paramPlaneList[i].lny;
							p.pos(2) = p.pos(2) + smallNumber*paramPlaneList[i].lnz;

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
							
	
							rotation(surf.cosThetaLeft, surf.sinThetaLeft, surf.cosPhiLeft, surf.sinPhiLeft,
                                                &p.rdata(FHD_realData::velx),&p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz));


							particle_tile.push_back(p);
						}
//                        ParallelDescriptor::ReduceIntSum(totalFluxInt);
//                        if(totalFluxInt != 0)
//                        {
//                            Print() << "Surface " << i << " generated " << totalFluxInt << " of species " << j << "\n";
//                        }
                        
					}
				}
			}
		}
	}
	Redistribute();
	SortParticles();
}

void FhdParticleContainer::SourcePhonons(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount) {
	int lev = 0;
	bool proc_enter = true;
	
	const Real* dx = Geom(lev).CellSize();
	Real smallNumber = dx[0];
	if(dx[1] < smallNumber){smallNumber = dx[1];}
	if(dx[2] < smallNumber){smallNumber = dx[2];}
	smallNumber = smallNumber*0.00000001;
	
	amrex::RandomEngine engine;
	
	for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi)
	{
		if(proc_enter)
		{
			proc_enter = false;//Make sure this runs only once incase of tiling

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
						Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
						Real fluxMean = density*area;

						Real totalFlux = dt*fluxMean;

						int totalFluxInt =  (int)floor(totalFlux);
						Real totalFluxLeftOver = totalFlux - totalFluxInt;

						if(amrex::Random() < totalFluxLeftOver)
						{
							totalFluxInt++;
						}

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

							//move the particle slightly off the surface so it doesn't intersect it when it moves
							p.pos(0) = p.pos(0) + smallNumber*paramPlaneList[i].lnx;
							p.pos(1) = p.pos(1) + smallNumber*paramPlaneList[i].lny;
							p.pos(2) = p.pos(2) + smallNumber*paramPlaneList[i].lnz;


							p.idata(FHD_intData::i) = -100;
							p.idata(FHD_intData::j) = -100;
							p.idata(FHD_intData::k) = -100;

							p.rdata(FHD_realData::timeFrac) = amrex::Random();


							const paramPlane surf = paramPlaneList[i];
							
							cosineLawHemisphere(surf.cosThetaLeft, surf.sinThetaLeft, surf.cosPhiLeft, surf.sinPhiLeft,
                                                &p.rdata(FHD_realData::velx),&p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz), phonon_sound_speed, engine);

                            p.rdata(FHD_realData::omega) = plankDist(surf.temperatureLeft, engine);
                            //I hope this is right?
                            p.rdata(FHD_realData::lambda) = phonon_sound_speed*2.0*M_PI/p.rdata(FHD_realData::omega);

                            //cout << "velx: " << p.rdata(FHD_realData::velx) << "\n";

							particle_tile.push_back(p);
						}
//                        ParallelDescriptor::ReduceIntSum(totalFluxInt);
//                        if(totalFluxInt != 0)
//                        {
//                            Print() << "Surface " << i << " generated " << totalFluxInt << " of species " << j << "\n";
//                        }
                        
					}
				}
				else if(paramPlaneList[i].sourceRight == 1)
				{
					for(int j=0; j< nspecies; j++)
					{
						Real density = paramPlaneList[i].densityRight[j];
						Real temp = paramPlaneList[i].temperatureRight;
						Real area = paramPlaneList[i].area/ParallelDescriptor::NProcs();
						Real fluxMean = density*area;

						Real totalFlux = dt*fluxMean;

						int totalFluxInt =  (int)floor(totalFlux);
						Real totalFluxLeftOver = totalFlux - totalFluxInt;

						if(amrex::Random() < totalFluxLeftOver)
						{
							totalFluxInt++;
						}


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

							p.pos(0) = p.pos(0) + smallNumber*paramPlaneList[i].rnx;
							p.pos(1) = p.pos(1) + smallNumber*paramPlaneList[i].rny;
							p.pos(2) = p.pos(2) + smallNumber*paramPlaneList[i].rnz;
							
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

							cosineLawHemisphere(surf.cosThetaRight, surf.sinThetaRight, surf.cosPhiRight, surf.sinPhiRight,
								&p.rdata(FHD_realData::velx), &p.rdata(FHD_realData::vely), &p.rdata(FHD_realData::velz),phonon_sound_speed, engine);

							particle_tile.push_back(p);
						}
//                        ParallelDescriptor::ReduceIntSum(totalFluxInt);
//                        if(totalFluxInt != 0)
//                        {
//                            Print() << "Surface " << i << " generated " << totalFluxInt << " of species " << j << "\n";
//                        }
					}
				}
			}
		}
	}
	Redistribute();
	SortParticles();
}

/*
void FhdParticleContainer::EvaluateStats(MultiFab& particleInstant, MultiFab& particleMeans,
                                         MultiFab& particleVars, const Real delt, int steps)
{
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
    
    for (FhdParIter pti(*this, lev); pti.isValid(); ++pti) 
    {

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

            for(int l=0;l<nspecies;l++)
            {

            }

        });

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            part_inst(i,j,k,1) = part_inst(i,j,k,1)*cellVolInv;
             
        });

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            part_mean(i,j,k,1)  = (part_mean(i,j,k,1)*stepsMinusOne + part_inst(i,j,k,1))*stepsInv;

            for(int l=0;l<nspecies;l++)
            {

            }
             
        });

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            Real del1 = part_inst(i,j,k,1) - part_mean(i,j,k,1);

            part_var(i,j,k,1)  = (part_var(i,j,k,1)*stepsMinusOne + del1*del1)*stepsInv;            
        });

    }

}
*/
