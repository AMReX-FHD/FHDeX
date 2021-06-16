#include "DsmcParticleContainer.H"
using namespace std;
void FhdParticleContainer::EvaluateStats(MultiFab& Cu,    
                                         MultiFab& CuMeans,
                                         MultiFab& CuVars,
													               MultiFab& CoVars,
                                         int steps) {
	BL_PROFILE_VAR("EvaluateStats()",EvaluateStats);

	const Real osteps = 1.0/steps;
	const int stepsMinusOne = steps-1;

   const int lev = 0;  
	for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {
		const int grid_id = pti.index();
		const int tile_id = pti.LocalTileIndex();
		const Box& tile_box  = pti.tilebox();
		auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
		auto& particles = particle_tile.GetArrayOfStructs();
		IntVect smallEnd = tile_box.smallEnd();
		IntVect bigEnd = tile_box.bigEnd();
		
		// Instantaneous Values
		/*
		0 	- n - number density
		1	- rho - mass density
		2 	- Jx - x-mom density
		3 	- Jy - y-mom density
		4 	- Jz - z-mom density
		5 	- tau_xx - xx shear stress
		6 	- tau_xy - xy shear stress
		7 	- tau_xz - xz shear stress
		8 	- tau_yy - yy shear stress
		9  - tau_yz - yz shear stress
		10 - tau_zz - zz shear stress
		11 - E - energy
		12 - T_g - granular temperature
		13 - X_i - mole fraction for spec. i
		14 - Y_i - mass fraction for spec. i
		15 - u_i - x-vel for spec. i
		16 - v_i - y-vel for spec. i
		17 - w_i - z-vel for spec. i
		18 - uu_i 
		19 - uv_i
		20 - uw_i
		21 - vv_i
		22 - vw_i
		23 - ww_i
		24 - E_i
		25 - T_g_i
		... (repeat for each add. species)
		*/
		int nstats = 13;
		Array4<Real> cu = Cu[pti].array();
		Array4<Real> cumean = CuMeans[pti].array();
		Array4<Real> cuvars = CuVars[pti].array();
		Array4<Real> covars = CoVars[pti].array();
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);
			int cnt = nstats;
			for (int l=0; l<nspecies; l++) {
				const long np_spec = m_cell_vectors[l][grid_id][imap].size();
				Real moV  = properties[l].mass*ocollisionCellVol*properties[l].Neff;     // density
				Real mass = properties[l].mass*properties[l].Neff;
				cu(i,j,k,cnt)   = np_spec*ocollisionCellVol;		   // number of species l
				cu(i,j,k,0)    += np_spec*ocollisionCellVol; 		// total number
				cu(i,j,k,cnt+1) = np_spec*moV;         		    	// mass density of species l
				cu(i,j,k,1)    += np_spec*moV;							// total mass density
				
				// Read particle data
				for (int m=0; m<np_spec; m++) {
					int pind = m_cell_vectors[l][grid_id][imap][m];
					ParticleType ptemp = particles[pind];
					ParticleType & p = ptemp;
					// ParticleType & p = particles[pind]; 
										
					cu(i,j,k,cnt+2)  += p.rdata(FHD_realData::velx); // x-mom of species l
					cu(i,j,k,cnt+3)  += p.rdata(FHD_realData::vely); // y-mom of species l
					cu(i,j,k,cnt+4)  += p.rdata(FHD_realData::velz); // z-mom of species l
					
					Real u = p.rdata(FHD_realData::velx);
					Real v = p.rdata(FHD_realData::vely);
					Real w = p.rdata(FHD_realData::velz);
					cu(i,j,k,cnt+5)  += pow(u,2);							   	    // tau_xx of spec l
					cu(i,j,k,cnt+6)  += u*v; 										    // tau_xy of spec l
					cu(i,j,k,cnt+7)  += u*w; 							 			    // tau_xz of spec l
					cu(i,j,k,cnt+8)  += pow(v,2); 									 // tau_yy of spec l
					cu(i,j,k,cnt+9)  += v*w; 							 	          // tau_yz of spec l
					cu(i,j,k,cnt+10) += pow(w,2); 								    // tau_zz of spec l
					cu(i,j,k,cnt+11) += (pow(u,2)+pow(v,2)+pow(w,2))*mass/3;  // Energy of spec l
					cu(i,j,k,cnt+12) += pow(u,2)+pow(v,2)+pow(w,2);           // Gran. Temp.
				}

				// Convert velocities to momentum densities
				cu(i,j,k,2) += cu(i,j,k,cnt+2)*moV;	  // total x-mom
				cu(i,j,k,3) += cu(i,j,k,cnt+3)*moV;   // total y-mom
				cu(i,j,k,4) += cu(i,j,k,cnt+4)*moV;   // total z-mom
				
				// Convert speeds squared to energy/stress densities
				// Stresses do not include cross terms atm
				cu(i,j,k,5)  += cu(i,j,k,cnt+5)*moV*0.5;  // total tau_xx
				cu(i,j,k,6)  += cu(i,j,k,cnt+6)*moV*0.5;  // total tau_xy
				cu(i,j,k,7)  += cu(i,j,k,cnt+7)*moV*0.5;  // total tau_xz
				cu(i,j,k,8)  += cu(i,j,k,cnt+8)*moV*0.5;  // total tau_yy
				cu(i,j,k,9)  += cu(i,j,k,cnt+9)*moV*0.5;  // total tau_yz
				cu(i,j,k,10) += cu(i,j,k,cnt+10)*moV*0.5; // total tau_zz
				cu(i,j,k,11) += cu(i,j,k,cnt+11)*moV*0.5; // total energy
				// Gran. Temp of species l x number density of species l
				cu(i,j,k,12) += cu(i,j,k,cnt+12)*cu(i,j,k,cnt); // total Gran. Temp.
				cnt += nstats;
			}
			
			// Mixture temperature (Divide by total number density)
			cu(i,j,k,12) /= cu(i,j,k,0);
			
			// Determine mole/mass fraction
			cnt = nstats;
			for (int l=0; l<nspecies; l++) {
				cu(i,j,k,cnt)   /= cu(i,j,k,0); // mole fraction of species l
				cu(i,j,k,cnt+1) /= cu(i,j,k,1); // mass fraction of species l
				cnt += nstats;
			}
		});
		
	   // Mean & Variances
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
      Vector<Real> delfluc((nspecies+1)*13, 0.0);
			for (int l=0; l<(nspecies+1)*13; l++) {
				cumean(i,j,k,l) = (cumean(i,j,k,l)*stepsMinusOne+cu(i,j,k,l))*osteps;
         	delfluc[l] = cu(i,j,k,l) - cumean(i,j,k,l);
         	cuvars(i,j,k,l) = (cuvars(i,j,k,l)*stepsMinusOne+delfluc[l]*delfluc[l])*osteps;
			}

     		// Covariances
      	/*
      		0  - drho.dJx
      		1  - drho.dJy
      		2  - drho.dJz
      		3  - drho.dE
      		4  - dJx.dJy
      		5  - dJx.dJz
      		6  - dJy.dJz
      		7 - dJx.dE
      		8 - dJy.dE
      		9 - dJz.dE
      		// Comment -- need to add some temperature-velocity covariances as well (see Balakrishnan)
      	*/
			covars(i,j,k,0)  = (covars(i,j,k,0)*stepsMinusOne+delfluc[1]*delfluc[2])*osteps;    // 0
			covars(i,j,k,1)  = (covars(i,j,k,1)*stepsMinusOne+delfluc[1]*delfluc[3])*osteps;    // 1
			covars(i,j,k,2)  = (covars(i,j,k,2)*stepsMinusOne+delfluc[1]*delfluc[4])*osteps;    // 2
			covars(i,j,k,3)  = (covars(i,j,k,3)*stepsMinusOne+delfluc[1]*delfluc[11])*osteps;   // 3
			covars(i,j,k,4)  = (covars(i,j,k,4)*stepsMinusOne+delfluc[2]*delfluc[3])*osteps;    // 4
			covars(i,j,k,5)  = (covars(i,j,k,5)*stepsMinusOne+delfluc[2]*delfluc[4])*osteps;    // 5
			covars(i,j,k,6)  = (covars(i,j,k,6)*stepsMinusOne+delfluc[3]*delfluc[4])*osteps;    // 6
			covars(i,j,k,7)  = (covars(i,j,k,7)*stepsMinusOne+delfluc[2]*delfluc[11])*osteps;   // 7
			covars(i,j,k,8)  = (covars(i,j,k,8)*stepsMinusOne+delfluc[3]*delfluc[11])*osteps;   // 8
			covars(i,j,k,9)  = (covars(i,j,k,9)*stepsMinusOne+delfluc[4]*delfluc[11])*osteps;  // 9
		});
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
