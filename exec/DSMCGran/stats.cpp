#include "DsmcParticleContainer.H"
using namespace std;
void FhdParticleContainer::EvaluateStats(MultiFab& mfcuInst,    MultiFab& mfcuMean,
													  MultiFab& mfcuDel, 	 MultiFab& mfcovar,
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
		nvars = 13;
		Array4<Real> cuInst = mfcuInst[pti].array();
		for (int i = smallEnd[0]; i <= bigEnd[0]; i++) {
		for (int j = smallEnd[1]; j <= bigEnd[1]; j++) {
		for (int k = smallEnd[2]; k <= bigEnd[2]; k++) {
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);
			int cnt = nvars;
			for (int l=0; l<nspecies; l++) {
				const long np_spec = m_cell_vectors[l][grid_id][imap].size();
				Real moV  = properties[l].mass*ocollisionCellVol;     // density
				Real mass = properties[l].mass;
				cuInst(i,j,k,cnt)   = np_spec*ocollisionCellVol;		// number of species l
				cuInst(i,j,k,0)    += np_spec*ocollisionCellVol; 		// total number
				cuInst(i,j,k,cnt+1) = np_spec*moV;         		    	// mass density of species l
				cuInst(i,j,k,1)    += np_spec*moV;							// total mass density
				
				// Read particle data
				for (int m=0; m<np_spec; m++) {
					int pind = m_cell_vectors[l][grid_id][imap][m];
					ParticleType & p = particles[pind]; 
										
					cuInst(i,j,k,cnt+2)  += p.rdata(FHD_realData::velx); // x-mom of species l
					cuInst(i,j,k,cnt+3)  += p.rdata(FHD_realData::vely); // y-mom of species l
					cuInst(i,j,k,cnt+4)  += p.rdata(FHD_realData::velz); // z-mom of species l
					
					Real u = p.rdata(FHD_realData::velx);
					Real v = p.rdata(FHD_realData::vely);
					Real w = p.rdata(FHD_realData::velz);
					cuInst(i,j,k,cnt+5)  += pow(u,2);							   	 // tau_xx of spec l
					cuInst(i,j,k,cnt+6)  += u*v; 										    // tau_xy of spec l
					cuInst(i,j,k,cnt+7)  += u*w; 							 			    // tau_xz of spec l
					cuInst(i,j,k,cnt+8)  += pow(v,2); 									 // tau_yy of spec l
					cuInst(i,j,k,cnt+9)  += v*w; 							 	          // tau_yz of spec l
					cuInst(i,j,k,cnt+10) += pow(w,2); 								    // tau_zz of spec l
					cuInst(i,j,k,cnt+11) += pow(u,2)+pow(v,2)+pow(w,2);          // Energy of spec l
					cuInst(i,j,k,cnt+12) += (pow(u,2)+pow(v,2)+pow(w,2))*mass/3; // Gran. Temp.
				}

				// Convert velocities to momentum densities
				cuInst(i,j,k,2) += cuInst(i,j,k,cnt+2)*moV;	 // total x-mom
				cuInst(i,j,k,3) += cuInst(i,j,k,cnt+3)*moV;   // total y-mom
				cuInst(i,j,k,4) += cuInst(i,j,k,cnt+4)*moV;   // total z-mom
				
				// Convert speeds squared to energy/stress densities
				// Stresses do not include cross terms atm
				cuInst(i,j,k,5)  += cuInst(i,j,k,cnt+5)*moV*0.5;  // total tau_xx
				cuInst(i,j,k,6)  += cuInst(i,j,k,cnt+6)*moV*0.5;  // total tau_xy
				cuInst(i,j,k,7)  += cuInst(i,j,k,cnt+7)*moV*0.5;  // total tau_xz
				cuInst(i,j,k,8)  += cuInst(i,j,k,cnt+8)*moV*0.5;  // total tau_yy
				cuInst(i,j,k,9)  += cuInst(i,j,k,cnt+9)*moV*0.5;  // total tau_yz
				cuInst(i,j,k,10) += cuInst(i,j,k,cnt+10)*moV*0.5; // total tau_zz
				cuInst(i,j,k,11) += cuInst(i,j,k,cnt+11)*moV*0.5; // total energy
				// Gran. Temp of species l x number density of species l
				cuInst(i,j,k,12) += cuInst(i,j,k,cnt+12)*cuInst(i,j,k,cnt); // total Gran. Temp.
				cnt += nvars;
			}
			
			// Mixture temperature (Divide by total number density)
			cuInst(i,j,k,12) /= cuInst(i,j,k,0);
			
			// Determine mole/mass fraction
			cnt = nvars;
			for (int l=0; l<nspecies; l++) {
				cuInst(i,j,k,cnt)   /= cuInst(i,j,k,0); // mole fraction of species l
				cuInst(i,j,k,cnt+1) /= cuInst(i,j,k,1); // mass fraction of species l
				cnt += nvars;
			}
		}
		}
		}
		
	   // Mean Values
		Array4<Real> cuMean = mfcuMean[pti].array();
		Array4<Real> cuDel = mfcuDel[pti].array();
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			for (int l=0; l<nvars; l++) {
				cuMean(i,j,k,l) = (cuMean(i,j,k,l)*stepsMinusOne+cuInst(i,j,k,l))*osteps;
				cuDel(i,j,k,l)  = cuInst(i,j,k,l) - cuMean(i,j,k,l);
			}

			int cnt = nvars;
			for (int l=0; l<nspecies; l++) {
				for (int m=0; m<nvars; m++) {
					cuMean(i,j,k,cnt+m) = (cuMean(i,j,k,cnt+m)*stepsMinusOne+cuInst(i,j,k,cnt+m))*osteps;
					cuDel(i,j,k,cnt+m)  = cuInst(i,j,k,cnt+m) - cuMean(i,j,k,cnt+m);
				}
				cnt += nvars;
			}
		});		

		// Co/variances
		/*
		0  - drho.drho
		1  - drho.dJx
		2  - drho.dJy
		3  - drho.dJz
		4  - drho.dE
		5  - dJx.dJx
		6  - dJx.dJy
		7  - dJx.dJz
		8  - dJy.dJy
		9  - dJy.dJz
		10 - dJz.dJz
		11 - dJx.dE
		12 - dJy.dE
		13 - dJz.dE
		14 - dE.dE
		*/
		Array4<Real> covar = mfcovar[pti].array();
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			covar(i,j,k,0)  = (covar(i,j,k,0)*stepsMinusOne+cuDel(i,j,k,1)*cuDel(i,j,k,1))*osteps;    // 0
			covar(i,j,k,1)  = (covar(i,j,k,1)*stepsMinusOne+cuDel(i,j,k,1)*cuDel(i,j,k,2))*osteps;    // 1
			covar(i,j,k,2)  = (covar(i,j,k,2)*stepsMinusOne+cuDel(i,j,k,1)*cuDel(i,j,k,3))*osteps;    // 2
			covar(i,j,k,3)  = (covar(i,j,k,3)*stepsMinusOne+cuDel(i,j,k,1)*cuDel(i,j,k,4))*osteps;    // 3
			covar(i,j,k,4)  = (covar(i,j,k,4)*stepsMinusOne+cuDel(i,j,k,1)*cuDel(i,j,k,11))*osteps;   // 4
			covar(i,j,k,5)  = (covar(i,j,k,5)*stepsMinusOne+cuDel(i,j,k,2)*cuDel(i,j,k,2))*osteps;    // 5
			covar(i,j,k,6)  = (covar(i,j,k,6)*stepsMinusOne+cuDel(i,j,k,2)*cuDel(i,j,k,3))*osteps;    // 6
			covar(i,j,k,7)  = (covar(i,j,k,7)*stepsMinusOne+cuDel(i,j,k,2)*cuDel(i,j,k,4))*osteps;    // 7
			covar(i,j,k,8)  = (covar(i,j,k,8)*stepsMinusOne+cuDel(i,j,k,2)*cuDel(i,j,k,11))*osteps;   // 8
			covar(i,j,k,9)  = (covar(i,j,k,9)*stepsMinusOne+cuDel(i,j,k,3)*cuDel(i,j,k,3))*osteps;    // 9
			covar(i,j,k,10) = (covar(i,j,k,10)*stepsMinusOne+cuDel(i,j,k,3)*cuDel(i,j,k,4))*osteps;   // 10
			covar(i,j,k,11) = (covar(i,j,k,11)*stepsMinusOne+cuDel(i,j,k,3)*cuDel(i,j,k,11))*osteps;  // 11
			covar(i,j,k,12) = (covar(i,j,k,12)*stepsMinusOne+cuDel(i,j,k,4)*cuDel(i,j,k,4))*osteps;   // 12
			covar(i,j,k,13) = (covar(i,j,k,13)*stepsMinusOne+cuDel(i,j,k,4)*cuDel(i,j,k,11))*osteps;  // 13
			covar(i,j,k,14) = (covar(i,j,k,14)*stepsMinusOne+cuDel(i,j,k,11)*cuDel(i,j,k,11))*osteps; // 14
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
