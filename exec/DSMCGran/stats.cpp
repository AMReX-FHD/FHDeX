#include "DsmcParticleContainer.H"
using namespace std;
void FhdParticleContainer::EvaluateStats(MultiFab& mfcuConInst,
   						 MultiFab& mfcuConMeans,
   						 MultiFab& mfcuPrimInst,
   						 MultiFab& mfcuPrimMeans,
   						 MultiFab& mfcuVars,
   						 MultiFab& mfcoVars,
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
		
		/*
			Conserved Vars:
			0  - rho = (1/V) += m
			1  - Jx  = (1/V) += mu
			2  - Jy  = (1/V) += mv
			3  - Jz  = (1/V) += mw
			4  - K   = (1/V) += m|v|^2
			... (repeat for each species)

		   Primitive Vars:
			0	- n   (X_ns)
			1  - rho (Y_ns)
			2  - u   (u_ns)
			3  - v   (v_ns)
			4  - w   (w_ns)
			5  - uu  (uu_ns)
			6  - uv  (uv_ns)
			7  - uw  (uw_ns)
			8  - vv  (vv_ns)
			9  - vw  (vw_ns)
			10  - ww  (ww_ns)
			11 - T   (P_ns)  = (1/3V) += m|v|^2
			12 - P   (T_ns)  = (1/3Nk_B) += m|v|^2
			13 - E   (E_ns)  = (1/2) += |v|^2 + c_v*T
			... (repeat for each species)
		*/
		
		int ncon  = (nspecies+1)*5;
		int nprim = (nspecies+1)*14;
		Array4<Real> cuConInst   = mfcuConInst[pti].array();
		Array4<Real> cuPrimInst  = mfcuPrimInst[pti].array();	
		
   	//////////////////////////////////////
   	// Instantaneous Values
   	//////////////////////////////////////
   	
   	// Primitive by species + Conserved by species/total
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);
			int icon = 5; int iprim = 14;
			Real cv      = 0.;
			for (int l=0; l<nspecies; l++) {
				const long np_spec = m_cell_vectors[l][grid_id][imap].size();
				Real mass = properties[l].mass*properties[l].Neff;
				Real moV  = mass*ocollisionCellVol;    						// density
				Real cv_l = 3.0*k_B*0.5/mass;                        		// assume ideal gas <- in general untrue for most granular cases
				cuPrimInst(i,j,k,iprim+0) = np_spec*ocollisionCellVol;	// number density of species l
				cuPrimInst(i,j,k,0)      += cuPrimInst(i,j,k,iprim+0);
				cuPrimInst(i,j,k,iprim+1) = np_spec*moV;
				cuPrimInst(i,j,k,1)      += cuPrimInst(i,j,k,iprim+1);
				cuConInst(i,j,k,icon+0)   = np_spec*moV;         		   // mass density of species l
				cv      					 	 += cv_l*cuConInst(i,j,k,icon);  // total cv
				
				// Read particle data
				for (int m=0; m<np_spec; m++) {
					int pind = m_cell_vectors[l][grid_id][imap][m];
					ParticleType ptemp = particles[pind];
					ParticleType & p = ptemp;
					// ParticleType & p = particles[pind];
					Real u = p.rdata(FHD_realData::velx);
					Real v = p.rdata(FHD_realData::vely);
					Real w = p.rdata(FHD_realData::velz);

					cuConInst(i,j,k,icon+1) += u;
					cuConInst(i,j,k,icon+2) += v;
					cuConInst(i,j,k,icon+3) += w;
					cuConInst(i,j,k,icon+4) += (pow(u,2)+pow(v,2)+pow(w,2));
				}
				
				// Momentum densities
				cuConInst(i,j,k,icon+1) = cuConInst(i,j,k,icon+1)*moV;  		// x-mom density
				cuConInst(i,j,k,icon+2) = cuConInst(i,j,k,icon+2)*moV;  		// y-mom density
				cuConInst(i,j,k,icon+3) = cuConInst(i,j,k,icon+3)*moV;  		// z-mom density

				// Energy density
				cuConInst(i,j,k,icon+4) = cuConInst(i,j,k,icon+4)*moV*0.5;

				// Total of cons. vars
				for (int m=0; m<5; m++) {cuConInst(i,j,k,m) += cuConInst(i,j,k,icon+m);}

				// Care not to use mean of instananeous primitives as mean (use conserved)
				cuPrimInst(i,j,k,iprim+ 2) = cuConInst(i,j,k,icon+1)/cuConInst(i,j,k,icon);					// u  of spec l
				cuPrimInst(i,j,k,iprim+ 3) = cuConInst(i,j,k,icon+2)/cuConInst(i,j,k,icon);					// v  of spec l
				cuPrimInst(i,j,k,iprim+ 4) = cuConInst(i,j,k,icon+3)/cuConInst(i,j,k,icon);					// w  of spec l
				cuPrimInst(i,j,k,iprim+ 5) = pow(cuPrimInst(i,j,k,iprim+2),2);										// uu of spec l
				cuPrimInst(i,j,k,iprim+ 6) = cuPrimInst(i,j,k,iprim+1)*cuPrimInst(i,j,k,iprim+2);			// uv of spec l
				cuPrimInst(i,j,k,iprim+ 7) = cuPrimInst(i,j,k,iprim+1)*cuPrimInst(i,j,k,iprim+3);			// uw of spec l
				cuPrimInst(i,j,k,iprim+ 8) = pow(cuPrimInst(i,j,k,iprim+3),2);										// vv of spec l
				cuPrimInst(i,j,k,iprim+ 9) = cuPrimInst(i,j,k,iprim+2)*cuPrimInst(i,j,k,iprim+3);			// vw of spec l
				cuPrimInst(i,j,k,iprim+10) = pow(cuPrimInst(i,j,k,iprim+4),2);										// ww of spec l
				
				Real vsqb = 0.5*(pow(cuConInst(i,j,k,icon+1),2) + 
								pow(cuConInst(i,j,k,icon+2),2) +
								pow(cuConInst(i,j,k,icon+3),2));
				
				cuPrimInst(i,j,k,iprim+11) = (cuConInst(i,j,k,icon+4)/cuConInst(i,j,k,icon)-vsqb*0.5)/cv_l;	// T  of spec l
				cuPrimInst(i,j,k,iprim+12) = cuPrimInst(i,j,k,iprim+11)*(k_B/mass)*cuConInst(i,j,k,icon);		// P  of spec l
				cuPrimInst(i,j,k,iprim+13) = vsqb*moV+cv_l*cuPrimInst(i,j,k,iprim+11)*cuConInst(i,j,k,icon);	// E  of spec l

				// Total of primitive vars
				// Handle temperature/energy seperately
				cuPrimInst(i,j,k,11) += cuPrimInst(i,j,k,iprim+11)*cuPrimInst(i,j,k,iprim+0);                // Mixture T = sum of (nk*Tk)/n
				cuPrimInst(i,j,k,12) += cuPrimInst(i,j,k,iprim+12);														// Mixture P = sum of partial P
				icon += 5; iprim += 14;
			}
			
			// Primitive total
			cuPrimInst(i,j,k,2) = cuConInst(i,j,k,1)/cuConInst(i,j,k,0);		// Bulk x-velocity
			cuPrimInst(i,j,k,3) = cuConInst(i,j,k,2)/cuConInst(i,j,k,0);		// Bulk y-velocity
			cuPrimInst(i,j,k,4) = cuConInst(i,j,k,3)/cuConInst(i,j,k,0);		// Bulk z-velocity
			
			cuPrimInst(i,j,k,5)  = pow(cuConInst(i,j,k,1),2)/pow(cuConInst(i,j,k,0),2);					// Bulk uu
			cuPrimInst(i,j,k,6)  = cuConInst(i,j,k,1)*cuConInst(i,j,k,2)/pow(cuConInst(i,j,k,0),2);	// Bulk uv
			cuPrimInst(i,j,k,6)  = cuConInst(i,j,k,1)*cuConInst(i,j,k,3)/pow(cuConInst(i,j,k,0),2);	// Bulk uw
			cuPrimInst(i,j,k,8)  = pow(cuConInst(i,j,k,2),2)/pow(cuConInst(i,j,k,0),2);					// Bulk vv
			cuPrimInst(i,j,k,6)  = cuConInst(i,j,k,2)*cuConInst(i,j,k,3)/pow(cuConInst(i,j,k,0),2);	// Bulk vw
			cuPrimInst(i,j,k,10) = pow(cuConInst(i,j,k,3),2)/pow(cuConInst(i,j,k,0),2);					// Bulk ww
			
			cuPrimInst(i,j,k,11) /= cuPrimInst(i,j,k,0);																// T
			
			cv /= cuConInst(i,j,k,0);
			// Energy Density
			cuPrimInst(i,j,k,13)  = pow(cuConInst(i,j,k,1),2)+pow(cuConInst(i,j,k,2),2)+pow(cuConInst(i,j,k,3),2);
			cuPrimInst(i,j,k,13) *= (0.5*cuConInst(i,j,k,0));   // Bulk energy
			cuPrimInst(i,j,k,13) += (cv*cuPrimInst(i,j,k,11));  // Total Particle KE
			
			// Convert n to Xk and rho to Yk
			for (int l=0; l<nspecies; l++) {
				cuPrimInst(i,j,k,14*l+0)  /= cuPrimInst(i,j,k,0); //X_k
				// cuPrimInst(i,j,k,14*l+1)  /= cuConInst(i,j,k,0);  //Y_k
			}
		});
		
	   // Mean & Variances
	   // Need to rewrite primitive means in terms of the conserved vars
	   Array4<Real> cuConMeans  = mfcuConMeans[pti].array();
	   Array4<Real> cuPrimMeans = mfcuPrimMeans[pti].array();
		Array4<Real> cuVars      = mfcuVars[pti].array();
		Array4<Real> coVars      = mfcoVars[pti].array();
		
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
    		Vector<Real> delCon(ncon, 0.0);
      	// Conserved Variances
			for (int l=0; l<ncon; l++) {
				cuConMeans(i,j,k,l) = (cuConMeans(i,j,k,l)*stepsMinusOne+cuConInst(i,j,k,l))*osteps;
         	delCon[l] = cuConInst(i,j,k,l) - cuConMeans(i,j,k,l);
         	cuVars(i,j,k,l) = (cuVars(i,j,k,l)*stepsMinusOne+delCon[l]*delCon[l])*osteps;
			}
			
			// Primitive Variances
			// Derive from conserved means
			// Append to end of cuvars
			int prevEnd = ncon;
			Vector<Real> delPrim(nprim, 0.0);
			for (int l=0; l<nprim; l++) {
				cuPrimMeans(i,j,k,l) = (cuPrimMeans(i,j,k,l)*stepsMinusOne+cuPrimInst(i,j,k,l))*osteps;
         	delPrim[l] = cuPrimInst(i,j,k,l) - cuPrimMeans(i,j,k,l);
         	cuVars(i,j,k,prevEnd+l) = (cuVars(i,j,k,prevEnd+l)*stepsMinusOne+delPrim[l]*delPrim[l])*osteps;
			}

			// Covariances
			/*
				0  - drho.dJx
				1  - drho.dJy
				2  - drho.dJz
				3  - drho.dT
				4  - drho.d(rho*E)
				5  - dJx.dJy
				6  - dJx.dJz
				7  - dJy.dJz
				8  - dJx.d(rho*E)
				9  - dJy.d(rho*E)
				10 - dJz.d(rho*E)
				11 - drho.du
				12 - drho.dv
				13 - drho.dw
				14 - du.dv
				15 - du.dw
				16 - dv.dw
				17 - drho.dT
				18 - du.dT
				19 - dv.dT
				20 - dw.dT
			*/
      	
			coVars(i,j,k, 0)  = (coVars(i,j,k, 0)*stepsMinusOne+delCon[0]*delCon[1])*osteps;
			coVars(i,j,k, 1)  = (coVars(i,j,k, 1)*stepsMinusOne+delCon[0]*delCon[2])*osteps;
			coVars(i,j,k, 2)  = (coVars(i,j,k, 2)*stepsMinusOne+delCon[0]*delCon[3])*osteps;
			coVars(i,j,k, 3)  = (coVars(i,j,k, 3)*stepsMinusOne+delCon[0]*delPrim[11])*osteps;
			coVars(i,j,k, 4)  = (coVars(i,j,k, 4)*stepsMinusOne+delCon[0]*delPrim[13])*osteps;
			coVars(i,j,k, 5)  = (coVars(i,j,k, 5)*stepsMinusOne+delCon[1]*delCon[2])*osteps;
			coVars(i,j,k, 6)  = (coVars(i,j,k, 6)*stepsMinusOne+delCon[1]*delCon[3])*osteps;
			coVars(i,j,k, 7)  = (coVars(i,j,k, 7)*stepsMinusOne+delCon[2]*delCon[3])*osteps;
			coVars(i,j,k, 8)  = (coVars(i,j,k, 8)*stepsMinusOne+delCon[1]*delPrim[13])*osteps;
			coVars(i,j,k, 9)  = (coVars(i,j,k, 9)*stepsMinusOne+delCon[2]*delPrim[13])*osteps;
			coVars(i,j,k,10)  = (coVars(i,j,k,10)*stepsMinusOne+delCon[3]*delPrim[13])*osteps;
			coVars(i,j,k,11)  = (coVars(i,j,k,11)*stepsMinusOne+delCon[0]*delPrim[1])*osteps;
			coVars(i,j,k,12)  = (coVars(i,j,k,12)*stepsMinusOne+delCon[0]*delPrim[2])*osteps;
			coVars(i,j,k,13)  = (coVars(i,j,k,13)*stepsMinusOne+delCon[0]*delPrim[3])*osteps;
			coVars(i,j,k,14)  = (coVars(i,j,k,14)*stepsMinusOne+delPrim[1]*delPrim[2])*osteps;
			coVars(i,j,k,15)  = (coVars(i,j,k,15)*stepsMinusOne+delPrim[1]*delPrim[3])*osteps;
			coVars(i,j,k,16)  = (coVars(i,j,k,16)*stepsMinusOne+delPrim[2]*delPrim[3])*osteps;
			coVars(i,j,k,17)  = (coVars(i,j,k,17)*stepsMinusOne+delCon[0]*delPrim[11])*osteps;
			coVars(i,j,k,18)  = (coVars(i,j,k,18)*stepsMinusOne+delPrim[1]*delPrim[11])*osteps;
			coVars(i,j,k,19)  = (coVars(i,j,k,19)*stepsMinusOne+delPrim[2]*delPrim[11])*osteps;
			coVars(i,j,k,20)  = (coVars(i,j,k,20)*stepsMinusOne+delPrim[3]*delPrim[11])*osteps;
		});
	}
}

/*
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
}*/
