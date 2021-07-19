#include "DsmcParticleContainer.H"
using namespace std;
void FhdParticleContainer::EvaluateStats(MultiFab& mfcuInst,
						MultiFab& mfcuMeans,
						MultiFab& mfcuVars,
						MultiFab& mfprimInst,
						MultiFab& mfprimMeans,
						MultiFab& mfprimVars,
						MultiFab& mfcoVars,
						int steps,
						Real time) {
	BL_PROFILE_VAR("EvaluateStats()",EvaluateStats);
	const Real osteps = 1.0/steps;
	const int stepsMinusOne = steps-1;

	// TODO: Add Heat Fluxes
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
		*/

		/*
		   Primitive Vars:
			0	 - n   (n_ns)
			1  - rho (Y_ns)
			2  - u   (u_ns)
			3  - v   (v_ns)
			4  - w   (w_ns)
			5  - G   (G_ns) = dot(u_mean,dJ)
			6  - T   (T_ns)
			7  - P   (P_ns)
			8  - E   (E_ns)
		*/


		int ncon  = (nspecies+1)*5;
		int nprim = (nspecies+1)*9;
		Array4<Real> cuInst     = mfcuInst[pti].array();
		Array4<Real> primInst   = mfprimInst[pti].array();
		Array4<Real> cuMeans    = mfcuMeans[pti].array();
		Array4<Real> primMeans  = mfprimMeans[pti].array();
		Array4<Real> cuVars     = mfcuVars[pti].array();
		Array4<Real> primVars   = mfprimVars[pti].array();
		Array4<Real> coVars     = mfcoVars[pti].array();

		//////////////////////////////////////
		// Primitve and Conserved Instantaneous Values
		//////////////////////////////////////

		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);
			int icon = 5; int iprim = 9;
			Real cv  = 0.;

			for (int l=0; l<nspecies; l++) {
				const long np_spec = m_cell_vectors[l][grid_id][imap].size();
				Real mass = properties[l].mass*properties[l].Neff;
				Real moV  = properties[l].mass*ocollisionCellVol;
				Real cv_l = 3.0*k_B*0.5/mass;
				primInst(i,j,k,iprim+0) = np_spec*ocollisionCellVol;
				primInst(i,j,k,0)      += np_spec*ocollisionCellVol;
				primInst(i,j,k,iprim+1) = np_spec*moV;
				primInst(i,j,k,1)	     += np_spec*moV;
				cuInst(i,j,k,icon+0)	  = np_spec*moV;
				cv			               += cv_l*np_spec*moV;

				// Read particle data
				for (int m=0; m<np_spec; m++) {
					int pind = m_cell_vectors[l][grid_id][imap][m];
					ParticleType ptemp = particles[pind];
					ParticleType & p = ptemp;
					// ParticleType & p = particles[pind];
					Real u = p.rdata(FHD_realData::velx);
					Real v = p.rdata(FHD_realData::vely);
					Real w = p.rdata(FHD_realData::velz);

					cuInst(i,j,k,icon+1) += u;
					cuInst(i,j,k,icon+2) += v;
					cuInst(i,j,k,icon+3) += w;

					Real spdsq = pow(u,2)+pow(v,2)+pow(w,2);
					cuInst(i,j,k,icon+4) += spdsq;
				}

				cuInst(i,j,k,icon+1) *= moV;       // x-mom density
				cuInst(i,j,k,icon+2) *= moV;       // y-mom density
				cuInst(i,j,k,icon+3) *= moV;       // z-mom density
				cuInst(i,j,k,icon+4) *= (moV*0.5); // K     density

				// Total Conserved Vars
				for (int m=0; m<5; m++) {cuInst(i,j,k,m) += cuInst(i,j,k,icon+m);}

				primInst(i,j,k,iprim+ 2) = cuInst(i,j,k,icon+1)/cuInst(i,j,k,icon+0);				// u_l
				primInst(i,j,k,iprim+ 3) = cuInst(i,j,k,icon+2)/cuInst(i,j,k,icon+0);				// v_l
				primInst(i,j,k,iprim+ 4) = cuInst(i,j,k,icon+3)/cuInst(i,j,k,icon+0);				// w_l

				Real vsqb = (pow(primInst(i,j,k,iprim+2),2)+pow(primInst(i,j,k,iprim+3),2) +
								     pow(primInst(i,j,k,iprim+4),2));

				primInst(i,j,k,iprim+ 5) = primInst(i,j,k,iprim+2)*cuInst(i,j,k,1) +                  //G_l
					primInst(i,j,k,iprim+3)*cuInst(i,j,k,2) +
					primInst(i,j,k,iprim+4)*cuInst(i,j,k,3);
				primInst(i,j,k,iprim+ 6) = (cuInst(i,j,k,icon+4)/cuInst(i,j,k,icon)-vsqb*0.5)/cv_l;   // T_l
				primInst(i,j,k,6) += primInst(i,j,k,iprim+6)*primInst(i,j,k,iprim);
				
				primInst(i,j,k,iprim+ 7) = primInst(i,j,k,iprim+6)*(k_B/mass)*cuInst(i,j,k,icon);     // P_l
				primInst(i,j,k,7) += primInst(i,j,k,iprim+7);
				
				primInst(i,j,k,iprim+ 8) = vsqb*moV+cv_l*primInst(i,j,k,iprim+6)*cuInst(i,j,k,icon);  // E_l
				primInst(i,j,k,8) += primInst(i,j,k,iprim+8);

				icon += 5; iprim += 9;
			}

			iprim = 9;
			for (int l=0; l<nspecies; l++) {
				primInst(i,j,k,iprim+0) /= primInst(i,j,k,0);
				iprim += 9;
			}

			primInst(i,j,k,2) = cuInst(i,j,k,1)/cuInst(i,j,k,0);  // Total x-velocity
			primInst(i,j,k,3) = cuInst(i,j,k,2)/cuInst(i,j,k,0);  // Total y-velocity
			primInst(i,j,k,4) = cuInst(i,j,k,3)/cuInst(i,j,k,0);  // Total z-velocity
			primInst(i,j,k,5) = primInst(i,j,k,2)*cuInst(i,j,k,1) + 
				primInst(i,j,k,3)*cuInst(i,j,k,2) +
				primInst(i,j,k,4)*cuInst(i,j,k,3);                  // G
			primInst(i,j,k,6) /= primInst(i,j,k,0);               // Mixture Temperature

			// Energy Density
			cv /= primInst(i,j,k,1);
			primInst(i,j,k,8)  = pow(primInst(i,j,k,2),2)+pow(primInst(i,j,k,3),2)+pow(primInst(i,j,k,4),2);
			primInst(i,j,k,8)  = 0.5*primInst(i,j,k,1)*primInst(i,j,k,8);								        // Bulk energy
			primInst(i,j,k,8)  = primInst(i,j,k,8) + (cv*primInst(i,j,k,6)*primInst(i,j,k,1));	// Total Particle KE

		});

		//////////////////////////////////////
		// Means
		//////////////////////////////////////
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			for (int l=0; l<ncon; l++) {
				cuMeans(i,j,k,l) = (cuMeans(i,j,k,l)*stepsMinusOne+cuInst(i,j,k,l))*osteps;
			}

			// Evaluate Primitive Means from Conserved Means
			primMeans(i,j,k,0)  = 0.;
			primMeans(i,j,k,1)  = cuMeans(i,j,k,0);
			primMeans(i,j,k,2)  = cuMeans(i,j,k,1)/cuMeans(i,j,k,0);
			primMeans(i,j,k,3)  = cuMeans(i,j,k,2)/cuMeans(i,j,k,0);
			primMeans(i,j,k,4)  = cuMeans(i,j,k,3)/cuMeans(i,j,k,0);

			// Zero out hydrodynamic means (derive from conserved means)
			primMeans(i,j,k,5) = 0.0;
			primMeans(i,j,k,6) = 0.0;
			primMeans(i,j,k,7) = 0.0;
			primMeans(i,j,k,8) = 0.0;
			int iprim = 9; int icon = 5;
			Real cv = 0.;
			for(int l=0; l<nspecies; l++) { 
				Real mass = properties[l].mass*properties[l].Neff;
				Real moV  = properties[l].mass*ocollisionCellVol;
				Real cv_l = 3.0*k_B*0.5/mass;

				primMeans(i,j,k,iprim+0) = cuMeans(i,j,k,icon+0)/mass;     // n_l
				primMeans(i,j,k,0)      += primMeans(i,j,k,iprim+0)/mass;  // n
				primMeans(i,j,k,iprim+1) = cuMeans(i,j,k,icon+0);          // rho_l

				primMeans(i,j,k,iprim+2) = cuMeans(i,j,k,icon+1)/cuMeans(i,j,k,icon+0); // u_l
				primMeans(i,j,k,iprim+3) = cuMeans(i,j,k,icon+2)/cuMeans(i,j,k,icon+0); // v_l
				primMeans(i,j,k,iprim+4) = cuMeans(i,j,k,icon+3)/cuMeans(i,j,k,icon+0); // w_l

				cv += cv_l*primMeans(i,j,k,iprim+1);
				Real vsqb = pow(primMeans(i,j,k,iprim+2),2)+pow(primMeans(i,j,k,iprim+3),2) +
								     pow(primMeans(i,j,k,iprim+4),2);

				primMeans(i,j,k,iprim+6) = (cuMeans(i,j,k,icon+4)/cuMeans(i,j,k,icon+0)-vsqb*0.5)/cv_l;
				primMeans(i,j,k,iprim+7) = primMeans(i,j,k,iprim+6)*(k_B/mass)*cuMeans(i,j,k,icon+0);
				primMeans(i,j,k,7)      += primMeans(i,j,k,iprim+7);
				primMeans(i,j,k,iprim+8) = vsqb*moV+cv_l*primMeans(i,j,k,iprim+6)*cuMeans(i,j,k,icon+0);
				
				iprim += 9; icon += 5;
			}

			primMeans(i,j,k,5) = primMeans(i,j,k,2)*cuMeans(i,j,k,1) +
				primMeans(i,j,k,3)*cuMeans(i,j,k,2) +
				primMeans(i,j,k,4)*cuMeans(i,j,k,3);
			cv /= primMeans(i,j,k,1);
			
			Real vsqb = pow(primMeans(i,j,k,2),2) + pow(primMeans(i,j,k,3),2) + pow(primMeans(i,j,k,4),2);
			primMeans(i,j,k,6) = (cuMeans(i,j,k,4)/cuMeans(i,j,k,0) - 0.5*vsqb)/cv;

			primMeans(i,j,k,8)  = pow(primMeans(i,j,k,2),2)+pow(primMeans(i,j,k,3),2)+pow(primMeans(i,j,k,4),2);
			primMeans(i,j,k,8)  = 0.5*primMeans(i,j,k,1)*primMeans(i,j,k,8);
			primMeans(i,j,k,8)  = primMeans(i,j,k,8) + (cv*primMeans(i,j,k,6)*primMeans(i,j,k,1));
		});

		//////////////////////////////////////
		// Variances
		//////////////////////////////////////
		// Covariances
		/*
			// Conserved
			0  - drho.dJx
			1  - drho.dJy
			2  - drho.dJz
			3  - drho.dK
			4  - dJx.dJy
			5  - dJx.dJz
			6  - dJx.dK
			7  - dJy.dJz
			8  - dJy.dK
			9  - dJz.dk
			
			// Energy
			10 - drho.dG
			11 - dJx.dG
			12 - dJy.dG
			13 - dJz.dG
			14 - dK.dG
			
			// Hydro
			15 - drho.du
			16 - drho.dv
			17 - drho.dw
			18 - du.dv
			19 - du.dw
			20 - dv.dw
			21 - drho.dT
			22 - du.dT
			23 - dv.dT
			24 - dw.dT
		*/
		
		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			// Conserved Variances
			Vector<Real> delCon(ncon, 0.0);
			for (int l=0; l<ncon; l++) {
				delCon[l]        = cuInst(i,j,k,l) - cuMeans(i,j,k,l);
				cuVars(i,j,k,l)  = (cuVars(i,j,k,l)*stepsMinusOne+delCon[l]*delCon[l])*osteps;
			}
			
			//Conserved Covariances
			coVars(i,j,k, 0)  = (coVars(i,j,k, 0)*stepsMinusOne+delCon[0]*delCon[1])*osteps; // drho.dJx
			coVars(i,j,k, 1)  = (coVars(i,j,k, 1)*stepsMinusOne+delCon[0]*delCon[2])*osteps; // drho.dJy
			coVars(i,j,k, 2)  = (coVars(i,j,k, 2)*stepsMinusOne+delCon[0]*delCon[3])*osteps; // drho.dJz
			coVars(i,j,k, 3)  = (coVars(i,j,k, 3)*stepsMinusOne+delCon[0]*delCon[4])*osteps; // drho.dK
			coVars(i,j,k, 4)  = (coVars(i,j,k, 4)*stepsMinusOne+delCon[1]*delCon[2])*osteps; // dJx.dJy
			coVars(i,j,k, 5)  = (coVars(i,j,k, 5)*stepsMinusOne+delCon[1]*delCon[3])*osteps; // dJx.dJz
			coVars(i,j,k, 6)  = (coVars(i,j,k, 6)*stepsMinusOne+delCon[1]*delCon[4])*osteps; // dJx.dK
			coVars(i,j,k, 7)  = (coVars(i,j,k, 7)*stepsMinusOne+delCon[2]*delCon[3])*osteps; // dJy.dJz
			coVars(i,j,k, 8)  = (coVars(i,j,k, 8)*stepsMinusOne+delCon[2]*delCon[4])*osteps; // dJy.dK
			coVars(i,j,k, 9)  = (coVars(i,j,k, 9)*stepsMinusOne+delCon[3]*delCon[4])*osteps; // dJz.dK

			// Primitive Variances
			Real odensity = 1.0/cuMeans(i,j,k,0);
			Real deln = primInst(i,j,k,0) - primMeans(i,j,k,0);
			primVars(i,j,k, 0) = (primVars(i,j,k,0)*stepsMinusOne+deln*deln)*osteps;         // dn.dn
			Real delrho = delCon[1];
			primVars(i,j,k, 1) = cuVars(i,j,k,0);                                            // drho.drho
			
			Real dudu = pow(odensity,2.0)*(cuVars(i,j,k,1) - 2.0*primMeans(i,j,k,2)*coVars(i,j,k,0)+pow(primMeans(i,j,k,2),2)*cuVars(i,j,k,0));
			primVars(i,j,k, 2) = (primVars(i,j,k,2)*stepsMinusOne+dudu)*osteps;
			Real dvdv = pow(odensity,2.0)*(cuVars(i,j,k,2) - 2.0*primMeans(i,j,k,3)*coVars(i,j,k,1)+pow(primMeans(i,j,k,3),2)*cuVars(i,j,k,0));
			primVars(i,j,k, 3) = (primVars(i,j,k,3)*stepsMinusOne+dvdv)*osteps;
			Real dwdw = pow(odensity,2.0)*(cuVars(i,j,k,3) - 2.0*primMeans(i,j,k,4)*coVars(i,j,k,2)+pow(primMeans(i,j,k,4),2)*cuVars(i,j,k,0));
			primVars(i,j,k, 4) = (primVars(i,j,k,4)*stepsMinusOne+dwdw)*osteps;
			
			Real dG = primMeans(i,j,k,2)*delCon[1]+primMeans(i,j,k,3)*delCon[2]+primMeans(i,j,k,4)*delCon[3];
			primVars(i,j,k, 5) = (primVars(i,j,k,5)*stepsMinusOne+dG*dG)*osteps;                 // dG.dG

			coVars(i,j,k,10)   = (coVars(i,j,k,10)*stepsMinusOne+delCon[0]*dG)*osteps;           // drho.dG
			coVars(i,j,k,11)   = (coVars(i,j,k,11)*stepsMinusOne+delCon[1]*dG)*osteps;           // dJx.dG
			coVars(i,j,k,12)   = (coVars(i,j,k,12)*stepsMinusOne+delCon[2]*dG)*osteps;           // dJy.dG
			coVars(i,j,k,13)   = (coVars(i,j,k,13)*stepsMinusOne+delCon[3]*dG)*osteps;           // dJz.dG
			coVars(i,j,k,14)   = (coVars(i,j,k,14)*stepsMinusOne+delCon[4]*dG)*osteps;           // dK.dG

			Real drhodu = odensity*(coVars(i,j,k,0) - primMeans(i,j,k,2)*cuVars(i,j,k,0));
			coVars(i,j,k,15)   = (coVars(i,j,k,15)*stepsMinusOne+drhodu)*osteps;                 // drho.du
			Real drhodv = odensity*(coVars(i,j,k,1) - primMeans(i,j,k,3)*cuVars(i,j,k,0));
			coVars(i,j,k,16)   = (coVars(i,j,k,16)*stepsMinusOne+drhodv)*osteps;                 // drho.dv
			Real drhodw = odensity*(coVars(i,j,k,2) - primMeans(i,j,k,4)*cuVars(i,j,k,0));
			coVars(i,j,k,17)   = (coVars(i,j,k,17)*stepsMinusOne+drhodw)*osteps;                 // drho.dw

			Real dudv = pow(odensity,2.0)*(coVars(i,j,k,4) - primMeans(i,j,k,2)*coVars(i,j,k,1)
				- primMeans(i,j,k,3)*coVars(i,j,k,0) + primMeans(i,j,k,2)*primMeans(i,j,k,3)*cuVars(i,j,k,0));
			coVars(i,j,k,18) = (coVars(i,j,k,18)*stepsMinusOne+dudv)*osteps;                     // du.dv
			Real dudw = pow(odensity,2.0)*(coVars(i,j,k,5) - primMeans(i,j,k,2)*coVars(i,j,k,2)
				- primMeans(i,j,k,4)*coVars(i,j,k,0) + primMeans(i,j,k,2)*primMeans(i,j,k,4)*cuVars(i,j,k,0));
			coVars(i,j,k,19) = (coVars(i,j,k,19)*stepsMinusOne+dudw)*osteps;                     // du.dw
			Real dvdw = pow(odensity,2.0)*(coVars(i,j,k,7) - primMeans(i,j,k,4)*coVars(i,j,k,2)
				- primMeans(i,j,k,3)*coVars(i,j,k,1) + primMeans(i,j,k,3)*primMeans(i,j,k,4)*cuVars(i,j,k,0));
			coVars(i,j,k,20) = (coVars(i,j,k,20)*stepsMinusOne+dvdw)*osteps;                     // dv.dw


			Real cv = 0.;
			int iprim = 9;
			for(int l=0; l<nspecies; l++) { 
				Real mass = properties[l].mass*properties[l].Neff;
				Real cv_l = 3.0*k_B*0.5/mass;
				cv       += cv_l*primMeans(i,j,k,iprim+1);
				iprim += 9;
			}
			cv /= primMeans(i,j,k,1);
			
			Real Qbar = cv*primMeans(i,j,k,6) - 0.5*(pow(primMeans(i,j,k,2),2)+pow(primMeans(i,j,k,3),2)+pow(primMeans(i,j,k,4),2));
			Real dTdT = pow(odensity/cv,2.0)*(cuVars(i,j,k,4) - primVars(i,j,k,5) + pow(Qbar,2)*cuVars(i,j,k,0)
				- 2.0*coVars(i,j,k,14) - 2.0*Qbar*coVars(i,j,k,3) + 2.0*Qbar*coVars(i,j,k,10));
			primVars(i,j,k, 6) = (primVars(i,j,k,6)*stepsMinusOne+dTdT)*osteps;                  // dT.dT
			
			Real drhodT = odensity/cv*(coVars(i,j,k,3) - coVars(i,j,k,10) - Qbar*cuVars(i,j,k,0));
			coVars(i,j,k,21) = (coVars(i,j,k,21)*stepsMinusOne+drhodT)*osteps;                   // drho.dT
			
			Real dudT = pow(odensity,2.0)/cv*(coVars(i,j,k,6) - primMeans(i,j,k,2)*coVars(i,j,k,3)
				- coVars(i,j,k,11) + primMeans(i,j,k,2)*coVars(i,j,k,10) - Qbar*coVars(i,j,k,0)
				+ primMeans(i,j,k,2)*Qbar*cuVars(i,j,k,0));
			coVars(i,j,k,22) = (coVars(i,j,k,22)*stepsMinusOne+dudT)*osteps;                     // du.dT
			
			Real dvdT = pow(odensity,2.0)/cv*(coVars(i,j,k,8) - primMeans(i,j,k,3)*coVars(i,j,k,3)
				- coVars(i,j,k,12) + primMeans(i,j,k,2)*coVars(i,j,k,10) - Qbar*coVars(i,j,k,1)
				+ primMeans(i,j,k,3)*Qbar*cuVars(i,j,k,0));
			coVars(i,j,k,23) = (coVars(i,j,k,23)*stepsMinusOne+dvdT)*osteps;                     // dv.dT
			
			Real dwdT = pow(odensity,2.0)/cv*(coVars(i,j,k,9) - primMeans(i,j,k,4)*coVars(i,j,k,3)
				- coVars(i,j,k,13) + primMeans(i,j,k,2)*coVars(i,j,k,10) - Qbar*coVars(i,j,k,2)
				+ primMeans(i,j,k,4)*Qbar*cuVars(i,j,k,0));
			coVars(i,j,k,24) = (coVars(i,j,k,24)*stepsMinusOne+dwdT)*osteps;                     // dw.dT
			
			// TODO: Add P and E variances
			// TODO: Add variances by species
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
