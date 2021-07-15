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
			10 - ww  (ww_ns)
			11 - T   (P_ns)  = (1/3V) += m|v|^2
			12 - P   (T_ns)  = (1/3Nk_B) += m|v|^2
			13 - E   (E_ns)  = (1/2) += |v|^2 + c_v*T
			14 - qx  (qx_ns)
			15 - qy  (qy_ns)
			16 - qz  (qz_ns)
			... (repeat for each species)
		*/

		int ncon  = (nspecies+1)*5;
		int nprim = (nspecies+1)*17;
		Array4<Real> cuInst     = mfcuInst[pti].array();
		Array4<Real> primInst   = mfprimInst[pti].array();
		Array4<Real> cuMeans    = mfcuMeans[pti].array();
		Array4<Real> primMeans  = mfprimMeans[pti].array();
		Array4<Real> cuVars     = mfcuVars[pti].array();
		Array4<Real> primVars   = mfprimVars[pti].array();
		Array4<Real> coVars     = mfcoVars[pti].array();

		//////////////////////////////////////
		// Instantaneous Values
		//////////////////////////////////////

		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			const IntVect& iv = {i,j,k};
			long imap = tile_box.index(iv);
			int icon = 5; int iprim = 17;
			Real cv  = 0.;
			for (int l=0; l<nspecies; l++) {
				const long np_spec = m_cell_vectors[l][grid_id][imap].size();
				Real mass = properties[l].mass*properties[l].Neff;
				Real moV  = properties[l].mass*ocollisionCellVol;    	 // density
				Real cv_l = 3.0*k_B*0.5/mass;                        	 // assume ideal gas <- in general untrue for most granular cases
				primInst(i,j,k,iprim+0) = np_spec*ocollisionCellVol;	 // number density of species l
				primInst(i,j,k,0)	+= np_spec*ocollisionCellVol;
				primInst(i,j,k,iprim+1) = np_spec*moV;         		    // mass density of species l
				primInst(i,j,k,1)	+= np_spec*moV;
				cuInst(i,j,k,icon+0)	= np_spec*moV;
				cv			+= cv_l*np_spec*moV;				 // total cv

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

					primInst(i,j,k,iprim+14) += spdsq*u;	// qx
					primInst(i,j,k,iprim+15) += spdsq*v;	// qy
					primInst(i,j,k,iprim+16) += spdsq*w;	// qz
				}

				// Momentum densities
				cuInst(i,j,k,icon+1) *= moV;  		// x-mom density
				cuInst(i,j,k,icon+2) *= moV;  		// y-mom density
				cuInst(i,j,k,icon+3) *= moV;  		// z-mom density

				// Energy density
				cuInst(i,j,k,icon+4) *= (moV*0.5);

				// Total of cons. vars
				for (int m=0; m<5; m++) {cuInst(i,j,k,m) += cuInst(i,j,k,icon+m);}

				// Care not to use mean of instananeous primitives as mean (use conserved)
				primInst(i,j,k,iprim+ 2) = cuInst(i,j,k,icon+1)/cuInst(i,j,k,icon+0);				// u  of spec l
				primInst(i,j,k,iprim+ 3) = cuInst(i,j,k,icon+2)/cuInst(i,j,k,icon+0);				// v  of spec l
				primInst(i,j,k,iprim+ 4) = cuInst(i,j,k,icon+3)/cuInst(i,j,k,icon+0);				// w  of spec l
				primInst(i,j,k,iprim+ 5) = pow(primInst(i,j,k,iprim+2),2);									// uu of spec l
				primInst(i,j,k,iprim+ 6) = primInst(i,j,k,iprim+2)*primInst(i,j,k,iprim+3);	// uv of spec l
				primInst(i,j,k,iprim+ 7) = primInst(i,j,k,iprim+2)*primInst(i,j,k,iprim+4);	// uw of spec l
				primInst(i,j,k,iprim+ 8) = pow(primInst(i,j,k,iprim+3),2);									// vv of spec l
				primInst(i,j,k,iprim+ 9) = primInst(i,j,k,iprim+3)*primInst(i,j,k,iprim+4);	// vw of spec l
				primInst(i,j,k,iprim+10) = pow(primInst(i,j,k,iprim+4),2);									// ww of spec l

				Real vsqb = 0.5*(pow(primInst(i,j,k,iprim+2),2)+pow(primInst(i,j,k,iprim+3),2) +
								     pow(primInst(i,j,k,iprim+4),2));

				primInst(i,j,k,iprim+11) = (cuInst(i,j,k,icon+4)/cuInst(i,j,k,icon)-vsqb*0.5)/cv_l;		// T  of spec l
				primInst(i,j,k,iprim+12) = primInst(i,j,k,iprim+11)*(k_B/mass)*cuInst(i,j,k,icon);		// P  of spec l
				primInst(i,j,k,iprim+13) = vsqb*moV+cv_l*primInst(i,j,k,iprim+11)*cuInst(i,j,k,icon);	// E  of spec l

				primInst(i,j,k,iprim+14) = primInst(i,j,k,iprim+ 2)*cuInst(i,j,k,icon+4);	// qx
				primInst(i,j,k,iprim+15) = 0.5*moV;	// qy
				primInst(i,j,k,iprim+16) = 0.5*moV;	// qz

				// Total of primitive vars
				// Handle temperature/energy seperately
				primInst(i,j,k,11) += primInst(i,j,k,iprim+11)*primInst(i,j,k,iprim);		// Mixture T = sum (nk*Tk)/n

				primInst(i,j,k,12) += primInst(i,j,k,iprim+12);													// Mixture P = sum Pk
				primInst(i,j,k,14) += primInst(i,j,k,iprim+14);
				primInst(i,j,k,15) += primInst(i,j,k,iprim+15);
				primInst(i,j,k,16) += primInst(i,j,k,iprim+16);
				icon += 5; iprim += 17;
			}

			// Primitive total
			primInst(i,j,k,2) = cuInst(i,j,k,1)/cuInst(i,j,k,0);		// Bulk x-velocity
			primInst(i,j,k,3) = cuInst(i,j,k,2)/cuInst(i,j,k,0);		// Bulk y-velocity
			primInst(i,j,k,4) = cuInst(i,j,k,3)/cuInst(i,j,k,0);		// Bulk z-velocity

			primInst(i,j,k,5)  = pow(primInst(i,j,k,2),2);					// Bulk uu
			primInst(i,j,k,6)  = primInst(i,j,k,2)*primInst(i,j,k,3);	// Bulk uv
			primInst(i,j,k,7)  = primInst(i,j,k,2)*primInst(i,j,k,4);	// Bulk uw
			primInst(i,j,k,8)  = pow(primInst(i,j,k,3),2);					// Bulk vv
			primInst(i,j,k,9)  = primInst(i,j,k,3)*primInst(i,j,k,4);	// Bulk vw
			primInst(i,j,k,10) = pow(primInst(i,j,k,4),2);					// Bulk ww

			// Mixture Temperature
			primInst(i,j,k,11) /= primInst(i,j,k,0);

			// Energy Density
			cv /= primInst(i,j,k,1);
			primInst(i,j,k,13)  = pow(primInst(i,j,k,2),2)+pow(primInst(i,j,k,3),2)+pow(primInst(i,j,k,4),2);
			primInst(i,j,k,13)  = 0.5*primInst(i,j,k,1)*primInst(i,j,k,13);								// Bulk energy
			primInst(i,j,k,13)  = primInst(i,j,k,13) + (cv*primInst(i,j,k,11)*primInst(i,j,k,1));	// Total Particle KE

			// Convert n to Xk and rho to Yk
			/*
			for (int l=0; l<nspecies; l++) {
				primInst(i,j,k,nprimvars*l+0)  /= primInst(i,j,k,0); //X_k
				primInst(i,j,k,nprimvars*l+1)  /= primInst(i,j,k,0);  //Y_k
			}*/
		});

		//////////////////////////////////////
		// Mean Values and Variances
		//////////////////////////////////////

		amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
			Vector<Real> delCon(ncon, 0.0);
			// Conserved Variances
			for (int l=0; l<ncon; l++) {
				cuMeans(i,j,k,l) = (cuMeans(i,j,k,l)*stepsMinusOne+cuInst(i,j,k,l))*osteps;
				delCon[l]        = cuInst(i,j,k,l) - cuMeans(i,j,k,l);
				cuVars(i,j,k,l)  = (cuVars(i,j,k,l)*stepsMinusOne+delCon[l]*delCon[l])*osteps;
			}

			// Primitive Variances
			Vector<Real> delPrim(nprim, 0.0);
			for (int l=0; l<nprim; l++) {
				primMeans(i,j,k,l) = (primMeans(i,j,k,l)*stepsMinusOne+primInst(i,j,k,l))*osteps;
				delPrim[l]         = primInst(i,j,k,l) - primMeans(i,j,k,l);
				primVars(i,j,k,l)  = (primVars(i,j,k,l)*stepsMinusOne+delPrim[l]*delPrim[l])*osteps;
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
			coVars(i,j,k,11)  = (coVars(i,j,k,11)*stepsMinusOne+delCon[0]*delPrim[2])*osteps;
			coVars(i,j,k,12)  = (coVars(i,j,k,12)*stepsMinusOne+delCon[0]*delPrim[3])*osteps;
			coVars(i,j,k,13)  = (coVars(i,j,k,13)*stepsMinusOne+delCon[0]*delPrim[4])*osteps;
			coVars(i,j,k,14)  = (coVars(i,j,k,14)*stepsMinusOne+delPrim[2]*delPrim[3])*osteps;
			coVars(i,j,k,15)  = (coVars(i,j,k,15)*stepsMinusOne+delPrim[2]*delPrim[4])*osteps;
			coVars(i,j,k,16)  = (coVars(i,j,k,16)*stepsMinusOne+delPrim[3]*delPrim[4])*osteps;
			coVars(i,j,k,17)  = (coVars(i,j,k,17)*stepsMinusOne+delCon[0]*delPrim[11])*osteps;
			coVars(i,j,k,18)  = (coVars(i,j,k,18)*stepsMinusOne+delPrim[2]*delPrim[11])*osteps;
			coVars(i,j,k,19)  = (coVars(i,j,k,19)*stepsMinusOne+delPrim[3]*delPrim[11])*osteps;
			coVars(i,j,k,20)  = (coVars(i,j,k,20)*stepsMinusOne+delPrim[4]*delPrim[11])*osteps;
		});

		// Global Granular Temperature
		/*
		Real Tgl[nspecies];
		Real npl[nspecies];
		for (int l=0; l<nspecies; l++) {Tgl[l] = 0.; npl[l] = 0.;}
		int np = particles.numParticles();

		for (int i = 0; i < np; ++i) {
			ParticleType & part = particles[i];

			int ispec = part.idata(FHD_intData::species);
			Real vsq = pow(part.rdata(FHD_realData::velx),2) +
				pow(part.rdata(FHD_realData::vely),2) +
				pow(part.rdata(FHD_realData::velz),2);
			vsq	    *= (properties[ispec].mass/3.0);
			Tgl[ispec] += vsq; npl[ispec] += 1;
		}

		// Gather from all proc
		for (int l=0; l<nspecies; l++) {
			Real tempTg = Tgl[l];
			Real tempnp = npl[l];
			ParallelDescriptor::ReduceRealSum(tempTg);
			ParallelDescriptor::ReduceRealSum(tempnp);
			npl[l] = tempnp;
			Tgl[l] = tempTg/tempnp;
		}

		// Print to files
		if (ParallelDescriptor::IOProcessor()) {
			ofstream fileTg, fileTgN;
			std::string Tgfname = "Tg.dat";
			std::string TgNfname = "TgN.dat";
			if(steps==1) {
			fileTg.open(Tgfname);	fileTgN.open(TgNfname);
		} else {
			fileTg.open(Tgfname, fstream::app);
			fileTgN.open(TgNfname, fstream::app);
		}
		fileTg << std::scientific << setprecision(8) << time << " ";
		fileTgN << std::scientific << setprecision(8) << time << " ";
		for(int l=0; l<nspecies; l++){
			if(steps==1) {Tg0[l] = Tgl[l];}
			fileTg << Tgl[l] << " ";
			fileTgN << Tgl[l]/Tg0[l] << " ";
		}
		fileTg << "\n"; fileTgN << "\n";
		fileTg.close(); fileTgN.close();
	}*/
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
