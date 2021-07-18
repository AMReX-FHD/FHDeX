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
			... (repeat for each species)
		*/

		int ncon  = 5;
		int nprim = 14;
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
			const long np = m_cell_vectors[0][grid_id][imap].size();
			// Read particle data
			for (int m=0; m<np; m++) {
				int pind = m_cell_vectors[0][grid_id][imap][m];
				ParticleType ptemp = particles[pind];
				ParticleType & p = ptemp;
				Real mass = p.rdata(FHD_realData::mass);
				Real u = p.rdata(FHD_realData::velx);
				Real v = p.rdata(FHD_realData::vely);
				Real w = p.rdata(FHD_realData::velz);

				cuInst(i,j,k,0) += mass;
				cuInst(i,j,k,1) += (mass*u);
				cuInst(i,j,k,2) += (mass*v);
				cuInst(i,j,k,3) += (mass*w);

				Real spdsq = pow(u,2)+pow(v,2)+pow(w,2);
				cuInst(i,j,k,4) += (mass*spdsq);
			}
			
			primInst(i,j,k,0) = np;
			cuInst(i,j,k,0) *= ocollisionCellVol;
			cuInst(i,j,k,1) *= ocollisionCellVol;
			cuInst(i,j,k,2) *= ocollisionCellVol;
			cuInst(i,j,k,3) *= ocollisionCellVol;
			cuInst(i,j,k,4) *= (0.5*ocollisionCellVol);
			
			primInst(i,j,k,0) *= ocollisionCellVol;
			primInst(i,j,k,1) = cuInst(i,j,k,0);
			
			// u,v,w
			primInst(i,j,k,2) = cuInst(i,j,k,1)/cuInst(i,j,k,0);
			primInst(i,j,k,3) = cuInst(i,j,k,2)/cuInst(i,j,k,0);
			primInst(i,j,k,4) = cuInst(i,j,k,3)/cuInst(i,j,k,0);
			
			// Second velocity moments
			primInst(i,j,k,5) = pow(primInst(i,j,k,2),2);	
			primInst(i,j,k,6) = primInst(i,j,k,2)*primInst(i,j,k,3);
			primInst(i,j,k,7) = primInst(i,j,k,2)*primInst(i,j,k,4);
			primInst(i,j,k,8) = pow(primInst(i,j,k,3),2);
			primInst(i,j,k,9) = primInst(i,j,k,3)*primInst(i,j,k,4);
			primInst(i,j,k,10) = pow(primInst(i,j,k,4),2);
			
			Real vsqb = (pow(primInst(i,j,k,2),2)+pow(primInst(i,j,k,3),2) +
				pow(primInst(i,j,k,4),2));
				
			Real mMass = cuInst(i,j,k,0)/primInst(i,j,k,0);
			Real mMassoV = mMass*ocollisionCellVol;
			Real R = k_B/mMass;
			// Naively calculate cv from average mass
			Real cv = 3.0*R*0.5;

			primInst(i,j,k,11) = (cuInst(i,j,k,4)/cuInst(i,j,k,0)-vsqb*0.5)/cv;
			primInst(i,j,k,12) = primInst(i,j,k,11)*R*cuInst(i,j,k,0);
			primInst(i,j,k,13) = vsqb*mMassoV+cv*primInst(i,j,k,11)*cuInst(i,j,k,0);

			//////////////////////////////////////
			// Means and Variances
			//////////////////////////////////////
			Vector<Real> delCon(ncon, 0.0);
			for (int l=0; l<ncon; l++) {
				cuMeans(i,j,k,l) = (cuMeans(i,j,k,l)*stepsMinusOne+cuInst(i,j,k,l))*osteps;
				delCon[l]        = cuInst(i,j,k,l) - cuMeans(i,j,k,l);
				cuVars(i,j,k,l)  = (cuVars(i,j,k,l)*stepsMinusOne+delCon[l]*delCon[l])*osteps;
			}

			primMeans(i,j,k,0)  = 0.;
			primMeans(i,j,k,1)  = cuMeans(i,j,k,0);
			primMeans(i,j,k,2)  = cuMeans(i,j,k,1)/cuMeans(i,j,k,0);
			primMeans(i,j,k,3)  = cuMeans(i,j,k,2)/cuMeans(i,j,k,0);
			primMeans(i,j,k,4)  = cuMeans(i,j,k,3)/cuMeans(i,j,k,0);
			primMeans(i,j,k,5)  = pow(primMeans(i,j,k,2),2);
			primMeans(i,j,k,6)  = primMeans(i,j,k,2)*primMeans(i,j,k,3);
			primMeans(i,j,k,7)  = primMeans(i,j,k,2)*primMeans(i,j,k,4);
			primMeans(i,j,k,8)  = pow(primMeans(i,j,k,3),2);
			primMeans(i,j,k,9)  = primMeans(i,j,k,3)*primMeans(i,j,k,4);
			primMeans(i,j,k,10) = pow(primMeans(i,j,k,4),2);

			vsqb = pow(primMeans(i,j,k,2),2)+pow(primMeans(i,j,k,3),2) +
				pow(primMeans(i,j,k,4),2);

			primMeans(i,j,k,11) = (cuMeans(i,j,k,4)/cuMeans(i,j,k,0)-vsqb*0.5)/cv;
			primMeans(i,j,k,12) = primMeans(i,j,k,11)*R*cuMeans(i,j,k,0);
			primMeans(i,j,k,13)  = 0.5*primMeans(i,j,k,1)*vsqb;
			primMeans(i,j,k,13)  = primMeans(i,j,k,13) + (cv*primMeans(i,j,k,11)*primMeans(i,j,k,1));

			// Primitive Variances
			Vector<Real> delPrim(nprim, 0.0);
			for (int l=0; l<nprim; l++) {
				// primMeans(i,j,k,l) = (primMeans(i,j,k,l)*stepsMinusOne+primInst(i,j,k,l))*osteps;
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
		/*Real Tgl[nspecies];
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
		fileTg.close(); fileTgN.close();*/
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
