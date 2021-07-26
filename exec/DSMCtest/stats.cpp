#include "DsmcParticleContainer.H"
using namespace std;
void FhdParticleContainer::EvaluateStats(MultiFab& mfcuInst,
						MultiFab& mfcuMeans,
						MultiFab& mfcuVars,
						MultiFab& mfprimInst,
						MultiFab& mfprimMeans,
						MultiFab& mfprimVars,
						MultiFab& mfcoVars,
            MultiFab& spatialCross1D, 
						const int steps,
						Real time) {
    BL_PROFILE_VAR("EvaluateStats()",EvaluateStats);
    const Real osteps = 1.0/steps;
    const Real stepsMinusOne = steps-1.;

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

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept 
        {
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
        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept 
        {
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
        
        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
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


    // Spatial Correlations (works only for 1D currently: 1 cell each in y and z directions)
    // in this order: conserved variables [rho, K, jx, jy, jz] + primitive variables [vx, vy, vz, T]
    // total: 2*9
    // add species later
    if (plot_cross) {
        int nstats = 18;

        // Get all nstats at xcross and store in GpuVector
        amrex::Gpu::ManagedVector<Real> data_xcross_in(nstats, 0.0); // values at x*
        for ( MFIter mfi(mfcuInst); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.validbox();

            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);

            const Array4<const Real> cumeans   = mfcuMeans.array(mfi);
            const Array4<const Real> primmeans = mfprimMeans.array(mfi);
            const Array4<const Real> prim      = mfprimInst.array(mfi);
            const Array4<const Real> cu        = mfcuInst.array(mfi);

            Real* data_xcross = data_xcross_in.data();
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i==cross_cell) {
                    data_xcross[0]  = cu(i,j,k,0);        // rho-instant
                    data_xcross[1]  = cumeans(i,j,k,0);   // rho-mean
                    data_xcross[2]  = cu(i,j,k,4);        // energy-instant
                    data_xcross[3]  = cumeans(i,j,k,4);   // energy-mean
                    data_xcross[4]  = cu(i,j,k,1);        // jx-instant
                    data_xcross[5]  = cumeans(i,j,k,1);   // jx-mean
                    data_xcross[6]  = cu(i,j,k,2);        // jy-instant
                    data_xcross[7]  = cumeans(i,j,k,2);   // jy-mean
                    data_xcross[8]  = cu(i,j,k,3);        // jz-instant
                    data_xcross[9]  = cumeans(i,j,k,3);   // jz-mean
                    data_xcross[10] = prim(i,j,k,2);      // velx-instant
                    data_xcross[11] = primmeans(i,j,k,2); // velx-mean
                    data_xcross[12] = prim(i,j,k,3);      // vely-instant
                    data_xcross[13] = primmeans(i,j,k,3); // vely-mean
                    data_xcross[14] = prim(i,j,k,4);      // velz-instant
                    data_xcross[15] = primmeans(i,j,k,4); // velz-mean
                    data_xcross[16] = prim(i,j,k,6);      // T-instant
                    data_xcross[17] = primmeans(i,j,k,6); // T-mean
               }
            }); // end MFITer
        }

        // Reduce across MPI processes
        ParallelDescriptor::ReduceRealSum(data_xcross_in.data(),nstats);

        // Update Spatial Correlations
        for ( MFIter mfi(mfcuInst); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.validbox();

            const Array4<const Real> cumeans   = mfcuMeans.array(mfi);
            const Array4<const Real> primmeans = mfprimMeans.array(mfi);
            const Array4<const Real> prim      = mfprimInst.array(mfi);
            const Array4<const Real> cu        = mfcuInst.array(mfi);

            const Array4<      Real> spatialCross = spatialCross1D.array(mfi);
            
            Real* data_xcross = data_xcross_in.data();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {

                //////////////////////////////////////////////
                // Get fluctuations at xcross for this j and k
                //////////////////////////////////////////////
                Real meanrhocross = data_xcross[1];
                Real meanuxcross  = data_xcross[11];

                Real delrhocross = data_xcross[0] - data_xcross[1];
                Real delKcross   = data_xcross[2] - data_xcross[3];
                Real deljxcross  = data_xcross[4] - data_xcross[5];
                Real deljycross  = data_xcross[6] - data_xcross[7];
                Real deljzcross  = data_xcross[8] - data_xcross[9];

                Real delTcross = data_xcross[16] - data_xcross[17];
                Real delvxcross = data_xcross[10] - data_xcross[11];
                
                // [IS] Andrew: how do we cv here? We have precomputed previously
                // we should use it
                //Real cvcross = 0.;
                //for (int l=0; l<nspecies; ++l) {
                //    cvcross = cvcross + hcv[l]*data_xcross[18+4*l+1]/data_xcross[1];
                //}

                // Real cvinvcross = 1.0/cvcross;
                Real vxmeancross = data_xcross[11];
                Real vymeancross = data_xcross[13];
                Real vzmeancross = data_xcross[15];
                // Real qmeancross = cvcross*data_xcross[17] - 
                                 0.5*(vxmeancross*vxmeancross + vymeancross*vymeancross + vzmeancross*vzmeancross); 

                // delG = \vec{v}\cdot\vec{\deltaj}
                Real delGcross = vxmeancross*deljxcross + vymeancross*deljycross + vzmeancross*deljzcross;

                ////////////////////////////////////////
                // Get fluctuations at this cell (i,j,k)
                ////////////////////////////////////////
                
                Real meanrho = cumeans(i,j,k,0);

                Real delrho = cu(i,j,k,0) - cumeans(i,j,k,0);
                Real delK   = cu(i,j,k,4) - cumeans(i,j,k,4);
                Real deljx  = cu(i,j,k,1) - cumeans(i,j,k,1);
                Real deljy  = cu(i,j,k,2) - cumeans(i,j,k,2);
                Real deljz  = cu(i,j,k,3) - cumeans(i,j,k,3);

                Real delT   = prim(i,j,k,6) - primmeans(i,j,k,6);
                Real delvx  = prim(i,j,k,2) - primmeans(i,j,k,2);
                
                // [IS] Andrew: how do we cv here? We have precomputed previously
                // we should use it
                //Real cv = 0.;
                //for (int l=0; l<nspecies; ++l) {
                //    cv = cv + hcv[l]*cumeans(i,j,k,5+l)/cumeans(i,j,k,0);
                //}
                //Real cvinv = 1.0/cv;

                Real vxmean = primmeans(i,j,k,2);
                Real vymean = primmeans(i,j,k,3);
                Real vzmean = primmeans(i,j,k,4);

                // Real qmean = cv*primmeans(i,j,k,4) - 0.5*(vxmean*vxmean + vymean*vymean + vzmean*vzmean);

                // delG = \vec{v}\cdot\vec{\deltaj}
                Real delG = vxmean*deljx + vymean*deljy +vzmean*deljz;

                // Spatial Correlation Calculations
                
                spatialCross(i,j,k,0) = (spatialCross(i,j,k,0)*stepsMinusOne + delrhocross*delrho)*osteps; // <delrho(x*)delrho(x)>
                spatialCross(i,j,k,1) = (spatialCross(i,j,k,1)*stepsMinusOne + delKcross*delK)*osteps;     // <delK(x*)delK(x)>
                spatialCross(i,j,k,2) = (spatialCross(i,j,k,2)*stepsMinusOne + deljxcross*deljx)*osteps;   // <deljx(x*)deljx(x)>
                spatialCross(i,j,k,3) = (spatialCross(i,j,k,3)*stepsMinusOne + deljycross*deljy)*osteps;   // <deljy(x*)deljy(x)>
                spatialCross(i,j,k,4) = (spatialCross(i,j,k,4)*stepsMinusOne + deljzcross*deljz)*osteps;   // <deljz(x*)deljz(x)>
                spatialCross(i,j,k,5) = (spatialCross(i,j,k,5)*stepsMinusOne + deljxcross*delrho)*osteps;  // <deljx(x*)delrho(x)>

                // [IS] Species-dependent stuff -- commented now -- add later
                //spatialCross(i,j,k,6) = (spatialCross(i,j,k,6)*stepsMinusOne + deljxcross*delrhoYk[0])*osteps;  // <deljx(x*)delrhoYkL(x)>
                //spatialCross(i,j,k,7) = (spatialCross(i,j,k,7)*stepsMinusOne + deljxcross*delrhoYk[nspecies-1])*osteps;  // <deljx(x*)delrhoYkH(x)>
                //spatialCross(i,j,k,8) = (spatialCross(i,j,k,8)*stepsMinusOne + delrhocross*delrhoYk[0])*osteps; // <delrho(x*)delrhoYkL(x)>
                //spatialCross(i,j,k,9) = (spatialCross(i,j,k,9)*stepsMinusOne + delrhocross*delrhoYk[nspecies-1])*osteps; // <delrho(x*)delrhoYkH(x)>
                //spatialCross(i,j,k,10)= (spatialCross(i,j,k,10)*stepsMinusOne + delrhoYkcross[0]*delrho)*osteps; // <delrhoYkL(x*)delrho(x)>
                //spatialCross(i,j,k,11)= (spatialCross(i,j,k,11)*stepsMinusOne + delrhoYkcross[nspecies-1]*delrho)*osteps; // <delrhoYkH(x*)delrho(x)>
                
                spatialCross(i,j,k,12) = (spatialCross(i,j,k,12)*stepsMinusOne + delGcross*delG)*osteps; // <delG(x*)delG(x)>
                spatialCross(i,j,k,13) = (spatialCross(i,j,k,13)*stepsMinusOne + delGcross*delK)*osteps; // <delG(x*)delK(x)>
                spatialCross(i,j,k,14) = (spatialCross(i,j,k,14)*stepsMinusOne + delKcross*delG)*osteps; // <delK(x*)delG(x)>
                spatialCross(i,j,k,15) = (spatialCross(i,j,k,15)*stepsMinusOne + delrhocross*delK)*osteps; // <delrho(x*)delK(x)>
                spatialCross(i,j,k,16) = (spatialCross(i,j,k,16)*stepsMinusOne + delKcross*delrho)*osteps; // <delK(x*)delrho(x)>
                spatialCross(i,j,k,17) = (spatialCross(i,j,k,17)*stepsMinusOne + delrhocross*delG)*osteps; // <delrho(x*)delG(x)>
                spatialCross(i,j,k,18) = (spatialCross(i,j,k,18)*stepsMinusOne + delGcross*delrho)*osteps; // <delG(x*)delrho(x)>

                // [IS] the next two require cv and heat flux -- compute later

                // <delT(x*)delT(x)> = (1/cv*/cv/<rho(x)>/<rho(x*)>)(<delK*delK> + <delG*delG> - <delG*delK> - <delK*delG> 
                //                      + <Q><Q*><delrho*delrho> - <Q*><delrho*delK> - <Q><delK*delrho> + <Q*><delrho*delG> + <Q><delG*delrho>)
                //spatialCross(i,j,k,19) = (cvinvcross*cvinv/(meanrhocross*meanrho))*
                //                            (spatialCross(i,j,k,1) + spatialCross(i,j,k,12) - spatialCross(i,j,k,13) - spatialCross(i,j,k,14)
                //                             + qmean*qmeancross*spatialCross(i,j,k,0) - qmeancross*spatialCross(i,j,k,15) - qmean*spatialCross(i,j,k,16)
                //                             + qmeancross*spatialCross(i,j,k,17) + qmean*spatialCross(i,j,k,18));

                // <delT(x*)delrho(x)> = (1/cv/<rho(x*)>)*(<delK*delrho> - <delG*delrho> - <Q*><delrhodelrho*>)
                // spatialCross(i,j,k,20) = (cvinvcross*meanrhocross)*(spatialCross(i,j,k,16) - spatialCross(i,j,k,18) - qmeancross*spatialCross(i,j,k,0));

                // <delu(x*)delrho> = (1/<rho(x*)>)*(<deljx(x*)delrho(x)> - <u(x*)><<delrho(x*)delrho(x)>) 
                spatialCross(i,j,k,21) = (1.0/meanrhocross)*(spatialCross(i,j,k,5) - meanuxcross*spatialCross(i,j,k,0));

                // [IS] Species-dependent stuff -- commented now -- add later
                // <delu(x*)del(rhoYkL)> = (1/<rho(x*)>)*(<deljx(x*)del(rhoYkL)> - <u(x*)><delrho(x*)del(rhoYkL)>)
                // spatialCross(i,j,k,22) = (1.0/meanrhocross)*(spatialCross(i,j,k,6) - meanuxcross*spatialCross(i,j,k,8));  

                // <delu(x*)del(rhoYkH)> = (1/<rho(x*)>)*(<deljx(x*)del(rhoYkH)> - <u(x*)><delrho(x*)del(rhoYkH)>)
                // spatialCross(i,j,k,23) = (1.0/meanrhocross)*(spatialCross(i,j,k,7) - meanuxcross*spatialCross(i,j,k,9));  

                // <delu(x*)del(YkL)> = (1/<rho(x*)>/<rho(x)>)*(<deljx(x*)del(rhoYkL) - <u(x*)><delrho(x*)del(rhoYkL)> 
                //                      - <YkL(x)><deljx(x*)delrho(x)> + <u(x*)><YkL(x)><delrho(x*)delrho(x)>)
                // spatialCross(i,j,k,24) = (1.0/(meanrho*meanrhocross))*(spatialCross(i,j,k,6) - meanuxcross*spatialCross(i,j,k,8)
                //                                                         - meanYk[0]*spatialCross(i,j,k,5) + meanuxcross*meanYk[0]*spatialCross(i,j,k,0));

                // <delu(x*)del(YkH)> = (1/<rho(x*)>/<rho(x)>)*(<deljx(x*)del(rhoYkH) - <u(x*)><delrho(x*)del(rhoYkH)> 
                //                      - <YkH(x)><deljx(x*)delrho(x)> + <u(x*)><YkH(x)><delrho(x*)delrho(x)>)
                // spatialCross(i,j,k,25) = (1.0/(meanrho*meanrhocross))*(spatialCross(i,j,k,7) - meanuxcross*spatialCross(i,j,k,9)
                //                                                         - meanYk[nspecies-1]*spatialCross(i,j,k,5) + meanuxcross*meanYk[nspecies-1]*spatialCross(i,j,k,0));

                // Direct -- <delT(x*)delT(x)>
                spatialCross(i,j,k,26) = (spatialCross(i,j,k,26)*stepsMinusOne + delTcross*delT)*osteps;

                // Direct -- <delT(x*)delrho(x)>
                spatialCross(i,j,k,27) = (spatialCross(i,j,k,27)*stepsMinusOne + delTcross*delrho)*osteps;

                // Direct -- <delT(x*)delu(x)>
                spatialCross(i,j,k,28) = (spatialCross(i,j,k,28)*stepsMinusOne + delTcross*delvx)*osteps;

                // Direct -- <delu(x*)delrho>
                spatialCross(i,j,k,29) = (spatialCross(i,j,k,29)*stepsMinusOne + delvxcross*delrho)*osteps;

                // [IS] Species-dependent stuff -- commented now -- add later
                // Direct -- <delu(x*)del(rhoYkL)
                // spatialCross(i,j,k,30) = (spatialCross(i,j,k,30)*stepsMinusOne + delvxcross*delYk[0])*osteps;

                // Direct -- <delu(x*)del(rhoYkH)
                // spatialCross(i,j,k,31) = (spatialCross(i,j,k,31)*stepsMinusOne + delvxcross*delYk[nspecies-1])*osteps;

                // Direct -- <delu(x*)del(YkL)
                // spatialCross(i,j,k,32) = (spatialCross(i,j,k,32)*stepsMinusOne + delvxcross*delrhoYk[0])*osteps;

                // Direct -- <delu(x*)del(YkH)
                // spatialCross(i,j,k,33) = (spatialCross(i,j,k,33)*stepsMinusOne + delvxcross*delrhoYk[nspecies-1])*osteps;

                // Direct <delYkL(x*)delYkL(x)>
                // spatialCross(i,j,k,34) = (spatialCross(i,j,k,34)*stepsMinusOne + delYkcross[0]*delYk[0])*osteps;

                // Direct <delYkH(x*)delYkH(x)>
                // spatialCross(i,j,k,35) = (spatialCross(i,j,k,35)*stepsMinusOne + delYkcross[nspecies-1]*delYk[nspecies-1])*osteps;

                // Direct <delYkL(x*)delYkH(x)>
                // spatialCross(i,j,k,36) = (spatialCross(i,j,k,36)*stepsMinusOne + delYkcross[0]*delYk[nspecies-1])*osteps;

                // Last we rhoYk for species
                //for (int ns=0; ns<nspecies; ++ns) {
                //    spatialCross(i,j,k,37+ns) = (spatialCross(i,j,k,37+ns)*stepsMinusOne + delrhoYkcross[ns]*delrhoYk[ns])*osteps; // <delrhoYk(x*)delrhoYk(x)>
                //}

                // <delYkL(x*)delYkL(x)> = (1/<rho(x*)>/<rho(x)>)*(<delrhoYkL(x*)delrhoYkL> - <YkL(x*)><delrho(x*)delrhoYkL(x)>
                //                                                 - <YkL(x)><delrhoYkL(x*)delrho(x) + <YkL(x*)><YkL(x)><delrho(x*)delrho(x)>)
                //Real delrhoYkdelrhoYk = (spatialCross(i,j,k,37)*stepsMinusOne + delrhoYkcross[0]*delrhoYk[0])*osteps;
                //spatialCross(i,j,k,37+nspecies) = (1.0/(meanrho*meanrhocross))*(delrhoYkdelrhoYk - meanYkcross[0]*spatialCross(i,j,k,8)
                //                                                        - meanYk[0]*spatialCross(i,j,k,10) + meanYkcross[0]*meanYk[0]*spatialCross(i,j,k,0));

                // <delYkH(x*)delYkH(x)> = (1/<rho(x*)>/<rho(x)>)*(<delrhoYkH(x*)delrhoYkH> - <YkH(x*)><delrho(x*)delrhoYkH(x)>
                //                                                 - <YkH(x)><delrhoYkH(x*)delrho(x) + <YkH(x*)><YkH(x)><delrho(x*)delrho(x)>)
                //delrhoYkdelrhoYk = (spatialCross(i,j,k,37+nspecies-1)*stepsMinusOne + delrhoYkcross[nspecies-1]*delrhoYk[nspecies-1])*osteps;
                //spatialCross(i,j,k,37+nspecies+1) = (1.0/(meanrho*meanrhocross))*(delrhoYkdelrhoYk - meanYkcross[nspecies-1]*spatialCross(i,j,k,9)
                //                                        - meanYk[nspecies-1]*spatialCross(i,j,k,11) + meanYkcross[nspecies-1]*meanYk[nspecies-1]*spatialCross(i,j,k,0));

            });
        }
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
