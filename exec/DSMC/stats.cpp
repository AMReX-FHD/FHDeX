#include "DsmcParticleContainer.H"
using namespace std;
void FhdParticleContainer::EvaluateStats(MultiFab& mfcuInst,
    MultiFab& mfcuMeans,
    MultiFab& mfcuVars,
    MultiFab& mfprimInst,
    MultiFab& mfprimMeans,
    MultiFab& mfprimVars,
    MultiFab& mfcvlInst,
    MultiFab& mfcvlMeans,
    MultiFab& mfQMeans,
    MultiFab& mfcoVars,
    MultiFab& spatialCross1D,
    const int steps,
    Real time) {
    BL_PROFILE_VAR("EvaluateStats()",EvaluateStats);
    const Real osteps = 1.0/steps;
    const Real stepsMinusOne = steps-1.;

    // Zero out instantaneous values
    mfcuInst.setVal(0.);
    mfprimInst.setVal(0.);
    mfcvlInst.setVal(0.);
    mfcvlInst.setVal(0.);

    // Zero out heat flux and heat capacity
    mfQMeans.setVal(0.);
    mfcvlMeans.setVal(0.);

    // TODO: Add Heat Fluxes
    const int lev = 0;

    int ncon  = (nspecies+1)*5;
    int nprim = (nspecies+1)*10;
    int ncross = 38+nspecies*nspecies + nspecies;

    for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& aos = particle_tile.GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();

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
          0  - n   (n_ns)
          1  - rho (Y_ns)
          2  - u   (u_ns)
          3  - v   (v_ns)
          4  - w   (w_ns)
          5  - G   (G_ns) = dot(u_mean,dJ)
          6  - T   (T_ns)
          7  - P   (P_ns)
          8  - E   (E_ns)
          9  - c   (c_ns) (per species only)
        */


        Array4<Real> cuInst  = mfcuInst[pti].array();
        Array4<Real> primInst  = mfprimInst[pti].array();
        Array4<Real> cuMeans  = mfcuMeans[pti].array();
        Array4<Real> primMeans  = mfprimMeans[pti].array();
        Array4<Real> cuVars  = mfcuVars[pti].array();
        Array4<Real> primVars  = mfprimVars[pti].array();
        Array4<Real> coVars  = mfcoVars[pti].array();

        Array4<Real> cvlInst = mfcvlInst[pti].array();
        Array4<Real> cvlMeans = mfcvlMeans[pti].array();
        Array4<Real> QMeans  = mfQMeans[pti].array();

        //dsmcSpecies* propertiesPtr = properties;

        GpuArray<dsmcSpecies, MAX_SPECIES> propertiesTmp;
        for(int i=0;i<nspecies;i++)
        {
            propertiesTmp[i].mass = properties[i].mass;
            propertiesTmp[i].Neff = properties[i].Neff;
            propertiesTmp[i].R = properties[i].R;
            //Print() << "in: " << properties[i].mass << ", out: " << propertiesTmp[i].mass << endl;
        }

        Real ocollisionCellVolTmp = ocollisionCellVol;

        auto inds = m_bins.permutationPtr();
        auto offs = m_bins.offsetsPtr();
        //////////////////////////////////////
        // Primitve and Conserved Instantaneous Values
        //////////////////////////////////////

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const IntVect& iv = {i,j,k};
            long imap = tile_box.index(iv);
            int icon = 5; int iprim = 10; int icvl = 1;
            cvlInst(i,j,k,0) = 0;


            for (int l=0; l<nspecies; l++) {
                unsigned int np_spec = getBinSize(offs,iv,l,tile_box);
                unsigned int* cellList = getCellList(inds,offs,iv,l,tile_box);
                //long np_spec2 = m_cell_vectors[l][grid_id][imap].size();
//
                //                //cout << "old: " << np_spec2 << ", new: " << np_spec << endl;
                //

                Real mass = propertiesTmp[l].mass*propertiesTmp[l].Neff;
                Real moV  = propertiesTmp[l].mass*ocollisionCellVolTmp;
                //Real mass =1;
                //Real moV =1;
                primInst(i,j,k,iprim+0) = np_spec*ocollisionCellVolTmp;
                primInst(i,j,k,0)      += np_spec*ocollisionCellVolTmp;
                primInst(i,j,k,iprim+1) = np_spec*moV;
                primInst(i,j,k,1)      += np_spec*moV;
                cuInst(i,j,k,icon+0)    = np_spec*moV;

                Real rho = cuInst(i,j,k,icon+0);
                cvlInst(i,j,k,icvl) = 3.0*k_B*0.5/mass;
                cvlInst(i,j,k,0) += (cvlInst(i,j,k,icvl)*rho);

                // Read particle data
                for (int m=0; m<np_spec; m++) {
                    //int pind1 = m_cell_vectors[l][grid_id][imap][m];
                    int pind = cellList[m];

                    ParticleType & p = particles[pind];
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

                Real jx = cuInst(i,j,k,icon+1);
                Real jy = cuInst(i,j,k,icon+2);
                Real jz = cuInst(i,j,k,icon+3);
                Real K = cuInst(i,j,k,icon+4);
                Real n = primInst(i,j,k,iprim+0);

                // Total Conserved Vars
                for (int m=0; m<5; m++) {cuInst(i,j,k,m) += cuInst(i,j,k,icon+m);}

                if(rho != 0)
                {
                    primInst(i,j,k,iprim+2) = jx/rho;  // u_l
                    primInst(i,j,k,iprim+3) = jy/rho;  // v_l
                    primInst(i,j,k,iprim+4) = jz/rho;  // w_l
                }
                Real u = primInst(i,j,k,iprim+2);
                Real v = primInst(i,j,k,iprim+3);
                Real w = primInst(i,j,k,iprim+4);

                primInst(i,j,k,iprim+5) = u*jx+v*jy+w*jz;  // G_l

                Real vsqb = pow(u,2.)+pow(v,2.)+pow(w,2.);
                Real cv = cvlInst(i,j,k,icvl);

                primInst(i,j,k,iprim+6) = (K/rho - 0.5*vsqb)/cv;  // T_l


//                uTemp = uTemp/np_spec;
//                vTemp = vTemp/np_spec;
//                wTemp = wTemp/np_spec;
////
//                primInst(i,j,k,iprim+6) = 0;
////
//                for (int m=0; m<np_spec; m++) {
//                    int pind = m_cell_vectors[l][grid_id][imap][m];
//                    ParticleType & p = particles[pind];
//
//                    primInst(i,j,k,iprim+6) += (p.rdata(FHD_realData::velx) - uTemp)*(p.rdata(FHD_realData::velx) - uTemp) + (p.rdata(FHD_realData::vely) - vTemp)*(p.rdata(FHD_realData::vely) - vTemp)
//                                                + (p.rdata(FHD_realData::velz) - wTemp)*(p.rdata(FHD_realData::velz) - wTemp);
//                }
//                primInst(i,j,k,iprim+6) = primInst(i,j,k,iprim+6)*mass/(3.0*k_B*np_spec);
//
                Real T = primInst(i,j,k,iprim+6);
                //primInst(i,j,k,6) += T*n;

                primInst(i,j,k,iprim+7) = rho*(k_B/mass)*T;  // P_l
                primInst(i,j,k,7) += primInst(i,j,k,iprim+7);  // P

                // E not correct
                primInst(i,j,k,iprim+8) = 0.5*rho*vsqb+rho*cv*T;  // E_l
                primInst(i,j,k,8) += primInst(i,j,k,iprim+8);
                icon += 5; iprim += 10; icvl++;
            }


            Real uTemp = 0;
            Real vTemp = 0;
            Real wTemp = 0;
            Real massTotal = 0;

            int specTotal = 0;

            for (int l=nspecies-1; l>=0; l--) {
                //long np_spec = m_cell_vectors[l][grid_id][imap].size();
                unsigned int np_spec = getBinSize(offs,iv,l,tile_box);
                unsigned int* cellList = getCellList(inds,offs,iv,l,tile_box);
                for (int m=0; m<np_spec; m++) {
                    //int pind = m_cell_vectors[l][grid_id][imap][m];
                    int pind = cellList[m];

                    ParticleType & p = particles[pind];

                    uTemp += propertiesTmp[l].mass*p.rdata(FHD_realData::velx);
                    vTemp += propertiesTmp[l].mass*p.rdata(FHD_realData::vely);
                    wTemp += propertiesTmp[l].mass*p.rdata(FHD_realData::velz);
                    specTotal++;
                    massTotal +=propertiesTmp[l].mass;
                }

            }
            if(massTotal != 0)
            {
                uTemp /= massTotal;
                vTemp /= massTotal;
                wTemp /= massTotal;
            }


            //Total temperature
            primInst(i,j,k,6) = 0;

            for (int l=nspecies-1; l>=0; l--) {
                //long np_spec = m_cell_vectors[l][grid_id][imap].size();
                unsigned int np_spec = getBinSize(offs,iv,l,tile_box);
                unsigned int* cellList = getCellList(inds,offs,iv,l,tile_box);
                for (int m=0; m<np_spec; m++) {
                    //int pind = m_cell_vectors[l][grid_id][imap][m];
                    int pind = cellList[m];
                    ParticleType & p = particles[pind];

                    primInst(i,j,k,6) += propertiesTmp[l].mass*((p.rdata(FHD_realData::velx) - primMeans(i,j,k,2))*(p.rdata(FHD_realData::velx) - primMeans(i,j,k,2)) + (p.rdata(FHD_realData::vely) - primMeans(i,j,k,3))*(p.rdata(FHD_realData::vely) - primMeans(i,j,k,3))
                                                + (p.rdata(FHD_realData::velz) - primMeans(i,j,k,4))*(p.rdata(FHD_realData::velz) - primMeans(i,j,k,4)));
                }
            }


            primInst(i,j,k,6) = primInst(i,j,k,6)/(3.0*k_B*specTotal);


            Real rho = cuInst(i,j,k,0);
            Real jx = cuInst(i,j,k,1);
            Real jy = cuInst(i,j,k,2);
            Real jz = cuInst(i,j,k,3);
            Real K = cuInst(i,j,k,4);
            Real n = primInst(i,j,k,0);

            primInst(i,j,k,2) = jx/rho;  // u
            primInst(i,j,k,3) = jy/rho;  // v
            primInst(i,j,k,4) = jz/rho;  // w


            Real u = primInst(i,j,k,2);
            Real v = primInst(i,j,k,3);
            Real w = primInst(i,j,k,4);


            primInst(i,j,k,5) = u*jx+v*jy+w*jz;  // G

            //primInst(i,j,k,6) /= n;  // Mixture T

            // Energy Density
            cvlInst(i,j,k,0) /= rho;
            Real vsqb = pow(u,2.)+pow(v,2.)+pow(w,2.);
            Real cv = cvlInst(i,j,k,0);

            // Concentrations
            primInst(i,j,k,9) = 1;
            iprim = 10;
            for (int l=0; l<nspecies; l++) {

                primInst(i,j,k,iprim+9) = primInst(i,j,k,iprim+1)/primInst(i,j,k,1);
                iprim += 10;
            }


        });

        //////////////////////////////////////
        // Means
        //////////////////////////////////////

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const IntVect& iv = {i,j,k};
            long imap = tile_box.index(iv);

            primMeans(i,j,k,9) = 1;

            for (int l=0; l<ncon; l++) {
                cuMeans(i,j,k,l) = (cuMeans(i,j,k,l)*stepsMinusOne+cuInst(i,j,k,l))*osteps;
            }

            // Zero out hydrodynamic means that are sums of partials
            primMeans(i,j,k,0) = 0.0;
            //primMeans(i,j,k,6) = 0.0;
            primMeans(i,j,k,7) = 0.0;
            primMeans(i,j,k,8) = 0.0;
            int iprim = 10; int icon = 5; int icvl = 1;

            cvlMeans(i,j,k,0) = 0.;
            for(int l=0; l<nspecies; l++) {
                Real mass = propertiesTmp[l].mass*propertiesTmp[l].Neff;
                Real moV  = propertiesTmp[l].mass*ocollisionCellVolTmp;

                cvlMeans(i,j,k,icvl) = 3.0*k_B*0.5/mass;

                Real cv = cvlMeans(i,j,k,icvl);
                Real rho = cuMeans(i,j,k,icon+0);
                Real jx = cuMeans(i,j,k,icon+1);
                Real jy = cuMeans(i,j,k,icon+2);
                Real jz = cuMeans(i,j,k,icon+3);
                Real K = cuMeans(i,j,k,icon+4);

                cvlMeans(i,j,k,0) += cv*rho;

                primMeans(i,j,k,iprim+0) = rho/mass;  // n_l
                Real n = primMeans(i,j,k,iprim+0);

                primMeans(i,j,k,0)      += n;  // n
                primMeans(i,j,k,iprim+1) = rho;  // rho_l

                primMeans(i,j,k,iprim+2) = jx/rho;  // u_l
                primMeans(i,j,k,iprim+3) = jy/rho;  // v_l
                primMeans(i,j,k,iprim+4) = jz/rho;  // w_l

                Real u = primMeans(i,j,k,iprim+2);
                Real v = primMeans(i,j,k,iprim+3);
                Real w = primMeans(i,j,k,iprim+4);

                Real vsqb = pow(u,2)+pow(v,2)+pow(w,2);

                //primMeans(i,j,k,iprim+6) = (K/rho-vsqb*0.5)/cv;



                Real T = primMeans(i,j,k,iprim+6);
                //primMeans(i,j,k,6) += T*n;

                primMeans(i,j,k,iprim+7) = rho*(k_B/mass)*T;
                primMeans(i,j,k,7) += primMeans(i,j,k,iprim+7);
                primMeans(i,j,k,iprim+8) = vsqb*rho+cv*rho*T;
                primMeans(i,j,k,8) += primMeans(i,j,k,iprim+8);

                primMeans(i,j,k,l) = (primMeans(i,j,k,iprim+9)*stepsMinusOne+primInst(i,j,k,iprim+9))*osteps;

                primMeans(i,j,k,iprim+9) = rho/cuMeans(i,j,k,0);

                iprim += 10; icon += 5; icvl++;
            }

            // Evaluate Primitive Means from Conserved Means
            Real rho = cuMeans(i,j,k,0);
            Real jx = cuMeans(i,j,k,1);
            Real jy = cuMeans(i,j,k,2);
            Real jz = cuMeans(i,j,k,3);
            Real K = cuMeans(i,j,k,4);
            Real n = primMeans(i,j,k,0);

            primMeans(i,j,k,1)  = rho;
            primMeans(i,j,k,2)  = jx/rho; // u
            primMeans(i,j,k,3)  = jy/rho; // v
            primMeans(i,j,k,4)  = jz/rho; // w

            Real u = primMeans(i,j,k,2);
            Real v = primMeans(i,j,k,3);
            Real w = primMeans(i,j,k,4);

            primMeans(i,j,k,5) = u*jx+v*jy+w*jz; // G
            cvlMeans(i,j,k,0) /= rho;
            Real cv = cvlMeans(i,j,k,0);

            Real vsqb = pow(u,2.)+pow(v,2.)+pow(w,2.);
            //primMeans(i,j,k,6) /= n; // T

            Real tTemp = 0;
            int specTotal = 0;
            for (int l=nspecies-1; l>=0; l--) {
                //long np_spec = m_cell_vectors[l][grid_id][imap].size();
                unsigned int np_spec = getBinSize(offs,iv,l,tile_box);
                unsigned int* cellList = getCellList(inds,offs,iv,l,tile_box);
                for (int m=0; m<np_spec; m++) {
                    //int pind = m_cell_vectors[l][grid_id][imap][m];
                    int pind = cellList[m];
                    //int pind = 1;
                    ParticleType & p = particles[pind];

                    tTemp += propertiesTmp[l].mass*((p.rdata(FHD_realData::velx) - u)*(p.rdata(FHD_realData::velx) - u) + (p.rdata(FHD_realData::vely) - v)*(p.rdata(FHD_realData::vely) - v)
                                                + (p.rdata(FHD_realData::velz) - w)*(p.rdata(FHD_realData::velz) - w));
                    specTotal++;
                }
            }

            tTemp = tTemp/(3.0*k_B*specTotal);

            primMeans(i,j,k,6)  = (primMeans(i,j,k,6)*stepsMinusOne+tTemp)*osteps;

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
          17 - drho0.du
          18 - drho0.dv
          19 - drho0.du0
          20 - drho0.dv0
          21 - drho.dT
          22 - du.dT
          23 - dv.dT
          24 - dw.dT
        */

        Gpu::DeviceVector<Real> delConExt(ncon, 0.);
        Real* delCon = delConExt.data();

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Conserved Variances
            //Vector<Real> delCon(ncon, 0.0);
            //Real delCon[ncon];

            for (int l=0; l<ncon; l++) {
                delCon[l]        = cuInst(i,j,k,l) - cuMeans(i,j,k,l);
                cuVars(i,j,k,l)  = (cuVars(i,j,k,l)*stepsMinusOne+delCon[l]*delCon[l])*osteps;
            }

            Real drho = delCon[0];
            Real djx = delCon[1];
            Real djy = delCon[2];
            Real djz = delCon[3];
            Real dK = delCon[4];
            Real drho0 = delCon[5];
            Real djx0 = delCon[6];
            Real drho1 = delCon[10];
            Real djx1 = delCon[11];

            Real du = primInst(i,j,k,2)-primMeans(i,j,k,2);
            Real du0 = primInst(i,j,k,12)-primMeans(i,j,k,12);
            Real du1 = primInst(i,j,k,22)-primMeans(i,j,k,22);
            Real dv = primInst(i,j,k,3)-primMeans(i,j,k,3);
            Real dv0 = primInst(i,j,k,13)-primMeans(i,j,k,13);


            //Conserved Covariances
            coVars(i,j,k,0)  = (coVars(i,j,k, 0)*stepsMinusOne+drho*djx)*osteps; // drho.dJx
            coVars(i,j,k,1)  = (coVars(i,j,k, 1)*stepsMinusOne+drho*djy)*osteps; // drho.dJy
            coVars(i,j,k,2)  = (coVars(i,j,k, 2)*stepsMinusOne+drho*djz)*osteps; // drho.dJz
            coVars(i,j,k,3)  = (coVars(i,j,k, 3)*stepsMinusOne+drho*dK)*osteps;  // drho.dK
            coVars(i,j,k,4)  = (coVars(i,j,k, 4)*stepsMinusOne+djx*djy)*osteps;  // dJx.dJy
            coVars(i,j,k,5)  = (coVars(i,j,k, 5)*stepsMinusOne+djx*djz)*osteps;  // dJx.dJz
            coVars(i,j,k,6)  = (coVars(i,j,k, 6)*stepsMinusOne+djx*dK)*osteps;   // dJx.dK
            coVars(i,j,k,7)  = (coVars(i,j,k, 7)*stepsMinusOne+djy*djz)*osteps;  // dJy.dJz
            coVars(i,j,k,8)  = (coVars(i,j,k, 8)*stepsMinusOne+djy*dK)*osteps;   // dJy.dK
            coVars(i,j,k,9)  = (coVars(i,j,k, 9)*stepsMinusOne+djz*dK)*osteps;   // dJz.dK

            // Primitive Variances
            Real orhomean = 1.0/cuMeans(i,j,k,0);
            Real orhomean0 = 1.0/cuMeans(i,j,k,5);
            Real dn = primInst(i,j,k,0) - primMeans(i,j,k,0);
            primVars(i,j,k,0) = (primVars(i,j,k,0)*stepsMinusOne+dn*dn)*osteps; // dn.dn
            primVars(i,j,k,1) = drho*drho; // drho.drho

            Real umean = primMeans(i,j,k,2);
            Real umean0 = primMeans(i,j,k,12);
            Real vmean = primMeans(i,j,k,3);
            Real vmean0 = primMeans(i,j,k,13);
            Real wmean = primMeans(i,j,k,4);

            primVars(i,j,k,2) = // du.du
                pow(orhomean,2.0)*(cuVars(i,j,k,1)-2.0*umean*coVars(i,j,k,0)+pow(umean,2)*cuVars(i,j,k,0));
            primVars(i,j,k,3) = // dv.dv
                pow(orhomean,2.0)*(cuVars(i,j,k,2)-2.0*vmean*coVars(i,j,k,1)+pow(vmean,2)*cuVars(i,j,k,0));
            primVars(i,j,k,4) = // dw.dw
                pow(orhomean,2.0)*(cuVars(i,j,k,3)-2.0*wmean*coVars(i,j,k,2)+pow(wmean,2)*cuVars(i,j,k,0));

            Real dG = umean*djx+vmean*djy+wmean*djz;
            primVars(i,j,k,5) = // dG.dG <---- Is this correct? [Ask IS]
                (primVars(i,j,k,5)*stepsMinusOne+dG*dG)*osteps;

            coVars(i,j,k,10)   = (coVars(i,j,k,10)*stepsMinusOne+drho0*djx0)*osteps; // drho.dG
            coVars(i,j,k,11)   = (coVars(i,j,k,11)*stepsMinusOne+drho1*djx1)*osteps;  // dJx.dG
            coVars(i,j,k,12)   = (coVars(i,j,k,12)*stepsMinusOne+drho0*djx1)*osteps;  // dJy.dG
            coVars(i,j,k,13)   = (coVars(i,j,k,13)*stepsMinusOne+drho0*du0)*osteps;  // dJz.dG
            coVars(i,j,k,14)   = (coVars(i,j,k,14)*stepsMinusOne+drho1*du1)*osteps;   // dK.dG

            //coVars(i,j,k,15) = orhomean*(coVars(i,j,k,0) - umean*cuVars(i,j,k,0)); // drho.du
            coVars(i,j,k,15)  = (coVars(i,j,k,15)*stepsMinusOne+drho*du)*osteps; // drho.du
            coVars(i,j,k,16)  = (coVars(i,j,k,16)*stepsMinusOne+drho*dv)*osteps; // drho.dv
            coVars(i,j,k,17)  = (coVars(i,j,k,17)*stepsMinusOne+drho0*du)*osteps; // drho.du
            coVars(i,j,k,18)  = (coVars(i,j,k,18)*stepsMinusOne+drho1*du)*osteps; // drho.du
            coVars(i,j,k,19)  = (coVars(i,j,k,19)*stepsMinusOne+drho0*djx)*osteps; // drho.du
            coVars(i,j,k,20)  = (coVars(i,j,k,20)*stepsMinusOne+drho1*djx)*osteps; // drho.du

           // coVars(i,j,k,16) = orhomean*(coVars(i,j,k,1) - vmean*cuVars(i,j,k,0)); // drho.dv

           // coVars(i,j,k,15) = orhomean*(coVars(i,j,k,0) - umean*cuVars(i,j,k,0)); // drho.du
           // coVars(i,j,k,16) = orhomean*(coVars(i,j,k,1) - vmean*cuVars(i,j,k,0)); // drho.dv

          //  coVars(i,j,k,17) = orhomean*(coVars(i,j,k,2) - wmean*cuVars(i,j,k,0)); // drho.dw

          //  coVars(i,j,k,18) = // du.dv
          //      pow(orhomean,2.0)*(coVars(i,j,k,4)-umean*coVars(i,j,k,1)-vmean*coVars(i,j,k,0)
          //      +umean*vmean*cuVars(i,j,k,0));
         //   coVars(i,j,k,19) = // du.dw
         //       pow(orhomean,2.0)*(coVars(i,j,k,5)-umean*coVars(i,j,k,2)-wmean*coVars(i,j,k,0)
         //       +umean*wmean*cuVars(i,j,k,0));
         //   coVars(i,j,k,20) = // dv.dw
         //       pow(orhomean,2.0)*(coVars(i,j,k,7)-vmean*coVars(i,j,k,2)-wmean*coVars(i,j,k,1)
         //       +vmean*wmean*cuVars(i,j,k,0));

            Real vsqb = pow(umean,2)+pow(vmean,2)+pow(wmean,2);
            Real cv = cvlMeans(i,j,k,0);
            QMeans(i,j,k,0) = cv*primMeans(i,j,k,6)-0.5*vsqb;
            Real Qbar = QMeans(i,j,k,0);

            primVars(i,j,k,6) = pow(orhomean/cv,2.0)* // dT.dT
                (cuVars(i,j,k,4)+primVars(i,j,k,5)+pow(Qbar,2.)*cuVars(i,j,k,0)
                 -2.0*coVars(i,j,k,14)-2.0*Qbar*coVars(i,j,k,3)+2.0*Qbar*coVars(i,j,k,10));

            coVars(i,j,k,21) = orhomean/cv*(coVars(i,j,k,3)-coVars(i,j,k,10)-Qbar*cuVars(i,j,k,0)); // drho.dT

            coVars(i,j,k,22) = // du.dT
                pow(orhomean,2.0)/cv*(coVars(i,j,k,6)-umean*coVars(i,j,k,3)-coVars(i,j,k,11)
                +umean*coVars(i,j,k,10)-Qbar*coVars(i,j,k,0)+umean*Qbar*cuVars(i,j,k,0));
            coVars(i,j,k,23) = // dv.dT
                pow(orhomean,2.0)/cv*(coVars(i,j,k,8)-vmean*coVars(i,j,k,3)-coVars(i,j,k,12)
                +vmean*coVars(i,j,k,10)-Qbar*coVars(i,j,k,1)+vmean*Qbar*cuVars(i,j,k,0));
            coVars(i,j,k,24) = // dw.dT
                pow(orhomean,2.0)/cv*(coVars(i,j,k,9)-wmean*coVars(i,j,k,3)-coVars(i,j,k,13)
                +wmean*coVars(i,j,k,10)-Qbar*coVars(i,j,k,2)+wmean*Qbar*cuVars(i,j,k,0));

            //primVars(i,j,k,7) = // Add (R_m/c_m) dP.dP
                            //(cuvars(i,j,k,8)
            // TODO: Add P and E variances
            // TODO: Add variances by species
        });
    }


    // Spatial Correlations (works only for 1D currently: 1 cell each in y and z directions)
    // in this order: conserved variables [rho, K, jx, jy, jz] + primitive variables [vx, vy, vz, T]
    // total: 2*9
    // add species later - I added a couple - Daniel
    if (plot_cross) {
        int nstats = 50;

        // Get all nstats at xcross and store in GpuVector
        amrex::Gpu::ManagedVector<Real> data_xcross_in(nstats, 0.0); // values at x*
        amrex::Gpu::ManagedVector<Real> cvq_xcross_in(2, 0.0); // values at x*
        for ( MFIter mfi(mfcuInst); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.validbox();

            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);

            const Array4<const Real> cumeans   = mfcuMeans.array(mfi);
            const Array4<const Real> primmeans = mfprimMeans.array(mfi);
            const Array4<const Real> prim      = mfprimInst.array(mfi);
            const Array4<const Real> cu        = mfcuInst.array(mfi);

            const Array4<const Real> cvlMeans  = mfcvlMeans.array(mfi);
            const Array4<const Real> QMeans    = mfQMeans.array(mfi);

            Real* data_xcross = data_xcross_in.data();
            Real* cvq_xcross = cvq_xcross_in.data();
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

                    data_xcross[17] = primmeans(i,j,k,6); // T-mean


                    cvq_xcross[0]   = cvlMeans(i,j,k,0);  // cv-mean
                    cvq_xcross[1]   = QMeans(i,j,k,0);    // Q-mean


                    for(int m = 0; m<nspecies; m++)
                    {

                        data_xcross[19 + 2*m] = cu(i,j,k,(m+1)*5); // rho_-instant
                        data_xcross[19 + 2*m +1] = cumeans(i,j,k,(m+1)*5); // rho_-mean
                    }
                    for(int m = 0; m<nspecies; m++)
                    {

                        data_xcross[19 + 2*nspecies + 2*m] = cu(i,j,k,(m+1)*5 + 1); // jx_-instant
                        data_xcross[19 + 2*nspecies + 2*m +1] = cumeans(i,j,k,(m+1)*5 + 1); // jx_-mean
                    }
                    for(int m = 0; m<nspecies; m++)
                    {

                        data_xcross[19 + 4*nspecies + 2*m] = prim(i,j,k,(m+1)*10 + 2); // ux_-instant
                        data_xcross[19 + 4*nspecies + 2*m +1] = primmeans(i,j,k,(m+1)*10 + 2); // ux_-mean
                    }

               }
            });
        }  // end MFITer

        // Reduce across MPI processes
        ParallelDescriptor::ReduceRealSum(data_xcross_in.data(),nstats);
        ParallelDescriptor::ReduceRealSum(cvq_xcross_in.data(),2);

        // Update Spatial Correlations
        for ( MFIter mfi(mfcuInst); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.validbox();

            const Array4<const Real> cumeans   = mfcuMeans.array(mfi);
            const Array4<const Real> primmeans = mfprimMeans.array(mfi);
            const Array4<const Real> prim      = mfprimInst.array(mfi);
            const Array4<const Real> cu        = mfcuInst.array(mfi);
            const Array4<const Real> cvlMeans  = mfcvlMeans.array(mfi);
            const Array4<const Real> QMeans    = mfQMeans.array(mfi);

            const Array4<      Real> spatialCross = spatialCross1D.array(mfi);

            Real* data_xcross = data_xcross_in.data();
            Real* cvq_xcross = cvq_xcross_in.data();

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {

                //////////////////////////////////////////////
                // Get fluctuations at xcross for this j and k
                //////////////////////////////////////////////
                Real meanrhocross = data_xcross[1];
                Real meanrho0cross = data_xcross[20];
                Real meanjx0cross = data_xcross[24];
                Real meanux0cross = data_xcross[28];
                Real meanuxcross  = data_xcross[11];

                Real delrhocross = data_xcross[0] - data_xcross[1];
                Real delKcross   = data_xcross[2] - data_xcross[3];
                Real deljxcross  = data_xcross[4] - data_xcross[5];
                Real deljycross  = data_xcross[6] - data_xcross[7];
                Real deljzcross  = data_xcross[8] - data_xcross[9];

                Real delTcross = data_xcross[16] - data_xcross[17];
                Real delvxcross = data_xcross[10] - data_xcross[11];

                Real cvinvcross = 1.0/cvq_xcross[0];
                Real qmeancross = cvq_xcross[1];
                Real vxmeancross = data_xcross[11];
                Real vymeancross = data_xcross[13];
                Real vzmeancross = data_xcross[15];

                //interspecies
                Real delrho0cross = data_xcross[19] - data_xcross[20];
                Real deljx0cross = data_xcross[23] - data_xcross[24];
                Real delux0cross = data_xcross[27] - data_xcross[28];

                Real deljx1cross = data_xcross[25] - data_xcross[26];
                Real delrho1cross;
                if(nspecies >1)
                {
                    delrho1cross = data_xcross[21] - data_xcross[22];
                }

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

                Real vxmean = primmeans(i,j,k,2);
                Real vymean = primmeans(i,j,k,3);
                Real vzmean = primmeans(i,j,k,4);

                Real delrho0 = cu(i,j,k,5+0) - cumeans(i,j,k,5+0);
                Real delrho1;
                if(nspecies > 1)
                {
                    delrho1 = cu(i,j,k,2*5+0) - cumeans(i,j,k,2*5+0);
                }

                Real cvinv = 1.0/cvlMeans(i,j,k,0);
                Real qmean = QMeans(i,j,k,0);

                // delG = \vec{v}\cdot\vec{\deltaj}
                Real delG = vxmean*deljx+vymean*deljy+vzmean*deljz;

                // Spatial Correlation Calculations

                spatialCross(i,j,k,0) = (spatialCross(i,j,k,0)*stepsMinusOne + delrhocross*delrho)*osteps;  // <delrho(x*)delrho(x)>
                spatialCross(i,j,k,1) = (spatialCross(i,j,k,1)*stepsMinusOne + delKcross*delK)*osteps;      // <delK(x*)delK(x)>
                //                spatialCross(i,j,k,1) = (spatialCross(i,j,k,1)*stepsMinusOne + delrhocross*delrho0)*osteps;      // <delK(x*)delK(x)>
                //                spatialCross(i,j,k,2) = (spatialCross(i,j,k,1)*stepsMinusOne + delrho0cross*delrho0)*osteps;      // <delK(x*)delK(x)>
                spatialCross(i,j,k,2) = (spatialCross(i,j,k,2)*stepsMinusOne + deljxcross*deljx)*osteps;    // <deljx(x*)deljx(x)>
                spatialCross(i,j,k,3) = (spatialCross(i,j,k,3)*stepsMinusOne + deljycross*deljy)*osteps;    // <deljy(x*)deljy(x)>
                spatialCross(i,j,k,4) = (spatialCross(i,j,k,4)*stepsMinusOne + deljzcross*deljz)*osteps;    // <deljz(x*)deljz(x)>

                spatialCross(i,j,k,5) = (spatialCross(i,j,k,5)*stepsMinusOne + deljxcross*delrho)*osteps;   // <deljx(x*)delrho(x)>

                spatialCross(i,j,k,6) = (spatialCross(i,j,k,6)*stepsMinusOne + deljxcross*delrho0)*osteps;   // <deljy(x*)delrho(x)>
                spatialCross(i,j,k,7) = (spatialCross(i,j,k,7)*stepsMinusOne + deljx0cross*delrho0)*osteps;   // <deljz(x*)delrho(x)>
                spatialCross(i,j,k,8) = (spatialCross(i,j,k,8)*stepsMinusOne + deljxcross*delrho1)*osteps;    // <delK(x*)delrho(x)>
                spatialCross(i,j,k,9) = (spatialCross(i,j,k,9)*stepsMinusOne + deljx1cross*delrho1)*osteps;   // <delrho(x*)deljx(x)>
                spatialCross(i,j,k,10) = (spatialCross(i,j,k,10)*stepsMinusOne + deljx0cross*delrho1)*osteps;  // <deljy(x*)deljx(x)>
                spatialCross(i,j,k,11) = (spatialCross(i,j,k,11)*stepsMinusOne + deljx1cross*delrho0)*osteps;  // <deljz(x*)deljx(x)>

                spatialCross(i,j,k,12) = (spatialCross(i,j,k,12)*stepsMinusOne + deljx0cross*delrho)*osteps;   // <delK(x*)deljx(x)>

                spatialCross(i,j,k,13) = (spatialCross(i,j,k,13)*stepsMinusOne + deljx1cross*delrho)*osteps; // <delrho(x*)deljy(x)>
                spatialCross(i,j,k,14) = (spatialCross(i,j,k,14)*stepsMinusOne + deljxcross*deljy)*osteps;  // <deljx(x*)deljy(x)>
                spatialCross(i,j,k,15) = (spatialCross(i,j,k,15)*stepsMinusOne + deljzcross*deljy)*osteps;  // <deljz(x*)deljy(x)>
                spatialCross(i,j,k,16) = (spatialCross(i,j,k,16)*stepsMinusOne + delKcross*deljy)*osteps;   // <delK(x*)deljy(x)>

                spatialCross(i,j,k,17) = (spatialCross(i,j,k,17)*stepsMinusOne + delrhocross*deljz)*osteps; // <delrho(x*)delK(x)>
                spatialCross(i,j,k,18) = (spatialCross(i,j,k,18)*stepsMinusOne + deljxcross*deljz)*osteps;  // <deljx(x*)delK(x)>
                spatialCross(i,j,k,19) = (spatialCross(i,j,k,19)*stepsMinusOne + deljycross*deljz)*osteps;  // <deljy(x*)delK(x)>
                spatialCross(i,j,k,20) = (spatialCross(i,j,k,20)*stepsMinusOne + delKcross*deljz)*osteps;   // <deljz(x*)delK(x)>

                // [IS] Species-dependent stuff -- commented now -- add later
                //spatialCross(i,j,k,6) = (spatialCross(i,j,k,6)*stepsMinusOne + deljxcross*delrhoYk[0])*osteps;  // <deljx(x*)delrhoYkL(x)>
                //spatialCross(i,j,k,7) = (spatialCross(i,j,k,7)*stepsMinusOne + deljxcross*delrhoYk[nspecies-1])*osteps;  // <deljx(x*)delrhoYkH(x)>
                //spatialCross(i,j,k,8) = (spatialCross(i,j,k,8)*stepsMinusOne + delrhocross*delrhoYk[0])*osteps; // <delrho(x*)delrhoYkL(x)>
                //spatialCross(i,j,k,9) = (spatialCross(i,j,k,9)*stepsMinusOne + delrhocross*delrhoYk[nspecies-1])*osteps; // <delrho(x*)delrhoYkH(x)>
                //spatialCross(i,j,k,10)= (spatialCross(i,j,k,10)*stepsMinusOne + delrhoYkcross[0]*delrho)*osteps; // <delrhoYkL(x*)delrho(x)>
                //spatialCross(i,j,k,11)= (spatialCross(i,j,k,11)*stepsMinusOne + delrhoYkcross[nspecies-1]*delrho)*osteps; // <delrhoYkH(x*)delrho(x)>

                spatialCross(i,j,k,21) = (spatialCross(i,j,k,21)*stepsMinusOne + delGcross*delG)*osteps;   // <delG(x*)delG(x)>
                spatialCross(i,j,k,22) = (spatialCross(i,j,k,22)*stepsMinusOne + delGcross*delK)*osteps;   // <delG(x*)delK(x)>
                spatialCross(i,j,k,23) = (spatialCross(i,j,k,23)*stepsMinusOne + delKcross*delG)*osteps;   // <delK(x*)delG(x)>
                spatialCross(i,j,k,24) = (spatialCross(i,j,k,24)*stepsMinusOne + delrhocross*delG)*osteps; // <delrho(x*)delG(x)>
                spatialCross(i,j,k,25) = (spatialCross(i,j,k,25)*stepsMinusOne + delGcross*delrho)*osteps; // <delG(x*)delrho(x)>
                spatialCross(i,j,k,26) = (spatialCross(i,j,k,26)*stepsMinusOne + deljxcross*delG)*osteps;  // <deljx(x*)delG(x)>
                spatialCross(i,j,k,27) = (spatialCross(i,j,k,27)*stepsMinusOne + delGcross*deljx)*osteps;  // <delG(x*)deljx(x)>
                spatialCross(i,j,k,28) = (spatialCross(i,j,k,28)*stepsMinusOne + deljycross*delG)*osteps;  // <deljy(x*)delG(x)>
                spatialCross(i,j,k,29) = (spatialCross(i,j,k,29)*stepsMinusOne + delGcross*deljy)*osteps;  // <delG(x*)deljy(x)>
                spatialCross(i,j,k,30) = (spatialCross(i,j,k,30)*stepsMinusOne + deljzcross*delG)*osteps;  // <deljz(x*)delG(x)>
                spatialCross(i,j,k,31) = (spatialCross(i,j,k,31)*stepsMinusOne + delGcross*deljz)*osteps;  // <delG(x*)deljz(x)>
                spatialCross(i,j,k,32) = (spatialCross(i,j,k,32)*stepsMinusOne + delKcross*delG)*osteps;   // <delK(x*)delG(x)>
                spatialCross(i,j,k,33) = (spatialCross(i,j,k,33)*stepsMinusOne + delGcross*delK)*osteps;   // <delG(x*)delK(x)>
                spatialCross(i,j,k,34) = (spatialCross(i,j,k,34)*stepsMinusOne + delTcross*delT)*osteps;   // <delT(x*)delT(x)>

                // <delT(x*)delT(x)> = (1/cv*/cv/<rho(x)>/<rho(x*)>)(<delK*delK> + <delG*delG> - <delG*delK> - <delK*delG>
                //                      + <Q><Q*><delrho*delrho> - <Q*><delrho*delK> - <Q><delK*delrho>
                //                                            + <Q*><delrho*delG> + <Q><delG*delrho>)
//                spatialCross(i,j,k,34) = (cvinvcross*cvinv/(meanrhocross*meanrho))*
//                                    (spatialCross(i,j,k,1) + spatialCross(i,j,k,21) - spatialCross(i,j,k,22) - spatialCross(i,j,k,23)
//                                    + qmean*qmeancross*spatialCross(i,j,k,0) - qmeancross*spatialCross(i,j,k,17) - qmean*spatialCross(i,j,k,8)
//                                    + qmeancross*spatialCross(i,j,k,24) + qmean*spatialCross(i,j,k,25));

                // <delT(x*)delrho(x)> = (1/cv/<rho(x*)>)*(<delK*delrho> - <delG*delrho> - <Q*><delrhodelrho*>)
                spatialCross(i,j,k,35) = (cvinvcross/meanrhocross)*
                    (spatialCross(i,j,k,8) - spatialCross(i,j,k,25) - qmeancross*spatialCross(i,j,k,0));

                // <delu(x*)delrho> = (1/<rho(x*)>)*(<deljx(x*)delrho(x)> - <u(x*)><<delrho(x*)delrho(x)>)
                spatialCross(i,j,k,36) = (1.0/meanrhocross)*(spatialCross(i,j,k,5) - meanuxcross*spatialCross(i,j,k,0));

                // <delT(x*)delu> = (1/cv'/<rho(x*)>/<rho(x)>)*(<deljx(x*)delK(x)> - <u(x)><<delrho(x)delK(x*)>
                //        - <deljx(x)delG(x*)> + <u(x)><delrho(x)delG(x*)> - Qbar(x*)<deljx(x)delrho(x*)> +
                //        <u(x)>Qbar(x*)<delrho(x)delrho(x*)>)
                spatialCross(i,j,k,37) = (cvinvcross/(meanrhocross*meanrho))*(spatialCross(i,j,k,18) - vxmean*spatialCross(i,j,k,8)
                    - spatialCross(i,j,k,27) + vxmean*spatialCross(i,j,k,25) - qmeancross*spatialCross(i,j,k,9)
                    + vxmean*qmeancross*spatialCross(i,j,k,0));


                //added a few species specific terms here
                spatialCross(i,j,k,38) = (spatialCross(i,j,k,38)*stepsMinusOne + delrho0cross*delrho0)*osteps;   // <delRho0(x*)delRho0(x)>
                int cnt = 39;
                if(nspecies > 1)
                {
                    //spatialCross(i,j,k,cnt) = (spatialCross(i,j,k,cnt)*stepsMinusOne + delrho0cross*delrho1)*osteps;   // <delRho1(x*)delRho1(x)>
                    //cnt++;
                    //spatialCross(i,j,k,cnt) = (spatialCross(i,j,k,cnt)*stepsMinusOne + delrho1cross*delrho0)*osteps;   // <delRho1(x*)delRho0(x)>
                    //cnt++;
                    //spatialCross(i,j,k,cnt) = (spatialCross(i,j,k,cnt)*stepsMinusOne + delrho1cross*delrho1)*osteps;   // <delRho0(x*)delRho1(x)>
                    //cnt++;
                }
                spatialCross(i,j,k,cnt) = (spatialCross(i,j,k,cnt)*stepsMinusOne + delrho1cross*delrho1)*osteps;   // <delRho0(x*)delRho0(x)>
                cnt++;
                spatialCross(i,j,k,cnt) = (spatialCross(i,j,k,cnt)*stepsMinusOne + delux0cross*delrho0)*osteps;
                cnt++;
                spatialCross(i,j,k,cnt) = (spatialCross(i,j,k,cnt)*stepsMinusOne + delvxcross*delrho0)*osteps;
                cnt++;
                spatialCross(i,j,k,cnt) = (spatialCross(i,j,k,cnt)*stepsMinusOne + delrho1cross*delrho0)*osteps;
                cnt++;
                spatialCross(i,j,k,cnt) = (spatialCross(i,j,k,cnt)*stepsMinusOne + delvxcross*delrho1)*osteps;
                cnt++;
                //spatialCross(i,j,k,cnt) = (1.0/meanrho0cross)*(-meanux0cross*spatialCross(i,j,k,2));
                //cnt++;
//                if(nspecies > 1)
//                {
//
//                    spatialCross(i,j,k,cnt) = (spatialCross(i,j,k,cnt)*stepsMinusOne + delvxcross*delrho1)*osteps;
//                    cnt++;
//                }

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
        });
    }
}
