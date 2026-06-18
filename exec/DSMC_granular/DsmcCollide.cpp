#include "DsmcParticleContainer.H"
#include <math.h>

int FhdParticleContainer::getSpeciesIndex(int species1, int species2)
{
    if(species1<species2)
    {
        return species2+nspecies*species1;
    } else
    {
        return species1+nspecies*species2;
    }
}

Real FhdParticleContainer::g0_Ma_Ahmadi(int ispec, int jspec, Real iphi, Real jphi)
{
    const Real C1 = 2.5, C2 = 4.5904;
    const Real C3 = 4.515439, CPOW = 0.67802;
    Real chi;
    if(ispec==jspec)
    {
        // Ma Ahmadi 1990s
        Real numer = 1 + C1*iphi + C2*pow(iphi,2) + C3*pow(iphi,3);
        numer = numer * 4.0*iphi;
        Real phiRatio = iphi/phi_max;
        Real denom = pow(1.0 - pow(phiRatio,3),CPOW);
        chi = 1+numer/denom;
        return std::min(chi,chi_max);
    }
    else
    {
        // Mansoori et al (1971) Equ. Thermo. Prop. of the Mix. of HS
        // Agreement looks good up to 0.50 solid fraction from paper
        const Real phiTotal = iphi + jphi;
        const Real irad = properties[ispec].radius;
        const Real jrad = properties[jspec].radius;
        Real a1 = 1.0/(1.0-phiTotal);
        Real a2 = 3.0*irad*jrad/(irad+jrad);
        Real eta2 = (iphi/irad) + (jphi/jrad);
        eta2 = eta2/pow(1.0-phiTotal,2);
        a2 = a2 * eta2;
        Real a3 = 2.0 * pow((irad*jrad)/(irad+jrad),2);
        a3 = a3 * pow(eta2,2) / pow(1.0-phiTotal,3);
        chi = a1 + a2 + a3;
        return std::min(chi,chi_max);
    }
}

void FhdParticleContainer::InitCollisionCells()
{
    BL_PROFILE_VAR("InitCollisionCells()",InitCollisionCells);
    int indx, ij_spec;
    int cnt = 0;
    Real totalPhi = 0;
    Real mindiam = 1000;
    for(int i_spec=0;i_spec<nspecies;i_spec++)
    {
        totalPhi += phi_domain[i_spec];
        mindiam = std::min(mindiam,2.0*properties[i_spec].radius);
        for(int j_spec=0;j_spec<nspecies;j_spec++)
        {
            ij_spec = getSpeciesIndex(i_spec,j_spec);
            interproperties[ij_spec].alpha = alpha_pp[cnt];
            interproperties[ij_spec].csx = pow(properties[i_spec].radius+properties[j_spec].radius,2)*pi_usr;
            cnt++;
        }
        countedCollisions[i_spec] = 0;
        expectedCollisions[i_spec] = 0;
    }

    const int lev = 0;
    for (FhdParIter pti(* this, lev); pti.isValid(); ++pti)
    {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();

        const Array4<Real> & arrvrmax = mfvrmax.array(pti);
        const Array4<Real> & arrphi = mfphi.array(pti);
        const Array4<Real> & arrselect = mfselect.array(pti);

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            int ij_spec;
            for (int i_spec=0; i_spec<nspecies; i_spec++) {
                for (int j_spec = i_spec; j_spec < nspecies; j_spec++) {
                    ij_spec = getSpeciesIndex(i_spec,j_spec);
                    arrselect(i,j,k,ij_spec) = 0.0;
                }
            }

            const IntVect& iv = {i,j,k};
            long imap = tile_box.index(iv);

            for (int i_spec=0; i_spec<nspecies; i_spec++)
            {
                arrphi(i,j,k,i_spec) = m_cell_vectors[i_spec][grid_id][imap].size()
                    *properties[i_spec].part2cellVol*properties[i_spec].Neff;
            }
        });
    }
}

// Compute selections here
void FhdParticleContainer::CalcSelections(Real dt)
{
    BL_PROFILE_VAR("CalcSelections()",CalcSelections);
    int lev = 0;
    mfselect.setVal(0.0);
    for(MFIter mfi(mfvrmax); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();

        const Array4<Real> & arrvrmax = mfvrmax.array(mfi);
        const Array4<Real> & arrphi = mfphi.array(mfi);
        const Array4<Real> & arrselect = mfselect.array(mfi);

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            int ij_spec;
            long np_i, np_j;

            const IntVect& iv = {i,j,k};
            long imap = tile_box.index(iv);

            Real vrmax;
            Real NSel;
            Real phi1, phi2, chi0;
            Real crossSection;

            for (int i_spec=0; i_spec<nspecies; i_spec++)
            {
                arrphi(i,j,k,i_spec) = m_cell_vectors[i_spec][grid_id][imap].size()*
                    properties[i_spec].Neff*properties[i_spec].part2cellVol;
            }

            for (int i_spec = 0; i_spec < nspecies; i_spec++)
            {
                for (int j_spec = i_spec; j_spec < nspecies; j_spec++) {
                    ij_spec = getSpeciesIndex(i_spec,j_spec);
                    np_i = m_cell_vectors[i_spec][grid_id][imap].size();
                    np_j = m_cell_vectors[j_spec][grid_id][imap].size();
                    phi1 = arrphi(i,j,k,i_spec);
                    phi2 = arrphi(i,j,k,j_spec);
                    // comment out if expecting dilute
                    chi0 = 1.0;
                    //chi0 = g0_Ma_Ahmadi(i_spec,j_spec, phi1, phi2);
                    vrmax = arrvrmax(i,j,k,i_spec);
                    crossSection = interproperties[ij_spec].csx;
                    NSel = 4.0*particle_neff*np_i*np_j*crossSection*vrmax*ocollisionCellVol*chi0*dt;
                    if(i_spec==j_spec) {NSel = NSel*0.5; np_j = np_i-1;}
                    arrselect(i,j,k,ij_spec) = std::floor(NSel + amrex::Random());
                }
            }
        });
    }
}

void FhdParticleContainer::CollideParticles(Real dt)
{
    BL_PROFILE_VAR("CollideParticles()",CollideParticles);
    int lev = 0;
    for(MFIter mfi(mfvrmax); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();

        const Array4<Real> & arrvrmax = mfvrmax.array(mfi);
        const Array4<Real> & arrselect = mfselect.array(mfi);

        const long np = particles.numParticles();
        //amrex::ParallelForRNG(tile_box,
        //  [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept {
        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();

        for (int i = smallEnd[0]; i <= bigEnd[0]; i++) {
        for (int j = smallEnd[1]; j <= bigEnd[1]; j++) {
        for (int k = smallEnd[2]; k <= bigEnd[2]; k++) {
            const IntVect& iv = {i,j,k};
            long imap = tile_box.index(iv);

            Real NSel[nspecies*nspecies], totalSel, selrun;
            int np[nspecies];
            int pindxi, pindxj; // index of randomly sampled particles
            int ij_spec;

            RealVect eij, vreij;
            RealVect vi, vj, vij;
            RealVect vijpost, boost;
            //Real oboostmag;
            Real massi, massj, massij;
            Real vrmag, vrmax, vreijmag;

            totalSel = 0;
            for (int i_spec = 0; i_spec<nspecies; i_spec++)
            {
                np[i_spec] = m_cell_vectors[i_spec][grid_id][imap].size();
                for (int j_spec = i_spec; j_spec < nspecies; j_spec++)
                {
                    ij_spec = getSpeciesIndex(i_spec,j_spec);
                    NSel[ij_spec] = (int)arrselect(i,j,k,ij_spec);
                    totalSel += NSel[ij_spec];
                }
            }

            int speci, specj, specij;
            while (totalSel>0)
            {
                Real RR = amrex::Random();
                bool spec_select = false;
                selrun = 0;
                speci = -1; specj = -1; specij = -1;
                for(int i_spec=0;i_spec<nspecies;i_spec++)
                {
                    for(int j_spec=i_spec;j_spec<nspecies;j_spec++)
                    {
                        int ij_spec = getSpeciesIndex(i_spec,j_spec);
                        selrun += NSel[ij_spec];
                        if(selrun/totalSel>RR && !spec_select)
                        {
                            spec_select = true;
                            specij = ij_spec;
                            speci = i_spec; specj = j_spec;
                            NSel[ij_spec] -= 1;
                        }
                    }
                }
                totalSel--;
                massi = properties[speci].mass;
                massj = properties[specj].mass;
                massij = properties[speci].mass + properties[specj].mass;
                vrmax = arrvrmax(i,j,k,specij);
                pindxi = floor(amrex::Random()*np[speci]);
                pindxj = floor(amrex::Random()*np[specj]);
                pindxi = m_cell_vectors[speci][grid_id][imap][pindxi];
                pindxj = m_cell_vectors[specj][grid_id][imap][pindxj];
                ParticleType &    parti = particles[pindxi];
                ParticleType & partj = particles[pindxj];

                vi[0] = parti.rdata(FHD_realData::velx);
                vi[1] = parti.rdata(FHD_realData::vely);
                vi[2] = parti.rdata(FHD_realData::velz);

                vj[0] = partj.rdata(FHD_realData::velx);
                vj[1] = partj.rdata(FHD_realData::vely);
                vj[2] = partj.rdata(FHD_realData::velz);

                vij[0] = vi[0]-vj[0]; vij[1] = vi[1]-vj[1]; vij[2] = vi[2]-vj[2];
                vrmag = sqrt(pow(vij[0],2)+pow(vij[1],2)+pow(vij[2],2));
                if(vrmag>vrmax) {vrmax = vrmag; arrvrmax(i,j,k,ij_spec) = vrmax;}

                Real eijmag = 0;
                Real theta = 2.0*pi_usr*amrex::Random();
                Real phi = std::acos(1.0-2.0*amrex::Random());
                eij[0] = std::sin(phi)*std::cos(theta);
                eij[1] = std::sin(phi)*std::sin(theta);
                eij[2] = std::cos(phi);
                eijmag = pow(eij[0],2)+pow(eij[1],2)+pow(eij[2],2);
                eijmag = pow(eijmag,0.5);
                for(int idim=0; idim<3; idim++)
                {
                eij[idim] /= eijmag;
                }

                vreijmag = vij[0]*eij[0]+vij[1]*eij[1]+vij[2]*eij[2];
                if(vreijmag>vrmax*amrex::Random())
                {
                    countedCollisions[speci] += 1;
                    countedCollisions[specj] += 1;

                    vreijmag = vreijmag*(1.0+interproperties[specij].alpha)/massij;
                    vreij[0] = vreijmag*eij[0];
                    vreij[1] = vreijmag*eij[1];
                    vreij[2] = vreijmag*eij[2];

                    parti.rdata(FHD_realData::velx) = vi[0] - vreij[0]*massj;
                    parti.rdata(FHD_realData::vely) = vi[1] - vreij[1]*massj;
                    parti.rdata(FHD_realData::velz) = vi[2] - vreij[2]*massj;
                    partj.rdata(FHD_realData::velx) = vj[0] + vreij[0]*massi;
                    partj.rdata(FHD_realData::vely) = vj[1] + vreij[1]*massi;
                    partj.rdata(FHD_realData::velz) = vj[2] + vreij[2]*massi;

                    /*
                    // Boosted Velocities
                    vijpost[0] = vi[0]-vj[0];
                    vijpost[1] = vi[1]-vj[1];
                    vijpost[2] = vi[2]-vj[2];

                    boost[0] = vijpost[0]-vij[0];
                    boost[1] = vijpost[1]-vij[1];
                    boost[2] = vijpost[2]-vij[2];
                    oboostmag = (properties[speci].radius+properties[specj].radius)/
                        (sqrt(pow(boost[0],2)+pow(boost[1],2)+pow(boost[2],2))*dt);
                    boost[0] = boost[0]*oboostmag;
                    boost[1] = boost[1]*oboostmag;
                    boost[2] = boost[2]*oboostmag;

                    parti.rdata(FHD_realData::boostx) = boost[0];
                    parti.rdata(FHD_realData::boosty) = boost[1];
                    parti.rdata(FHD_realData::boostz) = boost[2];
                    partj.rdata(FHD_realData::boostx) = -boost[0];
                    partj.rdata(FHD_realData::boosty) = -boost[1];
                    partj.rdata(FHD_realData::boostz) = -boost[2];
                    */
                }
            }
        }
        }
        }
    }
}