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

void FhdParticleContainer::InitCollisionCells()
{
    BL_PROFILE_VAR("InitCollisionCells()",InitCollisionCells);
    int indx, ij_spec;
    int cnt = 0;
    Real totalPhi = 0;
    for(int i_spec=0;i_spec<nspecies;i_spec++)
    {
        totalPhi += phi_domain[i_spec];
        for(int j_spec=0;j_spec<nspecies;j_spec++)
        {
            ij_spec = getSpeciesIndex(i_spec,j_spec);
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
        const Array4<Real> & arrselect = mfselect.array(pti);

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            int ij_spec;
            for (int i_spec=0; i_spec<nspecies; i_spec++) {
                for (int j_spec = i_spec; j_spec < nspecies; j_spec++) {
                    ij_spec = getSpeciesIndex(i_spec,j_spec);
                    arrselect(i,j,k,ij_spec) = 0.0;
                }
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
        const Array4<Real> & arrselect = mfselect.array(mfi);

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            int ij_spec;
            long np_i, np_j;

            const IntVect& iv = {i,j,k};
            long imap = tile_box.index(iv);

            Real vrmax;
            Real NSel;
            Real crossSection;

            for (int i_spec = 0; i_spec < nspecies; i_spec++)
            {
                for (int j_spec = i_spec; j_spec < nspecies; j_spec++) {
                    ij_spec = getSpeciesIndex(i_spec,j_spec);
                    np_i = m_cell_vectors[i_spec][grid_id][imap].size();
                    np_j = m_cell_vectors[j_spec][grid_id][imap].size();
                    vrmax = arrvrmax(i,j,k,i_spec);
                    crossSection = interproperties[ij_spec].csx;
                    NSel = 4.0*particle_neff*np_i*np_j*crossSection*vrmax*ocollisionCellVol*dt;
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
            int totalCol = 0;
            const IntVect& iv = {i,j,k};
            long imap = tile_box.index(iv);

            Real NSel[nspecies*nspecies], totalSel, selrun;
            int np[nspecies];
            int pindxi, pindxj;
            int ij_spec;

            RealVect eij, vreij;
            RealVect vi, vj, vij;
            RealVect vijpost, boost;
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
                    vreijmag = vreijmag*2.0/massij;
                    vreij[0] = vreijmag*eij[0];
                    vreij[1] = vreijmag*eij[1];
                    vreij[2] = vreijmag*eij[2];

                    parti.rdata(FHD_realData::velx) = vi[0] - vreij[0]*massj;
                    parti.rdata(FHD_realData::vely) = vi[1] - vreij[1]*massj;
                    parti.rdata(FHD_realData::velz) = vi[2] - vreij[2]*massj;
                    partj.rdata(FHD_realData::velx) = vj[0] + vreij[0]*massi;
                    partj.rdata(FHD_realData::vely) = vj[1] + vreij[1]*massi;
                    partj.rdata(FHD_realData::velz) = vj[2] + vreij[2]*massi;
                    totalCol++;
                }
            }
            //Print() << "total cols: " << totalCol << "\n";
        }
        }
        }
    }
}