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

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void getSpeciesIndexRet(int species1, int species2, int* ret)
{
    if(species1<species2)
    {
        ret[0] = species2+nspecies*species1;
    } else
    {
        ret[0] = species1+nspecies*species2;
    }
}

void FhdParticleContainer::InitCollisionCells()
{
    BL_PROFILE_VAR("InitCollisionCells()",InitCollisionCells);
    int indx;
    int cnt = 0;
    Real totalPhi = 0;
    for(int i_spec=0;i_spec<nspecies;i_spec++)
    {
        int ij_spec;
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
                    //ij_spec = getSpeciesIndex(i_spec,j_spec);
                    getSpeciesIndexRet(i_spec,j_spec, &ij_spec);
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

        auto inds = m_bins.permutationPtr();
        auto offs = m_bins.offsetsPtr();

        Real ocollisionCellVolTmp = ocollisionCellVol;
        Real particle_neff_tmp = particle_neff;

        Real csx[MAX_SPECIES*MAX_SPECIES];

        for(int i=0;i<(nspecies*nspecies);i++)
        {
            csx[i] = interproperties[i].csx;
        }

        amrex::ParallelForRNG(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept {
            int ij_spec;
            long np_i, np_j;

            const IntVect& iv = {i,j,k};
            long imap = tile_box.index(iv);

            Real vrmax;
            Real NSel;
            Real crossSection;

            for (int i_spec = 0; i_spec<nspecies; i_spec++)
            {
                for (int j_spec = i_spec; j_spec < nspecies; j_spec++) {
                    //ij_spec = getSpeciesIndex(i_spec,j_spec);
                    getSpeciesIndexRet(i_spec,j_spec,&ij_spec);
                    //np_i = m_cell_vectors[i_spec][grid_id][imap].size();
                    //np_j = m_cell_vectors[j_spec][grid_id][imap].size();
                    np_i = getBinSize(offs,iv,i_spec,tile_box);
                    np_j = getBinSize(offs,iv,j_spec,tile_box);
                    vrmax = arrvrmax(i,j,k,ij_spec);
                    crossSection = csx[ij_spec];
                    //crossSection = 0;
                    if(i_spec==j_spec) {np_j = np_i-1;}
                    NSel = particle_neff_tmp*np_i*np_j*crossSection*vrmax*ocollisionCellVolTmp*dt*2;
                    if(i_spec==j_spec) {NSel = NSel*0.5;}
                    arrselect(i,j,k,ij_spec) = std::floor(NSel + amrex::Random(engine));

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
        auto& aos = particle_tile.GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();
        //const long np = particle_tile.numParticles();
        //auto& particles = particle_tile.GetArrayOfStructs();

        const Array4<Real> & arrvrmax = mfvrmax.array(mfi);
        const Array4<Real> & arrselect = mfselect.array(mfi);

        auto inds = m_bins.permutationPtr();
        auto offs = m_bins.offsetsPtr();


        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();

        dsmcSpecies* propertiesTmp = properties.data();

        Real mass[MAX_SPECIES];

        for(int i=0;i<(nspecies);i++)
        {
            mass[i] = properties[i].mass;
        }

        //for (int i = smallEnd[0]; i <= bigEnd[0]; i++) {
        //for (int j = smallEnd[1]; j <= bigEnd[1]; j++) {
        //for (int k = smallEnd[2]; k <= bigEnd[2]; k++) {
        amrex::ParallelForRNG(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept {
            int totalCol = 0;
            const IntVect& iv = {i,j,k};
            long imap = tile_box.index(iv);

            Real NSel[MAX_SPECIES*MAX_SPECIES], totalSel, selrun;
            int np[MAX_SPECIES];
            int pindxi, pindxj, bini, binj;
            int ij_spec;


            RealVect eij, vreij;
            RealVect vi, vj, vij;
            RealVect vijpost, boost;
            Real massi, massj, massij;
            Real vrmag, vrmax, vreijmag;

            unsigned int* specLists[MAX_SPECIES];

            for(int i=0;i<nspecies;i++)
            {
                specLists[i] = getCellList(inds,offs,iv,i,tile_box);
            }

            totalSel = 0;
            //for (int i_spec = nspecies-1; i_spec>=0; i_spec--)
            for(int i_spec = 0; i_spec<nspecies; i_spec++)
            {
                //np[i_spec] = m_cell_vectors[i_spec][grid_id][imap].size();
                np[i_spec] = getBinSize(offs,iv,i_spec,tile_box);
                for (int j_spec = i_spec; j_spec < nspecies; j_spec++)
                {
                    //ij_spec = getSpeciesIndex(i_spec,j_spec);
                    getSpeciesIndexRet(i_spec,j_spec, &ij_spec);
                    NSel[ij_spec] = (int)arrselect(i,j,k,ij_spec);
                    totalSel += NSel[ij_spec];

                }

            }

            int totalSelect2 = totalSel;
            int speci, specj, specij;


            while (totalSel>0)
            {
                Real RR = amrex::Random(engine);
                bool spec_select = false;
                selrun = 0;
                speci = -1; specj = -1; specij = -1;


                //for (int i_spec = nspecies-1; i_spec>=0; i_spec--)
                for(int i_spec = 0; i_spec<nspecies; i_spec++)
                {
                    for(int j_spec=i_spec;j_spec<nspecies;j_spec++)
                    {
                        //int j_spec = i_spec;
                        getSpeciesIndexRet(i_spec,j_spec,&ij_spec);
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
                massi = mass[speci];
                massj = mass[specj];
                massij = mass[speci] + mass[specj];
                vrmax = arrvrmax(i,j,k,specij);
                //pindxi = (int)floor(amrex::Random()*m_cell_vectors[speci][grid_id][imap].size());
                //pindxj = (int)floor(amrex::Random()*m_cell_vectors[specj][grid_id][imap].size());
                pindxi = (int)floor(amrex::Random(engine)*getBinSize(offs,iv,speci,tile_box));
                pindxj = (int)floor(amrex::Random(engine)*getBinSize(offs,iv,specj,tile_box));
                pindxi = specLists[speci][pindxi];
                pindxj = specLists[specj][pindxj];
                //pindxi = m_cell_vectors[speci][grid_id][imap][pindxi];
                //pindxj = m_cell_vectors[specj][grid_id][imap][pindxj];

                //if(i==1){Print() << "spec " << speci << " part " << pindxi << " of "  << m_cell_vectors[speci][grid_id][imap].size() << endl;}
                //if(i==1){Print() << "spec " << specj << " part " << pindxj << " of "  << m_cell_vectors[specj][grid_id][imap].size() << endl;}

                ParticleType & parti = particles[pindxi];
                ParticleType & partj = particles[pindxj];

                vi[0] = parti.rdata(FHD_realData::velx);
                vi[1] = parti.rdata(FHD_realData::vely);
                vi[2] = parti.rdata(FHD_realData::velz);

                vj[0] = partj.rdata(FHD_realData::velx);
                vj[1] = partj.rdata(FHD_realData::vely);
                vj[2] = partj.rdata(FHD_realData::velz);

                //if(i==1){Print() << "vel1 " << vi[0] << " vel2 " << vj[0] << endl;}

                vij[0] = vi[0]-vj[0]; vij[1] = vi[1]-vj[1]; vij[2] = vi[2]-vj[2];
                vrmag = sqrt(pow(vij[0],2)+pow(vij[1],2)+pow(vij[2],2));
                if(vrmag>vrmax) {vrmax = vrmag; arrvrmax(i,j,k,specij) = 1.1*vrmax;}

                Real eijmag = 0;
                Real theta = 2.0*M_PI*amrex::Random(engine);
                Real phi = std::acos(1.0-2.0*amrex::Random(engine));
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
                if(vrmag>vrmax*amrex::Random(engine))
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

            //        }
            //        }
            //        }
        });
    }
}

void FhdParticleContainer::CollideParticles2(Real dt)
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

        //const long np = particles.numParticles();
        //amrex::ParallelForRNG(tile_box,
        //  [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept {
        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();

        for (int i = smallEnd[0]; i <= bigEnd[0]; i++) {
        for (int j = smallEnd[1]; j <= bigEnd[1]; j++) {
        for (int k = smallEnd[2]; k <= bigEnd[2]; k++) {

            const IntVect& iv = {i,j,k};
            long imap = tile_box.index(iv);
            int Nselect[nspecies*nspecies];
            Real Rselect[nspecies*nspecies];

            for(int k=0; k<nspecies*nspecies; k++)
            {
                Nselect[k] = 0;
            }

            for(int i_spec = 0; i_spec<nspecies; i_spec++)
            {
                for(int j_spec = 0; j_spec<nspecies; j_spec++)
                {
                   int ij_index = getSpeciesIndex(i_spec,j_spec);
                   int np_i = m_cell_vectors[i_spec][grid_id][imap].size();
                   int np_j = m_cell_vectors[j_spec][grid_id][imap].size();
                   if(i_spec == j_spec){np_i--;}
                   Real select = np_i*np_j*interproperties[ij_index].csx*particle_neff*arrvrmax(i,j,k,ij_index)*ocollisionCellVol*dt;
                   if(i_spec != j_spec){select = select*0.5;}
                   Rselect[ij_index] += select;
                }
            }

            for(int i_spec = 0; i_spec<nspecies; i_spec++)
            {
                for(int j_spec = 0; j_spec<nspecies; j_spec++)
                {
                   int ij_index = getSpeciesIndex(i_spec,j_spec);
                   Nselect[ij_index] = floor(Rselect[ij_index] + amrex::Random());
                }
            }

            int totalCol = 0;
            for(int k=0; k<nspecies*nspecies; k++)
            {
                totalCol += Nselect[k];
            }

            while(totalCol > 0)
            {
                Real RR = amrex::Random();
                bool spec_select = false;
                int selrun = 0;
                int speci = -1;
                int specj = -1;
                int specij = -1;

                for(int i_spec = 0; i_spec<nspecies; i_spec++)
                {
                    for(int j_spec=i_spec;j_spec<nspecies;j_spec++)
                    {

                        int ij_spec = getSpeciesIndex(i_spec,j_spec);
                        selrun += Nselect[ij_spec];

                        if(selrun/totalCol>RR && !spec_select)
                        {
                            spec_select = true;
                            specij = ij_spec;
                            speci = i_spec; specj = j_spec;


                            Nselect[ij_spec] -= 1;
                        }
                    }
                }
                totalCol--;

                Real massi = properties[speci].mass;
                Real massj = properties[specj].mass;
                Real massij = properties[speci].mass + properties[specj].mass;
                //vrmax = arrvrmax(i,j,k,specij);
                int pindxi = (int)floor(amrex::Random()*m_cell_vectors[speci][grid_id][imap].size());
                int pindxj = (int)floor(amrex::Random()*m_cell_vectors[specj][grid_id][imap].size());
                pindxi = m_cell_vectors[speci][grid_id][imap][pindxi];
                pindxj = m_cell_vectors[specj][grid_id][imap][pindxj];

                ParticleType & parti = particles[pindxi];
                ParticleType & partj = particles[pindxj];

                Real vi[3];
                Real vj[3];
                Real vij[3];
                Real eij[3];

                vi[0] = parti.rdata(FHD_realData::velx);
                vi[1] = parti.rdata(FHD_realData::vely);
                vi[2] = parti.rdata(FHD_realData::velz);

                vj[0] = partj.rdata(FHD_realData::velx);
                vj[1] = partj.rdata(FHD_realData::vely);
                vj[2] = partj.rdata(FHD_realData::velz);

                vij[0] = vi[0]-vj[0]; vij[1] = vi[1]-vj[1]; vij[2] = vi[2]-vj[2];
                Real vrmag = sqrt(pow(vij[0],2)+pow(vij[1],2)+pow(vij[2],2));
                if(vrmag>arrvrmax(i,j,k,specij)) {arrvrmax(i,j,k,specij) = 1.1*vrmag;}

                Real theta = 2.0*pi_usr*amrex::Random();
                Real phi = std::acos(1.0-2.0*amrex::Random());
                eij[0] = std::sin(phi)*std::cos(theta);
                eij[1] = std::sin(phi)*std::sin(theta);
                eij[2] = std::cos(phi);
                Real eijmag = pow(eij[0],2)+pow(eij[1],2)+pow(eij[2],2);
                eijmag = pow(eijmag,0.5);
                for(int idim=0; idim<3; idim++)
                {
                    eij[idim] /= eijmag;
                }

                Real vreijmag = vij[0]*eij[0]+vij[1]*eij[1]+vij[2]*eij[2];
                if(vrmag>(arrvrmax(i,j,k,specij)*amrex::Random()))
                {
                    Real vreij[3];
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


        }
        }
        }
    }
}