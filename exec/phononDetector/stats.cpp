#include "DsmcParticleContainer.H"
using namespace std;
void FhdParticleContainer::EvaluateStatsPhonon(MultiFab& mfcuInst,
    MultiFab& mfcuMeans,
    MultiFab& mfcuVars,
    const int steps,
    Real time) {
    BL_PROFILE_VAR("EvaluateStats()",EvaluateStats);
    const Real osteps = 1.0/steps;
    const Real stepsMinusOne = steps-1.;

    // Zero out instantaneous values
    mfcuInst.setVal(0.);
    const int lev = 0;

    int ncon  = 5;

    const Real* dx = Geom(lev).CellSize();
    const Real dxInv = 1.0/dx[0];
    const Real cellVolInv = 1.0/(dx[0]*dx[0]*dx[0]);

    for (FhdParIter pti(* this, lev); pti.isValid(); ++pti) {
        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();
        const Box& tile_box  = pti.tilebox();

        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& aos = particle_tile.GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();

        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();


        Array4<Real> cuInst  = mfcuInst[pti].array();
        Array4<Real> cuMeans  = mfcuMeans[pti].array();
        Array4<Real> cuVars  = mfcuVars[pti].array();

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

            for (int l=0; l<nspecies; l++) {
               // long np_spec = m_cell_vectors[l][grid_id][imap].size();
                unsigned int np_spec = getBinSize(offs,iv,l,tile_box);
                unsigned int* cellList = getCellList(inds,offs,iv,l,tile_box);

                cuInst(i,j,k,0) += np_spec;

                // Read particle data
                for (int m=0; m<np_spec; m++) {
                    //int pind = m_cell_vectors[l][grid_id][imap][m];
                    int pind = cellList[m];
                    //int pind = 1;
                    ParticleType & p = particles[pind];

                    Real u = p.rdata(FHD_realData::velx);
                    Real v = p.rdata(FHD_realData::vely);
                    Real w = p.rdata(FHD_realData::velz);

                    cuInst(i,j,k,1) += u*h_bar*p.rdata(FHD_realData::omega);
                    cuInst(i,j,k,2) += v*h_bar*p.rdata(FHD_realData::omega);
                    cuInst(i,j,k,3) += w*h_bar*p.rdata(FHD_realData::omega);
                    cuInst(i,j,k,4) += h_bar*p.rdata(FHD_realData::omega);
                }


            }
            for (int l=0; l<ncon; l++) {
                cuInst(i,j,k,l) *= cellVolInv;
            }
        });

        //////////////////////////////////////
        // Means
        //////////////////////////////////////
        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int l=0; l<ncon; l++) {
                cuMeans(i,j,k,l) = (cuMeans(i,j,k,l)*stepsMinusOne+cuInst(i,j,k,l))*osteps;
            }
        });

        amrex::ParallelFor(tile_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Conserved Variances
            //Vector<Real> delCon(ncon, 0.0);
            for (int l=0; l<ncon; l++) {
                Real delCon        = cuInst(i,j,k,l) - cuMeans(i,j,k,l);
                cuVars(i,j,k,l)  = (cuVars(i,j,k,l)*stepsMinusOne+delCon*delCon)*osteps;
            }

        });
    }
}