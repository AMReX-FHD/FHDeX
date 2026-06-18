#include "DsmcParticleContainer.H"
using namespace std;
void globalVars(Real time) {
    BL_PROFILE_VAR("globalVars()",EvaluateStats);

    const int lev = 0;
    for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi) {
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        const Box& tile_box  = mfi.tilebox();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id,tile_id)];
        auto& particles = particle_tile.GetArrayOfStructs();
        IntVect smallEnd = tile_box.smallEnd();
        IntVect bigEnd = tile_box.bigEnd();


        //////////////////////////////////////
        // Global Temperature
        //////////////////////////////////////

        Real Tgl[nspecies], npl[nspecies];
        for (int l=0; l<nspecies; l++) {Tgl[l] = 0.; npl[l] = 0.;}
        int np = particles.numParticles();

        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];

            int ispec = part.idata(FHD_intData::species);
            Real vsq = pow(part.rdata(FHD_realData::velx),2) +
                pow(part.rdata(FHD_realData::vely),2) +
                pow(part.rdata(FHD_realData::velz),2);
            Tgl[ispec] += vsq; npl[ispec] += 1;
        }
        for (int l=0; i<nspecies; l++) {
            Tgl[ispec] *= (properties[ispec].mass/3.0);
        }

        // Gather from all proc
        // Likely better way to do this
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
            fileTg.open(Tgfname, fstream::app);
            fileTgN.open(TgNfname, fstream::app);
            fileTg << std::scientific << setprecision(8) << time << " ";
            fileTgN << std::scientific << setprecision(8) << time << " ";
            for(int l=0; l<nspecies; l++){
                if(steps==1) {Tg0[l] = Tgl[l];}
                fileTg << Tgl[l] << " ";
                fileTgN << Tgl[l]/Tg0[l] << " ";
            }
            fileTg << "\n"; fileTgN << "\n";
            fileTg.close(); fileTgN.close();
        }
    }
}