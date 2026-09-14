#include "AMReX_PlotFileUtil.H"
#include "AMReX_PlotFileDataImpl.H"
#include "DsmcParticleContainer.H"
#include <sys/stat.h>

#include "common_functions.H"
#include "common_namespace.H"
#include "chrono"

#include "INS_functions.H"
#include "common_functions.H"
#include <sstream>
#include <string>
#include <fstream>

using namespace std::chrono;
using namespace common;

namespace {
    void GotoNextLine (std::istream& is) {
        constexpr std::streamsize bl_ignore_max { 100000 };
        is.ignore(bl_ignore_max, '\n');
    }
}

void WriteCheckPoint(int step,
                     const Real time,
                     const Real dt,
                     int statsCount,
                     const amrex::MultiFab& cuInst,
                     const amrex::MultiFab& cuMeans,
                     const amrex::MultiFab& cuVars,
                     const amrex::MultiFab& primInst,
                     const amrex::MultiFab& primMeans,
                     const amrex::MultiFab& primVars,
                     const amrex::MultiFab& coVars,
                     const FhdParticleContainer& particles,
                     const amrex::MultiFab& spatialCross1D,
                     const int ncross) {

    // timer for profiling
    BL_PROFILE_VAR("WriteCheckPoint()",WriteCheckPoint);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname =
        amrex::Concatenate(chk_base_name,step,12);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    BoxArray ba = cuInst.boxArray();

    // single level problem
    int nlevels = 1;

    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);
    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    //////////////////////////////////////
    // Header File
    //////////////////////////////////////
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        std::string HeaderFileName(checkpointname + "/Header");
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                          std::ofstream::trunc |
                          std::ofstream::binary);

        if( !HeaderFile.good()) { amrex::FileOpenFailed(HeaderFileName); }

        HeaderFile.precision(17);
        HeaderFile << "Checkpoint file for FHDeX/dsmc\n"; // write out title line
        HeaderFile << step << "\n"; // write out the time step number
        HeaderFile << time << "\n"; // write out time
        HeaderFile << dt << "\n"; // write out time
        HeaderFile << statsCount << "\n"; // write out statsCount
        ba.writeOn(HeaderFile); // write the BoxArray (fluid)
        HeaderFile << '\n';
    }

    int comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    int n_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

    //////////////////////////////////////
    // Save RNG State
    //////////////////////////////////////
    // don't write out all the rng states at once (overload filesystem)
    // one at a time write out the rng states to different files, one for each MPI rank
    for (int rank=0; rank<n_ranks; ++rank) {
        if (comm_rank == rank) {
            std::ofstream rngFile;
            rngFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

            // create filename, e.g. chk0000005/rng0000002
            const std::string& rngFileNameBase = (checkpointname + "/rng");
            const std::string& rngFileName = amrex::Concatenate(rngFileNameBase,comm_rank,7);

            rngFile.open(rngFileName.c_str(), std::ofstream::out   |
              std::ofstream::trunc |
              std::ofstream::binary);

            if( !rngFile.good()) {
                amrex::FileOpenFailed(rngFileName);
            }
            amrex::SaveRandomState(rngFile);
        }
        ParallelDescriptor::Barrier();
    }

    //////////////////////////////////////
    // Record MF and Particles
    //////////////////////////////////////

    // Stat MFs
    VisMF::Write(cuInst,
        amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuInst"));
    VisMF::Write(cuMeans,
        amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuMeans"));
    VisMF::Write(cuVars,
        amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuVars"));
    VisMF::Write(primInst,
        amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primInst"));
    VisMF::Write(primMeans,
        amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primMeans"));
    VisMF::Write(primVars,
        amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primVars"));
    VisMF::Write(coVars,
        amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "coVars"));
    VisMF::Write(spatialCross1D,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "spatialCross1D"));

    // checkpoint particles
    particles.Checkpoint(checkpointname,"particle");
}

void ReadCheckPoint(int& step,
                    Real& time,
                    Real& dt,
                    int& statsCount,
                    amrex::MultiFab& cuInst,
                    amrex::MultiFab& cuMeans,
                    amrex::MultiFab& cuVars,
                    amrex::MultiFab& primInst,
                    amrex::MultiFab& primMeans,
                    amrex::MultiFab& primVars,
                    amrex::MultiFab& coVars,
                    amrex::MultiFab& spatialCross1D,
                    const int ncross) {

    // timer for profiling
    BL_PROFILE_VAR("ReadCheckPoint()",ReadCheckPoint);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname =
        amrex::Concatenate(chk_base_name,restart,12);

    amrex::Print() << "Restart from checkpoint " << checkpointname << "\n";

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    std::string line, word;

    // Header
    {
        std::string File(checkpointname + "/Header");
        Vector<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream is(fileCharPtrString, std::istringstream::in);

        // read in title line
        std::getline(is, line);

        // read in time step number
        is >> step;
        ++step;
        GotoNextLine(is);

        // read in time
        is >> time;
        GotoNextLine(is);

        // read in dt
        is >> dt;
        GotoNextLine(is);

        // read in statsCount
        is >> statsCount;
        GotoNextLine(is);

        // read in BoxArray (fluid) from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // MFs
        // MAKE SURE TO UPDATE IF INCREASED NUMBER OF VARS IN MF
        int cnvars  = (nspecies+1)*5;
        int pnvars  = (nspecies+1)*9;
        int convars = 25;
        cuInst.define(ba,dm,cnvars,0);
        cuMeans.define(ba,dm,cnvars,0);
        cuVars.define(ba,dm,cnvars,0);
        primInst.define(ba,dm,pnvars,0);
        primMeans.define(ba,dm,pnvars,0);
        primVars.define(ba,dm,pnvars,0);
        coVars.define(ba,dm,convars,0);
        spatialCross1D.define(ba,dm,ncross,0);
    }
    // C++ random number engine
    // each MPI process reads in its own file
    int comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    int n_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

    if (seed < 0) {

#ifdef AMREX_USE_CUDA
        Abort("Restart with negative seed not supported on GPU");
#endif

        // read in rng state from checkpoint
        // don't read in all the rng states at once (overload filesystem)
        // one at a time read the rng states to different files, one for each MPI rank
        for (int rank=0; rank<n_ranks; ++rank) {

            if (comm_rank == rank) {

                if (seed < 0) {
                    // create filename, e.g. chk0000005/rng0000002
                    std::string FileBase(checkpointname + "/rng");
                    std::string File = amrex::Concatenate(FileBase,comm_rank,7);

                    // read in contents
                    Vector<char> fileCharPtr;
                    ReadFile(File, fileCharPtr);
                    std::string fileCharPtrString(fileCharPtr.dataPtr());
                    std::istringstream is(fileCharPtrString, std::istringstream::in);

                    // restore random state
                    amrex::RestoreRandomState(is, 1, 0);
                }
            }

            ParallelDescriptor::Barrier();

        }

    }
    else if (seed == 0) {
        // initializes seed for random number calls from clock
        auto now = time_point_cast<nanoseconds>(system_clock::now());
        int randSeed = now.time_since_epoch().count();
        // broadcast the same root seed to all processors
        ParallelDescriptor::Bcast(&randSeed,1,ParallelDescriptor::IOProcessorNumber());
        InitRandom(randSeed+ParallelDescriptor::MyProc(),
                   ParallelDescriptor::NProcs(),
                   randSeed+ParallelDescriptor::MyProc());
    } else {
        // initializes the seed for C++ random number calls
        InitRandom(seed+ParallelDescriptor::MyProc(),
                   ParallelDescriptor::NProcs(),
                   seed+ParallelDescriptor::MyProc());
    }

    // read in the MultiFab data
    VisMF::Read(cuInst,
        amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuInst"));
    VisMF::Read(primInst,
        amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primInst"));

    // Set all stats to zero if reset stats, else read
    if (reset_stats == 1) {
        cuMeans.setVal(0.0);
        cuVars.setVal(0.0);
        primMeans.setVal(0.0);
        primVars.setVal(0.0);
        coVars.setVal(0.0);
        spatialCross1D.setVal(0.0);
    }
    else {
        VisMF::Read(cuMeans,amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuMeans"));
        VisMF::Read(cuVars,amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuVars"));
        VisMF::Read(primMeans,amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primMeans"));
        VisMF::Read(primVars,amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primVars"));
        VisMF::Read(coVars,amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "coVars"));
        VisMF::Read(spatialCross1D,amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "spatialCross1D"));
    }
}

void ReadCheckPointParticles(FhdParticleContainer& particles) {

    // timer for profiling
    BL_PROFILE_VAR("ReadCheckPointParticles()",ReadCheckPointParticles);

    std::string checkpointname;

    if (particle_restart > 0) {
        checkpointname = amrex::Concatenate(chk_base_name,particle_restart,12);
    }
    else {
        checkpointname = amrex::Concatenate(chk_base_name,restart,12);
    }

    amrex::Print() << "Restart particles from checkpoint " << checkpointname << "\n";

    std::string line, word;

    int temp;
    std::string File(checkpointname + "/Header");
    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    // read in title line
    std::getline(is, line);

    // read in time step number
    is >> temp;
    GotoNextLine(is);

    // read in time
    is >> temp;
    GotoNextLine(is);

    // read in dt
    is >> temp;
    GotoNextLine(is);

    // read in statsCount
    is >> temp;
    GotoNextLine(is);

    // read in BoxArray (fluid) from Header
    BoxArray ba;
    ba.readFrom(is);
    GotoNextLine(is);

    // create a distribution mapping
    DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

    RealBox realDomain({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                       {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    //cout << "Restoring!\n";
    particles.Restart(checkpointname,"particle");
    //cout << "Restored!\n";
    int np = particles.TotalNumberOfParticles();
    //particlesTemp.Checkpoint("testcheck","particle");
    Print() << "Checkpoint contains " << np << " particles." <<std::endl;
    particles.ReInitParticles();
}

void ReadFile(const std::string& filename, Vector<char>& charBuf,
              bool bExitOnError) {

    enum { IO_Buffer_Size = 262144 * 8 };

#ifdef BL_SETBUF_SIGNED_CHAR
    typedef signed char Setbuf_Char_Type;
#else
    typedef char Setbuf_Char_Type;
#endif

    Vector<Setbuf_Char_Type> io_buffer(IO_Buffer_Size);

    Long fileLength(0), fileLengthPadded(0);

    std::ifstream iss;

    iss.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    iss.open(filename.c_str(), std::ios::in);
    if ( ! iss.good()) {
        if(bExitOnError) {
            amrex::FileOpenFailed(filename);
        } else {
            fileLength = -1;
        }
    } else {
        iss.seekg(0, std::ios::end);
        fileLength = static_cast<std::streamoff>(iss.tellg());
        iss.seekg(0, std::ios::beg);
    }

    if(fileLength == -1) {
        return;
    }

    fileLengthPadded = fileLength + 1;
    //    fileLengthPadded += fileLengthPadded % 8;
    charBuf.resize(fileLengthPadded);

    iss.read(charBuf.dataPtr(), fileLength);
    iss.close();

    charBuf[fileLength] = '\0';
}