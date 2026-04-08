//#include "Checkpoint.H"
#include "INS_functions.H"

#include "AMReX_PlotFileUtil.H"
#include "AMReX_PlotFileDataImpl.H"

#include <sys/stat.h>
#include <chrono>

using namespace std::chrono;

namespace {
    void GotoNextLine (std::istream& is)
    {
        constexpr std::streamsize bl_ignore_max { 100000 };
        is.ignore(bl_ignore_max, '\n');
    }
}

void WriteCheckPoint(int step,
                     amrex::Real time,
                     int statsCount,
                     const std::array< MultiFab, AMREX_SPACEDIM >& umac,
                     const std::array< MultiFab, AMREX_SPACEDIM >& umacM,
                     const MultiFab& pres,
                     const FhdParticleContainer& particles,
                     const MultiFab& particleMeans,
                     const MultiFab& particleVars,
                     const MultiFab& chargeM,
                     const MultiFab& potential,
                     const MultiFab& potentialM,
                     const MultiFab& struct_cc_numdens0_real,
                     const MultiFab& struct_cc_numdens0_imag)
{
    // timer for profiling
    BL_PROFILE_VAR("WriteCheckPoint()",WriteCheckPoint);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate(chk_base_name,step,9);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    BoxArray ba = pres.boxArray();
    BoxArray bc = particleMeans.boxArray();
    BoxArray bp = potential.boxArray();

    // single level problem
    int nlevels = 1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    // write Header file
    if (ParallelDescriptor::IOProcessor()) {

        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        std::string HeaderFileName(checkpointname + "/Header");
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                        std::ofstream::trunc |
                        std::ofstream::binary);

        if( !HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);

        // write out title line
        HeaderFile << "Checkpoint file for FHDeX/hydro\n";

        // write out the time step number
        HeaderFile << step << "\n";

        // write out time
        HeaderFile << time << "\n";

        // write out statsCount
        HeaderFile << statsCount << "\n";

        // write the BoxArray (fluid)
        ba.writeOn(HeaderFile);
        HeaderFile << '\n';

        // write the BoxArray (particle)
        bc.writeOn(HeaderFile);
        HeaderFile << '\n';

        // write the BoxArray (electric potential)
        bp.writeOn(HeaderFile);
        HeaderFile << '\n';
    }

    // C++ random number engine
    // have each MPI process write its random number state to a different file
    int comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    int n_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

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

    // write the MultiFab data to, e.g., chk00010/Level_0/

    // umac
    VisMF::Write(umac[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "umac"));
    VisMF::Write(umac[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "vmac"));
#if (AMREX_SPACEDIM == 3)
    VisMF::Write(umac[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "wmac"));
#endif

    // umacM (mean)
    VisMF::Write(umacM[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "umacM"));
    VisMF::Write(umacM[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "vmacM"));
#if (AMREX_SPACEDIM == 3)
    VisMF::Write(umacM[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "wmacM"));
#endif

    // pressure
    VisMF::Write(pres,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "pressure"));

    // particle mean and variance
    VisMF::Write(particleMeans,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "particleMeans"));
    //VisMF::Write(particleVars,
    //             amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "particleVars"));

    // charge
    VisMF::Write(chargeM,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "chargeM"));

    // electrostatic potential
    VisMF::Write(potential,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "potential"));
    VisMF::Write(potentialM,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "potentialM"));

    // initial dsf state
    VisMF::Write(struct_cc_numdens0_real,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "struct_cc_numdens0_real"));
    VisMF::Write(struct_cc_numdens0_imag,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "struct_cc_numdens0_imag"));

    // checkpoint particles
    particles.Checkpoint(checkpointname,"particle");
}

void ReadCheckPoint(int& step,
                    amrex::Real& time,
                    int& statsCount,
                    std::array< MultiFab, AMREX_SPACEDIM >& umac,
                    std::array< MultiFab, AMREX_SPACEDIM >& umacM,
                    MultiFab& pres,
                    MultiFab& particleMeans,
                    MultiFab& particleVars,
                    MultiFab& chargeM,
                    MultiFab& potential,
                    MultiFab& potentialM,
                    MultiFab& struct_cc_numdens0_real,
                    MultiFab& struct_cc_numdens0_imag)
{
    // timer for profiling
    BL_PROFILE_VAR("ReadCheckPoint()",ReadCheckPoint);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate(chk_base_name,restart,9);

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

        // read in statsCount
        is >> statsCount;
        GotoNextLine(is);

        // read in BoxArray (fluid) from Header
        BoxArray ba;
        ba.readFrom(is);
        GotoNextLine(is);

        // read in BoxArray (particle) from Header
        BoxArray bc;
        bc.readFrom(is);
        GotoNextLine(is);

        BoxArray bp;
        bp.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        //set number of ghost cells to fit whole peskin kernel
        int ang = 1;
        if(*(std::max_element(pkernel_fluid.begin(),pkernel_fluid.begin()+nspecies)) == 3) {
            ang = 2;
        }
        else if(*(std::max_element(pkernel_fluid.begin(),pkernel_fluid.begin()+nspecies)) == 4) {
            ang = 3;
        }
        else if(*(std::max_element(pkernel_fluid.begin(),pkernel_fluid.begin()+nspecies)) == 6) {
            ang = 4;
        }
        else if (*(std::max_element(eskernel_fluid.begin(),eskernel_fluid.begin()+nspecies)) > 0) {
            ang = static_cast<int>(floor(*(std::max_element(eskernel_fluid.begin(),eskernel_fluid.begin()+nspecies)))/2+1);
        }

        // build MultiFab data

        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            umac [d].define(convert(ba,nodal_flag_dir[d]), dm, 1, ang);
            umacM[d].define(convert(ba,nodal_flag_dir[d]), dm, 1, 1);
        }

        // pressure
        pres.define(ba,dm,1,1);

        // particle means and variances
        particleMeans.define(bc,dm,14,0);
        particleVars .define(bc,dm,18,0);

        // cell centred es potential
        int ngp = 1;
        // using maximum number of peskin kernel points to determine the ghost cells for the whole grid.
        //     not sure if it will cause problem for BCs.
        if (*(std::max_element(pkernel_es.begin(),pkernel_es.begin()+nspecies)) == 3) {
            ngp = 2;
        }
        else if (*(std::max_element(pkernel_es.begin(),pkernel_es.begin()+nspecies)) == 4) {
            ngp = 3;
        }
        else if (*(std::max_element(pkernel_es.begin(),pkernel_es.begin()+nspecies)) == 6) {
            ngp = 4;
        }
        else if (*(std::max_element(eskernel_fluid.begin(),eskernel_fluid.begin()+nspecies)) > 0) {
            ngp = static_cast<int>(floor(*(std::max_element(eskernel_fluid.begin(),eskernel_fluid.begin()+nspecies)))/2+1);
        }
        //// TODO: need a better way to determine ghost cells for bonds
        //if (bond_tog != 0) {
        //    ngp = std::max(ngp, 6);
        //}

        chargeM.define(bp,dm,1,1);
        potential.define(bp,dm,1,ngp);
        potentialM.define(bp,dm,1,1);
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

        // initializes the seed for C++ random number calls based on the clock
        auto now = time_point_cast<nanoseconds>(system_clock::now());
        int randSeed = now.time_since_epoch().count();
        // broadcast the same root seed to all processors
        ParallelDescriptor::Bcast(&randSeed,1,ParallelDescriptor::IOProcessorNumber());

        InitRandom(randSeed+ParallelDescriptor::MyProc(),
                   ParallelDescriptor::NProcs(),
                   randSeed+ParallelDescriptor::MyProc());
    }
    else {
        // initializes the seed for C++ random number calls
        InitRandom(seed+ParallelDescriptor::MyProc(),
                   ParallelDescriptor::NProcs(),
                   seed+ParallelDescriptor::MyProc());
    }

    // read in the MultiFab data

    // umac
    VisMF::Read(umac[0],
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "umac"));
    VisMF::Read(umac[1],
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "vmac"));
#if (AMREX_SPACEDIM == 3)
    VisMF::Read(umac[2],
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "wmac"));
#endif

    // umacM
    VisMF::Read(umacM[0],
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "umacM"));
    VisMF::Read(umacM[1],
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "vmacM"));
#if (AMREX_SPACEDIM == 3)
    VisMF::Read(umacM[2],
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "wmacM"));
#endif

    // pressure
    VisMF::Read(pres,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "pressure"));

    // particle means and variances
    VisMF::Read(particleMeans,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "particleMeans"));
    //VisMF::Read(particleVars,
    //            amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "particleVars"));

    // charge
    VisMF::Read(chargeM,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "chargeM"));

    // electrostatic potential
    VisMF::Read(potential,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "potential"));
    VisMF::Read(potentialM,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "potentialM"));

    // initial dsf state
    VisMF::Read(struct_cc_numdens0_real,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "struct_cc_numdens0_real"));
    VisMF::Read(struct_cc_numdens0_imag,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "struct_cc_numdens0_imag"));
}

void ReadCheckPointParticles(FhdParticleContainer& particles, species* particleInfo, const Real* dxp) {

    // timer for profiling
    BL_PROFILE_VAR("ReadCheckPointParticles()",ReadCheckPointParticles);

    std::string checkpointname;

    // checkpoint file name, e.g., chk0000010
    if (particle_restart > 0) {
        checkpointname = amrex::Concatenate(chk_base_name,particle_restart,9);
    }
    else {
        checkpointname = amrex::Concatenate(chk_base_name,restart,9);
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

    // read in statsCount
    is >> temp;
    GotoNextLine(is);

    // read in BoxArray (fluid) from Header
    BoxArray ba;
    ba.readFrom(is);
    GotoNextLine(is);

    // read in BoxArray (particle) from Header
    BoxArray bc;
    bc.readFrom(is);
    GotoNextLine(is);

    BoxArray bp;
    bp.readFrom(is);
    GotoNextLine(is);

    // create a distribution mapping
    DistributionMapping dm { bc, ParallelDescriptor::NProcs() };

    //set number of ghost cells to fit whole peskin kernel
    int ang = 1;
    if(*(std::max_element(pkernel_fluid.begin(),pkernel_fluid.begin()+nspecies)) == 3) {
        ang = 2;
    }
    else if(*(std::max_element(pkernel_fluid.begin(),pkernel_fluid.begin()+nspecies)) == 4) {
        ang = 3;
    }
    else if(*(std::max_element(pkernel_fluid.begin(),pkernel_fluid.begin()+nspecies)) == 6) {
        ang = 4;
    }
    else if (*(std::max_element(eskernel_fluid.begin(),eskernel_fluid.begin()+nspecies)) > 0) {
        ang = static_cast<int>(floor(*(std::max_element(eskernel_fluid.begin(),eskernel_fluid.begin()+nspecies)))/2+1);
    }

    Box minBox = bc.minimalBox();

//    //IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
//    //IntVect dom_hi(AMREX_D_DECL(n_cells[0]/2-1, n_cells[1]/2-1, n_cells[2]/2-1));

//    IntVect dom_lo
//    IntVect dom_hi(AMREX_D_DECL(n_cells[0]/2-1, n_cells[1]/2-1, n_cells[2]/2-1));


//    Box domain(dom_lo, dom_hi);

    RealBox realDomain({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                       {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    Vector<int> is_periodic_c(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_vel_lo[i] == -1 && bc_vel_hi[i] == -1) {
            is_periodic_c[i] = 1;
        }
    }

    Geometry geomC(minBox,&realDomain,CoordSys::cartesian,is_periodic_c.data());

//    Print() <<  "domain: " << domain << std::endl;
//    Print() <<  "geom: " << geomC << std::endl;
//    Print() <<  "Box Array: " << bc << std::endl;
//    Print() <<  "Dist Map: " << dm << std::endl;

    // restore particles

    //cout << "Restoring!\n";
    particles.Restart(checkpointname,"particle");
    //cout << "Restored!\n";
    int np = particles.TotalNumberOfParticles();
    //particlesTemp.Checkpoint("testcheck","particle");
    Print() << "Checkpoint contains " << np << " particles." <<std::endl;

    particles.ReInitParticles();

    particles.PostRestart();
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
