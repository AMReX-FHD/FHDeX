#include "AMReX_PlotFileUtil.H"
#include "AMReX_PlotFileDataImpl.H"

#include <sys/stat.h>

#include "common_functions.H"

#include "compressible_functions_stag.H"

#include "common_namespace.H"

#include "chrono"

using namespace std::chrono;

using namespace common;

namespace {
    void GotoNextLine (std::istream& is)
    {
        constexpr std::streamsize bl_ignore_max { 100000 };
        is.ignore(bl_ignore_max, '\n');
    }
}

void WriteCheckPoint3D(int step,
                       const amrex::Real time,
                       int statsCount,
                       const amrex::Geometry geom,
                       const amrex::MultiFab& cu,
                       const amrex::MultiFab& cuMeans,
                       const amrex::MultiFab& cuVars,
                       const amrex::MultiFab& prim,
                       const amrex::MultiFab& primMeans,
                       const amrex::MultiFab& primVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velVars,
                       const amrex::MultiFab& coVars,
                       const Vector<Real>& spatialCross, int ncross)
{
    // timer for profiling
    BL_PROFILE_VAR("WriteCheckPoint3D()",WriteCheckPoint3D);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate("chk",step,9);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    BoxArray ba = cu.boxArray();

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
        HeaderFile << "Checkpoint file for FHDeX/compressStag\n";

        // write out the time step number
        HeaderFile << step << "\n";

        // write out time
        HeaderFile << time << "\n";

        // write out statsCount
        HeaderFile << statsCount << "\n";
        
        // write the BoxArray (fluid)
        ba.writeOn(HeaderFile);
        HeaderFile << '\n';

        // Write all the vectors associated with cross averages into the Header file
        if (plot_cross) {

            // spatialCross
            for (int i=0; i<n_cells[0]*ncross; i++) {
                HeaderFile << std::setprecision(16) << spatialCross[i] << "\n";
            }
        }
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

    // cu, cuMeans and cuVars
    VisMF::Write(cu,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cu"));
    VisMF::Write(cuMeans,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuMeans"));
    VisMF::Write(cuVars,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuVars"));

    // prim, primMeans and primVars
    VisMF::Write(prim,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "prim"));
    VisMF::Write(primMeans,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primMeans"));
    VisMF::Write(primVars,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primVars"));

    // velocity and momentum (instantaneous, means, variances)
    VisMF::Write(vel[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velx"));
    VisMF::Write(vel[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "vely"));
    VisMF::Write(vel[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velz"));
    VisMF::Write(velMeans[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velmeanx"));
    VisMF::Write(velMeans[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velmeany"));
    VisMF::Write(velMeans[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velmeanz"));
    VisMF::Write(velVars[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velvarx"));
    VisMF::Write(velVars[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velvary"));
    VisMF::Write(velVars[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velvarz"));

    VisMF::Write(cumom[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomx"));
    VisMF::Write(cumom[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomy"));
    VisMF::Write(cumom[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomz"));
    VisMF::Write(cumomMeans[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumommeanx"));
    VisMF::Write(cumomMeans[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumommeany"));
    VisMF::Write(cumomMeans[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumommeanz"));
    VisMF::Write(cumomVars[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomvarx"));
    VisMF::Write(cumomVars[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomvary"));
    VisMF::Write(cumomVars[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomvarz"));

    // coVars
    VisMF::Write(coVars,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "coVars"));
}

void WriteCheckPoint1D(int step,
                       const amrex::Real time,
                       int statsCount,
                       const amrex::Geometry geom,
                       const amrex::MultiFab& cu,
                       const amrex::MultiFab& cuMeans,
                       const amrex::MultiFab& cuVars,
                       const amrex::MultiFab& prim,
                       const amrex::MultiFab& primMeans,
                       const amrex::MultiFab& primVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velVars,
                       const amrex::MultiFab& coVars,
                       const amrex::MultiFab& spatialCross, int ncross)
{
    // timer for profiling
    BL_PROFILE_VAR("WriteCheckPoint1D()",WriteCheckPoint1D);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate("chk",step,9);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    BoxArray ba = cu.boxArray();

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
        HeaderFile << "Checkpoint file for FHDeX/compressStag\n";

        // write out the time step number
        HeaderFile << step << "\n";

        // write out time
        HeaderFile << time << "\n";

        // write out statsCount
        HeaderFile << statsCount << "\n";
        
        // write the BoxArray (fluid)
        ba.writeOn(HeaderFile);
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

    // cu, cuMeans and cuVars
    VisMF::Write(cu,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cu"));
    VisMF::Write(cuMeans,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuMeans"));
    VisMF::Write(cuVars,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuVars"));

    // prim, primMeans and primVars
    VisMF::Write(prim,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "prim"));
    VisMF::Write(primMeans,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primMeans"));
    VisMF::Write(primVars,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primVars"));

    // velocity and momentum (instantaneous, means, variances)
    VisMF::Write(vel[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velx"));
    VisMF::Write(vel[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "vely"));
    VisMF::Write(vel[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velz"));
    VisMF::Write(velMeans[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velmeanx"));
    VisMF::Write(velMeans[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velmeany"));
    VisMF::Write(velMeans[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velmeanz"));
    VisMF::Write(velVars[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velvarx"));
    VisMF::Write(velVars[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velvary"));
    VisMF::Write(velVars[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velvarz"));

    VisMF::Write(cumom[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomx"));
    VisMF::Write(cumom[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomy"));
    VisMF::Write(cumom[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomz"));
    VisMF::Write(cumomMeans[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumommeanx"));
    VisMF::Write(cumomMeans[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumommeany"));
    VisMF::Write(cumomMeans[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumommeanz"));
    VisMF::Write(cumomVars[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomvarx"));
    VisMF::Write(cumomVars[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomvary"));
    VisMF::Write(cumomVars[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomvarz"));

    // coVars
    VisMF::Write(coVars,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "coVars"));

    // spatialCross
    VisMF::Write(spatialCross,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "spatialCross"));
}

void ReadCheckPoint3D(int& step,
                     amrex::Real& time,
                     int& statsCount,
                     amrex::Geometry geom,
                     const amrex::Box& domain,
                     amrex::MultiFab& cu,
                     amrex::MultiFab& cuMeans,
                     amrex::MultiFab& cuVars,
                     amrex::MultiFab& prim,
                     amrex::MultiFab& primMeans,
                     amrex::MultiFab& primVars,
                     std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                     std::array<MultiFab, AMREX_SPACEDIM>& cumomMeans,
                     std::array<MultiFab, AMREX_SPACEDIM>& cumomVars,
                     std::array<MultiFab, AMREX_SPACEDIM>& vel,
                     std::array<MultiFab, AMREX_SPACEDIM>& velMeans,
                     std::array<MultiFab, AMREX_SPACEDIM>& velVars,
                     amrex::MultiFab& coVars,
                     Vector<Real>& spatialCross,
                     int ncross,
                     BoxArray& ba, DistributionMapping& dmap)
{
    // timer for profiling
    BL_PROFILE_VAR("ReadCheckPoint3D()",ReadCheckPoint3D);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate("chk",restart,9);

    amrex::Print() << "Restart from checkpoint " << checkpointname << "\n";

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    std::string line, word;

    // read in old boxarray, and create old distribution map (this is to read in MFabs)
    BoxArray ba_old;
    DistributionMapping dmap_old;

    // initialize new boxarray
    ba.define(domain);
    ba.maxSize(IntVect(max_grid_size));
    dmap.define(ba, ParallelDescriptor::NProcs());

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
        ba_old.readFrom(is);
        GotoNextLine(is);

        // Read all the vectors associated with cross averages from the Header file
        if (plot_cross) {
            if (reset_stats == 1) {
                spatialCross.assign(spatialCross.size(), 0.0);
            }
            else {
                Real val;
                // spatialCross
                for (int i=0; i<n_cells[0]*ncross; i++) {
                    is >> val;
                    GotoNextLine(is);
                    spatialCross[i] = val;
                }
            }
        }

        // create old distribution mapping
        dmap_old.define(ba_old, ParallelDescriptor::NProcs());

        // Define these multifabs using new ba and dmap
        // cu, cuMeans, cuVars
        cu.define(ba,dmap,nvars,ngc);
        cuMeans.define(ba,dmap,nvars,ngc);
        cuVars.define(ba,dmap,nvars,ngc);

        // prim, primMeans, primVars
        prim.define(ba,dmap,nprimvars,ngc);
        primMeans.define(ba,dmap,nprimvars+3,ngc);
        primVars.define(ba,dmap,nprimvars+5,ngc);

        // velocity and momentum (instantaneous, means, variances)
        for (int d=0; d<AMREX_SPACEDIM; d++) {
            vel[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ngc);
            cumom[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ngc);
            velMeans[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
            cumomMeans[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
            velVars[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
            cumomVars[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
        }

        // coVars
        coVars.define(ba,dmap,26,0);

    }

    // C++ random number engine
    // each MPI process reads in its own file
    int comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    int n_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
    
    // Need to add guard for restarting with more MPI ranks than the previous checkpointing run
    if (seed < 0) {

#ifdef AMREX_USE_CUDA
        Abort("Restart with negative seed not supported on GPU");
#endif
        
        // read in rng state from checkpoint
        // don't read in all the rng states at once (overload filesystem)
        // one at a time write out the rng states to different files, one for each MPI rank        
        for (int rank=0; rank<n_ranks; ++rank) {

            if (comm_rank == rank) {
    
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
    Read_Copy_MF_Checkpoint(cu,"cu",checkpointname,ba_old,dmap_old,nvars,1);
    Read_Copy_MF_Checkpoint(prim,"prim",checkpointname,ba_old,dmap_old,nprimvars,1);

    Read_Copy_MF_Checkpoint(vel[0],"velx",checkpointname,ba_old,dmap_old,1,1,0);
    Read_Copy_MF_Checkpoint(vel[1],"vely",checkpointname,ba_old,dmap_old,1,1,1);
    Read_Copy_MF_Checkpoint(vel[2],"velz",checkpointname,ba_old,dmap_old,1,1,2);

    Read_Copy_MF_Checkpoint(cumom[0],"cumomx",checkpointname,ba_old,dmap_old,1,1,0);
    Read_Copy_MF_Checkpoint(cumom[1],"cumomy",checkpointname,ba_old,dmap_old,1,1,1);
    Read_Copy_MF_Checkpoint(cumom[2],"cumomz",checkpointname,ba_old,dmap_old,1,1,2);

    // Set all stats to zero if reset stats, else read
    if (reset_stats == 1) {
        cuMeans.setVal(0.0);
        cuVars.setVal(0.0);
        primMeans.setVal(0.0);
        primVars.setVal(0.0);
        for (int d=0; d<AMREX_SPACEDIM; d++) {
            velMeans[d].setVal(0.);
            velVars[d].setVal(0.);
            cumomMeans[d].setVal(0.);
            cumomVars[d].setVal(0.);
        }
        coVars.setVal(0.0);
    }
    else {
        Read_Copy_MF_Checkpoint(cuMeans,"cuMeans",checkpointname,ba_old,dmap_old,nvars,1);
        Read_Copy_MF_Checkpoint(cuVars,"cuVars",checkpointname,ba_old,dmap_old,nvars,1);

        Read_Copy_MF_Checkpoint(primMeans,"primMeans",checkpointname,ba_old,dmap_old,nprimvars+3,1);
        Read_Copy_MF_Checkpoint(primVars,"primVars",checkpointname,ba_old,dmap_old,nprimvars+5,1);

        Read_Copy_MF_Checkpoint(coVars,"coVars",checkpointname,ba_old,dmap_old,26,0);

        Read_Copy_MF_Checkpoint(velMeans[0],"velmeanx",checkpointname,ba_old,dmap_old,1,0,0);
        Read_Copy_MF_Checkpoint(velMeans[1],"velmeany",checkpointname,ba_old,dmap_old,1,0,1);
        Read_Copy_MF_Checkpoint(velMeans[2],"velmeanz",checkpointname,ba_old,dmap_old,1,0,2);
        Read_Copy_MF_Checkpoint(velVars[0],"velvarx",checkpointname,ba_old,dmap_old,1,0,0);
        Read_Copy_MF_Checkpoint(velVars[1],"velvary",checkpointname,ba_old,dmap_old,1,0,1);
        Read_Copy_MF_Checkpoint(velVars[2],"velvarz",checkpointname,ba_old,dmap_old,1,0,2);

        Read_Copy_MF_Checkpoint(cumomMeans[0],"cumommeanx",checkpointname,ba_old,dmap_old,1,0,0);
        Read_Copy_MF_Checkpoint(cumomMeans[1],"cumommeany",checkpointname,ba_old,dmap_old,1,0,1);
        Read_Copy_MF_Checkpoint(cumomMeans[2],"cumommeanz",checkpointname,ba_old,dmap_old,1,0,2);
        Read_Copy_MF_Checkpoint(cumomVars[0],"cumomvarx",checkpointname,ba_old,dmap_old,1,0,0);
        Read_Copy_MF_Checkpoint(cumomVars[1],"cumomvary",checkpointname,ba_old,dmap_old,1,0,1);
        Read_Copy_MF_Checkpoint(cumomVars[2],"cumomvarz",checkpointname,ba_old,dmap_old,1,0,2);
    }

    // FillBoundaries
    cu.FillBoundary(geom.periodicity());
    cuMeans.FillBoundary(geom.periodicity());
    cuVars.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());
    primMeans.FillBoundary(geom.periodicity());
    primVars.FillBoundary(geom.periodicity());
    vel[0].FillBoundary(geom.periodicity());
    vel[1].FillBoundary(geom.periodicity());
    vel[2].FillBoundary(geom.periodicity());
    cumom[0].FillBoundary(geom.periodicity());
    cumom[1].FillBoundary(geom.periodicity());
    cumom[2].FillBoundary(geom.periodicity());
}

void ReadCheckPoint1D(int& step,
                     amrex::Real& time,
                     int& statsCount,
                     amrex::Geometry geom,
                     const amrex::Box& domain,
                     amrex::MultiFab& cu,
                     amrex::MultiFab& cuMeans,
                     amrex::MultiFab& cuVars,
                     amrex::MultiFab& prim,
                     amrex::MultiFab& primMeans,
                     amrex::MultiFab& primVars,
                     std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                     std::array<MultiFab, AMREX_SPACEDIM>& cumomMeans,
                     std::array<MultiFab, AMREX_SPACEDIM>& cumomVars,
                     std::array<MultiFab, AMREX_SPACEDIM>& vel,
                     std::array<MultiFab, AMREX_SPACEDIM>& velMeans,
                     std::array<MultiFab, AMREX_SPACEDIM>& velVars,
                     amrex::MultiFab& coVars,
                     amrex::MultiFab& spatialCross,
                     int ncross,
                     BoxArray& ba, DistributionMapping& dmap)
{
    // timer for profiling
    BL_PROFILE_VAR("ReadCheckPoint1D()",ReadCheckPoint1D);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate("chk",restart,9);

    amrex::Print() << "Restart from checkpoint " << checkpointname << "\n";

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    std::string line, word;

    // read in old boxarray, and create old distribution map (this is to read in MFabs)
    BoxArray ba_old;
    DistributionMapping dmap_old;

    // initialize new boxarray
    ba.define(domain);
    ba.maxSize(IntVect(max_grid_size));
    dmap.define(ba, ParallelDescriptor::NProcs());

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
        ba_old.readFrom(is);
        GotoNextLine(is);

        // create old distribution mapping
        dmap_old.define(ba_old, ParallelDescriptor::NProcs());

        // Define these multifabs using new ba and dmap
        // cu, cuMeans, cuVars
        cu.define(ba,dmap,nvars,ngc);
        cuMeans.define(ba,dmap,nvars,ngc);
        cuVars.define(ba,dmap,nvars,ngc);

        // prim, primMeans, primVars
        prim.define(ba,dmap,nprimvars,ngc);
        primMeans.define(ba,dmap,nprimvars+3,ngc);
        primVars.define(ba,dmap,nprimvars+5,ngc);

        // velocity and momentum (instantaneous, means, variances)
        for (int d=0; d<AMREX_SPACEDIM; d++) {
            vel[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ngc);
            cumom[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, ngc);
            velMeans[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
            cumomMeans[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
            velVars[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
            cumomVars[d].define(convert(ba,nodal_flag_dir[d]), dmap, 1, 0);
        }

        // coVars
        coVars.define(ba,dmap,26,0);

        // spatialCross
        spatialCross.define(ba,dmap,ncross,0);

    }

    // C++ random number engine
    // each MPI process reads in its own file
    int comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    int n_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

    // don't read in all the rng states at once (overload filesystem)
    // one at a time write out the rng states to different files, one for each MPI rank
    // Need to add guard for restarting with more MPI ranks than the previous checkpointing run
    if (seed == -1) {
        for (int rank=0; rank<n_ranks; ++rank) {

            if (comm_rank == rank) {
    
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
    Read_Copy_MF_Checkpoint(cu,"cu",checkpointname,ba_old,dmap_old,nvars,1);
    Read_Copy_MF_Checkpoint(prim,"prim",checkpointname,ba_old,dmap_old,nprimvars,1);

    Read_Copy_MF_Checkpoint(vel[0],"velx",checkpointname,ba_old,dmap_old,1,1,0);
    Read_Copy_MF_Checkpoint(vel[1],"vely",checkpointname,ba_old,dmap_old,1,1,1);
    Read_Copy_MF_Checkpoint(vel[2],"velz",checkpointname,ba_old,dmap_old,1,1,2);

    Read_Copy_MF_Checkpoint(cumom[0],"cumomx",checkpointname,ba_old,dmap_old,1,1,0);
    Read_Copy_MF_Checkpoint(cumom[1],"cumomy",checkpointname,ba_old,dmap_old,1,1,1);
    Read_Copy_MF_Checkpoint(cumom[2],"cumomz",checkpointname,ba_old,dmap_old,1,1,2);

    // Set all stats to zero if reset stats, else read
    if (reset_stats == 1) {
        cuMeans.setVal(0.0);
        cuVars.setVal(0.0);
        primMeans.setVal(0.0);
        primVars.setVal(0.0);
        for (int d=0; d<AMREX_SPACEDIM; d++) {
            velMeans[d].setVal(0.);
            velVars[d].setVal(0.);
            cumomMeans[d].setVal(0.);
            cumomVars[d].setVal(0.);
        }
        coVars.setVal(0.0);
        spatialCross.setVal(0.0);
    }
    else {
        Read_Copy_MF_Checkpoint(cuMeans,"cuMeans",checkpointname,ba_old,dmap_old,nvars,1);
        Read_Copy_MF_Checkpoint(cuVars,"cuVars",checkpointname,ba_old,dmap_old,nvars,1);

        Read_Copy_MF_Checkpoint(primMeans,"primMeans",checkpointname,ba_old,dmap_old,nprimvars+3,1);
        Read_Copy_MF_Checkpoint(primVars,"primVars",checkpointname,ba_old,dmap_old,nprimvars+5,1);

        Read_Copy_MF_Checkpoint(coVars,"coVars",checkpointname,ba_old,dmap_old,26,0);

        Read_Copy_MF_Checkpoint(spatialCross,"spatialCross",checkpointname,ba_old,dmap_old,ncross,0);

        Read_Copy_MF_Checkpoint(velMeans[0],"velmeanx",checkpointname,ba_old,dmap_old,1,0,0);
        Read_Copy_MF_Checkpoint(velMeans[1],"velmeany",checkpointname,ba_old,dmap_old,1,0,1);
        Read_Copy_MF_Checkpoint(velMeans[2],"velmeanz",checkpointname,ba_old,dmap_old,1,0,2);
        Read_Copy_MF_Checkpoint(velVars[0],"velvarx",checkpointname,ba_old,dmap_old,1,0,0);
        Read_Copy_MF_Checkpoint(velVars[1],"velvary",checkpointname,ba_old,dmap_old,1,0,1);
        Read_Copy_MF_Checkpoint(velVars[2],"velvarz",checkpointname,ba_old,dmap_old,1,0,2);

        Read_Copy_MF_Checkpoint(cumomMeans[0],"cumommeanx",checkpointname,ba_old,dmap_old,1,0,0);
        Read_Copy_MF_Checkpoint(cumomMeans[1],"cumommeany",checkpointname,ba_old,dmap_old,1,0,1);
        Read_Copy_MF_Checkpoint(cumomMeans[2],"cumommeanz",checkpointname,ba_old,dmap_old,1,0,2);
        Read_Copy_MF_Checkpoint(cumomVars[0],"cumomvarx",checkpointname,ba_old,dmap_old,1,0,0);
        Read_Copy_MF_Checkpoint(cumomVars[1],"cumomvary",checkpointname,ba_old,dmap_old,1,0,1);
        Read_Copy_MF_Checkpoint(cumomVars[2],"cumomvarz",checkpointname,ba_old,dmap_old,1,0,2);
    }

    // FillBoundaries
    cu.FillBoundary(geom.periodicity());
    cuMeans.FillBoundary(geom.periodicity());
    cuVars.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());
    primMeans.FillBoundary(geom.periodicity());
    primVars.FillBoundary(geom.periodicity());
    vel[0].FillBoundary(geom.periodicity());
    vel[1].FillBoundary(geom.periodicity());
    vel[2].FillBoundary(geom.periodicity());
    cumom[0].FillBoundary(geom.periodicity());
    cumom[1].FillBoundary(geom.periodicity());
    cumom[2].FillBoundary(geom.periodicity());
}


void ReadFile (const std::string& filename, Vector<char>& charBuf,
          bool bExitOnError)
{
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

void Read_Copy_MF_Checkpoint(amrex::MultiFab& mf, std::string mf_name, const std::string& checkpointname, 
                             BoxArray& ba_old, DistributionMapping& dmap_old,
                             int NVARS, int ghost, int nodal_flag)
{
    // define temporary MF
    MultiFab mf_temp;
    if (nodal_flag < 0) {
        if (ghost) {
            mf_temp.define(ba_old,dmap_old,NVARS,ngc);
        }
        else {
            mf_temp.define(ba_old,dmap_old,NVARS,0);
        }

    }
    else {
        if (ghost) {
            mf_temp.define(convert(ba_old,nodal_flag_dir[nodal_flag]),dmap_old,NVARS,ngc);
        }
        else {
            mf_temp.define(convert(ba_old,nodal_flag_dir[nodal_flag]),dmap_old,NVARS,0);
        }

    }
    
    // Read into temporary MF from file
    VisMF::Read(mf_temp,amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", mf_name));

    // Copy temporary MF into the new MF
    if (ghost) {
        mf.ParallelCopy(mf_temp, 0, 0, NVARS, ngc, ngc);
    }
    else {
        mf.ParallelCopy(mf_temp, 0, 0, NVARS, 0, 0);
    }
}
