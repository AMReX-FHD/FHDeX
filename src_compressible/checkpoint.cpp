#include "AMReX_PlotFileUtil.H"
#include "AMReX_PlotFileDataImpl.H"

#include <sys/stat.h>

#include "common_functions.H"

#include "compressible_functions.H"

#include "common_namespace.H"

using namespace common;

namespace {
    void GotoNextLine (std::istream& is)
    {
        constexpr std::streamsize bl_ignore_max { 100000 };
        is.ignore(bl_ignore_max, '\n');
    }
}

void WriteCheckPoint(int step,
                     const amrex::Real time,
                     int statsCount,
                     const amrex::Geometry& /*geom*/,
                     const amrex::MultiFab& cu,
                     const amrex::MultiFab& cuMeans,
                     const amrex::MultiFab& cuVars,
                     const amrex::MultiFab& prim,
                     const amrex::MultiFab& primMeans,
                     const amrex::MultiFab& primVars,
                     const amrex::MultiFab& spatialCross,
                     const amrex::MultiFab& miscStats,
                     const amrex::MultiFab& eta,
                     const amrex::MultiFab& kappa)
{
    // timer for profiling
    BL_PROFILE_VAR("WriteCheckPoint()",WriteCheckPoint);

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
        HeaderFile << "Checkpoint file for FHDeX/compress\n";

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

    // spatialCross
    VisMF::Write(spatialCross,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "spatialCross"));

    // miscStats
    VisMF::Write(miscStats,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "miscStats"));

    // eta
    VisMF::Write(eta,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "eta"));

    // kappa
    VisMF::Write(kappa,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "kappa"));
}

void ReadCheckPoint(int& step,
                     amrex::Real& time,
                     int& statsCount,
                     amrex::Geometry& /*geom*/,
                     amrex::MultiFab& cu,
                     amrex::MultiFab& cuMeans,
                     amrex::MultiFab& cuVars,
                     amrex::MultiFab& prim,
                     amrex::MultiFab& primMeans,
                     amrex::MultiFab& primVars,
                     amrex::MultiFab& spatialCross,
                     amrex::MultiFab& miscStats,
                     amrex::MultiFab& eta,
                     amrex::MultiFab& kappa)
{
    // timer for profiling
    BL_PROFILE_VAR("ReadCheckPoint()",ReadCheckPoint);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate("chk",restart,9);

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

        // create a distribution mapping
        DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

        // cu, cuMeans, cuVars
        cu.define(ba,dm,nvars,ngc);
        cuMeans.define(ba,dm,nvars,ngc);
        cuVars.define(ba,dm,nvars,ngc);

        // prim, primMeans, primVars
        prim.define(ba,dm,nprimvars,ngc);
        primMeans.define(ba,dm,nprimvars,ngc);
        primVars.define(ba,dm,nprimvars + 5,ngc);

        // spatialCross
        spatialCross.define(ba,dm,8,ngc);

        // miscStats
        miscStats.define(ba,dm,10,ngc);

        //eta and kappa
        eta.define(ba,dm,1,ngc);
        kappa.define(ba,dm,1,ngc);
    }

    // read in the MultiFab data
    // cu
    VisMF::Read(cu,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cu"));

    // prim
    VisMF::Read(prim,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "prim"));

    // eta and kappa
    VisMF::Read(eta,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "eta"));
    VisMF::Read(kappa,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "kappa"));

    if (reset_stats == 0) {
        // cuMeans, cuVars
        VisMF::Read(cuMeans,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuMeans"));
        VisMF::Read(cuVars,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuVars"));

        // primMeans, primVars
        VisMF::Read(primMeans,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primMeans"));
        VisMF::Read(primVars,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primVars"));

        // spatialCross
        VisMF::Read(spatialCross,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "spatialCross"));

        // miscStats
        VisMF::Read(miscStats,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "miscStats"));
    }
    else if (reset_stats == 1) {
        cuMeans.setVal(0.0);
        cuVars.setVal(0.0);
        primMeans.setVal(0.0);
        primVars.setVal(0.0);
        spatialCross.setVal(0.0);
        miscStats.setVal(0.0);
        statsCount = 1;
    }
}


