#include "AMReX_PlotFileUtil.H"
#include "AMReX_PlotFileDataImpl.H"

#include <sys/stat.h>

#include "common_functions.H"

#include "rng_functions_F.H"

#include "compressible_functions_stag.H"

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
                     const Vector<Real>& spatialCross)
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

    int ncross = 16+nspecies;

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

    int check;
    char str[80];
    
    strcpy (str,checkpointname.c_str());
    strcat (str,"/rng_eng_fhd");
    check = mkdir(str,0777);
    
    strcpy (str,checkpointname.c_str());
    strcat (str,"/rng_eng_particle");
    check = mkdir(str,0777);
    
    strcpy (str,checkpointname.c_str());
    strcat (str,"/rng_eng_select");
    check = mkdir(str,0777);
    
    strcpy (str,checkpointname.c_str());
    strcat (str,"/rng_eng_scatter_theta");
    check = mkdir(str,0777);
    
    strcpy (str,checkpointname.c_str());
    strcat (str,"/rng_eng_scatter_phi");
    check = mkdir(str,0777);
    
    strcpy (str,checkpointname.c_str());
    strcat (str,"/rng_eng_general");
    check = mkdir(str,0777);
    
    // random number engines (fortran interface)
    int n_digits = 7;
    rng_checkpoint(& step, & n_digits);
    
}

void ReadCheckPoint(int& step,
                     amrex::Real& time,
                     int& statsCount,
                     amrex::Geometry geom,
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
                     BoxArray& ba, DistributionMapping& dmap)
{
    // timer for profiling
    BL_PROFILE_VAR("ReadCheckPoint()",ReadCheckPoint);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate("chk",restart,9);

    amrex::Print() << "Restart from checkpoint " << checkpointname << "\n";

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    std::string line, word;

    int ncross = 16+nspecies;

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
        ba.readFrom(is);
        GotoNextLine(is);

        // Read all the vectors associated with cross averages from the Header file
        if (plot_cross) {

            Real val;
            // spatialCross
            for (int i=0; i<n_cells[0]*ncross; i++) {
                is >> val;
                GotoNextLine(is);
                spatialCross[i] = val;
            }
        }

        // create a distribution mapping
        dmap.define(ba, ParallelDescriptor::NProcs());
        
        // cu, cuMeans, cuVars
        cu.define(ba,dmap,nvars,ngc);
        cuMeans.define(ba,dmap,nvars,ngc);
        cuVars.define(ba,dmap,nvars,ngc);

        // prim, primMeans, primVars
        prim.define(ba,dmap,nprimvars,ngc);
        primMeans.define(ba,dmap,nprimvars,ngc);
        primVars.define(ba,dmap,nprimvars + 5,ngc);

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
        coVars.define(ba,dmap,21,0);

    }

    // C++ random number engine
    // each MPI process reads in its own file
    int comm_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

    int n_ranks;
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

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

    // read in the MultiFab data
    // cu, cuMeans, cuVars
    VisMF::Read(cu,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cu"));
    VisMF::Read(cuMeans,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuMeans"));
    VisMF::Read(cuVars,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cuVars"));

    // prim, primMeans, primVars
    VisMF::Read(prim,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "prim"));
    VisMF::Read(primMeans,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primMeans"));
    VisMF::Read(primVars,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "primVars"));

    // velocity and momentum (instantaneous, means, variances)
    VisMF::Read(vel[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velx"));
    VisMF::Read(vel[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "vely"));
    VisMF::Read(vel[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velz"));
    VisMF::Read(velMeans[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velmeanx"));
    VisMF::Read(velMeans[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velmeany"));
    VisMF::Read(velMeans[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velmeanz"));
    VisMF::Read(velVars[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velvarx"));
    VisMF::Read(velVars[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velvary"));
    VisMF::Read(velVars[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "velvarz"));

    VisMF::Read(cumom[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomx"));
    VisMF::Read(cumom[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomy"));
    VisMF::Read(cumom[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomz"));
    VisMF::Read(cumomMeans[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumommeanx"));
    VisMF::Read(cumomMeans[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumommeany"));
    VisMF::Read(cumomMeans[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumommeanz"));
    VisMF::Read(cumomVars[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomvarx"));
    VisMF::Read(cumomVars[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomvary"));
    VisMF::Read(cumomVars[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "cumomvarz"));

    // coVars
    VisMF::Read(coVars,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "coVars"));

    // random number engines
    int digits = 7;
    rng_restart(&restart,&digits);
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
