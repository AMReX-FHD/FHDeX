#include "hydro_test_functions.H"

#include "rng_functions_F.H"

#include "AMReX_PlotFileUtil.H"
#include "AMReX_PlotFileDataImpl.H"

#include <sys/stat.h> 
using namespace amrex;

namespace {
    void GotoNextLine (std::istream& is)
    {
        constexpr std::streamsize bl_ignore_max { 100000 };
        is.ignore(bl_ignore_max, '\n');
    }
}

void WriteCheckPoint(int step,
                     const amrex::Real time,
                     std::array< MultiFab, AMREX_SPACEDIM >& umac,
                     const MultiFab& tracer, TurbForcing& turbforce)
{
    // timer for profiling
    BL_PROFILE_VAR("WriteCheckPoint()",WriteCheckPoint);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate(chk_base_name,step,7);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    BoxArray ba = tracer.boxArray();

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

        // write the BoxArray
        ba.writeOn(HeaderFile);
        HeaderFile << '\n';

	if (turbForcing == 1) {
	  // write turbulent forcing U's
	  for (int i=0; i<132; ++i) {
            HeaderFile << turbforce.getU(i) << '\n';
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
    VisMF::Write(umac[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "umac"));
    VisMF::Write(umac[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "vmac"));
#if (AMREX_SPACEDIM == 3)
    VisMF::Write(umac[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "wmac"));
#endif
    VisMF::Write(tracer,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "tracer"));

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
                    std::array< MultiFab, AMREX_SPACEDIM >& umac,
                    MultiFab& tracer, TurbForcing& turbforce,
                    BoxArray& ba, DistributionMapping& dmap)
{
    // timer for profiling
    BL_PROFILE_VAR("ReadCheckPoint()",ReadCheckPoint);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate(chk_base_name,restart,7);

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
        GotoNextLine(is);
        ++step;

        // read in time
        is >> time;
        GotoNextLine(is);

        // read in level 'lev' BoxArray from Header
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        dmap.define(ba, ParallelDescriptor::NProcs());

	if (turbForcing == 1) {
	  turbforce.define(ba,dmap,turb_a,turb_b);
	}

	if (turbForcing == 1) {
	  // read in turbulent forcing U's
	  Real utemp;
	  for (int i=0; i<132; ++i) {
            is >> utemp;
            turbforce.setU(i,utemp);
	  }        
	}

        // build MultiFab data
        umac[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);
        umac[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);
#if (AMREX_SPACEDIM == 3)
        umac[2].define(convert(ba,nodal_flag_z), dmap, 1, 1);
#endif
        tracer.define(ba, dmap, 1, 1);
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
    VisMF::Read(umac[0],
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "umac"));
    VisMF::Read(umac[1],
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "vmac"));
#if (AMREX_SPACEDIM == 3)
    VisMF::Read(umac[2],
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "wmac"));
#endif
    VisMF::Read(tracer,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "tracer"));
    
    // random number engines (fortran interface)
    int digits = 7;
    rng_restart(&restart,&digits);
}

void
ReadFile (const std::string& filename, Vector<char>& charBuf,
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
