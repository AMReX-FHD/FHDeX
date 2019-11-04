#include "INS_functions.H"

#include "AMReX_PlotFileUtil.H"
#include "AMReX_PlotFileDataImpl.H"

#include <sys/stat.h>

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
                     const std::array< MultiFab, AMREX_SPACEDIM >& umacV,
                     const MultiFab& pres,
                     const FhdParticleContainer& particles,
                     const MultiFab& particleMeans,
                     const MultiFab& particleVars,
                     const MultiFab& potential)
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

    // umacV (variance)
    VisMF::Write(umacV[0],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "umacV"));
    VisMF::Write(umacV[1],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "vmacV"));
#if (AMREX_SPACEDIM == 3)
    VisMF::Write(umacV[2],
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "wmacV"));
#endif

    // pressure
    VisMF::Write(pres,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "pressure"));

    // particle mean and variance
    VisMF::Write(particleMeans,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "particleMeans"));
    VisMF::Write(particleVars,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "particleVars"));

    // electrostatic potential
    VisMF::Write(potential,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "potential"));

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
    
    // random number engines
    int digits = 9;
    rng_checkpoint(&step,&digits);
    
    // checkpoint particles
    particles.Checkpoint(checkpointname,"particle");
}

void ReadCheckPoint(int& step,
                    amrex::Real& time,
                    int& statsCount,
                    std::array< MultiFab, AMREX_SPACEDIM >& umac,
                    std::array< MultiFab, AMREX_SPACEDIM >& umacM,
                    std::array< MultiFab, AMREX_SPACEDIM >& umacV,
                    MultiFab& pres,
                    MultiFab& particleMeans,
                    MultiFab& particleVars,
                    MultiFab& potential)
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
        if(pkernel_fluid == 3) {
            ang = 2;
        }
        else if(pkernel_fluid == 4) {
            ang = 3;
        }
        else if(pkernel_fluid == 6) {
            ang = 4;
        }
    
        // build MultiFab data

        // umac
        umac[0].define(convert(ba,nodal_flag_x), dm, 1, ang);
        umac[1].define(convert(ba,nodal_flag_y), dm, 1, ang);
#if (AMREX_SPACEDIM == 3)
        umac[2].define(convert(ba,nodal_flag_z), dm, 1, ang);
#endif

        // umacM
        umacM[0].define(convert(ba,nodal_flag_x), dm, 1, ang);
        umacM[1].define(convert(ba,nodal_flag_y), dm, 1, ang);
#if (AMREX_SPACEDIM == 3)
        umacM[2].define(convert(ba,nodal_flag_z), dm, 1, ang);
#endif

        // umacV
        umacV[0].define(convert(ba,nodal_flag_x), dm, 1, ang);
        umacV[1].define(convert(ba,nodal_flag_y), dm, 1, ang);
#if (AMREX_SPACEDIM == 3)
        umacV[2].define(convert(ba,nodal_flag_z), dm, 1, ang);
#endif

        // pressure
        pres.define(ba,dm,1,1);

        // particle means and variances
        particleMeans.define(bc,dm,14,0);
        particleVars .define(bc,dm,18,0);
        
        // cell centred es potential
        int ngp = 1;
        if (pkernel_es == 3) {
            ngp = 2;
        }
        else if (pkernel_es == 4) {
            ngp = 3;
        }
        else if (pkernel_es == 6) {
            ngp = 4;
        }

        potential.define(bp,dm,1,ngp);
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

    // umacV
    VisMF::Read(umacV[0],
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "umacV"));
    VisMF::Read(umacV[1],
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "vmacV"));
#if (AMREX_SPACEDIM == 3)
    VisMF::Read(umacV[2],
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "wmacV"));
#endif

    // pressure
    VisMF::Read(pres,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "pressure"));
        
    // particle means and variances
    VisMF::Read(particleMeans,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "particleMeans"));
    VisMF::Read(particleVars,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "particleVars"));

    // electrostatic potential
    VisMF::Read(potential,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "potential"));
    
    // random number engines
    int digits = 9;
    rng_restart(&restart,&digits);
}

void ReadCheckPointParticles(FhdParticleContainer& particles) {
    
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
    
    // restore particles
    particles.Restart(checkpointname,"particle");

    particles.PostRestart();

}
