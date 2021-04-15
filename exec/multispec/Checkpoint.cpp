#include "multispec_test_functions.H"
#include "multispec_functions.H"

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
                     const amrex::Real dt,
                     const MultiFab& rho,
                     const MultiFab& rhotot,
                     const MultiFab& pi,
                     std::array< MultiFab, AMREX_SPACEDIM >& umac,
                     const MultiFab& Epot,
                     std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot)
{
    // timer for profiling
    BL_PROFILE_VAR("WriteCheckPoint()",WriteCheckPoint);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate(chk_base_name,step,7);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    BoxArray ba = rho.boxArray();

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
        HeaderFile << "Checkpoint file for FHDeX/multispecies\n";

        // write out the time step number
        HeaderFile << step << "\n";

        // write out time
        HeaderFile << time << "\n";

        // write out dt
        HeaderFile << dt << "\n";

        // C++ random number engine
        amrex::SaveRandomState(HeaderFile);

        // write the BoxArray
        ba.writeOn(HeaderFile);
        HeaderFile << '\n';
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

    if (use_charged_fluid) {
        VisMF::Write(grad_Epot[0],
                     amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "grad_Epotx"));
        VisMF::Write(grad_Epot[1],
                     amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "grad_Epoty"));
#if (AMREX_SPACEDIM == 3)
        VisMF::Write(grad_Epot[2],
                     amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "grad_Epotz"));
#endif
    }
    
    VisMF::Write(rho,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "rho"));
    VisMF::Write(rhotot,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "rhotot"));
    VisMF::Write(pi,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "pi"));

    if (use_charged_fluid) {
        VisMF::Write(Epot,
                     amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "Epot"));
    }

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
                    amrex::Real& dt,
                    MultiFab& rho,
                    MultiFab& rhotot,
                    MultiFab& pi,
                    std::array< MultiFab, AMREX_SPACEDIM >& umac,
                    MultiFab& Epot,
                    std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot,
                    BoxArray& ba,
                    DistributionMapping& dmap)
{
    // timer for profiling
    BL_PROFILE_VAR("ReadCheckPoint()",ReadCheckPoint);

    // checkpoint file name, e.g., chk0000010
    const std::string& checkpointname = amrex::Concatenate(chk_base_name,restart,7);

    amrex::Print() << "Restart from checkpoint " << checkpointname << "\n";

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    std::string line, word;

    int ng_s; // ghost cells for density MultiFabs
    if (advection_type == 0) {
        ng_s = 2; // centered advection
    }
    else if (advection_type <= 3) {
        ng_s = 3; // bilinear limited, biliniear unlimited, or unlimited quad bds
    }
    else if (advection_type == 4) {
        ng_s = 4; // limited quad bds
    }
    
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

        // read in dt
        is >> dt;
        GotoNextLine(is);

        // C++ random number engine
        amrex::RestoreRandomState(is, 1, 0);

        // read in level 'lev' BoxArray from Header
        ba.readFrom(is);
        GotoNextLine(is);

        // create a distribution mapping
        dmap.define(ba, ParallelDescriptor::NProcs());

        // build MultiFab data
        umac[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);
        umac[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);
#if (AMREX_SPACEDIM == 3)
        umac[2].define(convert(ba,nodal_flag_z), dmap, 1, 1);
#endif

        if (use_charged_fluid) {
            grad_Epot[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);
            grad_Epot[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);
#if (AMREX_SPACEDIM == 3)
            grad_Epot[2].define(convert(ba,nodal_flag_z), dmap, 1, 1);
#endif
        }
        
        rho   .define(ba, dmap, nspecies, ng_s);
        rhotot.define(ba, dmap,        1, ng_s);
        pi    .define(ba, dmap,        1, 1);

        if (use_charged_fluid) {
            Epot.define(ba, dmap, 1, 1);
        }
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

    if (use_charged_fluid) {
        VisMF::Read(grad_Epot[0],
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "grad_Epotx"));
        VisMF::Read(grad_Epot[1],
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "grad_Epoty"));
#if (AMREX_SPACEDIM == 3)
        VisMF::Read(grad_Epot[2],
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "grad_Epotz"));
#endif
    }
    
    VisMF::Read(rho,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "rho"));
    VisMF::Read(rhotot,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "rhotot"));
    VisMF::Read(pi,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "pi"));

    if (use_charged_fluid) {
        VisMF::Read(Epot,
                    amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "Epot"));
    }
    
    // random number engines (fortran interface)
    int digits = 7;
    rng_restart(&restart,&digits);
}


