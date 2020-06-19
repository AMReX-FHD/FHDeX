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
                     const MultiFab& chargeM,
                     const MultiFab& chargeV,
                     const MultiFab& potential,
                     const MultiFab& potentialM,
                     const MultiFab& potentialV)
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

    // charge
    VisMF::Write(chargeM,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "chargeM"));
    VisMF::Write(chargeV,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "chargeV"));
    
    // electrostatic potential
    VisMF::Write(potential,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "potential"));
    VisMF::Write(potentialM,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "potentialM"));
    VisMF::Write(potentialV,
                 amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "potentialV"));

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
                    MultiFab& chargeM,
                    MultiFab& chargeV,
                    MultiFab& potential,
                    MultiFab& potentialM,
                    MultiFab& potentialV)
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
        
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            umac [d].define(convert(ba,nodal_flag_dir[d]), dm, 1, ang);
            umacM[d].define(convert(ba,nodal_flag_dir[d]), dm, 1, 1);
            umacV[d].define(convert(ba,nodal_flag_dir[d]), dm, 1, 1);
        }

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

        chargeM.define(bp,dm,1,1);
        chargeV.define(bp,dm,1,1);
        potential.define(bp,dm,1,ngp);
        potentialM.define(bp,dm,1,1);
        potentialV.define(bp,dm,1,1);
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

    // charge
    VisMF::Read(chargeM,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "chargeM"));
    VisMF::Read(chargeV,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "chargeV"));
    
    // electrostatic potential
    VisMF::Read(potential,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "potential"));
    VisMF::Read(potentialM,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "potentialM"));
    VisMF::Read(potentialV,
                amrex::MultiFabFileFullPrefix(0, checkpointname, "Level_", "potentialV"));
    
    // random number engines
    int digits = 9;
    rng_restart(&restart,&digits);
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
        if(pkernel_fluid == 3) {
            ang = 2;
        }
        else if(pkernel_fluid == 4) {
            ang = 3;
        }
        else if(pkernel_fluid == 6) {
            ang = 4;
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

    FhdParticleContainer particlesTemp(geomC, dm, bc, 1);
    
    // restore particles

    //cout << "Restoring!\n";
    particlesTemp.Restart(checkpointname,"particle");
    //cout << "Restored!\n";
    int np = particlesTemp.TotalNumberOfParticles();
    //particlesTemp.Checkpoint("testcheck","particle");
    Print() << "Checkpoint contains " << np << " particles." <<std::endl;

    Real posx[np];
    Real posy[np];
    Real posz[np];
    Real charge[np];

    Real sigma[np];
    Real epsilon[np];
    int speciesV[np];

    Real diffwet[np];
    Real diffdry[np];
    Real difftot[np];

    //std::cout << "Proc " << ParallelDescriptor::MyProc() << " pull down 1 started." <<std::endl;

    particlesTemp.PullDown(0, posx, -1, np);
    //std::cout << "Proc " << ParallelDescriptor::MyProc() << " pull down 1 finished." <<std::endl;

    //std::cout << "Proc " << ParallelDescriptor::MyProc() << " particle 5 xPos: " << posx[4] << std:: endl;

    //ParallelDescriptor::Barrier();

    //std::cout << "Proc " << ParallelDescriptor::MyProc() << " through barrier." << std:: endl;

    //std::cout << "Proc " << ParallelDescriptor::MyProc() << " pull down 2 started." <<std::endl;
    particlesTemp.PullDown(0, posy, -2, np);
    //std::cout << "Proc " << ParallelDescriptor::MyProc() << " pull down 2 finished." <<std::endl;

    particlesTemp.PullDown(0, posz, -3, np);

    particlesTemp.PullDown(0, charge, 27, np);

    particlesTemp.PullDown(0, sigma, 43, np);
    particlesTemp.PullDown(0, epsilon, 44, np);

    particlesTemp.PullDown(0, diffdry, 40, np);
    particlesTemp.PullDown(0, diffwet, 41, np);
    particlesTemp.PullDown(0, difftot, 42, np);

    particlesTemp.PullDownInt(0, speciesV, 4, np);



    particles.ReInitParticles(particleInfo, dxp, posx, posy, posz, charge, sigma, epsilon, speciesV, diffdry, diffwet, difftot);


    //particles.PostRestart();
}















