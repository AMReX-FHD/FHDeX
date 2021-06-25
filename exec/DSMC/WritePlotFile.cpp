#include "INS_functions.H"
#include "common_functions.H"

#include "AMReX_PlotFileUtil.H"

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry cgeom,
                   const MultiFab& particleInstant,
                   const MultiFab& particleMeans,
                   const MultiFab& particleVars,
                   FhdParticleContainer& particles) 
{
    BL_PROFILE_VAR("WritePlotFile()",WritePlotFile);

    std::string cplotfilename = Concatenate("cplt",step,9);

    BoxArray cba = particleInstant.boxArray();
    DistributionMapping cdmap = particleInstant.DistributionMap();

    int cnPlot = 16+2*nspecies;

    MultiFab cplotfile(cba, cdmap, cnPlot, 0);

    Vector<std::string> cvarNames(cnPlot);

    amrex::MultiFab::Copy(cplotfile,particleInstant,0,0 ,8+nspecies,0);
    amrex::MultiFab::Copy(cplotfile,particleMeans  ,0,8+nspecies,8+nspecies,0);

    cvarNames[0] = "membersInstant";
    cvarNames[1] = "densityInstant";
    cvarNames[2] = "velxInstant";
    cvarNames[3] = "velyInstant";
    cvarNames[4] = "velzInstant";
    cvarNames[5] = "pxInstant";
    cvarNames[6] = "pyInstant";
    cvarNames[7] = "pzInstant";

    int ccount = 8;

    for(int i=0;i<nspecies;i++)
    {
        std::string specname = Concatenate("densityInstantSpecies",i);
        cvarNames[ccount+i]= specname;
    }

    cvarNames[8+nspecies] = "membersMean";
    cvarNames[9+nspecies] = "densityMean";
    cvarNames[10+nspecies] = "velxMean";
    cvarNames[11+nspecies] = "velyMean";
    cvarNames[12+nspecies] = "velzMean";
    cvarNames[13+nspecies] = "pxMean";
    cvarNames[14+nspecies] = "pyMean";
    cvarNames[15+nspecies] = "pzMean";

    ccount = 16+nspecies;

    for(int i=0;i<nspecies;i++)
    {
        std::string specname = Concatenate("densityMeanSpecies",i);
        cvarNames[ccount+i]= specname;
    }

    // timer
    Real t1 = ParallelDescriptor::second();
    
    WriteSingleLevelPlotfile(cplotfilename,cplotfile,cvarNames,cgeom,time,step);
    
    Real t2 = ParallelDescriptor::second() - t1;
    ParallelDescriptor::ReduceRealMax(t2);
    amrex::Print() << "Time spent writing plotfile cplt " << t2 << std::endl;
    
    // particle ASCII
    if(plot_ascii == 1)
    {

        std::string asciiName = Concatenate("ascii_means",step,9);
        outputMFAscii(particleMeans, asciiName);

    }

    // particle in cplt file
    Vector<std::string> real_comp_names;
    Vector<std::string>  int_comp_names;



    Vector<int> write_real_comp = {
        1, // radius
        1, // velx
        1, // vely
        1, // velz
        0, // boostx
        0, // boosty
        0, // boostz
        0, // R
        0 // timeFrac
    };

    Vector<int> write_int_comp = {
        0, // sorted
        0, // i
        0, // j
        0, // k
        1  // species
    };

    t1 = ParallelDescriptor::second();
    
    particles.WritePlotFile(cplotfilename, "particles",
                            write_real_comp, write_int_comp, FHD_realData::names(), FHD_intData::names());
    
    t2 = ParallelDescriptor::second() - t1;
    ParallelDescriptor::ReduceRealMax(t2);
    amrex::Print() << "Time spent writing particle plotfile " << t2 << std::endl;

}
