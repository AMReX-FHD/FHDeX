#include "INS_functions.H"
#include "common_functions.H"

#include "AMReX_PlotFileUtil.H"

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const amrex::Geometry cgeom,
                   const amrex::Geometry egeom,
                   const MultiFab& particleInstant,
                   const MultiFab& particleMeans,
                   const MultiFab& particleVars,
                   FhdParticleContainer& particles,
                   const MultiFab& charge,
                   const MultiFab& potential,
                   const std::array< MultiFab, AMREX_SPACEDIM >& efield,
                   const MultiFab& mobility) 
{

    std::string cplotfilename = Concatenate("cplt",step,9);
    std::string eplotfilename = Concatenate("eplt",step,9);
    std::string pplotfilename = Concatenate("parplt",step,9);
    std::string mplotfilename = Concatenate("mplt",step,9);

    BoxArray cba = particleInstant.boxArray();
    DistributionMapping cdmap = particleInstant.DistributionMap();

    BoxArray eba = charge.boxArray();
    DistributionMapping edmap = charge.DistributionMap();

 
//    int cnPlot = 40;
    int cnPlot = 46;

    int enPlot = 2+AMREX_SPACEDIM;

    int fnPlot = nspecies*AMREX_SPACEDIM;

    MultiFab cplotfile(cba, cdmap, cnPlot, 0);

    MultiFab eplotfile(eba, edmap, enPlot, 0);

    Vector<std::string> cvarNames(cnPlot);
    Vector<std::string> evarNames(enPlot);

    Vector<std::string> fvarNames(fnPlot);

    amrex::MultiFab::Copy(eplotfile,charge,0,0,1,0);
    amrex::MultiFab::Copy(eplotfile,potential,0,1,1,0);

    amrex::MultiFab::Copy(eplotfile,efield[0],0,2,1,0);
    amrex::MultiFab::Copy(eplotfile,efield[1],0,3,1,0);
    amrex::MultiFab::Copy(eplotfile,efield[2],0,4,1,0);


    amrex::MultiFab::Copy(cplotfile,particleInstant,0,0,14,0);
    amrex::MultiFab::Copy(cplotfile,particleMeans,0,14,14,0);
    amrex::MultiFab::Copy(cplotfile,particleVars,0,28,18,0);

    for (int l=0; l<nspecies; ++l)
    {
        for (int d=0; d<AMREX_SPACEDIM; ++d)
        {
            fvarNames[l*AMREX_SPACEDIM + d] = Concatenate("spec ", l, 0) + Concatenate(" dim ", d, 0);

          //  Print() << fvarNames[l*AMREX_SPACEDIM + d] << "\n";
        }
    }

    evarNames[0] = "chargeInstant";
    evarNames[1] = "potentialInstant";
    evarNames[2] = "ExInstant";
    evarNames[3] = "EyInstant";

#if (AMREX_SPACEDIM==3)
    evarNames[4] = "EzInstant";
#endif

    cvarNames[0] = "membersInstant";
    cvarNames[1] = "densityInstant";
    cvarNames[2] = "velxInstant";
    cvarNames[3] = "velyInstant";
    cvarNames[4] = "velzInstant";
    cvarNames[5] = "temperatureInstant";
    cvarNames[6] = "jxInstant";
    cvarNames[7] = "jyInstant";
    cvarNames[8] = "jzInstant";
    cvarNames[9] = "energyInstant";
    cvarNames[10] = "pressureInstant";
    cvarNames[11] = "ixInstant";
    cvarNames[12] = "iyInstant";
    cvarNames[13] = "izInstant";

    cvarNames[14] = "membersMean";
    cvarNames[15] = "densityMean";
    cvarNames[16] = "velxMean";
    cvarNames[17] = "velyMean";
    cvarNames[18] = "velzMean";
    cvarNames[19] = "temperatureMean";
    cvarNames[20] = "jxMean";
    cvarNames[21] = "jyMean";
    cvarNames[22] = "jzMean";
    cvarNames[23] = "energyMean";
    cvarNames[24] = "pressureMean";
    cvarNames[25] = "ixMean";
    cvarNames[26] = "iyMean";
    cvarNames[27] = "izMean";

    cvarNames[28] = "membersVar";
    cvarNames[29] = "densityVar";
    cvarNames[30] = "velxVar";
    cvarNames[31] = "velyVar";
    cvarNames[32] = "velzVar";
    cvarNames[33] = "temperatureVar";
    cvarNames[34] = "jxVar";
    cvarNames[35] = "jyVar";
    cvarNames[36] = "jzVar";
    cvarNames[37] = "energyVar";
    cvarNames[38] = "pressureVar";
    cvarNames[39] = "GVar";
    cvarNames[40] = "KGCross";
    cvarNames[41] = "KRhoCross";
    cvarNames[42] = "RhoGCross";
    cvarNames[43] = "ixVar";
    cvarNames[44] = "iyVar";
    cvarNames[45] = "izVar";


    WriteSingleLevelPlotfile(cplotfilename,cplotfile,cvarNames,cgeom,time,step);

    //WriteSingleLevelPlotfile(mplotfilename,mobility,fvarNames,geom,time,step);

    WriteSingleLevelPlotfile(eplotfilename,eplotfile,evarNames,egeom,time,step);

    particles.Checkpoint(pplotfilename, "particle0");

    if(plot_ascii == 1)
    {
        std::string asciiName = Concatenate("asciiMeans",step,9);
        outputMFAscii(particleMeans, asciiName);

        std::string asciiPName = Concatenate("asciiParticles",step,9);
        particles.WriteParticlesAscii(asciiPName);
    }

}
