#include "INS_functions.H"
#include "common_functions.H"

#include "AMReX_PlotFileUtil.H"

void particleGen(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const MultiFab& particleInstant,
                   const MultiFab& particleMeans,
                   const MultiFab& particleVars,
                   FhdParticleContainer& particles)
{

    std::string cplotfilename = Concatenate("cplt",step,9);
    std::string pplotfilename = Concatenate("pplt",step,9);

    BoxArray cba = particleInstant.boxArray();
    DistributionMapping cdmap = particleInstant.DistributionMap();


//    int cnPlot = 40;
    int cnPlot = 46;

    MultiFab cplotfile(cba, cdmap, cnPlot, 0);

    Vector<std::string> cvarNames(cnPlot);

    amrex::MultiFab::Copy(cplotfile,particleInstant,0,0,14,0);
    amrex::MultiFab::Copy(cplotfile,particleMeans,0,14,14,0);
    amrex::MultiFab::Copy(cplotfile,particleVars,0,28,18,0);

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


    WriteSingleLevelPlotfile(cplotfilename,cplotfile,cvarNames,geom,time,step);

    particles.Checkpoint(pplotfilename, "particle0");

}
