#include "INS_functions.H"
#include "common_functions.H"

#include "AMReX_PlotFileUtil.H"

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const amrex::Geometry cgeom,
                   const MultiFab& particleInstant,
                   const MultiFab& particleMeans,
                   const MultiFab& particleVars,
                   const MultiFab& particleMembraneFlux,
                   FhdParticleContainer& particles) 
{

    std::string cplotfilename = Concatenate("cplt",step,9);
    std::string pplotfilename = Concatenate("pplt",step,9);

    BoxArray cba = particleInstant.boxArray();
    DistributionMapping cdmap = particleInstant.DistributionMap();

 
//    int cnPlot = 40;
    int cnPlot = 44;

    MultiFab cplotfile(cba, cdmap, cnPlot, 0);

    Vector<std::string> cvarNames(cnPlot);

    amrex::MultiFab::Copy(cplotfile,particleInstant,0,0,11,0);
    amrex::MultiFab::Copy(cplotfile,particleMeans,0,11,12,0);
    amrex::MultiFab::Copy(cplotfile,particleVars,0,23,21,0);

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

    cvarNames[11] = "membersMean";
    cvarNames[12] = "densityMean";
    cvarNames[13] = "velxMean";
    cvarNames[14] = "velyMean";
    cvarNames[15] = "velzMean";
    cvarNames[16] = "temperatureMean";
    cvarNames[17] = "jxMean";
    cvarNames[18] = "jyMean";
    cvarNames[19] = "jzMean";
    cvarNames[20] = "energyMean";
    cvarNames[21] = "pressureMean";
    cvarNames[22] = "speedMean";

    cvarNames[23] = "membersVar";
    cvarNames[24] = "densityVar";
    cvarNames[25] = "velxVar";
    cvarNames[26] = "velyVar";
    cvarNames[27] = "velzVar";
    cvarNames[28] = "temperatureVar";
    cvarNames[29] = "jxVar";
    cvarNames[30] = "jyVar";
    cvarNames[31] = "jzVar";
    cvarNames[32] = "energyVar";
    cvarNames[33] = "pressureVar";
    cvarNames[34] = "GVar";
    cvarNames[35] = "KGCross";
    cvarNames[36] = "KRhoCross";
    cvarNames[37] = "RhoGCross";
    cvarNames[38] = "Energy-densityCross";
    cvarNames[39] = "Energy-energyCross";
    cvarNames[40] = "Momentum-densityCross";
    cvarNames[41] = "Temperature-densityCross";
    cvarNames[42] = "Temperature-temperatueCross";
    cvarNames[43] = "velocity-densityCross";

    cplotfile.mult(0.001,2,1);    //cgs coords density
    cplotfile.mult(0.001,12,1);
    cplotfile.mult(0.000001,24,1);   

    cplotfile.mult(100,2,1);  //cgs coords velocity
    cplotfile.mult(100,3,1);
    cplotfile.mult(100,4,1);
    cplotfile.mult(100,13,1);
    cplotfile.mult(100,14,1);
    cplotfile.mult(100,15,1);
    cplotfile.mult(10000,25,1);
    cplotfile.mult(10000,26,1);
    cplotfile.mult(10000,27,1);

    cplotfile.mult(0.1,6,1);  //cgs coords momentum density
    cplotfile.mult(0.1,7,1);
    cplotfile.mult(0.1,8,1);
    cplotfile.mult(0.1,17,1);
    cplotfile.mult(0.1,18,1);
    cplotfile.mult(0.1,19,1);
    cplotfile.mult(0.01,29,1);
    cplotfile.mult(0.01,30,1);
    cplotfile.mult(0.01,31,1);


    cplotfile.mult(10,9,1); //cgs coords energy density
    cplotfile.mult(10,20,1);
    cplotfile.mult(100,32,1);

    cplotfile.mult(0.1,10,1); //cgs coords pressure
    cplotfile.mult(0.1,21,1);
    cplotfile.mult(0.01,33,1);

    cplotfile.mult(10*0.001,38,1); //cgscoords energy/density cross
    cplotfile.mult(10*10,39,1); //cgscoords energy/energy cross
    cplotfile.mult(0.1*0.001,40,1); //cgscoords momentum/density cross

    cplotfile.mult(1*0.001,41,1); //cgscoords temperature/density cross
    cplotfile.mult(1,42,1); //cgscoords temperature/temperature cross
    cplotfile.mult(100*0.001,43,1); //cgscoords velocity/density cross

    WriteSingleLevelPlotfile(cplotfilename,cplotfile,cvarNames,cgeom,time,step);


    particles.Checkpoint(pplotfilename, "particle0");

}
