#include "INS_functions.H"
#include "common_functions.H"

#include "AMReX_PlotFileUtil.H"

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const amrex::Geometry cgeom,
                   const MultiFab& rhotot,
                   const std::array< MultiFab, AMREX_SPACEDIM >& umac,
	               const MultiFab& div,
                   const MultiFab& particleMembers,
                   const MultiFab& particleDensity,
                   const std::array< MultiFab, 3 >& particleVelocity,
                   const MultiFab& particleTemperature,
                   const MultiFab& particlePressure,
                   const MultiFab& particleEnergy,
                   const MultiFab& spatialCorrelation,
                   const MultiFab& densityVar,
                   const MultiFab& density,
                   const MultiFab& particleMembraneFlux,
                   FhdParticleContainer& particles) 
{

    std::string plotfilename = Concatenate("plt",step,8);
    std::string cplotfilename = Concatenate("cplt",step,8);

    BoxArray ba = rhotot.boxArray();
    BoxArray cba = particleMembers.boxArray();

    DistributionMapping dmap = rhotot.DistributionMap();
    DistributionMapping cdmap = particleMembers.DistributionMap();

    int nPlot = 2+AMREX_SPACEDIM;
    int cnPlot = 12;

    MultiFab plotfile(ba, dmap, nPlot, 0);
    MultiFab cplotfile(cba, cdmap, cnPlot, 0);

    Vector<std::string> varNames(nPlot);
    Vector<std::string> cvarNames(cnPlot);

    // keep a counter for plotfile variables
    int cnt = 0;

    varNames[cnt++] = "rhotot";

    cvarNames[0] = "particles";
    cvarNames[1] = "particle density";
    cvarNames[2] = "particle xVel";
    cvarNames[3] = "particle yVel";
    cvarNames[4] = "particle zVel";
    cvarNames[5] = "particle temperature";
    cvarNames[6] = "particle pressure";
    cvarNames[7] = "particle energy";
    cvarNames[8] = "spatial correlation";
    cvarNames[9] = "density variance";
    cvarNames[10] = "density instant";
    cvarNames[11] = "diff flux";

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        std::string x = "vel";
        x += (120+i);
        varNames[cnt++] = x;
    }

	varNames[cnt++] = "div";

    // reset plotfile variable counter
    cnt = 0;


    // copy rhotot into plotfile
    plotfile.copy(rhotot,0,cnt,1);
    cnt++;

    cplotfile.copy(particleMembers,0,0,1); 
    cplotfile.copy(particleDensity,0,1,1);
    cplotfile.copy(particleVelocity[0],0,2,1);
    cplotfile.copy(particleVelocity[1],0,3,1);
    cplotfile.copy(particleVelocity[2],0,4,1);
    cplotfile.copy(particleTemperature,0,5,1);
    cplotfile.copy(particlePressure,0,6,1);
    cplotfile.copy(particleEnergy,0,7,1);
    cplotfile.copy(spatialCorrelation,0,8,1);
    cplotfile.copy(densityVar,0,9,1);
    cplotfile.copy(density,0,10,1);
    cplotfile.copy(particleMembraneFlux,0,11,1);

    // average staggered velocities to cell-centers and copy into plotfile
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        AverageFaceToCC(umac[i],0,plotfile,cnt,1);
        cnt++;
    }

    plotfile.copy(div,0,cnt,1);

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

    cplotfile.mult(0.001,1,1);    
    cplotfile.mult(10,7,1);
    cplotfile.mult(0.000001,9,1);   
    cplotfile.mult(0.001,10,1);

    WriteSingleLevelPlotfile(cplotfilename,cplotfile,cvarNames,cgeom,time,step);


    particles.Checkpoint(plotfilename, "particle0");

}
