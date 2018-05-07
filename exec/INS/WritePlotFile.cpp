#include "INS_functions.H"
#include "common_functions.H"

#include "AMReX_PlotFileUtil.H"

void WritePlotFile(const int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const MultiFab& rhotot,
                   const std::array< MultiFab, AMREX_SPACEDIM >& umac) {

    std::string plotfilename = Concatenate("plt",step,7);

    BoxArray ba = rhotot.boxArray();
    DistributionMapping dmap = rhotot.DistributionMap();

    // rhotot plus velocity
    int nPlot = 1+AMREX_SPACEDIM;

    MultiFab plotfile(ba, dmap, nPlot, 0);

    Vector<std::string> varNames(nPlot);

    // keep a counter for plotfile variables
    int cnt = 0;

    varNames[cnt++] = "rhotot";

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        std::string x = "vel";
        x += (120+i);
        varNames[cnt++] = x;
    }

    // reset plotfile variable counter
    cnt = 0;

    // copy rhotot into plotfile
    plotfile.copy(rhotot,0,cnt,1);
    cnt++;

    // average staggered velocities to cell-centers and copy into plotfile
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        AverageFaceToCC(umac[i],0,plotfile,cnt,1);
        cnt++;
    }

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

}
