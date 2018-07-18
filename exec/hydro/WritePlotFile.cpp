#include "hydro_test_functions.H"

#include "AMReX_PlotFileUtil.H"

#include "AMReX_MultiFab.H"

#include "common_functions.H"

#include "common_namespace.H"

using namespace common;

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const std::array< MultiFab, AMREX_SPACEDIM >& umac,
		   const MultiFab& pres)
{
    const std::string plotfilename = Concatenate("plt",step,7);

    BoxArray ba = umac[0].boxArray();
    DistributionMapping dmap = umac[0].DistributionMap();

    // plot all the velocity variables
    int nPlot = AMREX_SPACEDIM+1;

    MultiFab plotfile(ba, dmap, nPlot, 0);

    Vector<std::string> varNames(nPlot);

    // keep a counter for plotfile variables
    int cnt = 0;

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        std::string x = "vel";
        x += (120+i);
        varNames[cnt++] = x;
    }

    varNames[nPlot-1] = "pres";

    // reset plotfile variable counter
    cnt = 0;

    // average staggered velocities to cell-centers and copy into plotfile
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        AverageFaceToCC(umac[i],0,plotfile,cnt,1);
        cnt++;
    }

    MultiFab::Copy(plotfile, pres, 0, nPlot-1, 1, 0);

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

    // staggered velocity
    const std::string plotfilenamex = Concatenate("stagx",step,7);
    const std::string plotfilenamey = Concatenate("stagy",step,7);
    const std::string plotfilenamez = Concatenate("stagz",step,7);

    WriteSingleLevelPlotfile(plotfilenamex,umac[0],{"umac"},geom,time,step);
    WriteSingleLevelPlotfile(plotfilenamey,umac[1],{"vmac"},geom,time,step);
#if (AMREX_SPACEDIM == 3)
    WriteSingleLevelPlotfile(plotfilenamez,umac[2],{"wmac"},geom,time,step);
#endif

}
