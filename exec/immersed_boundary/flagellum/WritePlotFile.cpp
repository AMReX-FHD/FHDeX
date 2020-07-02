#include "main_driver.H"

#include "AMReX_PlotFileUtil.H"

#include "AMReX_MultiFab.H"

#include "common_functions.H"



void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM > & umac,
                   std::array< MultiFab, AMREX_SPACEDIM > & umac_avg,
                   std::array< MultiFab, AMREX_SPACEDIM > & force_ib,
                   const MultiFab& pres,
                   const IBMarkerContainer & ib_pc)
{

    BL_PROFILE_VAR("WritePlotFile()", WritePlotFile);

    const std::string plotfilename = Concatenate(plot_base_name,step,7);

    BoxArray ba = pres.boxArray();
    DistributionMapping dmap = pres.DistributionMap();

    // plot all the velocity variables (cell-center average)
    // plot all the velocity variables (cc + time-averaged)
    // plot all spread forces (cell-centere average)
    // plot pressure
    // plot divergence
    int nPlot = 3*AMREX_SPACEDIM + 2;

    MultiFab plotfile(ba, dmap, nPlot, 0);

    Vector<std::string> varNames(nPlot);

    // keep a counter for plotfile variables
    int cnt = 0;

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        std::string x = "cc_vel";
        x += (120+i); // 120 is x in ASCII => 120+i = x, y, z
        varNames[cnt++] = x;
    }

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        std::string x = "avg_vel";
        x += (120+i); // 120 is x in ASCII => 120+i = x, y, z
        varNames[cnt++] = x;
    }

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        std::string x = "cc_force";
        x += (120+i); // 120 is x in ASCII => 120+i = x, y, z
        varNames[cnt++] = x;
    }


    varNames[cnt++] = "pres";
    varNames[cnt++] = "divergence";

    // reset plotfile variable counter
    cnt = 0;

    // average staggered velocities to cell-centers and copy into plotfile
    AverageFaceToCC(umac, plotfile, cnt);
    cnt+=AMREX_SPACEDIM;

    // average staggered average (in time) velocities to cell-centers and copy
    // into plotfile
    AverageFaceToCC(umac_avg, plotfile, cnt);
    cnt+=AMREX_SPACEDIM;

    AverageFaceToCC(force_ib, plotfile, cnt);
    cnt+=AMREX_SPACEDIM;

    // copy pressure into plotfile
    MultiFab::Copy(plotfile, pres, 0, cnt, 1, 0);
    cnt++;

    // compute divergence and store result in plotfile
    ComputeDiv(plotfile, umac, 0, cnt, 1, geom, 0);

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

    // add immersed boundary markers data to plot file
    ib_pc.WritePlotFile(plotfilename, "immbdy_markers",
                        IBMReal::names(), IBMInt::names());
}
