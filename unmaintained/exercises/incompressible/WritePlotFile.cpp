#include "multispec_test_functions.H"

#include "AMReX_PlotFileUtil.H"

#include "AMReX_MultiFab.H"

#include "common_functions.H"



void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM >& umac,
                   const MultiFab& rho,
                   const MultiFab& rhoxav,
                   const MultiFab& tracer,
                   const MultiFab& pres)
{

    BL_PROFILE_VAR("WritePlotFile()",WritePlotFile);

    const std::string plotfilename = Concatenate(plot_base_name,step,7);

    BoxArray ba = pres.boxArray();
    DistributionMapping dmap = pres.DistributionMap();

    int nspecies = rho.nComp();

    // plot all the velocity variables (averaged)
    // plot all the velocity variables (shifted)
    // plot pressure
    // plot tracer
    // plot divergence
    int nPlot = 2*AMREX_SPACEDIM + nspecies + 3 + nspecies;

    MultiFab plotfile(ba, dmap, nPlot, 0);

    Vector<std::string> varNames(nPlot);

    // keep a counter for plotfile variables
    int cnt = 0;

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        std::string x = "averaged_vel";
        x += (120+i);
        varNames[cnt++] = x;
    }

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        std::string x = "shifted_vel";
        x += (120+i);
        varNames[cnt++] = x;
    }

    for (int i=0; i<nspecies; ++i) {
        std::string x = "rho";
        x += (48+i);
        varNames[cnt++] = x;
    }

    varNames[cnt++] = "tracer";
    varNames[cnt++] = "pres";
    varNames[cnt++] = "divergence";

    for (int i=0; i<nspecies; ++i) {
        std::string x = "rhomean";
        x += (48+i);
        varNames[cnt++] = x;
    }

    // reset plotfile variable counter
    cnt = 0;

    // average staggered velocities to cell-centers and copy into plotfile
    AverageFaceToCC(umac,plotfile,cnt);
    cnt+=AMREX_SPACEDIM;

    // shift staggered velocities to cell-centers and copy into plotfile
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        ShiftFaceToCC(umac[i],0,plotfile,cnt,1);
        cnt++;
    }

    // copy species concentration density into plotfile
    for (int i=0; i<nspecies; ++i) {
        MultiFab::Copy(plotfile, rho, i, cnt, 1, 0);
        cnt++;
    }

    // copy tracer into plotfile
    MultiFab::Copy(plotfile, tracer, 0, cnt, 1, 0);
    cnt++;

    // copy pressure into plotfile
    MultiFab::Copy(plotfile, pres, 0, cnt, 1, 0);
    cnt++;

    // compute divergence and store result in plotfile
    ComputeDiv(plotfile, umac, 0, cnt, 1, geom, 0);
    cnt++;

    for (int i=0; i<nspecies; ++i) {
        MultiFab::Copy(plotfile, rhoxav, i, cnt, 1, 0);
        cnt++;
    }

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

    // staggered velocity
    if (plot_stag == 1) {
        const std::string plotfilenamex = Concatenate("stagx",step,7);
        const std::string plotfilenamey = Concatenate("stagy",step,7);
        const std::string plotfilenamez = Concatenate("stagz",step,7);

        WriteSingleLevelPlotfile(plotfilenamex,umac[0],{"umac"},geom,time,step);
        WriteSingleLevelPlotfile(plotfilenamey,umac[1],{"vmac"},geom,time,step);
#if (AMREX_SPACEDIM == 3)
        WriteSingleLevelPlotfile(plotfilenamez,umac[2],{"wmac"},geom,time,step);
#endif
    }

}