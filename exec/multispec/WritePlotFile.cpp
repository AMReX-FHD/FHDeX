#include "multispec_test_functions.H"

#include "AMReX_PlotFileUtil.H"

#include "common_functions.H"

#include "multispec_namespace.H"

using namespace multispec;

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM >& umac,
                   const MultiFab& rhotot,
                   const MultiFab& rho,
                   const MultiFab& pres,
                   const MultiFab& charge,
                   const MultiFab& Epot)
{

    BL_PROFILE_VAR("WritePlotFile()",WritePlotFile);

    const std::string plotfilename = Concatenate(plot_base_name,step,7);

    BoxArray ba = pres.boxArray();
    DistributionMapping dmap = pres.DistributionMap();

    // rhotot        1
    // rho           nspecies
    // c             nspecies
    // averaged vel  AMREX_SPACEDIM
    // shifted  vel  AMREX_SPACEDIM
    // pres          1
    int nPlot = 2*AMREX_SPACEDIM + 2*nspecies + 2;

    if (use_charged_fluid) {
        // charge
        // Epot
        nPlot += 2;
    }

    MultiFab plotfile(ba, dmap, nPlot, 0);

    Vector<std::string> varNames(nPlot);

    // keep a counter for plotfile variables
    int cnt = 0;

    varNames[cnt++] = "rho";

    for (int i=0; i<nspecies; ++i) {
        std::string x = "rho";
        x += (49+i);
        varNames[cnt++] = x;
    }

    for (int i=0; i<nspecies; ++i) {
        std::string x = "c";
        x += (49+i);
        varNames[cnt++] = x;
    }

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

    varNames[cnt++] = "pres";

    if (use_charged_fluid) {
        varNames[cnt++] = "charge";
        varNames[cnt++] = "Epot";
    }

    // reset plotfile variable counter
    cnt = 0;

    // copy rhotot into plotfile
    MultiFab::Copy(plotfile, rhotot, 0, cnt, 1, 0);
    cnt++;

    // copy densities into plotfile
    for (int i=0; i<nspecies; ++i) {
        MultiFab::Copy(plotfile, rho, i, cnt, 1, 0);
        cnt++;
    }

    // copy densities and convert to concentrations
    for (int i=0; i<nspecies; ++i) {
        MultiFab::Copy(plotfile, rho, i, cnt, 1, 0);
        MultiFab::Divide(plotfile, rhotot, 0, cnt, 1, 0);
        cnt++;
    }

    // average staggered velocities to cell-centers and copy into plotfile
    AverageFaceToCC(umac,plotfile,cnt);
    cnt+=AMREX_SPACEDIM;

    // shift staggered velocities to cell-centers and copy into plotfile
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        ShiftFaceToCC(umac[i],0,plotfile,cnt,1);
        cnt++;
    }

    // copy pressure into plotfile
    MultiFab::Copy(plotfile, pres, 0, cnt, 1, 0);
    cnt++;

    // copy charge and Epot into plotfile
    if (use_charged_fluid) {
        MultiFab::Copy(plotfile, charge, 0, cnt, 1, 0);
        cnt++;
        MultiFab::Copy(plotfile, Epot, 0, cnt, 1, 0);
        cnt++;
    }

    // write a plotfile
    // timer
    Real t1 = ParallelDescriptor::second();

    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

    Real t2 = ParallelDescriptor::second() - t1;
    ParallelDescriptor::ReduceRealMax(t2);
    amrex::Print() << "Time spent writing plotfile " << t2 << std::endl;

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