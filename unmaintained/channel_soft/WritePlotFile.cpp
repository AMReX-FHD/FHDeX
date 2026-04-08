#include "main_driver.H"

#include "AMReX_PlotFileUtil.H"

#include "AMReX_MultiFab.H"

#include "common_functions.H"


#include "AmrCoreAdv.H"


void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const std::array< MultiFab, AMREX_SPACEDIM > & umac,
                   const MultiFab & tracer,
                   const MultiFab & pres,
                   const std::array< MultiFab, AMREX_SPACEDIM> & force_ibm,
                   const IBParticleContainer & ib_pc,
                   AmrCoreAdv & amr_core_adv,
                   int lev
                    )
{

    Vector< std::unique_ptr<MultiFab> > conNew;
    Vector< std::unique_ptr<MultiFab> > DconNew_x;
    Vector< std::unique_ptr<MultiFab> > DconNew_y;
    Vector< std::unique_ptr<MultiFab> > DconNew_z;
    Vector< std::unique_ptr<MultiFab> > magDC;
    Vector< std::unique_ptr<MultiFab> > DconNewc_x;
    Vector< std::unique_ptr<MultiFab> > DconNewc_y;
    Vector< std::unique_ptr<MultiFab> > DconNewc_z;

    Vector< std::unique_ptr<MultiFab> > mf(lev+1);

    BL_PROFILE_VAR("WritePlotFile()",WritePlotFile);

    const std::string plotfilename = Concatenate(plot_base_name,step,7);

    BoxArray ba = pres.boxArray();
    DistributionMapping dmap = pres.DistributionMap();

    // plot all the velocity variables (averaged)
    // plot all the velocity variables (shifted)
    // plot pressure
    // plot tracer
    // plot divergence
    // plot concentration and it's gradient tangent to the level set
    int nPlot = 5*AMREX_SPACEDIM+8*(lev+1);

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

    varNames[cnt++] = "tracer";
    varNames[cnt++] = "pres";
    varNames[cnt++] = "divergence";
    varNames[cnt++] = "force_ibm_x";
    varNames[cnt++] = "force_ibm_y";
    varNames[cnt++] = "force_ibm_z";
    varNames[cnt++] = "shifted_force_ibm_x";
    varNames[cnt++] = "shifted_force_ibm_y";
    varNames[cnt++] = "shifted_force_ibm_z";
    varNames[cnt++] = "C";
    varNames[cnt++] = "dCdx";
    varNames[cnt++] = "dCdy";
    varNames[cnt++] = "dCdz";
    varNames[cnt++] = "|Dc|";
    varNames[cnt++] = "dCdx_cen";
    varNames[cnt++] = "dCdy_cen";
    varNames[cnt++] = "dCdz_cen";

   // for (int i=0; i<nPlot; ++i) {
   //     std::cout<<" i= "<< varNames[i]<<std::endl ;
   // }


    // reset plotfile variable counter
    cnt = 0;

    // average staggered velocities to cell-centers and copy into plotfile
    AverageFaceToCC(umac, plotfile, cnt);
    cnt+=AMREX_SPACEDIM;

    // shift staggered velocities to cell-centers and copy into plotfile
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        ShiftFaceToCC(umac[d], 0, plotfile, cnt, 1);
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

    AverageFaceToCC(force_ibm, plotfile, cnt);
    cnt+=AMREX_SPACEDIM;

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        ShiftFaceToCC(force_ibm[d], 0, plotfile, cnt, 1);
        cnt++;
    }
    // copy concentration in to the plot file
    amr_core_adv.con_new_copy(lev, mf, 0);
    MultiFab::Copy(plotfile, * mf[lev], 0, cnt, 1, 0);
    cnt++;
    // copy x component of concentration gradient tangent to the level set in to the plot file, face centered
    amr_core_adv.con_new_copy(lev, mf,1);
    ShiftFaceToCC(*mf[lev], 0, plotfile, cnt, 1);
    cnt++;
    // copy y component of concentration gradient tangent to the level set in to the plot file, face centered

    amr_core_adv.con_new_copy(lev, mf, 2);
    ShiftFaceToCC(*mf[lev], 0, plotfile, cnt, 1);
    cnt++;
    // copy z component of concentration gradient tangent to the level set in to the plot file, face centered

    amr_core_adv.con_new_copy(lev, mf, 3);
    ShiftFaceToCC(*mf[lev], 0, plotfile, cnt, 1);
    cnt++;
    // copy magnitude of concentration gradient tangent to the level set in to the plot file, cell centered

    amr_core_adv.con_new_copy(lev, mf, 4);
    MultiFab::Copy(plotfile, * mf[lev], 0, cnt, 1, 0);
    cnt++;
    // copy x component of concentration gradient  tangent to the level set in to the plot file, cell centered

    amr_core_adv.con_new_copy(lev, mf, 5);
    MultiFab::Copy(plotfile, * mf[lev], 0, cnt, 1, 0);
    cnt++;
    // copy y component of concentration gradient tangent to the level set in to the plot file, cell centered
    amr_core_adv.con_new_copy(lev, mf, 6);
    MultiFab::Copy(plotfile, * mf[lev], 0, cnt, 1, 0);
    cnt++;
    // copy z component of concentration gradient tangent to the level set in to the plot file, cell centered
    amr_core_adv.con_new_copy(lev, mf, 7);
    MultiFab::Copy(plotfile, * mf[lev], 0, cnt, 1, 0);
    cnt++;


    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename, plotfile, varNames, geom, time, step);

    // add particle data to plot file
    ib_pc.WritePlotFile(plotfilename, "immersed_boundaries",
                        IBP_realData::names(), IBP_intData::names());


    // staggered velocity
    if (plot_stag == 1) {
      const std::string plotfilenamex = Concatenate("stagx", step, 7);
      const std::string plotfilenamey = Concatenate("stagy", step, 7);
      const std::string plotfilenamez = Concatenate("stagz", step, 7);

      WriteSingleLevelPlotfile(plotfilenamex, umac[0], {"umac"}, geom, time, step);
      WriteSingleLevelPlotfile(plotfilenamey, umac[1], {"vmac"}, geom, time, step);
#if (AMREX_SPACEDIM == 3)
      WriteSingleLevelPlotfile(plotfilenamez, umac[2], {"wmac"}, geom, time, step);
#endif
    }

}
