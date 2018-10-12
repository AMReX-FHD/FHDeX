#include "AMReX_PlotFileUtil.H"

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
	           const amrex::MultiFab& cu,
	           const amrex::MultiFab& cuMeans,
	           const amrex::MultiFab& cuVars,
	           const amrex::MultiFab& prim,
	           const amrex::MultiFab& primMeans,
	           const amrex::MultiFab& primVars) 
{


    amrex::BoxArray ba = cuMeans.boxArray();
    amrex::DistributionMapping dmap = cuMeans.DistributionMap();

    amrex::MultiFab plotfile(ba, dmap, 16, 0);

    amrex::MultiFab::Copy(plotfile,cuMeans,0,0,5,0);
    amrex::MultiFab::Copy(plotfile,primMeans,0,5,5,0);
    amrex::MultiFab::Copy(plotfile,prim,0,10,5,0);

    amrex::MultiFab::Copy(plotfile,cuVars,0,15,1,0);

    std::string plotfilename = amrex::Concatenate("plt",step,8);

    amrex::Vector<std::string> varNames(16);

    varNames[0] = "rhoMean";
    varNames[1] = "jxMean";
    varNames[2] = "jyMean";
    varNames[3] = "jzMean";
    varNames[4] = "eMean";

    varNames[5] = "rhoMean";
    varNames[6] = "uxMean";
    varNames[7] = "uyMean";
    varNames[8] = "uzMean";
    varNames[9] = "tMean";

    varNames[10] = "rhoInstant";
    varNames[11] = "uxInstant";
    varNames[12] = "uyInstant";
    varNames[13] = "uzInstant";
    varNames[14] = "tInstant";

    varNames[15] = "rhoVar";

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);
}
