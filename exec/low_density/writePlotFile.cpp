#include "AMReX_PlotFileUtil.H"

#include "common_functions.H"


void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
	           const amrex::MultiFab& cu,
	           const amrex::MultiFab& cuMeans,
	           const amrex::MultiFab& cuVars) 
{

    int cnt, numvars, i = 0;
    int nplot = 3*5 + 4*6 + 2*1 + 3*nspecies;

    amrex::BoxArray ba = cuMeans.boxArray();
    amrex::DistributionMapping dmap = cuMeans.DistributionMap();

    amrex::MultiFab plotfile(ba, dmap, nplot, 0);

    std::string x;
    std::string plotfilename = amrex::Concatenate("plt",step,9);
    amrex::Vector<std::string> varNames(nplot);

    // Load into plotfile MF
    
    cnt = 0;

    numvars = 1;
    amrex::MultiFab::Copy(plotfile,cu,0,cnt,numvars,0);
    cnt+=numvars;

    numvars = 1;
    amrex::MultiFab::Copy(plotfile,cuMeans,0,cnt,numvars,0);
    cnt+=numvars;

    numvars = 1;
    amrex::MultiFab::Copy(plotfile,cuVars,0,cnt,numvars,0);


    // Set variable names

    cnt = 0;

    varNames[cnt++] = "rhoInstant";
    varNames[cnt++] = "rhoMean";
    varNames[cnt++] = "rhoVar";

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

}
