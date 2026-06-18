#include "AMReX_PlotFileUtil.H"

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const amrex::MultiFab& phi,
                   amrex::Real& umbrella, amrex::Real& phi0)
{

    amrex::BoxArray ba = phi.boxArray();
    amrex::DistributionMapping dmap = phi.DistributionMap();

    amrex::MultiFab plotfile(ba, dmap, 1, 0);

    amrex::MultiFab::Copy(plotfile,phi,0,0,1,0);

    std::string plotfilename = amrex::Concatenate("plt",step,9);

    amrex::Vector<std::string> varNames(1);

    varNames[0] = "Conc";

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

    std::string scalar_param_loc=plotfilename;
    scalar_param_loc=scalar_param_loc+ "/";
    scalar_param_loc=scalar_param_loc+"params.txt";

    std::ofstream file_scalar(scalar_param_loc.c_str());
    file_scalar << umbrella  << "\n" ;
    file_scalar << phi0 << "\n" << std::endl;

}