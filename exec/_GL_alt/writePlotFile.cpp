#include "AMReX_PlotFileUtil.H"

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
	           const amrex::MultiFab& phi)
{



    amrex::BoxArray ba = phi.boxArray();
    amrex::DistributionMapping dmap = phi.DistributionMap();

    amrex::MultiFab plotfile(ba, dmap, 1, 0);

    amrex::MultiFab::Copy(plotfile,phi,0,0,1,0);


    std::string plotfilename = amrex::Concatenate("plt",step,9);

    amrex::Vector<std::string> varNames(1);

    varNames[0] = "phi";

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

}
