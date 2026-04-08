#include "AMReX_PlotFileUtil.H"
#include "reactDiff_functions.H"

void WritePlotFile(int step,
    const amrex::Real time,
    const amrex::Geometry& geom,
    const MultiFab& n_in)
{

    BL_PROFILE_VAR("WritePlotFile()",WritePlotFile);

    std::string plotfilename = Concatenate(plot_base_name,step,7);

    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    BoxArray ba = n_in.boxArray();
    DistributionMapping dmap = n_in.DistributionMap();

    Vector<std::string> varNames(nspecies);

    // keep a counter for plotfile variables
    int cnt = 0;

    for (int i=0; i<nspecies; ++i) {
        std::string x = "n";
        x += char(49+i);
        varNames[cnt++] = x;
    }

    WriteSingleLevelPlotfile(plotfilename,n_in,varNames,geom,time,step);

}