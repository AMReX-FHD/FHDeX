#include "AMReX_PlotFileUtil.H"

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
	           const amrex::MultiFab& cu,
	           const amrex::MultiFab& cuMeans,
	           const amrex::MultiFab& cuVars,
	           const amrex::MultiFab& prim,
	           const amrex::MultiFab& primMeans,
	           const amrex::MultiFab& primVars,
	           const amrex::MultiFab& spatialCross) 
{

    amrex::Print() << "Start plot: " << step << "\n"; 

    amrex::BoxArray ba = cuMeans.boxArray();
    amrex::DistributionMapping dmap = cuMeans.DistributionMap();

    amrex::MultiFab plotfile(ba, dmap, 34, 0);

    amrex::Print() << "Start copy 1: " << step << "\n";

    amrex::MultiFab::Copy(plotfile,cuMeans,0,0,5,0);

    amrex::Print() << "Start copy 2: " << step << "\n";
    amrex::MultiFab::Copy(plotfile,primMeans,0,5,6,0);
    amrex::Print() << "Start copy 3: " << step << "\n";
    amrex::MultiFab::Copy(plotfile,cu,0,11,5,0);
    amrex::Print() << "Start copy 4: " << step << "\n";
    amrex::MultiFab::Copy(plotfile,prim,0,16,5,0);
    amrex::Print() << "Start copy 5: " << step << "\n";
    amrex::MultiFab::Copy(plotfile,cuVars,0,21,5,0);
    amrex::Print() << "Start copy 6: " << step << "\n";
    amrex::MultiFab::Copy(plotfile,primVars,0,26,5,0);
    amrex::Print() << "Start copy 7: " << step << "\n";
    amrex::MultiFab::Copy(plotfile,spatialCross,0,31,3,0);

    amrex::Print() << "End copy: " << step << "\n";

    std::string plotfilename = amrex::Concatenate("plt",step,9);

    amrex::Vector<std::string> varNames(34);

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
    varNames[10] = "pMean";

    amrex::Print() << "varN: " << varNames[10] << "\n";

    varNames[11] = "rhoInstant";
    varNames[12] = "jxInstant";
    varNames[13] = "jyInstant";
    varNames[14] = "jzInstant";
    varNames[15] = "eInstant";

    amrex::Print() << "varN: " << varNames[15] << "\n";

    varNames[16] = "rhoInstant";
    varNames[17] = "uxInstant";
    varNames[18] = "uyInstant";
    varNames[19] = "uzInstant";
    varNames[20] = "tInstant";

    amrex::Print() << "varN: " << varNames[20] << "\n";

    varNames[21] = "rhoVar";
    varNames[22] = "jxVar";
    varNames[23] = "jyVar";
    varNames[24] = "jzVar";
    varNames[25] = "eVar";

    amrex::Print() << "varN: " << varNames[25] << "\n";

    varNames[26] = "rhoVar";
    varNames[27] = "uxVar";
    varNames[28] = "uyVar";
    varNames[29] = "uzVar";
    varNames[30] = "tVar";

    amrex::Print() << "varN: " << varNames[30] << "\n";

    varNames[31] = "spatialCross1";
    varNames[32] = "spatialCross2";

    amrex::Print() << "varN: " << varNames[32] << "\n";

    varNames[33] = "spatialCross3";

    amrex::Print() << "varN: " << varNames[33] << "\n";

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

    amrex::Print() << "End plot: " << step << "\n";
}
