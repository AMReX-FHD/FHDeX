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


    amrex::BoxArray ba = cuMeans.boxArray();
    amrex::DistributionMapping dmap = cuMeans.DistributionMap();

    amrex::MultiFab plotfile(ba, dmap, 34, 0);

    amrex::MultiFab::Copy(plotfile,cuMeans,0,0,5,0);
    amrex::MultiFab::Copy(plotfile,primMeans,0,5,6,0);
    amrex::MultiFab::Copy(plotfile,cu,0,11,5,0);
    amrex::MultiFab::Copy(plotfile,prim,0,16,5,0);
    amrex::MultiFab::Copy(plotfile,cuVars,0,21,5,0);
    amrex::MultiFab::Copy(plotfile,primVars,0,26,5,0);
    amrex::MultiFab::Copy(plotfile,spatialCross,0,31,3,0);

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

    varNames[11] = "rhoInstant";
    varNames[12] = "jxInstant";
    varNames[13] = "jyInstant";
    varNames[14] = "jzInstant";
    varNames[15] = "eInstant";

    varNames[16] = "rhoInstant";
    varNames[17] = "uxInstant";
    varNames[18] = "uyInstant";
    varNames[19] = "uzInstant";
    varNames[20] = "tInstant";

    varNames[21] = "rhoVar";
    varNames[22] = "jxVar";
    varNames[23] = "jyVar";
    varNames[24] = "jzVar";
    varNames[25] = "eVar";

    varNames[26] = "rhoVar";
    varNames[27] = "uxVar";
    varNames[28] = "uyVar";
    varNames[29] = "uzVar";
    varNames[30] = "tVar";

    varNames[31] = "spatialCross1";
    varNames[32] = "spatialCross2";
    varNames[33] = "spatialCross3";

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);
}
