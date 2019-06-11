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
	           const amrex::MultiFab& spatialCross, const amrex::MultiFab& etaMeanAv, const amrex::MultiFab& kappaMeanAv) 
{



    amrex::BoxArray ba = cuMeans.boxArray();
    amrex::DistributionMapping dmap = cuMeans.DistributionMap();

    amrex::MultiFab plotfile(ba, dmap, 40, 0);

    amrex::MultiFab::Copy(plotfile,cuMeans,0,0,5,0);
    amrex::MultiFab::Copy(plotfile,primMeans,0,5,6,0);

    amrex::MultiFab::Copy(plotfile,cu,0,11,5,0);

    amrex::MultiFab::Copy(plotfile,prim,0,16,6,0);

    amrex::MultiFab::Copy(plotfile,cuVars,0,22,5,0);

    amrex::MultiFab::Copy(plotfile,primVars,0,27,5,0);

    amrex::MultiFab::Copy(plotfile,spatialCross,0,32,6,0);

    amrex::MultiFab::Copy(plotfile,etaMeanAv,0,38,1,0);
    amrex::MultiFab::Copy(plotfile,kappaMeanAv,0,39,1,0);


    std::string plotfilename = amrex::Concatenate("plt",step,9);

    amrex::Vector<std::string> varNames(40);

    int cnt = 0;

    varNames[cnt++] = "rhoMean";
    varNames[cnt++] = "jxMean";
    varNames[cnt++] = "jyMean";
    varNames[cnt++] = "jzMean";
    varNames[cnt++] = "eMean";

    varNames[cnt++] = "rhoMean";
    varNames[cnt++] = "uxMean";
    varNames[cnt++] = "uyMean";
    varNames[cnt++] = "uzMean";
    varNames[cnt++] = "tMean";
    varNames[cnt++] = "pMean";

    varNames[cnt++] = "rhoInstant";
    varNames[cnt++] = "jxInstant";
    varNames[cnt++] = "jyInstant";
    varNames[cnt++] = "jzInstant";
    varNames[cnt++] = "eInstant";

    varNames[cnt++] = "rhoInstant";
    varNames[cnt++] = "uxInstant";
    varNames[cnt++] = "uyInstant";
    varNames[cnt++] = "uzInstant";
    varNames[cnt++] = "tInstant";
    varNames[cnt++] = "pInstant";

    varNames[cnt++] = "rhoVar";
    varNames[cnt++] = "jxVar";
    varNames[cnt++] = "jyVar";
    varNames[cnt++] = "jzVar";
    varNames[cnt++] = "eVar";

    varNames[cnt++] = "rhoVar";
    varNames[cnt++] = "uxVar";
    varNames[cnt++] = "uyVar";
    varNames[cnt++] = "uzVar";
    varNames[cnt++] = "tVar";

    varNames[cnt++] = "Energy-densityCross";
    varNames[cnt++] = "Energy-energyCross";
    varNames[cnt++] = "Momentum-densityCross";

    varNames[cnt++] = "Temperature-temperatureCross";
    varNames[cnt++] = "Temperature-densityCross";
    varNames[cnt++] = "Velelocity-densityCross";

    varNames[cnt++] = "eta";
    varNames[cnt++] = "kappa";


    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

}
