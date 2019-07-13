#include "AMReX_PlotFileUtil.H"

#include "common_functions.H"
#include "common_namespace.H"

using namespace common;

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
	           const amrex::MultiFab& cu,
	           const amrex::MultiFab& cuMeans,
	           const amrex::MultiFab& cuVars,
	           const amrex::MultiFab& prim,
	           const amrex::MultiFab& primMeans,
	           const amrex::MultiFab& primVars,
		   const amrex::MultiFab& etaMeanAv, 
		   const amrex::MultiFab& kappaMeanAv) 
{

    int cnt, numvars, i = 0;
    int nplot = 3*5 + 3*6 + 2*1 + 3*nspecies;

    amrex::BoxArray ba = cuMeans.boxArray();
    amrex::DistributionMapping dmap = cuMeans.DistributionMap();

    amrex::MultiFab plotfile(ba, dmap, nplot, 0);

    std::string x;
    std::string plotfilename = amrex::Concatenate("plt",step,9);
    amrex::Vector<std::string> varNames(nplot);

    // Load into plotfile MF
    
    cnt = 0;

    numvars = 5+nspecies;
    amrex::MultiFab::Copy(plotfile,cu,0,cnt,numvars,0);
    cnt+=numvars;

    numvars = 6+2*nspecies;
    amrex::MultiFab::Copy(plotfile,prim,0,cnt,numvars,0);
    cnt+=numvars;

    numvars = 5;
    amrex::MultiFab::Copy(plotfile,cuMeans,0,cnt,numvars,0);
    cnt+=numvars;

    numvars = 6;
    amrex::MultiFab::Copy(plotfile,primMeans,0,cnt,numvars,0);
    cnt+=numvars;

    numvars = 5;
    amrex::MultiFab::Copy(plotfile,cuVars,0,cnt,numvars,0);
    cnt+=numvars;

    numvars = 6;
    amrex::MultiFab::Copy(plotfile,primVars,0,cnt,numvars,0);
    cnt+=numvars;

    numvars = 1;
    amrex::MultiFab::Copy(plotfile,etaMeanAv,0,cnt,numvars,0);
    cnt+=numvars;

    numvars = 1;
    amrex::MultiFab::Copy(plotfile,kappaMeanAv,0,cnt,numvars,0);
    cnt+=numvars;

    // Set variable names

    cnt = 0;

    varNames[cnt++] = "rhoInstant";
    varNames[cnt++] = "jxInstant";
    varNames[cnt++] = "jyInstant";
    varNames[cnt++] = "jzInstant";
    varNames[cnt++] = "eInstant";
    x = "rhoYkInstant_";
    for (i=0; i<nspecies; i++) {
        varNames[cnt] = x;
        varNames[cnt++] += 48+i;
    }

    varNames[cnt++] = "rhoInstant";
    varNames[cnt++] = "uxInstant";
    varNames[cnt++] = "uyInstant";
    varNames[cnt++] = "uzInstant";
    varNames[cnt++] = "tInstant";
    varNames[cnt++] = "pInstant";
    x = "YkInstant_";
    for (i=0; i<nspecies; i++) {
        varNames[cnt] = x;
        varNames[cnt++] += 48+i;
    }
    x = "XkInstant_";
    for (i=0; i<nspecies; i++) {
        varNames[cnt] = x;
        varNames[cnt++] += 48+i;
    }

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

    varNames[cnt++] = "eta";
    varNames[cnt++] = "kappa";

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

}
