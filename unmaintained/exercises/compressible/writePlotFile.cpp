#include "AMReX_PlotFileUtil.H"

#include "common_functions.H"



void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const amrex::MultiFab& cu,
                   const amrex::MultiFab& cuMeans,
                   const amrex::MultiFab& cuVars,
                   const amrex::MultiFab& prim,
                   const amrex::MultiFab& primMeans,
                   const amrex::MultiFab& primVars,
                   const amrex::MultiFab& eta,
                   const amrex::MultiFab& kappa)
{

    int cnt, numvars, i = 0;

    // instantaneous values
    // 5 + nspecies (conserved)
    // 6 + 2*nspecies (primitive)
    // 2 (eta and kappa)
    int nplot = (5+nspecies) + (6+2*nspecies) + 2;

    if (plot_means == 1) {
        nplot += 11;
    }

    if (plot_vars == 1) {
        nplot += 10;
    }

    amrex::BoxArray ba = cuMeans.boxArray();
    amrex::DistributionMapping dmap = cuMeans.DistributionMap();

    amrex::MultiFab plotfile(ba, dmap, nplot, 0);

    std::string x;
    std::string plotfilename = amrex::Concatenate(plot_base_name,step,9);
    amrex::Vector<std::string> varNames(nplot);

    // Load into plotfile MF

    cnt = 0;

    // instantaneous values of conserved variables
    // rho, jx, jy, jz, e, rhoYk
    numvars = 5+nspecies;
    amrex::MultiFab::Copy(plotfile,cu,0,cnt,numvars,0);
    cnt+=numvars;

    // instantaneous values of primitive variables
    // rho, ux, uy, uz, temp, pres, Yk, Xk
    numvars = 6+2*nspecies;
    amrex::MultiFab::Copy(plotfile,prim,0,cnt,numvars,0);
    cnt+=numvars;

    if (plot_means == 1) {

        // mean values of conserved variables
        // rho, jx, jy, jz, e
        numvars = 5;
        amrex::MultiFab::Copy(plotfile,cuMeans,0,cnt,numvars,0);
        cnt+=numvars;

        // mean values of primitive variables
        // rho, ux, uy, uz, temp, pres
        numvars = 6;
        amrex::MultiFab::Copy(plotfile,primMeans,0,cnt,numvars,0);
        cnt+=numvars;
    }

    if (plot_vars == 1) {

        // variance of conserved variables
        // rho, jx, jy, jz, e
        numvars = 5;
        amrex::MultiFab::Copy(plotfile,cuVars,0,cnt,numvars,0);
        cnt+=numvars;

        // variances of primitive variables
        // rho, ux, uy, uz, temp
        numvars = 5;
        amrex::MultiFab::Copy(plotfile,primVars,0,cnt,numvars,0);
        cnt+=numvars;
    }

    // eta
    numvars = 1;
    amrex::MultiFab::Copy(plotfile,eta,0,cnt,numvars,0);
    cnt+=numvars;

    // kappa
    numvars = 1;
    amrex::MultiFab::Copy(plotfile,kappa,0,cnt,numvars,0);
    cnt+=numvars;

    // Set variable names
    cnt = 0;

    varNames[cnt++] = "rhoInstant";
    varNames[cnt++] = "jxInstant";
    varNames[cnt++] = "jyInstant";
    varNames[cnt++] = "jzInstant";
    varNames[cnt++] = "rhoEInstant";
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

    if (plot_means == 1) {
        varNames[cnt++] = "rhoMean";
        varNames[cnt++] = "jxMean";
        varNames[cnt++] = "jyMean";
        varNames[cnt++] = "jzMean";
        varNames[cnt++] = "rhoEMean";

        varNames[cnt++] = "rhoMean";
        varNames[cnt++] = "uxMean";
        varNames[cnt++] = "uyMean";
        varNames[cnt++] = "uzMean";
        varNames[cnt++] = "tMean";
        varNames[cnt++] = "pMean";
    }

    if (plot_vars == 1) {
        varNames[cnt++] = "rhoVar";
        varNames[cnt++] = "jxVar";
        varNames[cnt++] = "jyVar";
        varNames[cnt++] = "jzVar";
        varNames[cnt++] = "rhoEVar";

        varNames[cnt++] = "rhoVar";
        varNames[cnt++] = "uxVar";
        varNames[cnt++] = "uyVar";
        varNames[cnt++] = "uzVar";
        varNames[cnt++] = "tVar";
    }

    varNames[cnt++] = "eta";
    varNames[cnt++] = "kappa";

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

}