#include "INS_functions.H"
#include "common_functions.H"
#include "species.H"
#include "DsmcParticleContainer.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;
using namespace std;

void FhdParticleContainer::writePlotFile(const MultiFab& mfcuInst,
                                         const MultiFab& mfcuMeans,
                                         const MultiFab& mfcuVars,
                                         const MultiFab& mfprimInst,
                                         const MultiFab& mfprimMeans,
                                         const MultiFab& mfprimVars,
                                         const MultiFab& mfcoVars,
                                         const MultiFab& mfspatialCorr1d,
                                         const Geometry& geom,
                                         Real time,
                                         const int ncross,
                                         int step) {
    BL_PROFILE_VAR("writePlotFile()",writePlotFile);

    int ncon    = (nspecies+1)*5;
    int nprim   = (nspecies+1)*9;
    int ncovar  = 25;
    int nvars   = ncovar + ncon + nprim; // covariances + prim. vars + cons. vars

    amrex::BoxArray ba = mfcuInst.boxArray();
    amrex::DistributionMapping dmap = mfcuInst.DistributionMap();

    std::string pltcu    = amrex::Concatenate("pltcu",step,12);
    std::string pltprim  = amrex::Concatenate("pltprim",step,12);
    std::string pltvar   = amrex::Concatenate("pltvar",step,12);

    amrex::MultiFab mfcuplt(ba, dmap, 2*ncon, 0);
    amrex::MultiFab mfprimplt(ba, dmap, 2*nprim, 0);
    amrex::MultiFab mfvarplt(ba, dmap, nvars, 0);

    amrex::Vector<std::string> cuNames(2*ncon);
    amrex::Vector<std::string> primNames(2*nprim);
    amrex::Vector<std::string> varNames(nvars);

    //////////////////////////////////////
    // Conserved Means and Instants
    //////////////////////////////////////

    /*
      Conserved Vars:
      0  - rho = (1/V) += m
      1  - Jx  = (1/V) += mu
      2  - Jy  = (1/V) += mv
      3  - Jz  = (1/V) += mw
      4  - K   = (1/V) += m|v|^2
      ... (repeat for each species)
    */

    int cnt = 0;
    // Instant values
    cuNames[cnt++] = "rhoInstant";
    cuNames[cnt++] = "JxInstant";
    cuNames[cnt++] = "JyInstant";
    cuNames[cnt++] = "JzInstant";
    cuNames[cnt++] = "KInstant";

    for(int ispec=0;ispec<nspecies;ispec++) {
        cuNames[cnt++] = amrex::Concatenate("rhoInstant_",ispec,2);
        cuNames[cnt++] = amrex::Concatenate("JxInstant_",ispec,2);
        cuNames[cnt++] = amrex::Concatenate("JyInstant_",ispec,2);
        cuNames[cnt++] = amrex::Concatenate("JzInstant_",ispec,2);
        cuNames[cnt++] = amrex::Concatenate("KInstant_",ispec,2);
    }
    MultiFab::Copy(mfcuplt, mfcuInst, 0, ncon*0, ncon, 0);

    // Mean values
    cuNames[cnt++] = "rhoMean";
    cuNames[cnt++] = "JxMean";
    cuNames[cnt++] = "JyMean";
    cuNames[cnt++] = "JzMean";
    cuNames[cnt++] = "KMean";

    for(int ispec=0;ispec<nspecies;ispec++) {
        cuNames[cnt++] = amrex::Concatenate("rhoMean_",ispec,2);
        cuNames[cnt++] = amrex::Concatenate("JxMean_",ispec,2);
        cuNames[cnt++] = amrex::Concatenate("JyMean_",ispec,2);
        cuNames[cnt++] = amrex::Concatenate("JzMean_",ispec,2);
        cuNames[cnt++] = amrex::Concatenate("KMean_",ispec,2);
    }
    MultiFab::Copy(mfcuplt, mfcuMeans, 0, ncon*1, ncon, 0);
    WriteSingleLevelPlotfile(pltcu, mfcuplt, cuNames, geom, time, step);

    //////////////////////////////////////
    // Primitive Means and Instants
    //////////////////////////////////////

    cnt = 0;
    // Instant values
    primNames[cnt++] = "nInstant";
    primNames[cnt++] = "rhoInstant";
    primNames[cnt++] = "uInstant";
    primNames[cnt++] = "vInstant";
    primNames[cnt++] = "wInstant";
    primNames[cnt++] = "GInstant";
    primNames[cnt++] = "TInstant";
    primNames[cnt++] = "PInstant";
    primNames[cnt++] = "EInstant";

    for(int ispec=0;ispec<nspecies;ispec++) {
        primNames[cnt++] = amrex::Concatenate("nInstant_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("rhoInstant_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("uInstant_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("vInstant_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("wInstant_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("GInstant_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("TInstant_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("PInstant_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("EInstant_",ispec,2);
    }
    MultiFab::Copy(mfprimplt, mfprimInst, 0, nprim*0, nprim, 0);

    // Mean values
    primNames[cnt++] = "nMean";
    primNames[cnt++] = "rhoMean";
    primNames[cnt++] = "uMean";
    primNames[cnt++] = "vMean";
    primNames[cnt++] = "wMean";
    primNames[cnt++] = "GMean";
    primNames[cnt++] = "TMean";
    primNames[cnt++] = "PMean";
    primNames[cnt++] = "EMean";

    for(int ispec=0;ispec<nspecies;ispec++) {
        primNames[cnt++] = amrex::Concatenate("nMean_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("rhoMean_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("uMean_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("vMean_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("wMean_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("GMean_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("TMean_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("PMean_",ispec,2);
        primNames[cnt++] = amrex::Concatenate("EMean_",ispec,2);
    }
    MultiFab::Copy(mfprimplt, mfprimMeans, 0, nprim*1, nprim, 0);
    WriteSingleLevelPlotfile(pltprim, mfprimplt, primNames, geom, time, step);

    //////////////////////////////////////
    // Variances
    //////////////////////////////////////

    cnt = 0;
    int istart = 0;

    // Conserved Variances
    varNames[cnt++] = "rhoVar";
    varNames[cnt++] = "JxVar";
    varNames[cnt++] = "JyVar";
    varNames[cnt++] = "JzVar";
    varNames[cnt++] = "KVar";

    for(int ispec=0;ispec<nspecies;ispec++) {
        varNames[cnt++] = amrex::Concatenate("rhoVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("JxVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("JyVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("JzVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("KVar_",ispec,2);
    }
    MultiFab::Copy(mfvarplt, mfcuVars, 0, istart, ncon, 0);
    istart += ncon;

    // Primitive Variances
    varNames[cnt++] = "nVar";
    varNames[cnt++] = "rhoVar";
    varNames[cnt++] = "uVar";
    varNames[cnt++] = "vVar";
    varNames[cnt++] = "wVar";
    varNames[cnt++] = "GVar";
    varNames[cnt++] = "TVar";
    varNames[cnt++] = "PVar";
    varNames[cnt++] = "EVar";

    for(int ispec=0;ispec<nspecies;ispec++) {
        varNames[cnt++] = amrex::Concatenate("nVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("rhoVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("uVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("vVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("wVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("GVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("TVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("PVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("EVar_",ispec,2);
    }
    MultiFab::Copy(mfvarplt, mfprimVars, 0, istart, nprim, 0);
    istart += nprim;

    varNames[cnt++] = "rho.Jx";
    varNames[cnt++] = "rho.Jy";
    varNames[cnt++] = "rho.Jz";
    varNames[cnt++] = "rho.K";
    varNames[cnt++] = "Jx.Jy";
    varNames[cnt++] = "Jx.Jz";
    varNames[cnt++] = "Jx.K";
    varNames[cnt++] = "Jy.Jz";
    varNames[cnt++] = "Jy.K";
    varNames[cnt++] = "Jz.K";
    varNames[cnt++] = "rho.G";
    varNames[cnt++] = "Jx.G";
    varNames[cnt++] = "Jy.G";
    varNames[cnt++] = "Jz.G";
    varNames[cnt++] = "K.G";
    varNames[cnt++] = "rho.u";
    varNames[cnt++] = "rho.v";
    varNames[cnt++] = "rho.w";
    varNames[cnt++] = "u.v";
    varNames[cnt++] = "u.w";
    varNames[cnt++] = "v.w";
    varNames[cnt++] = "rho.T";
    varNames[cnt++] = "u.T";
    varNames[cnt++] = "v.T";
    varNames[cnt++] = "w.T";

    MultiFab::Copy(mfvarplt, mfcoVars, 0, istart, ncovar, 0);
    WriteSingleLevelPlotfile(pltvar, mfvarplt, varNames, geom, time, step);

    if (plot_cross) {
        std::string file_prefix = "spatialCross1D_";
        WriteHorizontalAverage(mfspatialCorr1d,0,0,ncross,step,geom,file_prefix);
    }

}