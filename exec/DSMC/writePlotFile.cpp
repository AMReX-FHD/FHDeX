#include "LocalFunctions.H"
#include "common_functions.H"
#include "species.H"
#include "DsmcParticleContainer.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;
using namespace std;

void writePlotFile(const MultiFab& mfcuInst,
    const MultiFab& mfcuMeans,
    const MultiFab& mfcuVars,
    const MultiFab& mfprimInst,
    const MultiFab& mfprimMeans,
    const MultiFab& mfprimVars,
    const MultiFab& mfcoVars,
    const MultiFab& mfspatialCorr1d,
    FhdParticleContainer& particles,
    const Geometry& geom,
    Real time,
    const int ncross,
    int step) {
    BL_PROFILE_VAR("writePlotFile()",writePlotFile);

    int ncon    = (nspecies+1)*5;
    int nprim   = (nspecies+1)*10;
    int ncovar  = 25;
    int nvars   = ncovar + ncon + nprim + ncross; // covariances + prim. vars + cons. vars


    amrex::BoxArray ba = mfcuInst.boxArray();
    amrex::DistributionMapping dmap = mfcuInst.DistributionMap();

    std::string pltcu    = amrex::Concatenate("pltcu",step,12);
    std::string pltprim  = amrex::Concatenate("pltprim",step,12);
    std::string pltvar   = amrex::Concatenate("pltvar",step,12);
    std::string pltpart   = amrex::Concatenate("pltpart",step,12);

    amrex::MultiFab mfcuplt(ba, dmap, 2*ncon, 0);
    amrex::MultiFab mfprimplt(ba, dmap, 2*nprim, 0);
    amrex::MultiFab mfvarplt(ba, dmap, nvars, 0);

    amrex::MultiFab mfcrossav(ba, dmap, ncross, 0);

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
    primNames[cnt++] = "cInstant";

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
        primNames[cnt++] = amrex::Concatenate("cInstant_",ispec,2);
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
    primNames[cnt++] = "cMean";

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
        primNames[cnt++] = amrex::Concatenate("cMean_",ispec,2);
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
    varNames[cnt++] = "rhoVarP";
    varNames[cnt++] = "uVar";
    varNames[cnt++] = "vVar";
    varNames[cnt++] = "wVar";
    varNames[cnt++] = "GVar";
    varNames[cnt++] = "TVar";
    varNames[cnt++] = "PVar";
    varNames[cnt++] = "EVar";
    varNames[cnt++] = "cVar";

    for(int ispec=0;ispec<nspecies;ispec++) {
        varNames[cnt++] = amrex::Concatenate("nVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("rhoVarP_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("uVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("vVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("wVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("GVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("TVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("PVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("EVar_",ispec,2);
        varNames[cnt++] = amrex::Concatenate("cVar_",ispec,2);
    }
    MultiFab::Copy(mfvarplt, mfprimVars, 0, istart, nprim, 0);
    istart += nprim;

    //Covariances
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
    varNames[cnt++] = "rho0.Jx0";
    varNames[cnt++] = "rho1.Jx1";
    varNames[cnt++] = "rho0.Jx1";
    varNames[cnt++] = "rho0.u0";
    varNames[cnt++] = "rho1.u1";
    varNames[cnt++] = "rho.u";
    varNames[cnt++] = "rho.v";
    varNames[cnt++] = "rho0.u";
    varNames[cnt++] = "rho0.v";
    varNames[cnt++] = "rho0.Jx";
    varNames[cnt++] = "rho1.Jx";
    varNames[cnt++] = "rho.T";
    varNames[cnt++] = "u.T";
    varNames[cnt++] = "v.T";
    varNames[cnt++] = "w.T";

    MultiFab::Copy(mfvarplt, mfcoVars, 0, istart, ncovar, 0);
    istart += ncovar;

    //Cross correlations
    varNames[cnt++] = "rho*.rho";
    varNames[cnt++] = "K*.K";
    varNames[cnt++] = "Jx*.Jx";
    varNames[cnt++] = "Jy*.Jy";
    varNames[cnt++] = "Jz*.Jz";
    varNames[cnt++] = "Jx*.rho";
    varNames[cnt++] = "Jx*.rho_00";
    varNames[cnt++] = "Jx_00*.rho_00";
    varNames[cnt++] = "Jx*.rho_01";
    varNames[cnt++] = "Jx_01*.rho_01";
    varNames[cnt++] = "Jx_00*.rho_01";
    varNames[cnt++] = "Jx_01*.rho_00";
    varNames[cnt++] = "Jx_00*.rho";
    varNames[cnt++] = "Jx_01*.rho";
    varNames[cnt++] = "Jx*.Jy";
    varNames[cnt++] = "Jz*.Jy";
    varNames[cnt++] = "K*.Jy";
    varNames[cnt++] = "rho*.K";
    varNames[cnt++] = "Jx*.K";
    varNames[cnt++] = "Jy*.K";
    varNames[cnt++] = "Jz*.K";
    varNames[cnt++] = "G*.G";
    varNames[cnt++] = "G*.K";
    varNames[cnt++] = "K*.G";
    varNames[cnt++] = "rho*.G";
    varNames[cnt++] = "G*.rho";
    varNames[cnt++] = "Jx*.G";
    varNames[cnt++] = "G*.Jx";
    varNames[cnt++] = "Jy*.G";
    varNames[cnt++] = "G*.Jy";
    varNames[cnt++] = "Jz*.G";
    varNames[cnt++] = "G*.Jz";
    varNames[cnt++] = "K*.G";
    varNames[cnt++] = "G*.K";
    varNames[cnt++] = "T*.T";
    varNames[cnt++] = "T*.rho";
    varNames[cnt++] = "u*.rho";
    varNames[cnt++] = "T*.u";
//    for(int i=0;i<nspecies;i++)
//    {
//        for(int j=0;j<nspecies;j++)
//        {
//               std::string a = amrex::Concatenate("rho_",i,2);
//               std::string b = amrex::Concatenate("*.rho_",j,2);
//               varNames[cnt++] = a+b;
//        }
//    }

    varNames[cnt++] = "rho_00*.rho_00";
    //varNames[cnt++] = "rho_01*.rho_00";
    //varNames[cnt++] = "rho_00*.rho_01";
    varNames[cnt++] = "rho_01*.rho_01";
    varNames[cnt++] = "ux_00*.rho_00";
    varNames[cnt++] = "ux*.rho_00";
    varNames[cnt++] = "rho_01*.rho_00";
    varNames[cnt++] = "ux*.rho_01";
//    for(int i=0;i<nspecies;i++)
//    {
//        varNames[cnt++] = amrex::Concatenate("u*.rho_",i,2);
//
//    }

    //WriteHorizontalAverage(mfspatialCorr1d,mfcrossav,0,ncross);
    MultiFab::Copy(mfvarplt, mfspatialCorr1d, 0, istart, ncross, 0);

    WriteSingleLevelPlotfile(pltvar, mfvarplt, varNames, geom, time, step);


        // particle in cplt file
//    Vector<std::string> real_comp_names = FHD_realData::names();
//    Vector<std::string>  int_comp_names = FHD_intData::names();

//    Vector<int> write_real_comp(real_comp_names.size(),1);
//    Vector<int> write_int_comp(int_comp_names.size(),1);
//
//    particles.WritePlotFile(pltpart, "particles",
//                            write_real_comp, write_int_comp, real_comp_names, int_comp_names);

}