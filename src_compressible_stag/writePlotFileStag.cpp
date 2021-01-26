#include "AMReX_PlotFileUtil.H"
 
#include "common_functions.H"

#include "compressible_functions_stag.H"

void WritePlotFileStag(int step,
                       const amrex::Real time,
                       const amrex::Geometry geom,
                       const amrex::MultiFab& cu,
                       const amrex::MultiFab& cuMeans,
                       const amrex::MultiFab& cuVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                       const amrex::MultiFab& prim,
                       const amrex::MultiFab& primMeans,
                       const amrex::MultiFab& primVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                       const amrex::MultiFab& spatialCross,
                       const amrex::MultiFab& eta,
                       const amrex::MultiFab& kappa)
{
    BL_PROFILE_VAR("writePlotFile()",writePlotFile);
    
    int cnt, numvars, i = 0;

    // instantaneous values
    // 5 + nspecies (conserved) + 3 shifted momentum
    // 6 + 2*nspecies (primitive) + 3 shifted velocities
    // 2 (eta and kappa)
    int nplot = (5+nspecies+3) + (6+2*nspecies+3) + 2;

    if (plot_means == 1) {
        nplot += 11;
    }
    
    if (plot_vars == 1) {
        nplot += 10;
    }
   
    nplot += 6; //spatial correl

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

    // shifted momentum: jx_shifted, jy_shifted, jz_shifted
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        ShiftFaceToCC(cumom[d],0,plotfile,cnt,1);
        ++cnt;
    }
    

    // instantaneous values of primitive variables
    // rho, ux, uy, uz, temp, pres, Yk, Xk
    numvars = 6+2*nspecies;
    amrex::MultiFab::Copy(plotfile,prim,0,cnt,numvars,0);
    cnt+=numvars;

    // shifted momentum: velx_shifted, vely_shifted, velz_shifted
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        ShiftFaceToCC(vel[d],0,plotfile,cnt,1);
        ++cnt;
    }

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

    numvars = 6;
    amrex::MultiFab::Copy(plotfile,spatialCross,0,cnt,numvars,0);
    cnt+=numvars;

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
    varNames[cnt++] = "jx_avg_Instant"; 
    varNames[cnt++] = "jy_avg_Instant";
    varNames[cnt++] = "jz_avg_Instant";
    varNames[cnt++] = "rhoEInstant";
    x = "rhoYkInstant_";
    for (i=0; i<nspecies; i++) {
        varNames[cnt] = x;
        varNames[cnt++] += 48+i;
    }
    varNames[cnt++] = "jx_shifted_Instant";  
    varNames[cnt++] = "jy_shifted_Instant";  
    varNames[cnt++] = "jz_shifted_Instant";  

    varNames[cnt++] = "rhoInstant";
    varNames[cnt++] = "ux_avg_Instant";
    varNames[cnt++] = "uy_avg_Instant";
    varNames[cnt++] = "uz_avg_Instant";
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
    varNames[cnt++] = "ux_shifted_Instant";  
    varNames[cnt++] = "uy_shifted_Instant";  
    varNames[cnt++] = "uz_shifted_Instant";  

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

    varNames[cnt++] = "Energy-densityCross";
    varNames[cnt++] = "Energy-energyCross";
    varNames[cnt++] = "Momentum-densityCross";

    varNames[cnt++] = "Temperature-temperatureCross";
    varNames[cnt++] = "Temperature-densityCross";
    varNames[cnt++] = "Velocity-densityCross";

    varNames[cnt++] = "eta";
    varNames[cnt++] = "kappa";

    // write a plotfile
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);

}
