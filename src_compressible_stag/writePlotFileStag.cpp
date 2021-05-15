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
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& cumomVars,
                       const amrex::MultiFab& prim,
                       const amrex::MultiFab& primMeans,
                       const amrex::MultiFab& primVars,
                       const std::array<MultiFab, AMREX_SPACEDIM>& vel,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velMeans,
                       const std::array<MultiFab, AMREX_SPACEDIM>& velVars,
//                       const amrex::MultiFab& spatialCross,
                       const amrex::MultiFab& eta,
                       const amrex::MultiFab& kappa)
{
    BL_PROFILE_VAR("writePlotFileStag()",writePlotFileStag);
    
    int cnt, numvars, i = 0;

    // instantaneous values
    // 5 + nspecies (conserved) + 3 shifted momentum
    // 6 + 2*nspecies (primitive) + 3 shifted velocities
    // 2 (eta and kappa)
    int nplot = (5+nspecies+3) + (6+2*nspecies+3) + 2;

    if (plot_means == 1) {
        nplot += 17;
    }
    
    if (plot_vars == 1) {
        nplot += 16;
    }
   
    //nplot += 6; //spatial correl

    amrex::BoxArray ba = cuMeans.boxArray();
    amrex::DistributionMapping dmap = cuMeans.DistributionMap();

    amrex::MultiFab plotfile(ba, dmap, nplot, 0);

    std::string x;
    std::string plotfilename = amrex::Concatenate(plot_base_name,step,9);
    amrex::Vector<std::string> varNames(nplot);

    // Load into plotfile MF
    
    cnt = 0;

    // instantaneous values of conserved variables
    // rho, jx (avgd), jy (avgd), jz (avgd), e, rhoYk
    numvars = 5+nspecies;
    amrex::MultiFab::Copy(plotfile,cu,0,cnt,numvars,0);
    cnt+=numvars;

    // shifted momentum: 
    // jx , jy , jz
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        ShiftFaceToCC(cumom[d],0,plotfile,cnt,1);
        ++cnt;
    }

    // instantaneous values of primitive variables
    // rho, ux (avgd), uy (avgd), uz (avgd), temp, pres, Yk, Xk
    numvars = 6+2*nspecies;
    amrex::MultiFab::Copy(plotfile,prim,0,cnt,numvars,0);
    cnt+=numvars;

    // shifted velocities: 
    // velx, vely, velz
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        ShiftFaceToCC(vel[d],0,plotfile,cnt,1);
        ++cnt;
    }

    if (plot_means == 1) {
    
        // mean values of conserved variables
        // rho, jx (avgd), jy (avgd), jz (avgd), e
        numvars = 5;
        amrex::MultiFab::Copy(plotfile,cuMeans,0,cnt,numvars,0);
        cnt+=numvars;

        // mean shifted momentum: 
        // jx , jy , jz
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            ShiftFaceToCC(cumomMeans[d],0,plotfile,cnt,1);
            ++cnt;
        }
    
        // mean values of primitive variables
        // rho, ux (avgd), uy (avgd), uz (avgd), temp, pres
        numvars = 6;
        amrex::MultiFab::Copy(plotfile,primMeans,0,cnt,numvars,0);
        cnt+=numvars;

        // mean shifted velocities: 
        // velx, vely, velz
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            ShiftFaceToCC(velMeans[d],0,plotfile,cnt,1);
            ++cnt;
        }
    }

    if (plot_vars == 1) {
    
        // variance of conserved variables
        // rho, jx (avgd), jy (avgd), jz (avgd), e
        numvars = 5;
        amrex::MultiFab::Copy(plotfile,cuVars,0,cnt,numvars,0);
        cnt+=numvars;

        // variance of shifted momentum: 
        // jx , jy , jz
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            ShiftFaceToCC(cumomVars[d],0,plotfile,cnt,1);
            ++cnt;
        }
    
        // variances of primitive variables
        // rho, ux, uy, uz, temp
        numvars = 5;
        amrex::MultiFab::Copy(plotfile,primVars,0,cnt,numvars,0);
        cnt+=numvars;

        // variance of shifted velocities: 
        // velx, vely, velz
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            ShiftFaceToCC(velVars[d],0,plotfile,cnt,1);
            ++cnt;
        }
    }

    //numvars = 6;
    //amrex::MultiFab::Copy(plotfile,spatialCross,0,cnt,numvars,0);
    //cnt+=numvars;

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
    varNames[cnt++] = "jxInstantCC"; 
    varNames[cnt++] = "jyInstantCC";
    varNames[cnt++] = "jzInstantCC";
    varNames[cnt++] = "rhoEInstant";
    x = "rhoYkInstant_";
    for (i=0; i<nspecies; i++) {
        varNames[cnt] = x;
        varNames[cnt++] += 48+i;
    }
    varNames[cnt++] = "jxInstantFACE";  
    varNames[cnt++] = "jyInstantFACE";  
    varNames[cnt++] = "jzInstantFACE";  

    varNames[cnt++] = "rhoInstant";
    varNames[cnt++] = "uxInstantCC";
    varNames[cnt++] = "uyInstantCC";
    varNames[cnt++] = "uzInstantCC";
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
    varNames[cnt++] = "uxInstantFACE";  
    varNames[cnt++] = "uyInstantFACE";  
    varNames[cnt++] = "uzInstantFACE";  

    if (plot_means == 1) {
        varNames[cnt++] = "rhoMean";
        varNames[cnt++] = "jxMeanCC";
        varNames[cnt++] = "jyMeanCC";
        varNames[cnt++] = "jzMeanCC";
        varNames[cnt++] = "rhoEMean";
        varNames[cnt++] = "jxMeanFACE";
        varNames[cnt++] = "jyMeanFACE";
        varNames[cnt++] = "jzMeanFACE";

        varNames[cnt++] = "rhoMean";
        varNames[cnt++] = "uxMeanCC";
        varNames[cnt++] = "uyMeanCC";
        varNames[cnt++] = "uzMeanCC";
        varNames[cnt++] = "tMean";
        varNames[cnt++] = "pMean";
        varNames[cnt++] = "uxMeanFACE";
        varNames[cnt++] = "uyMeanFACE";
        varNames[cnt++] = "uzMeanFACE";
    }

    if (plot_vars == 1) {
        varNames[cnt++] = "rhoVar";
        varNames[cnt++] = "jxVarCC";
        varNames[cnt++] = "jyVarCC";
        varNames[cnt++] = "jzVarCC";
        varNames[cnt++] = "rhoEVar";
        varNames[cnt++] = "jxVarFACE";
        varNames[cnt++] = "jyVarFACE";
        varNames[cnt++] = "jzVarFACE";

        varNames[cnt++] = "rhoVar";
        varNames[cnt++] = "uxVarCC";
        varNames[cnt++] = "uyVarCC";
        varNames[cnt++] = "uzVarCC";
        varNames[cnt++] = "tVar";
        varNames[cnt++] = "uxVarFACE";
        varNames[cnt++] = "uyVarFACE";
        varNames[cnt++] = "uzVarFACE";
    }

    //varNames[cnt++] = "Energy-densityCross";
    //varNames[cnt++] = "Energy-energyCross";
    //varNames[cnt++] = "Momentum-densityCross";

    //varNames[cnt++] = "Temperature-temperatureCross";
    //varNames[cnt++] = "Temperature-densityCross";
    //varNames[cnt++] = "Velocity-densityCross";

    varNames[cnt++] = "eta";
    varNames[cnt++] = "kappa";

    // write a plotfile
    // timer
    Real t1 = ParallelDescriptor::second();
    
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);
    
    Real t2 = ParallelDescriptor::second() - t1;
    ParallelDescriptor::ReduceRealMax(t2);
    amrex::Print() << "Time spent writing plotfile " << t2 << std::endl;
}
