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
                       const amrex::MultiFab& coVars,
                       const amrex::MultiFab& eta,
                       const amrex::MultiFab& kappa)
{
    BL_PROFILE_VAR("writePlotFileStag()",writePlotFileStag);
    
    int cnt, numvars, i = 0;

    int nplot = 0; 

    // instantaneous values
    // cu: [rho, jx, jy, jz, rhoE, rhoYk] -- nvars
    // shifted [jx, jy, jz] -- 3
    // prim: [vx, vy, vz, T, p, Yk, Xk] -- 5 + 2*nspecies
    // shifted [vx, vy, vz] -- 3
    // eta, kappa -- 2
    nplot += nvars;
    nplot += 3;
    nplot += 5 + 2*nspecies;
    nplot += 3;
    nplot += 2;

    // mean values
    // cu: [rho, jx, jy, jz, rhoE, rhoYk] -- nvars
    // shifted [jx, jy, jz] -- 3
    // prim: [vx, vy, vz, T, p, Yk, Xk] -- 5 + 2*nspecies
    // shifted [vx, vy, vz] -- 3
    // center of mass [ux, uy, uz] -- 3
    if (plot_means == 1) {
        nplot += nvars;
        nplot += 3;
        nplot += 5 + 2*nspecies;
        nplot += 3;
        nplot += 3;
    }
    
    // variances
    // cu: [rho, jx, jy, jz, rhoE, rhoYk] -- nvars
    // shifted [jx, jy, jz] -- 3
    // prim: [vx, vy, vz, Tdirect, Yk, gvar, kgcross, krcorss, rgcross, T] -- 9 + nspecies
    // shifted: [vx, vy, vz] -- 3
    if (plot_vars == 1) {
        nplot += nvars;
        nplot += 3;
        nplot += 9 + nspecies;
        nplot += 3;
    }
   
    // co-variances -- see the list in main_driver
    if (plot_covars == 1) {
        nplot += 26;
    }
   
    amrex::BoxArray ba = cuMeans.boxArray();
    amrex::DistributionMapping dmap = cuMeans.DistributionMap();

    amrex::MultiFab plotfile(ba, dmap, nplot, 0);

    std::string x;
    std::string plotfilename = amrex::Concatenate(plot_base_name,step,9);
    amrex::Vector<std::string> varNames(nplot);

    // Load into plotfile MF
    cnt = 0;

    // instantaneous 
    // cu: [rho, jx, jy, jz, rhoE, rhoYk] -- nvars
    numvars = nvars;
    amrex::MultiFab::Copy(plotfile,cu,0,cnt,numvars,0);
    cnt+=numvars;

    // instantaneous 
    // shifted [jx, jy, jz] -- 3
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        ShiftFaceToCC(cumom[d],0,plotfile,cnt,1);
        ++cnt;
    }

    // instantaneous 
    // prim: [vx, vy, vz, T, p, Yk, Xk] -- 5 + 2*nspecies
    numvars = 5+2*nspecies;
    amrex::MultiFab::Copy(plotfile,prim,1,cnt,numvars,0);
    cnt+=numvars;

    // instantaneous 
    // shifted [vx, vy, vz] -- 3
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        ShiftFaceToCC(vel[d],0,plotfile,cnt,1);
        ++cnt;
    }

    // instantaneous 
    // eta -- 1
    numvars = 1;
    amrex::MultiFab::Copy(plotfile,eta,0,cnt,numvars,0);
    cnt+=numvars;

    // instantaneous 
    // kappa -- 1
    numvars = 1;
    amrex::MultiFab::Copy(plotfile,kappa,0,cnt,numvars,0);
    cnt+=numvars;

    if (plot_means == 1) {
    
        // cu: [rho, jx, jy, jz, rhoE, rhoYk] -- nvars
        numvars = nvars;
        amrex::MultiFab::Copy(plotfile,cuMeans,0,cnt,numvars,0);
        cnt+=numvars;

        // shifted [jx, jy, jz] -- 3
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            ShiftFaceToCC(cumomMeans[d],0,plotfile,cnt,1);
            ++cnt;
        }
    
        // prim: [vx, vy, vz, T, p, Yk, Xk] -- 5 + nspecies
        numvars = 5 + 2*nspecies;
        amrex::MultiFab::Copy(plotfile,primMeans,1,cnt,numvars,0);
        cnt+=numvars;

        // shifted [vx, vy, vz] -- 3
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            ShiftFaceToCC(velMeans[d],0,plotfile,cnt,1);
            ++cnt;
        }

        // prim: center of mass velocity [ux, uy, uz] -- 3
        numvars = 3;
        amrex::MultiFab::Copy(plotfile,primMeans,nprimvars,cnt,numvars,0);
        cnt+=numvars;
    }

    if (plot_vars == 1) {
    
        // cu: [rho, jx, jy, jz, rhoE, rhoYk] -- nvars
        numvars = nvars;
        amrex::MultiFab::Copy(plotfile,cuVars,0,cnt,numvars,0);
        cnt+=numvars;

        // shifted [jx, jy, jz] -- 3
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            ShiftFaceToCC(cumomVars[d],0,plotfile,cnt,1);
            ++cnt;
        }
    
        // prim: [vx, vy, vz, T] -- 4 
        numvars = 4;
        amrex::MultiFab::Copy(plotfile,primVars,1,cnt,numvars,0);
        cnt+=numvars;

        // prim: Yk -- nspecies
        numvars = nspecies;
        amrex::MultiFab::Copy(plotfile,primVars,6,cnt,numvars,0);
        cnt+=numvars;

        // <delg delg> 
        amrex::MultiFab::Copy(plotfile,primVars,nprimvars+0,cnt,1,0);
        ++cnt;
        // <delg delenergy>
        amrex::MultiFab::Copy(plotfile,primVars,nprimvars+1,cnt,1,0);
        ++cnt;
        // <delrho delenergy>
        amrex::MultiFab::Copy(plotfile,primVars,nprimvars+2,cnt,1,0);
        ++cnt;
        // <delrho delg>
        amrex::MultiFab::Copy(plotfile,primVars,nprimvars+3,cnt,1,0);
        ++cnt;
        // <delT delT> -- direct computation
        amrex::MultiFab::Copy(plotfile,primVars,nprimvars+4,cnt,1,0);
        ++cnt;

        // shifted: [vx, vy, vz] -- 3
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            ShiftFaceToCC(velVars[d],0,plotfile,cnt,1);
            ++cnt;
        }
    }

    if (plot_covars == 1) {
        // co-variances -- see the list in main_driver
        numvars = 26;
        amrex::MultiFab::Copy(plotfile,coVars,0,cnt,numvars,0);
        cnt+=numvars;
    }

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

    varNames[cnt++] = "eta";
    varNames[cnt++] = "kappa";

    if (plot_means == 1) {
        varNames[cnt++] = "rhoMean";
        varNames[cnt++] = "jxMeanCC";
        varNames[cnt++] = "jyMeanCC";
        varNames[cnt++] = "jzMeanCC";
        varNames[cnt++] = "rhoEMean";
        x = "rhoYkMean_";
        for (i=0; i<nspecies; i++) {
            varNames[cnt] = x;
            varNames[cnt++] += 48+i;
        }

        varNames[cnt++] = "jxMeanFACE";
        varNames[cnt++] = "jyMeanFACE";
        varNames[cnt++] = "jzMeanFACE";

        varNames[cnt++] = "uxMeanCC";
        varNames[cnt++] = "uyMeanCC";
        varNames[cnt++] = "uzMeanCC";
        varNames[cnt++] = "tMean";
        varNames[cnt++] = "pMean";

        x = "YkMean_";
        for (i=0; i<nspecies; i++) {
            varNames[cnt] = x;
            varNames[cnt++] += 48+i;
        }

        x = "XkMean_";
        for (i=0; i<nspecies; i++) {
            varNames[cnt] = x;
            varNames[cnt++] += 48+i;
        }

        varNames[cnt++] = "uxMeanFACE";
        varNames[cnt++] = "uyMeanFACE";
        varNames[cnt++] = "uzMeanFACE";

        varNames[cnt++] = "uxMeanCOM";
        varNames[cnt++] = "uyMeanCOM";
        varNames[cnt++] = "uzMeanCOM";
    }

    if (plot_vars == 1) {
        varNames[cnt++] = "rhoVar";
        varNames[cnt++] = "jxVarCC";
        varNames[cnt++] = "jyVarCC";
        varNames[cnt++] = "jzVarCC";
        varNames[cnt++] = "rhoEVar";
        x = "rhoYkVar_";
        for (i=0; i<nspecies; i++) {
            varNames[cnt] = x;
            varNames[cnt++] += 48+i;
        }

        varNames[cnt++] = "jxVarFACE";
        varNames[cnt++] = "jyVarFACE";
        varNames[cnt++] = "jzVarFACE";

        varNames[cnt++] = "uxVarCC";
        varNames[cnt++] = "uyVarCC";
        varNames[cnt++] = "uzVarCC";
        varNames[cnt++] = "TVarDirect";

        x = "YkVar_";
        for (i=0; i<nspecies; i++) {
            varNames[cnt] = x;
            varNames[cnt++] += 48+i;
        }

        varNames[cnt++] = "g-g";
        varNames[cnt++] = "g-energy";
        varNames[cnt++] = "rho-energy";
        varNames[cnt++] = "rho-g";
        varNames[cnt++] = "TVar";

        varNames[cnt++] = "uxVarFACE";
        varNames[cnt++] = "uyVarFACE";
        varNames[cnt++] = "uzVarFACE";
    }

    if (plot_covars == 1) {
        varNames[cnt++] = "rho-jx";
        varNames[cnt++] = "rho-jy";
        varNames[cnt++] = "rho-jz";
        varNames[cnt++] = "jx-jy";
        varNames[cnt++] = "jy-jz";
        varNames[cnt++] = "jx-jz";
        varNames[cnt++] = "rho-rhoE";
        varNames[cnt++] = "rhoE-jx";
        varNames[cnt++] = "rhoE-jy";
        varNames[cnt++] = "rhoE-jz";
        varNames[cnt++] = "rhoL-rhoH";
        varNames[cnt++] = "rho-vx";
        varNames[cnt++] = "rho-vy";
        varNames[cnt++] = "rho-vz";
        varNames[cnt++] = "vx-vy";
        varNames[cnt++] = "vy-vz";
        varNames[cnt++] = "vx-vz";
        varNames[cnt++] = "rho-T";
        varNames[cnt++] = "vx-T";
        varNames[cnt++] = "vy-T";
        varNames[cnt++] = "vz-T";
        varNames[cnt++] = "YkL-YkH";
        varNames[cnt++] = "YkL-vx";
        varNames[cnt++] = "YkH-vx";
        varNames[cnt++] = "rhoYkL-vx";
        varNames[cnt++] = "rhoYkH-vx";
    }

    // write a plotfile
    // timer
    Real t1 = ParallelDescriptor::second();
    
    WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);
    
    Real t2 = ParallelDescriptor::second() - t1;
    ParallelDescriptor::ReduceRealMax(t2);
    amrex::Print() << "Time spent writing plotfile " << t2 << std::endl;
}

void WriteSpatialCross3D(const Vector<Real>& spatialCross, int step, const Geometry& geom, const int ncross) 
{
    if (ParallelDescriptor::IOProcessor()) {

        // write out spatial correlation
        std::string filename = amrex::Concatenate("spatialCross",step,9);
        std::ofstream outfile;
        outfile.open(filename);

        // cell size
        Real h = geom.CellSize(0);
    
        for (auto i=0; i<n_cells[0]; ++i) {
            outfile << prob_lo[0] + (i+0.5)*h << " "; 
            for (auto n=0; n<ncross; ++n) {
                outfile << std::setprecision(16) << spatialCross[i*ncross+n] << " ";
            }
            outfile << std::endl;
        }
        outfile.close();
    }
}

void WriteSpatialCross1D(const amrex::MultiFab& spatialCross, int step, const Geometry& geom, const int ncross) 
{
    if (all_correl == 0) { // single spatial correlation file
        std::string file_prefix = "spatialCross1D_";
        WriteHorizontalAverage(spatialCross,0,0,ncross,step,geom,file_prefix);
    }
    else { // five spatial correlation files
        for (int i=0; i<5; ++i) {
            std::string file_prefix = std::to_string(i) + "_spatialCross1D_";
            WriteHorizontalAverage(spatialCross,0,i*ncross,ncross,step,geom,file_prefix);
        }
    }
}

void WritePlotFilesSF_2D(const amrex::MultiFab& mag, const amrex::MultiFab& realimag, const amrex::Geometry& geom,
                         const int step, const Real time, const amrex::Vector< std::string >& names, std::string plotfile_base) {

      // Magnitude of the Structure Factor
      std::string name = plotfile_base;
      name += "_mag";
      const std::string plotfilename1 = amrex::Concatenate(name,step,9);
      
      Real dx0 = geom.CellSize(0);
      Real dx1 = geom.CellSize(1);
      Real dx2 = geom.CellSize(2);
      Real pi = 3.1415926535897932;
      Box domain = geom.Domain();
      RealBox real_box({AMREX_D_DECL(-pi/dx0,-pi/dx1,-pi/dx2)},
                       {AMREX_D_DECL( pi/dx0, pi/dx1, pi/dx2)});
      Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
      for (int i=0; i<AMREX_SPACEDIM; ++i) {
          is_periodic[i] = geom.isPeriodic(i);
      }
      Geometry geom2;
      geom2.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());

      Vector<std::string> varNames;
      varNames.resize(names.size());
      for (int n=0; n<names.size(); n++) {
          varNames[n] = names[n];
      }

      WriteSingleLevelPlotfile(plotfilename1,mag,varNames,geom2,time,step);

      // Components of the Structure Factor
      name = plotfile_base;
      name += "_real_imag";
      const std::string plotfilename2 = amrex::Concatenate(name,step,9);

      varNames.resize(2*names.size());
      int cnt = 0; // keep a counter for plotfile variables
      for (int n=0; n<names.size(); n++) {
          varNames[cnt] = names[cnt];
          varNames[cnt] += "_real";
          cnt++;
      }

      int index = 0;
      for (int n=0; n<names.size(); n++) {
          varNames[cnt] = names[index];
          varNames[cnt] += "_imag";
          index++;
          cnt++;
      }

      WriteSingleLevelPlotfile(plotfilename2,realimag,varNames,geom2,time,step);
}
