#include "AMReX_PlotFileUtil.H"
 
#include "common_functions.H"

#include "compressible_functions_stag.H"

void WritePlotFileStag(int step,
                       const amrex::Real time,
                       const amrex::Geometry& geom,
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
                       const amrex::MultiFab& surfcov,
                       const amrex::MultiFab& surfcovMeans,
                       const amrex::MultiFab& surfcovVars,
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

    if (nspec_surfcov>0) nplot += nspec_surfcov;

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

        if (nspec_surfcov>0) nplot += nspec_surfcov;
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

        if (nspec_surfcov>0) nplot += nspec_surfcov;
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

    // instantaneous
    // surfcov -- nspec_surfcov
    if (nspec_surfcov>0) {
        numvars = nspec_surfcov;
        amrex::MultiFab::Copy(plotfile,surfcov,0,cnt,numvars,0);
        cnt+=numvars;
    }

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

        // surfcov -- nspec_surfcov
        if (nspec_surfcov>0) {
            numvars = nspec_surfcov;
            amrex::MultiFab::Copy(plotfile,surfcovMeans,0,cnt,numvars,0);
            cnt+=numvars;
        }
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

        // surfcov -- nspec_surfcov
        if (nspec_surfcov>0) {
            numvars = nspec_surfcov;
            amrex::MultiFab::Copy(plotfile,surfcovVars,0,cnt,numvars,0);
            cnt+=numvars;
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

    if (nspec_surfcov>0) {
        x = "surfcov_";
        for (i=0; i<nspec_surfcov; i++) {
            varNames[cnt] = x;
            varNames[cnt++] += 48+i;
        }
    }

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

        if (nspec_surfcov>0) {
            x = "surfcovMean_";
            for (i=0; i<nspec_surfcov; i++) {
                varNames[cnt] = x;
                varNames[cnt++] += 48+i;
            }
        }
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

        if (nspec_surfcov>0) {
            x = "surfcovVar_";
            for (i=0; i<nspec_surfcov; i++) {
                varNames[cnt] = x;
                varNames[cnt++] += 48+i;
            }
        }
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

    AMREX_ASSERT(cnt==nplot);

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

void EvaluateWritePlotFileVelGrad(int step,
                                  const amrex::Real time,
                                  const amrex::Geometry& geom,
                                  const std::array<MultiFab, AMREX_SPACEDIM>& vel)
{
    BL_PROFILE_VAR("EvaluateWritePlotFileVelGrad()",EvaluateWritePlotFileVelGrad);

    // Evaluate velocity gradient components and divergence and vorticity
    MultiFab vel_grad;
    
    // Cell-Centered Velocity Gradient Stats (1,2,3 are directions)
    // 0: u_1,1
    // 1: u_2,2
    // 2: u_3,3
    // 3: u_1,2
    // 4: u_1,3
    // 5: u_2,3
    // 6: divergence = u_1,1 + u_2,2 + u_3,3
    // 7: vorticity = sqrt(wx + wy + wz)
    vel_grad.define(convert(vel[0].boxArray(),IntVect(AMREX_D_DECL(0,0,0))), vel[0].DistributionMap(), 8, 0);
    vel_grad.setVal(0.0);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(vel_grad,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.tilebox();
        
        const Array4<Real> & vgrad = vel_grad.array(mfi);
        
        AMREX_D_TERM(Array4<Real const> const& velx = vel[0].array(mfi);,
                     Array4<Real const> const& vely = vel[1].array(mfi);,
                     Array4<Real const> const& velz = vel[2].array(mfi););

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // u_1,1
            vgrad(i,j,k,0) = (velx(i+1,j,k) - velx(i,j,k))/dx[0];
            
            // u_2,2
            vgrad(i,j,k,1) = (vely(i,j+1,k) - vely(i,j,k))/dx[1];
            
            // u_3,3
            vgrad(i,j,k,2) = (velz(i,j,k+1) - velz(i,j,k))/dx[2];

            // divergence
            vgrad(i,j,k,6) = vgrad(i,j,k,0) + vgrad(i,j,k,1) + vgrad(i,j,k,2);

            // on edges: u_1,2 and u_2,1 and curl w1 = u_2,1 - u_1,2
            Real u12_mm = (velx(i,j,k) - velx(i,j-1,k))/dx[1];
            Real u21_mm = (vely(i,j,k) - vely(i-1,j,k))/dx[0];
            Real w1_mm  = u21_mm - u12_mm;
            Real u12_mp = (velx(i,j+1,k) - velx(i,j,k))/dx[1];
            Real u21_mp = (vely(i,j+1,k) - vely(i-1,j+1,k))/dx[0];
            Real w1_mp  = u21_mp - u12_mp;
            Real u12_pm = (velx(i+1,j,k) - velx(i+1,j-1,k))/dx[1];
            Real u21_pm = (vely(i+1,j,k) - vely(i,j,k))/dx[0];
            Real w1_pm  = u21_pm - u12_pm;
            Real u12_pp = (velx(i+1,j+1,k) - velx(i+1,j,k))/dx[1];
            Real u21_pp = (vely(i+1,j+1,k) - vely(i,j+1,k))/dx[0];
            Real w1_pp  = u21_pp - u12_pp;

            // u_1,2
            vgrad(i,j,k,3) = 0.25*(u12_mm + u12_mp + u12_pm + u12_pp);

            // on edges: u_1,3 and u_3,1 and curl w2 = u_1,3 - u_3,1
            Real u13_mm = (velx(i,j,k) - velx(i,j,k-1))/dx[2];
            Real u31_mm = (velz(i,j,k) - velz(i-1,j,k))/dx[0];
            Real w2_mm  = u13_mm - u31_mm;
            Real u13_mp = (velx(i,j,k+1) - velx(i,j,k))/dx[2];
            Real u31_mp = (velz(i,j,k+1) - velz(i-1,j,k+1))/dx[0];
            Real w2_mp  = u13_mp - u31_mp;
            Real u13_pm = (velx(i+1,j,k) - velx(i+1,j,k-1))/dx[2];
            Real u31_pm = (velz(i+1,j,k) - velz(i,j,k))/dx[0];
            Real w2_pm  = u13_pm - u31_pm;
            Real u13_pp = (velx(i+1,j,k+1) - velx(i+1,j,k))/dx[2];
            Real u31_pp = (velz(i+1,j,k+1) - velz(i,j,k+1))/dx[0];
            Real w2_pp  = u13_pp - u31_pp;

            // u_1,3
            vgrad(i,j,k,4) = 0.25*(u13_mm + u13_mp + u13_pm + u13_pp);

            // on edges: u_2,3 and u_3,2 and curl w2 = u_3,2 - u_2,3
            Real u23_mm = (vely(i,j,k) - vely(i,j,k-1))/dx[2];
            Real u32_mm = (velz(i,j,k) - velz(i,j-1,k))/dx[1];
            Real w3_mm  = u32_mm - u23_mm;
            Real u23_mp = (vely(i,j,k+1) - vely(i,j,k))/dx[2];
            Real u32_mp = (velz(i,j,k+1) - velz(i,j-1,k+1))/dx[1];
            Real w3_mp  = u32_mp - u23_mp;
            Real u23_pm = (vely(i,j+1,k) - vely(i,j+1,k-1))/dx[2];
            Real u32_pm = (velz(i,j+1,k) - velz(i,j,k))/dx[1];
            Real w3_pm  = u32_pm - u23_pm;
            Real u23_pp = (vely(i,j+1,k+1) - vely(i,j+1,k))/dx[2];
            Real u32_pp = (velz(i,j+1,k+1) - velz(i,j,k+1))/dx[1];
            Real w3_pp  = u32_pp - u23_pp;

            // u_2,3
            vgrad(i,j,k,5) = 0.25*(u23_mm + u23_mp + u23_pm + u23_pp);

            // vorticity magnitude: sqrt(w1*w1 + w2*w2 + w3*w3)
            vgrad(i,j,k,7) = sqrt(0.25*(w1_mm*w1_mm + w1_mp*w1_mp + w1_pm*w1_pm + w1_pp*w1_pp +
                                        w2_mm*w2_mm + w2_mp*w2_mp + w2_pm*w2_pm + w2_pp*w2_pp +
                                        w3_mm*w3_mm + w3_mp*w3_mp + w3_pm*w3_pm + w3_pp*w3_pp));
        });
    }

    // Write on a plotfile
    std::string plotfilename = amrex::Concatenate("vel_grad",step,9);
    amrex::Vector<std::string> varNames(8);
    varNames[0] = "ux_x";
    varNames[1] = "uy_y";
    varNames[2] = "uz_z";
    varNames[3] = "ux_y";
    varNames[4] = "ux_z";
    varNames[5] = "uy_z";
    varNames[6] = "div";
    varNames[7] = "vort";
    WriteSingleLevelPlotfile(plotfilename,vel_grad,varNames,geom,time,step);
}


