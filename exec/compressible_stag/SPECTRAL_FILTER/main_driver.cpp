#include "spectral_functions.H"
#include <AMReX_Vector.H>
#include <AMReX_MPMD.H>
#include "AMReX_ParmParse.H"

#include "chrono"

using namespace std::chrono;
using namespace amrex;

// argv contains the name of the inputs file entered at the command line
void main_driver(const char* argv)
{
    BL_PROFILE_VAR("main_driver()",main_driver);
    
    amrex::Vector<amrex::IntVect> nodal_flag_dir;
    amrex::IntVect                nodal_flag_x;
    amrex::IntVect                nodal_flag_y;
    amrex::IntVect                nodal_flag_z;
    nodal_flag_dir.resize(3);

    for (int i=0; i<3; ++i) {
        nodal_flag_x[i] = int(i==0);
        nodal_flag_y[i] = int(i==1);
        nodal_flag_z[i] = int(i==2);
        AMREX_D_TERM(nodal_flag_dir[0][i] = nodal_flag_x[i];,
                     nodal_flag_dir[1][i] = nodal_flag_y[i];,
                     nodal_flag_dir[2][i] = nodal_flag_z[i];);
    }

    // timer
    Real ts1 = ParallelDescriptor::second();
    
    std::string inputs_file = argv;

    ParmParse pp;
    amrex::Vector<amrex::Real> temp_real(3,0.);
    amrex::Vector<int>         temp_int (3,0 );

    amrex::Vector<int>   max_grid_size(3,1 );
    amrex::Vector<int>   n_cells(3,0 );
    amrex::Vector<Real>  prob_lo(3,0 );
    amrex::Vector<Real>  prob_hi(3,0 );
    
    if (pp.queryarr("n_cells",temp_int,0,3)) {
        for (int i=0; i<3; ++i) {
            n_cells[i] = temp_int[i];
        }
    }
    int npts = n_cells[0]*n_cells[1]*n_cells[2];
    if (pp.queryarr("prob_lo",temp_real,0,3)) {
        for (int i=0; i<3; ++i) {
            prob_lo[i] = temp_real[i];
        }
    }
    if (pp.queryarr("prob_hi",temp_real,0,3)) {
        for (int i=0; i<3; ++i) {
            prob_hi[i] = temp_real[i];
        }
    }
    pp.queryarr("max_grid_size",max_grid_size,0,3);

    int restart;
    pp.query("restart",restart);

    int nprimvars;
    pp.query("nprimvars",nprimvars);

    amrex::IntVect ngc;
    for (int i=0; i<3; ++i) {
        ngc[i] = 1;           // number of ghost cells
    }
    if (pp.queryarr("ngc",temp_int,0,3)) {
        for (int i=0; i<3; ++i) {
            ngc[i] = temp_int[i];
        }
    }

    amrex::Real kmin;
    pp.query("kmin",kmin);

    amrex::Real kmax;
    pp.query("kmax",kmax);
    
    std::array< MultiFab, 3 > vel;
    MultiFab prim;
    
    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    DistributionMapping dmap;

    IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
    IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
    Box domain(dom_lo, dom_hi);

    // This defines the physical box, [-1,1] in each direction.
    RealBox real_box({AMREX_D_DECL(prob_lo[0],prob_lo[1],prob_lo[2])},
                     {AMREX_D_DECL(prob_hi[0],prob_hi[1],prob_hi[2])});

    // This defines a Geometry object
    Vector<int> is_periodic(3,1);  // force to be periodic -- can change later
    geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());

    const Real* dx = geom.CellSize();
    const RealBox& realDomain = geom.ProbDomain();

    SpectralReadCheckPoint(geom, domain, prim, vel, ba, dmap, n_cells, nprimvars, max_grid_size, ngc, restart);
    
    MultiFab MFTurbScalar;
    MultiFab MFTurbVel;
    MultiFab vel_decomp_filter;
    MultiFab scalar_filter;
    MFTurbVel.define(ba, dmap, 3, 0);
    MFTurbScalar.define(ba, dmap, 1, 0);
    vel_decomp_filter.define(ba, dmap, 9, 0);
    scalar_filter.define(ba, dmap, 1, 0);
    vel_decomp_filter.setVal(0.0);
    scalar_filter.setVal(0.0);

    // Set BC: 1) fill boundary 2) physical
    for (int d=0; d<3; d++) {
        vel[d].FillBoundary(geom.periodicity());
    }
    prim.FillBoundary(geom.periodicity());

    for(int d=0; d<3; d++) {
        ShiftFaceToCC(vel[d], 0, MFTurbVel, d, 1);
    }
    MultiFab::Copy(MFTurbScalar, prim, 0, 0, 1, 0);
    
    SpectralVelDecomp(MFTurbVel, vel_decomp_filter, kmin, kmax, geom, n_cells);
    SpectralScalarDecomp(MFTurbScalar, scalar_filter, kmin, kmax, geom, n_cells);

    SpectralWritePlotFile(restart, kmin, kmax, geom, vel_decomp_filter, scalar_filter);

    // Turbulence Diagnostics
    Real u_rms, u_rms_s, u_rms_d, delta_u_rms;
    Real taylor_len, taylor_Re_eta;
    Real skew, skew_s, skew_d, kurt, kurt_s, kurt_d;
    {
      vel_decomp_filter.FillBoundary(geom.periodicity());
    
      Vector<Real> dProb(3);
      dProb[0] = 1.0/((n_cells[0]+1)*n_cells[1]*n_cells[2]);
      dProb[1] = 1.0/((n_cells[1]+1)*n_cells[2]*n_cells[0]);
      dProb[2] = 1.0/((n_cells[2]+1)*n_cells[0]*n_cells[1]);

      // Setup temp MultiFabs
      std::array< MultiFab, AMREX_SPACEDIM > gradU;
      std::array< MultiFab, AMREX_SPACEDIM > faceTemp;
      MultiFab sound_speed;
      MultiFab ccTemp;
      MultiFab ccTempA;
      AMREX_D_TERM(gradU[0].define(convert(prim.boxArray(),nodal_flag_x), prim.DistributionMap(), 6, 0);,
                   gradU[1].define(convert(prim.boxArray(),nodal_flag_y), prim.DistributionMap(), 6, 0);,
                   gradU[2].define(convert(prim.boxArray(),nodal_flag_z), prim.DistributionMap(), 6, 0););   
      AMREX_D_TERM(faceTemp[0].define(convert(prim.boxArray(),nodal_flag_x), prim.DistributionMap(), 1, 0);,
                   faceTemp[1].define(convert(prim.boxArray(),nodal_flag_y), prim.DistributionMap(), 1, 0);,
                   faceTemp[2].define(convert(prim.boxArray(),nodal_flag_z), prim.DistributionMap(), 1, 0););   
      sound_speed.define(prim.boxArray(),prim.DistributionMap(),1,0);
      ccTemp.define(prim.boxArray(),prim.DistributionMap(),1,0);
      ccTempA.define(prim.boxArray(),prim.DistributionMap(),1,0);

      // Setup temp variables
      Vector<Real> gradU2_temp(3);
      Vector<Real> gradU2(3);
      Vector<Real> gradU3(3);
      Vector<Real> gradU4(3);
      Vector<Real> gradU2_s(3);
      Vector<Real> gradU3_s(3);
      Vector<Real> gradU4_s(3);
      Vector<Real> gradU2_d(3);
      Vector<Real> gradU3_d(3);
      Vector<Real> gradU4_d(3);
      
      Vector<int> comps   {0,1,2};
      Vector<int> comps_s{3,4,5};
      Vector<int> comps_d{6,7,8};

      // turbulent kinetic energy (total)
      ccTemp.setVal(0.0);
      MultiFab::AddProduct(ccTemp,vel_decomp_filter,0,vel_decomp_filter,0,0,1,0); //uu
      MultiFab::AddProduct(ccTemp,vel_decomp_filter,1,vel_decomp_filter,1,0,1,0); //vv
      MultiFab::AddProduct(ccTemp,vel_decomp_filter,2,vel_decomp_filter,2,0,1,0); //ww
      u_rms = ccTemp.sum(0)/npts;
      u_rms = sqrt(u_rms/3.0);
      MultiFab::Multiply(ccTemp,prim,0,0,1,0); // rho*(uu+vv+ww)

      // turbulent kinetic energy (solenoidal)
      ccTemp.setVal(0.0);
      MultiFab::AddProduct(ccTemp,vel_decomp_filter,3,vel_decomp_filter,0,0,1,0); //uu
      MultiFab::AddProduct(ccTemp,vel_decomp_filter,4,vel_decomp_filter,1,0,1,0); //vv
      MultiFab::AddProduct(ccTemp,vel_decomp_filter,5,vel_decomp_filter,2,0,1,0); //ww
      u_rms_s = ccTemp.sum(0)/npts;
      u_rms_s = sqrt(u_rms_s/3.0);
      MultiFab::Multiply(ccTemp,prim,0,0,1,0); // rho*(uu+vv+ww)

      // turbulent kinetic energy (dilatational)
      ccTemp.setVal(0.0);
      MultiFab::AddProduct(ccTemp,vel_decomp_filter,6,vel_decomp_filter,0,0,1,0); //uu
      MultiFab::AddProduct(ccTemp,vel_decomp_filter,7,vel_decomp_filter,1,0,1,0); //vv
      MultiFab::AddProduct(ccTemp,vel_decomp_filter,8,vel_decomp_filter,2,0,1,0); //ww
      u_rms_d = ccTemp.sum(0)/npts;
      u_rms_d = sqrt(u_rms_d/3.0);
      MultiFab::Multiply(ccTemp,prim,0,0,1,0); // rho*(uu+vv+ww)
      
      // ratio of turbulent kinetic energies
      delta_u_rms  = u_rms_d/u_rms_s;

      // compute gradU = [du/dx dv/dy dw/dz] at cell-centers
      ComputeGrad(vel_decomp_filter,gradU,0,0,9,-1,geom,0);

      // Compute Velocity gradient moment sum
      // 2nd moment (total)
      FCMoments(gradU,comps,faceTemp,2,gradU2_temp);
      gradU2[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
      gradU2[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
      gradU2[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));
      ccTemp.setVal(0.0);
      ccTempA.setVal(0.0);
      ShiftFaceToCC(faceTemp[0],0,ccTempA,0,1);
      MultiFab::Add(ccTemp,ccTempA,0,0,1,0);
      ShiftFaceToCC(faceTemp[1],0,ccTempA,0,1);
      MultiFab::Add(ccTemp,ccTempA,0,0,1,0);
      ShiftFaceToCC(faceTemp[2],0,ccTempA,0,1);
      MultiFab::Add(ccTemp,ccTempA,0,0,1,0);
      Real avg_mom2 = ccTemp.sum(0)/npts;
      // 2nd moment (solenoidal)
      FCMoments(gradU,comps_s,faceTemp,2,gradU2_temp);
      gradU2_s[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
      gradU2_s[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
      gradU2_s[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));
      // 2nd moment (dilatational)
      FCMoments(gradU,comps_d,faceTemp,2,gradU2_temp);
      gradU2_d[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
      gradU2_d[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
      gradU2_d[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));

      // Taylor Mach
      //ComputeSoundSpeed(sound_speed,prim,2);
      //Real c_speed = sound_speed.sum(0)/npts;
      Real rho_avg = prim.sum(0)/npts;
      // Taylor Ma
      //taylor_Ma = sqrt(3.0)*u_rms/c_speed;
      // Taylor Microscale
      taylor_len = sqrt(3.0)*u_rms/sqrt(avg_mom2); // from Wang et al., JFM, 2012
      taylor_Re_eta = rho_avg*taylor_len*u_rms; // from from John, Donzis, Sreenivasan, PRL 2019

      // Compute Velocity gradient moment sum
      // 3rd moment (total)
      FCMoments(gradU,comps,faceTemp,3,gradU2_temp);
      gradU3[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
      gradU3[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
      gradU3[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));
      // 3rd moment (solenoidal)
      FCMoments(gradU,comps_s,faceTemp,3,gradU2_temp);
      gradU3_s[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
      gradU3_s[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
      gradU3_s[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));
      // 3rd moment (dilatational)
      FCMoments(gradU,comps_d,faceTemp,3,gradU2_temp);
      gradU3_d[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
      gradU3_d[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
      gradU3_d[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));

      // Compute Velocity gradient moment sum
      // 4th moment (total)
      FCMoments(gradU,comps,faceTemp,4,gradU2_temp);
      gradU4[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
      gradU4[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
      gradU4[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));
      // 4th moment (solenoidal)
      FCMoments(gradU,comps_s,faceTemp,4,gradU2_temp);
      gradU4_s[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
      gradU4_s[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
      gradU4_s[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));
      // 4th moment (dilatational)
      FCMoments(gradU,comps_d,faceTemp,4,gradU2_temp);
      gradU4_d[0] = dProb[0]*(faceTemp[0].sum_unique(0,false,geom.periodicity()));
      gradU4_d[1] = dProb[1]*(faceTemp[1].sum_unique(0,false,geom.periodicity()));
      gradU4_d[2] = dProb[2]*(faceTemp[2].sum_unique(0,false,geom.periodicity()));
              
      // Skewness
      // <\sum_i (du_i/dx_i)^3> / (\sum_i <(du_i/dx_i)^2>^1.5)
      skew   = (gradU3[0] + gradU3[1] + gradU3[2])/
               (pow(gradU2[0],1.5) + pow(gradU2[1],1.5) + pow(gradU2[2],1.5)); 
      skew_s = (gradU3_s[0] + gradU3_s[1] + gradU3_s[2])/
               (pow(gradU2_s[0],1.5) + pow(gradU2_s[1],1.5) + pow(gradU2_s[2],1.5)); 
      skew_d = (gradU3_d[0] + gradU3_d[1] + gradU3_d[2])/
               (pow(gradU2_d[0],1.5) + pow(gradU2_d[1],1.5) + pow(gradU2_d[2],1.5)); 
              
      // Kurtosis
      // <\sum_i (du_i/dx_i)^4> / (\sum_i <(du_i/dx_i)^2>^2)
      kurt   = (gradU4[0] + gradU4[1] + gradU4[2])/
               (pow(gradU2[0],2.0) + pow(gradU2[1],2.0) + pow(gradU2[2],2.0)); 
      kurt_s = (gradU4_s[0] + gradU4_s[1] + gradU4_s[2])/
               (pow(gradU2_s[0],2.0) + pow(gradU2_s[1],2.0) + pow(gradU2_s[2],2.0)); 
      kurt_d = (gradU4_d[0] + gradU4_d[1] + gradU4_d[2])/
               (pow(gradU2_d[0],2.0) + pow(gradU2_d[1],2.0) + pow(gradU2_d[2],2.0)); 
    }
    std::string turbfilename = "turbstats_";
    std::ostringstream os;
    os << std::setprecision(3) << kmin;
    turbfilename += os.str();;
    std::ostringstream oss;
    oss << std::setprecision(3) << kmax;
    turbfilename += oss.str();
    
    std::ofstream turboutfile;
    turboutfile.open(turbfilename, std::ios::app);
    turboutfile << "u_rms " << "u_rms_s " << "u_rms_d " << "delta_u_rms " 
                << "TaylorLen " << "TaylorRe*Eta "
                << "skew " << "skew_s " << "skew_d "
                << "kurt " << "kurt_s " << "kurt_d "
                << std::endl;
     turboutfile << u_rms << " ";
     turboutfile << u_rms_s << " ";
     turboutfile << u_rms_d << " ";
     turboutfile << delta_u_rms << " ";
     turboutfile << taylor_len << " ";
     turboutfile << taylor_Re_eta << " ";
     turboutfile << skew << " ";
     turboutfile << skew_s << " ";
     turboutfile << skew_d << " ";
     turboutfile << kurt << " ";
     turboutfile << kurt_s << " ";
     turboutfile << kurt_d << " ";
     turboutfile << std::endl;

    // timer
    Real ts2 = ParallelDescriptor::second() - ts1;
    ParallelDescriptor::ReduceRealMax(ts2,  ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "Time (spectral filtering) " << ts2 << " seconds\n";
    
    // MultiFab memory usage
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    amrex::Long min_fab_megabytes  = amrex::TotalBytesAllocatedInFabsHWM()/1048576;
    amrex::Long max_fab_megabytes  = min_fab_megabytes;

    ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
    ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

    amrex::Print() << "High-water FAB megabyte spread across MPI nodes: ["
                   << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";

    min_fab_megabytes  = amrex::TotalBytesAllocatedInFabs()/1048576;
    max_fab_megabytes  = min_fab_megabytes;

    ParallelDescriptor::ReduceLongMin(min_fab_megabytes, IOProc);
    ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);

    amrex::Print() << "Curent     FAB megabyte spread across MPI nodes: ["
                   << min_fab_megabytes << " ... " << max_fab_megabytes << "]\n";
    
    turboutfile.close();
}
