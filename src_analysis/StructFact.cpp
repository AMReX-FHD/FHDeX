
#include "common_functions.H"
#include "StructFact.H"

#include <AMReX_MultiFabUtil.H>
#include "AMReX_PlotFileUtil.H"
#include "AMReX_BoxArray.H"

#ifdef AMREX_USE_CUDA
std::string cufftErrorToString (const cufftResult& err)
{
    switch (err) {
    case CUFFT_SUCCESS:  return "CUFFT_SUCCESS";
    case CUFFT_INVALID_PLAN: return "CUFFT_INVALID_PLAN";
    case CUFFT_ALLOC_FAILED: return "CUFFT_ALLOC_FAILED";
    case CUFFT_INVALID_TYPE: return "CUFFT_INVALID_TYPE";
    case CUFFT_INVALID_VALUE: return "CUFFT_INVALID_VALUE";
    case CUFFT_INTERNAL_ERROR: return "CUFFT_INTERNAL_ERROR";
    case CUFFT_EXEC_FAILED: return "CUFFT_EXEC_FAILED";
    case CUFFT_SETUP_FAILED: return "CUFFT_SETUP_FAILED";
    case CUFFT_INVALID_SIZE: return "CUFFT_INVALID_SIZE";
    case CUFFT_UNALIGNED_DATA: return "CUFFT_UNALIGNED_DATA";
    default: return std::to_string(err) + " (unknown error code)";
    }
}
#endif

StructFact::StructFact()
{}

// var_names contains the names of all variables under consideration
// this constructor computes the covariances of all possible pairs of variables
// var_scaling must be sized to match the total number of pairs of variables
StructFact::StructFact(const BoxArray& ba_in, const DistributionMapping& dmap_in,
		       const Vector< std::string >& var_names,
		       const Vector< Real >& var_scaling_in,
		       const int& verbosity_in) {

    this->define(ba_in,dmap_in,var_names,var_scaling_in,verbosity_in);

}

// var_names contains the names of all variables under consideration
// this constructor compute the covariances of the pairs of variables defined in s_pairA/B_in
// var_scaling must be sized to match the total number of pairs of variables
StructFact::StructFact(const BoxArray& ba_in, const DistributionMapping& dmap_in,
		       const Vector< std::string >& var_names,
		       const Vector< Real >& var_scaling_in,
		       const Vector< int >& s_pairA_in,
		       const Vector< int >& s_pairB_in,
		       const int& verbosity_in) {

    this->define(ba_in,dmap_in,var_names,var_scaling_in,s_pairA_in,s_pairB_in,verbosity_in);

}


// this builds a list of all possible pairs of variables and calls define()
void StructFact::define(const BoxArray& ba_in, const DistributionMapping& dmap_in,
                        const Vector< std::string >& var_names,
                        const Vector< Real >& var_scaling_in,
                        const int& verbosity_in) {

    NVAR = var_names.size();

    Vector<int> s_pairA(NVAR*(NVAR+1)/2);
    Vector<int> s_pairB(NVAR*(NVAR+1)/2);

    int counter=0;
    for (int i=0; i<NVAR; ++i) {
        for (int j=i; j<NVAR; ++j) {
            s_pairA[counter] = j;
            s_pairB[counter] = i;
            ++counter;
        }
    }      

    define(ba_in, dmap_in, var_names, var_scaling_in, s_pairA, s_pairB, verbosity_in);

}

void StructFact::define(const BoxArray& ba_in, const DistributionMapping& dmap_in,
                        const Vector< std::string >& var_names,
                        const Vector< Real >& var_scaling_in,
                        const Vector< int >& s_pairA_in,
                        const Vector< int >& s_pairB_in,
                        const int& verbosity_in) {

  BL_PROFILE_VAR("StructFact::define()",StructFactDefine);

  verbosity = verbosity_in;
  
  if (s_pairA_in.size() != s_pairB_in.size())
        amrex::Error("StructFact::define() - Must have an equal number of components");

  NVAR = var_names.size();

  //////////////////////////////////////////////////////
  // Reorder pair selecting vectors & error checking

  NCOV = s_pairA_in.size();

  if ( NCOV != var_scaling_in.size() )
      amrex::Error("StructFact::define() - NCOV != var_scaling_in.size()");

  scaling.resize(NCOV);
  for (int n=0; n<NCOV; n++) {
      scaling[n] = 1.0/var_scaling_in[n];
  }
  
  s_pairA.resize(NCOV);
  s_pairB.resize(NCOV);
  
  // Set vectors identifying covariance pairs
  for (int n=0; n<NCOV; n++) {
    s_pairA[n] = s_pairA_in[n];
    s_pairB[n] = s_pairB_in[n];
  }
  
  // Create vector of unique variable indices to select which to take the FFT
  NVARU = 2*NCOV;    // temporary before selecting unique variables
  amrex::Vector< int > varu_temp (NVARU);
  
  // "Stack" on vectors A and B
  int indx = 0;
  for (int n=0; n<NCOV; n++) {
    varu_temp[indx] = s_pairA_in[n];
    indx++;
  }
  for (int n=0; n<NCOV; n++) {
    varu_temp[indx] = s_pairB_in[n];
    indx++;
  }

  // Reorder: ascending order using "bubble sort"
  int loop = 1;
  int tot_iters = 0;
  int x_temp;
  while (loop == 1) {
    loop = 0;
    for (int n=1; n<NVARU; n++) {
      if (varu_temp[n] < varu_temp[n-1]) {

	x_temp = varu_temp[n];
        varu_temp[n  ] = varu_temp[n-1];
        varu_temp[n-1] = x_temp;

	loop = 1;
      }
    }
    tot_iters++;
    if (tot_iters > 2*NVARU) {
      loop = 0;
      amrex::Error("StructFact::StructFact() - Bubble sort failed to converge");
    }
  }

  /*
  for (int n=0; n<NVARU; n++) {
    Print() << "HACK 1: vector (" << n << ") = " << varu_temp[n] << std::endl;
  }
  */

  // Identify number of repeats
  int N_dup = 0;
  for (int n=1; n<NVARU; n++) {
    if (varu_temp[n] == varu_temp[n-1]) {
      N_dup = N_dup+1;
    }
  }
  
  // Resize based on number of unique variables
  int N_u = NVARU - N_dup;
  var_u.resize(N_u);
  
  // Only select unique pairs
  indx = 0;
  var_u[indx] = varu_temp[0];
  for (int n=1; n<NVARU; n++) {
    if (varu_temp[n] != varu_temp[n-1]) {
      indx = indx+1;
      var_u[indx] = varu_temp[n];
    }
  }
  
  // Update number of variables according to unique vars
  NVARU = N_u;
  
  for (int n=0; n<NCOV; n++) {
    Print() << "SF pairs (" << n << ") = " << s_pairA[n] << ", " << s_pairB[n] << std::endl;
  }
  Print() << "SF numPairs = " << NCOV << std::endl;

  for (int n=0; n<NCOV; n++) {
    if(s_pairA[n]<0 || s_pairA[n]>=NVAR || s_pairB[n]<0 || s_pairB[n]>=NVAR)
       amrex::Error("StructFact::StructFact() - Invalid pair select values: must be between 0 and (num of varibles - 1)");
  }
  //////////////////////////////////////////////////////

  // Note that we are defining with NO ghost cells

  cov_real.define(ba_in, dmap_in, NCOV, 0);
  cov_imag.define(ba_in, dmap_in, NCOV, 0);
  cov_mag.define( ba_in, dmap_in, NCOV, 0);
  cov_real.setVal(0.0);
  cov_imag.setVal(0.0);
  cov_mag.setVal( 0.0);

  cov_names.resize(NCOV);
  std::string x;
  int cnt = 0;
  for (int n=0; n<NCOV; n++) {
    x = "struct_fact";
    x += '_';
    x += var_names[s_pairB[n]];
    x += '_';
    x += var_names[s_pairA[n]];
    cov_names[cnt] = x;
    cnt++;
  }
}

void StructFact::FortStructure(const MultiFab& variables, const Geometry& geom,
                               const int& reset) {

  BL_PROFILE_VAR("StructFact::FortStructure()",FortStructure);

  const BoxArray& ba = variables.boxArray();
  const DistributionMapping& dm = variables.DistributionMap();
  
  MultiFab variables_dft_real, variables_dft_imag;
  variables_dft_real.define(ba, dm, NVAR, 0);
  variables_dft_imag.define(ba, dm, NVAR, 0);

  ComputeFFT(variables, variables_dft_real, variables_dft_imag, geom);

  // temporary storage built on BoxArray and DistributionMapping of "variables"
  // One case where "variables" and "cov_real/imag/mag" may have different DistributionMappings
  // is for flattened MFs with one grid newly built flattened MFs may be on a different
  // processor than the flattened MF used to build cov_real/imag/mag
  // or in general, problems that are not perfectly load balanced
  MultiFab cov_temp;
  cov_temp.define(ba, dm, 1, 0);

  // temporary storage built on BoxArray and DistributionMapping of "cov_real/imag/mag"
  MultiFab cov_temp2;
  cov_temp2.define(cov_real.boxArray(), cov_real.DistributionMap(), 1, 0);
  
  int index = 0;
  for (int n=0; n<NCOV; n++) {
    int i = s_pairA[n];
    int j = s_pairB[n];
    
    // Compute temporary real and imaginary components of covariance

    // Real component of covariance
    cov_temp.setVal(0.0);
    MultiFab::AddProduct(cov_temp,variables_dft_real,i,variables_dft_real,j,0,1,0);
    MultiFab::AddProduct(cov_temp,variables_dft_imag,i,variables_dft_imag,j,0,1,0);

    // copy into a MF with same ba and dm as cov_real/imag/mag
    cov_temp2.ParallelCopy(cov_temp,0,0,1);
        
    if (reset == 1) {
        MultiFab::Copy(cov_real,cov_temp2,0,index,1,0);
    } else {        
        MultiFab::Add(cov_real,cov_temp2,0,index,1,0);
    }

    // Imaginary component of covariance
    cov_temp.setVal(0.0);
    MultiFab::AddProduct(cov_temp,variables_dft_imag,i,variables_dft_real,j,0,1,0);
    cov_temp.mult(-1.0,0);
    MultiFab::AddProduct(cov_temp,variables_dft_real,i,variables_dft_imag,j,0,1,0);

    // copy into a MF with same ba and dm as cov_real/imag/mag
    cov_temp2.ParallelCopy(cov_temp,0,0,1);
    
    if (reset == 1) {
        MultiFab::Copy(cov_imag,cov_temp2,0,index,1,0);
    } else {
        MultiFab::Add(cov_imag,cov_temp2,0,index,1,0);
    }

    index++;

  }

  bool write_data = false;
  if (write_data) {
    std::string plotname; 
    plotname = "a_VAR";
    VisMF::Write(variables,plotname);
    plotname = "a_COV_REAL"; 
    VisMF::Write(cov_real,plotname);
    plotname = "a_COV_IMAG"; 
    VisMF::Write(cov_imag,plotname);
  }

  if (reset == 1) {
      nsamples = 1;
  } else {
      nsamples++;
  }
  
}

void StructFact::Reset() {

    BL_PROFILE_VAR("StructFact::Reset()", StructFactReset);
  
    cov_real.setVal(0.);
    cov_imag.setVal(0.);
    nsamples = 0;
    
}

void StructFact::ComputeFFT(const MultiFab& variables,
			    MultiFab& variables_dft_real, 
			    MultiFab& variables_dft_imag,
			    const Geometry& geom) {

    BL_PROFILE_VAR("StructFact::ComputeFFT()", ComputeFFT);

#ifdef AMREX_USE_CUDA
    Print() << "Using cuFFT\n";
#else
    Print() << "Using FFTW\n";
#endif

    bool is_flattened = false;

    long npts;

    // Initialize the boxarray "ba_onegrid" from the single box "domain"
    BoxArray ba_onegrid;
    {
      Box domain = geom.Domain();
      ba_onegrid.define(domain);

      if (domain.bigEnd(AMREX_SPACEDIM-1) == 0) {
          is_flattened = true;
      }

#if (AMREX_SPACEDIM == 2)
      npts = (domain.length(0)*domain.length(1));
#elif (AMREX_SPACEDIM == 3)
      npts = (domain.length(0)*domain.length(1)*domain.length(2));
#endif

    }

    Real sqrtnpts = std::sqrt(npts);

    DistributionMapping dmap_onegrid(ba_onegrid);

    // we will take one FFT at a time and copy the answer into the
    // corresponding component
    MultiFab variables_onegrid;
    MultiFab variables_dft_real_onegrid;
    MultiFab variables_dft_imag_onegrid;
    variables_onegrid.define(ba_onegrid, dmap_onegrid, 1, 0);
    variables_dft_real_onegrid.define(ba_onegrid, dmap_onegrid, 1, 0);
    variables_dft_imag_onegrid.define(ba_onegrid, dmap_onegrid, 1, 0);

//    fftw_mpi_init();

#ifdef AMREX_USE_CUDA
    using FFTplan = cufftHandle;
    using FFTcomplex = cuDoubleComplex;
#else
    using FFTplan = fftw_plan;
    using FFTcomplex = fftw_complex;
#endif

    // contain to store FFT - note it is shrunk by "half" in x
    Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > > spectral_field;

    Vector<FFTplan> forward_plan;

    // for CUDA builds we only need to build the plan once; track whether we did
    bool built_plan = false;
    
    for (int comp=0; comp<NVAR; comp++) {

        bool comp_fft = false;
        for (int i=0; i<NVARU; i++) {
            if (comp == var_u[i]) {
                comp_fft = true;
                break;
            }
        }

	if (comp_fft == false) continue;

        variables_onegrid.ParallelCopy(variables,comp,0,1);

        if (!built_plan) {

            for (MFIter mfi(variables_onegrid); mfi.isValid(); ++mfi) {

                // grab a single box including ghost cell range
                Box realspace_bx = mfi.fabbox();

                // size of box including ghost cell range
                IntVect fft_size = realspace_bx.length(); // This will be different for hybrid FFT

                // this is the size of the box, except the 0th component is 'halved plus 1'
                IntVect spectral_bx_size = fft_size;
                spectral_bx_size[0] = fft_size[0]/2 + 1;

                // spectral box
                Box spectral_bx = Box(IntVect(0), spectral_bx_size - IntVect(1));

                spectral_field.emplace_back(new BaseFab<GpuComplex<Real> >(spectral_bx,1,
                                                                       The_Device_Arena()));
                spectral_field.back()->setVal<RunOn::Device>(0.0); // touch the memory

                FFTplan fplan;

#ifdef AMREX_USE_CUDA
                if (is_flattened) {
#if (AMREX_SPACEDIM == 2)
                    cufftResult result = cufftPlan1d(&fplan, fft_size[0], CUFFT_D2Z, 1);
                    if (result != CUFFT_SUCCESS) {
                        amrex::AllPrint() << " cufftplan1d forward failed! Error: "
                                          << cufftErrorToString(result) << "\n";
                    }
#elif (AMREX_SPACEDIM == 3)
                    cufftResult result = cufftPlan2d(&fplan, fft_size[1], fft_size[0], CUFFT_D2Z);
                    if (result != CUFFT_SUCCESS) {
                        amrex::AllPrint() << " cufftplan2d forward failed! Error: "
                                          << cufftErrorToString(result) << "\n";
                    }
#endif
                } else {
#if (AMREX_SPACEDIM == 2)
                    cufftResult result = cufftPlan2d(&fplan, fft_size[1], fft_size[0], CUFFT_D2Z);
                    if (result != CUFFT_SUCCESS) {
                        amrex::AllPrint() << " cufftplan2d forward failed! Error: "
                                          << cufftErrorToString(result) << "\n";
                    }
#elif (AMREX_SPACEDIM == 3)
                    cufftResult result = cufftPlan3d(&fplan, fft_size[2], fft_size[1], fft_size[0], CUFFT_D2Z);
                    if (result != CUFFT_SUCCESS) {
                        amrex::AllPrint() << " cufftplan3d forward failed! Error: "
                                          << cufftErrorToString(result) << "\n";
                    }
#endif
                }
#else // host

                if (is_flattened) {
#if (AMREX_SPACEDIM == 2)
                    fplan = fftw_plan_dft_r2c_1d(fft_size[0],
                                                 variables_onegrid[mfi].dataPtr(),
                                                 reinterpret_cast<FFTcomplex*>
                                                 (spectral_field.back()->dataPtr()),
                                                 FFTW_ESTIMATE);
#elif (AMREX_SPACEDIM == 3)
                    fplan = fftw_plan_dft_r2c_2d(fft_size[1], fft_size[0],
                                                 variables_onegrid[mfi].dataPtr(),
                                                 reinterpret_cast<FFTcomplex*>
                                                 (spectral_field.back()->dataPtr()),
                                                 FFTW_ESTIMATE);
#endif
                } else {
#if (AMREX_SPACEDIM == 2)
                    fplan = fftw_plan_dft_r2c_2d(fft_size[1], fft_size[0],
                                                 variables_onegrid[mfi].dataPtr(),
                                                 reinterpret_cast<FFTcomplex*>
                                                 (spectral_field.back()->dataPtr()),
                                                 FFTW_ESTIMATE);
#elif (AMREX_SPACEDIM == 3)
                    fplan = fftw_plan_dft_r2c_3d(fft_size[2], fft_size[1], fft_size[0],
                                                 variables_onegrid[mfi].dataPtr(),
                                                 reinterpret_cast<FFTcomplex*>
                                                 (spectral_field.back()->dataPtr()),
                                                 FFTW_ESTIMATE);
#endif
                }
#endif

                forward_plan.push_back(fplan);
            }

	    built_plan = true;
        
        }

        ParallelDescriptor::Barrier();

        // ForwardTransform
        for (MFIter mfi(variables_onegrid); mfi.isValid(); ++mfi) {
            int i = mfi.LocalIndex();
#ifdef AMREX_USE_CUDA
            cufftSetStream(forward_plan[i], amrex::Gpu::gpuStream());
            cufftResult result = cufftExecD2Z(forward_plan[i],
                                              variables_onegrid[mfi].dataPtr(),
                                              reinterpret_cast<FFTcomplex*>
                                                  (spectral_field[i]->dataPtr()));
            if (result != CUFFT_SUCCESS) {
	      amrex::AllPrint() << " forward transform using cufftExec failed! Error: "
				<< cufftErrorToString(result) << "\n";
	    }
#else
            fftw_execute(forward_plan[i]);
#endif
        }

        // copy data to a full-sized MultiFab
        // this involves copying the complex conjugate from the half-sized field
        // into the appropriate place in the full MultiFab
        for (MFIter mfi(variables_dft_real_onegrid); mfi.isValid(); ++mfi) {

            Array4< GpuComplex<Real> > spectral = (*spectral_field[0]).array();

            Array4<Real> const& realpart = variables_dft_real_onegrid.array(mfi);
            Array4<Real> const& imagpart = variables_dft_imag_onegrid.array(mfi);

            Box bx = mfi.fabbox();

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i <= bx.length(0)/2) {
                    // copy value
                    realpart(i,j,k) = spectral(i,j,k).real();
                    imagpart(i,j,k) = spectral(i,j,k).imag();
                } else {
                    // copy complex conjugate
                    int iloc = bx.length(0)-i;
                    int jloc, kloc;
                    if (is_flattened) {
#if (AMREX_SPACEDIM == 2)
                        jloc = 0;
#elif (AMREX_SPACEDIM == 3)
                        jloc = (j == 0) ? 0 : bx.length(1)-j;
#endif
                        kloc = 0;
                    } else {
                        jloc = (j == 0) ? 0 : bx.length(1)-j;
#if (AMREX_SPACEDIM == 2)
                        kloc = 0;
#elif (AMREX_SPACEDIM == 3)
                        kloc = (k == 0) ? 0 : bx.length(2)-k;
#endif
                    }

                    realpart(i,j,k) =  spectral(iloc,jloc,kloc).real();
                    imagpart(i,j,k) = -spectral(iloc,jloc,kloc).imag();
                }

                realpart(i,j,k) /= sqrtnpts;
                imagpart(i,j,k) /= sqrtnpts;
            });

            /*
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                std::cout << "HACKFFT " << i << " " << j << " " << k << " "
                          << realpart(i,j,k) << " + " << imagpart(i,j,k) << "i"
                          << std::endl;
            });
            */
        }

        variables_dft_real.ParallelCopy(variables_dft_real_onegrid,0,comp,1);
        variables_dft_imag.ParallelCopy(variables_dft_imag_onegrid,0,comp,1);

    }

    // destroy fft plan
    for (int i = 0; i < forward_plan.size(); ++i) {
#ifdef AMREX_USE_CUDA
        cufftDestroy(forward_plan[i]);
#else
        fftw_destroy_plan(forward_plan[i]);
#endif
    }

//    fftw_mpi_cleanup();

}

void StructFact::WritePlotFile(const int step, const Real time, const Geometry& geom,
                               std::string plotfile_base,
                               const int& zero_avg) {
  
  BL_PROFILE_VAR("StructFact::WritePlotFile()",StructFactWritePlotFile);

  MultiFab plotfile;
  Vector<std::string> varNames;
  int nPlot = 1;

  // Build temp real & imag components
  const BoxArray& ba = cov_mag.boxArray();
  const DistributionMapping& dm = cov_mag.DistributionMap();

  MultiFab cov_real_temp(ba, dm, NCOV, 0);
  MultiFab cov_imag_temp(ba, dm, NCOV, 0);
  MultiFab::Copy(cov_real_temp, cov_real, 0, 0, NCOV, 0);
  MultiFab::Copy(cov_imag_temp, cov_imag, 0, 0, NCOV, 0);

  // Finalize covariances - scale & compute magnitude
  Finalize(cov_real_temp, cov_imag_temp, geom, zero_avg);

  //////////////////////////////////////////////////////////////////////////////////
  // Write out structure factor magnitude to plot file
  //////////////////////////////////////////////////////////////////////////////////

  std::string name = plotfile_base;
  name += "_mag";
  
  const std::string plotfilename1 = amrex::Concatenate(name,step,9);
  nPlot = NCOV;
  plotfile.define(cov_mag.boxArray(), cov_mag.DistributionMap(), nPlot, 0);
  varNames.resize(nPlot);

  for (int n=0; n<NCOV; n++) {
      varNames[n] = cov_names[n];
  }
  
  MultiFab::Copy(plotfile, cov_mag, 0, 0, NCOV, 0); // copy structure factor into plotfile

  Real dx = geom.CellSize(0);
  Real pi = 3.1415926535897932;
  Box domain = geom.Domain();

  RealBox real_box({AMREX_D_DECL(-pi/dx,-pi/dx,-pi/dx)},
                   {AMREX_D_DECL( pi/dx, pi/dx, pi/dx)});
  
  // check bc_vel_lo/hi to determine the periodicity
  Vector<int> is_periodic(AMREX_SPACEDIM,0);  // set to 0 (not periodic) by default
  for (int i=0; i<AMREX_SPACEDIM; ++i) {
      is_periodic[i] = geom.isPeriodic(i);
  }

  Geometry geom2;
  geom2.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    
  // write a plotfile
  WriteSingleLevelPlotfile(plotfilename1,plotfile,varNames,geom2,time,step);
  
  //////////////////////////////////////////////////////////////////////////////////
  // Write out real and imaginary components of structure factor to plot file
  //////////////////////////////////////////////////////////////////////////////////

  name = plotfile_base;
  name += "_real_imag";
  
  const std::string plotfilename2 = amrex::Concatenate(name,step,9);
  nPlot = 2*NCOV;
  plotfile.define(cov_mag.boxArray(), cov_mag.DistributionMap(), nPlot, 0);
  varNames.resize(nPlot);

  int cnt = 0; // keep a counter for plotfile variables
  for (int n=0; n<NCOV; n++) {
      varNames[cnt] = cov_names[cnt];
      varNames[cnt] += "_real";
      cnt++;
  }

  int index = 0;
  for (int n=0; n<NCOV; n++) {
      varNames[cnt] = cov_names[index];
      varNames[cnt] += "_imag";
      index++;
      cnt++;
  }

  MultiFab::Copy(plotfile,cov_real_temp,0,0,   NCOV,0);
  MultiFab::Copy(plotfile,cov_imag_temp,0,NCOV,NCOV,0);

  // write a plotfile
  WriteSingleLevelPlotfile(plotfilename2,plotfile,varNames,geom2,time,step);
}

void StructFact::Finalize(MultiFab& cov_real_in, MultiFab& cov_imag_in,
                          const Geometry& geom, const int& zero_avg) {

  BL_PROFILE_VAR("StructFact::Finalize()",StructFactFinalize);
  
  Real nsamples_inv = 1.0/(Real)nsamples;
  
  ShiftFFT(cov_real_in,geom,zero_avg);
  ShiftFFT(cov_imag_in,geom,zero_avg);

  cov_real_in.mult(nsamples_inv);
  for (int d=0; d<NCOV; d++) {
      cov_real_in.mult(scaling[d],d,1);
  }

  cov_imag_in.mult(nsamples_inv);
  for (int d=0; d<NCOV; d++) {
      cov_imag_in.mult(scaling[d],d,1);
  }
  
  cov_mag.setVal(0.0);
  MultiFab::AddProduct(cov_mag,cov_real_in,0,cov_real_in,0,0,NCOV,0);
  MultiFab::AddProduct(cov_mag,cov_imag_in,0,cov_imag_in,0,0,NCOV,0);

  SqrtMF(cov_mag);

}

// Finalize covariances - scale & compute magnitude
void StructFact::CallFinalize( const Geometry& geom,
                               const int& zero_avg) {
  
  BL_PROFILE_VAR("CallFinalize()",CallFinalize);

  // Build temp real & imag components
  const BoxArray& ba = cov_mag.boxArray();
  const DistributionMapping& dm = cov_mag.DistributionMap();

  MultiFab cov_real_temp(ba, dm, NCOV, 0);
  MultiFab cov_imag_temp(ba, dm, NCOV, 0);
  MultiFab::Copy(cov_real_temp, cov_real, 0, 0, NCOV, 0);
  MultiFab::Copy(cov_imag_temp, cov_imag, 0, 0, NCOV, 0);

  // Finalize covariances - scale & compute magnitude
  Finalize(cov_real_temp, cov_imag_temp, geom, zero_avg);
}



void StructFact::ShiftFFT(MultiFab& dft_out, const Geometry& geom, const int& zero_avg) {

  BL_PROFILE_VAR("StructFact::ShiftFFT()",ShiftFFT);

  /*
    Shifting rules:

    For domains from (0,0,0) to (Nx-1,Ny-1,Nz-1)

    For any cells with i index >= Nx/2, these values are complex conjugates of the corresponding
    entry where (Nx-i,Ny-j,Nz-k) UNLESS that index is zero, in which case you use 0.

    e.g. for an 8^3 domain, any cell with i index 

    Cell (6,2,3) is complex conjugate of (2,6,5)

    Cell (4,1,0) is complex conjugate of (4,7,0)  (note that the FFT is computed for 0 <= i <= Nx/2)
  */

  BoxArray ba_onegrid;
  {
      Box domain = geom.Domain();
      
      // Initialize the boxarray "ba" from the single box "bx"
      ba_onegrid.define(domain);
  }

  DistributionMapping dmap_onegrid(ba_onegrid);

  MultiFab dft_onegrid;
  MultiFab dft_onegrid_temp;
  dft_onegrid     .define(ba_onegrid, dmap_onegrid, 1, 0);
  dft_onegrid_temp.define(ba_onegrid, dmap_onegrid, 1, 0);

  for (int d=0; d<NCOV; d++) {
    dft_onegrid_temp.ParallelCopy(dft_out, d, 0, 1);

    // Shift DFT by N/2+1 (pi)
    for (MFIter mfi(dft_onegrid); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real>& dft = dft_onegrid.array(mfi);
        const Array4<Real>& dft_temp = dft_onegrid_temp.array(mfi);

        if (zero_avg == 1) {
#if (AMREX_SPACEDIM == 2)
            dft_temp(bx.smallEnd(0),bx.smallEnd(1),0) = 0.;
#elif (AMREX_SPACEDIM == 3)
            dft_temp(bx.smallEnd(0),bx.smallEnd(1),bx.smallEnd(2)) = 0.;
#endif
        }

        int nx = bx.length(0);
        int nxh = (nx+1)/2;
        int ny = bx.length(1);
        int nyh = (ny+1)/2;
#if (AMREX_SPACEDIM == 3)
        int nz = bx.length(2);
        int nzh = (nz+1)/2;
#endif

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ip = (i+nxh)%nx;
            int jp = (j+nyh)%ny;
            int kp = 0;
#if (AMREX_SPACEDIM == 3)
            kp = (k+nzh)%nz;
#endif
            dft(ip,jp,kp) = dft_temp(i,j,k);
        });
    }

    dft_out.ParallelCopy(dft_onegrid, 0, d, 1);
  }

}

// integrate cov_mag over k shells
void StructFact::IntegratekShells(const int& step, const Geometry& geom) {

    BL_PROFILE_VAR("StructFact::IntegratekShells",IntegratekShells);

    GpuArray<int,AMREX_SPACEDIM> center;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        center[d] = n_cells[d]/2;
    }

    int npts = n_cells[0]/2-1;
    int npts_sq = npts*npts;

    Gpu::DeviceVector<Real> phisum_vect(npts);
    Gpu::DeviceVector<int>  phicnt_vect(npts);

//    Gpu::DeviceVector<Real> phisum_vect_large(npts_sq);
//    Gpu::DeviceVector<int>  phicnt_vect_large(npts_sq);

    for (int d=0; d<npts; ++d) {
        phisum_vect[d] = 0.;
        phicnt_vect[d] = 0;
    }

//    for (int d=0; d<npts_sq; ++d) {
//        phisum_vect_large[d] = 0.;
//        phicnt_vect_large[d] = 0;
//    }
    
    Real* phisum_gpu = phisum_vect.dataPtr();  // pointer to data
    int*  phicnt_gpu = phicnt_vect.dataPtr();  // pointer to data
    
  //  Real* phisum_large_gpu = phisum_vect_large.dataPtr();  // pointer to data
 //   int*  phicnt_large_gpu = phicnt_vect_large.dataPtr();  // pointer to data

    // only consider cells that are within 15k of the center point
    
    for ( MFIter mfi(cov_mag,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.tilebox();

        const Array4<Real> & cov = cov_mag.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ilen = amrex::Math::abs(i-center[0]);
            int jlen = amrex::Math::abs(j-center[1]);
            int klen = (AMREX_SPACEDIM == 3) ? amrex::Math::abs(k-center[2]) : 0;

            Real dist = (ilen*ilen + jlen*jlen + klen*klen);
 //           int idist = (ilen*ilen + jlen*jlen + klen*klen);
            dist = std::sqrt(dist);
            
            if ( dist <= center[0]-0.5) {
	        dist = dist+0.5;
                int cell = int(dist);
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
		  amrex::HostDevice::Atomic::Add(&(phisum_gpu[cell]), cov(i,j,k,d));
//		    phisum_large_gpu[idist]  += cov(i,j,k,d);
                }
		amrex::HostDevice::Atomic::Add(&(phicnt_gpu[cell]),1);
 //               ++phicnt_large_gpu[idist];
            }
        });
    }
        
    for (int d=1; d<npts; ++d) {
        ParallelDescriptor::ReduceRealSum(phisum_vect[d]);
        ParallelDescriptor::ReduceIntSum(phicnt_vect[d]);
    }
        
#if 0
    for (int d=1; d<npts_sq; ++d) {
        ParallelDescriptor::ReduceRealSum(phisum_vect_large[d]);
        ParallelDescriptor::ReduceIntSum(phicnt_vect_large[d]);
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream turb_disc;
        std::string turbNamedisc = "turb_disc";
        turbNamedisc += std::to_string(step);
        turbNamedisc += ".txt";
        
        turb_disc.open(turbNamedisc);
        for (int d=1; d<npts_sq; ++d) {
	    if(phicnt_vect_large[d]>0) {
		Real dreal = d;
                turb_disc << sqrt(dreal) << " " << 4.*M_PI*d*phisum_vect_large[d]/phicnt_vect_large[d] << std::endl;
	    }
        }
    }

    Real dk = 1.;
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream turb_alt;
        std::string turbNamealt = "turb_alt";
        turbNamealt += std::to_string(step);
        turbNamealt += ".txt";
        
        turb_alt.open(turbNamealt);
        for (int d=1; d<npts; ++d) {
            turb_alt << d << " " << phisum_vect[d] << std::endl;
        }
    }
#endif

    Real dk = 1.;
    
#if (AMREX_SPACEDIM == 2)
    for (int d=1; d<npts; ++d) {
      //  phisum_vect[d] *= 2.*M_PI*d*dk*dk/phicnt_vect[d];
        phisum_vect[d] *= 2.*M_PI*(d*dk+.5*dk*dk)/phicnt_vect[d];
    }
#else
    for (int d=1; d<npts; ++d) {
      //  phisum_vect[d] *= 4.*M_PI*(d*d)*dk*dk*dk/phicnt_vect[d];
      //  phisum_vect[d] *= 4.*M_PI*(d*d*dk+d*dk*dk+dk*dk*dk/3.)/phicnt_vect[d];
        phisum_vect[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_vect[d];
    }
#endif
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream turb;
        std::string turbBaseName = "turb";
        //turbName += std::to_string(step);
        std::string turbName = Concatenate(turbBaseName,step,7);
        turbName += ".txt";
        
        turb.open(turbName);
        for (int d=1; d<npts; ++d) {
            turb << d << " " << phisum_vect[d] << std::endl;
        }
    }
}

void StructFact::AddToExternal(MultiFab& x_mag, MultiFab& x_realimag, const Geometry& geom, const int& zero_avg) {

    BL_PROFILE_VAR("StructFact::AddToExternal",AddToExternal);

    MultiFab plotfile;
    int nPlot = 1;

    // Build temp real & imag components
    const BoxArray& ba = cov_mag.boxArray();
    const DistributionMapping& dm = cov_mag.DistributionMap();

    MultiFab cov_real_temp(ba, dm, NCOV, 0);
    MultiFab cov_imag_temp(ba, dm, NCOV, 0);
    MultiFab::Copy(cov_real_temp, cov_real, 0, 0, NCOV, 0);
    MultiFab::Copy(cov_imag_temp, cov_imag, 0, 0, NCOV, 0);

    // Finalize covariances - scale & compute magnitude
    Finalize(cov_real_temp, cov_imag_temp, geom, zero_avg);

    nPlot = NCOV;
    plotfile.define(cov_mag.boxArray(), cov_mag.DistributionMap(), nPlot, 0);
    MultiFab::Copy(plotfile, cov_mag, 0, 0, NCOV, 0); // copy structure factor into plotfile
    MultiFab::Add(x_mag,plotfile,0,0,NCOV,0);

    nPlot = 2*NCOV;
    plotfile.define(cov_mag.boxArray(), cov_mag.DistributionMap(), nPlot, 0);
    MultiFab::Copy(plotfile,cov_real_temp,0,0,   NCOV,0);
    MultiFab::Copy(plotfile,cov_imag_temp,0,NCOV,NCOV,0);
    MultiFab::Add(x_realimag,plotfile,0,0,2*NCOV,0);

}
