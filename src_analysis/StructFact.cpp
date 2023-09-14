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

#ifdef AMREX_USE_HIP
std::string rocfftErrorToString (const rocfft_status err)
{
    if              (err == rocfft_status_success) {
        return std::string("rocfft_status_success");
    } else if       (err == rocfft_status_failure) {
        return std::string("rocfft_status_failure");
    } else if       (err == rocfft_status_invalid_arg_value) {
        return std::string("rocfft_status_invalid_arg_value");
    } else if       (err == rocfft_status_invalid_dimensions) {
        return std::string("rocfft_status_invalid_dimensions");
    } else if       (err == rocfft_status_invalid_array_type) {
        return std::string("rocfft_status_invalid_array_type");
    } else if       (err == rocfft_status_invalid_strides) {
        return std::string("rocfft_status_invalid_strides");
    } else if       (err == rocfft_status_invalid_distance) {
        return std::string("rocfft_status_invalid_distance");
    } else if       (err == rocfft_status_invalid_offset) {
        return std::string("rocfft_status_invalid_offset");
    } else {
        return std::to_string(err) + " (unknown error code)";
    }
}

void assert_rocfft_status (std::string const& name, rocfft_status status)
{
    if (status != rocfft_status_success) {
        amrex::AllPrint() <<  name + " failed! Error: " + rocfftErrorToString(status) << "\n";;
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

    Vector<int> l_s_pairA(NVAR*(NVAR+1)/2);
    Vector<int> l_s_pairB(NVAR*(NVAR+1)/2);

    int counter=0;
    for (int i=0; i<NVAR; ++i) {
        for (int j=i; j<NVAR; ++j) {
            l_s_pairA[counter] = j;
            l_s_pairB[counter] = i;
            ++counter;
        }
    }      

    define(ba_in, dmap_in, var_names, var_scaling_in, l_s_pairA, l_s_pairB, verbosity_in);

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

void StructFact::defineDecomp(const amrex::BoxArray& ba_in,
                              const amrex::DistributionMapping& dmap_in,
                              const Vector< std::string >& /*var_names*/,
                              const amrex::Vector< amrex::Real >& var_scaling_in,
                              const Vector< int >& s_pairA_in,
                              const Vector< int >& s_pairB_in)
{

  BL_PROFILE_VAR("StructFact::defineDecomp()",StructFactDefineDecomp);

  decompose = true;
  
  if (s_pairA_in.size() != s_pairB_in.size())
        amrex::Error("StructFact::define() - Must have an equal number of components");

  NVAR = 3;
  NCOV = 6;
  scaling.resize(NCOV);
  for (int n=0; n<NCOV; n++) {
      scaling[n] = 1.0/var_scaling_in[n];
  }
  
  s_pairA.resize(3);
  s_pairB.resize(3);
  
  // Set vectors identifying covariance pairs
  for (int n=0; n<3; n++) {
    s_pairA[n] = s_pairA_in[n];
    s_pairB[n] = s_pairB_in[n];
  }
  
  // Create vector of unique variable indices to select which to take the FFT
  NVARU = NVAR;    // temporary before selecting unique variables
  var_u.resize(NVARU);
  for (int n=0; n<NVARU; n++) {
    var_u[n] = s_pairA[n];
  }

  //BoxArray ba_onegrid;
  //{
  //  Box domain = geom.Domain();
  //  ba_onegrid.define(domain);
  //}
  //DistributionMapping dmap_onegrid(ba_onegrid);

  vel_sol_real.define(ba_in, dmap_in, 3, 0);
  vel_sol_imag.define(ba_in, dmap_in, 3, 0);
  vel_dil_real.define(ba_in, dmap_in, 3, 0);
  vel_dil_imag.define(ba_in, dmap_in, 3, 0);

  cov_real.define(ba_in, dmap_in, 6, 0);
  cov_imag.define(ba_in, dmap_in, 6, 0);
  cov_mag.define( ba_in, dmap_in, 6, 0);
  cov_real.setVal(0.0);
  cov_imag.setVal(0.0);
  cov_mag.setVal (0.0);
  
  NCOV = 6;
  NVAR = 3;
  scaling.resize(NCOV);
  for (int n=0; n<NCOV; n++) {
      scaling[n] = 1.0/var_scaling_in[n];
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

void StructFact::FortStructureDecomp(const MultiFab& vel, const Geometry& geom,
                                     const int& reset) 
{


  BL_PROFILE_VAR("StructFact::FortStructureDecomp()",FortStructureDecomp);

  if (!decompose) amrex::Error("StructFact::FortStructureDecomp() is specific for vel decomposition in turbulence");

  const BoxArray& ba = vel.boxArray();
  const DistributionMapping& dm = vel.DistributionMap();
  
  MultiFab vel_dft_real, vel_dft_imag;
  //BoxArray ba_onegrid;
  //{
  //  Box domain = geom.Domain();
  //  ba_onegrid.define(domain);
  //}
  //DistributionMapping dmap_onegrid(ba_onegrid);
  vel_dft_real.define(ba, dm, 3, 0);
  vel_dft_imag.define(ba, dm, 3, 0);

  ComputeFFT(vel, vel_dft_real, vel_dft_imag, geom);
  
  DecomposeVelFourier(vel_dft_real, vel_dft_imag, geom);

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
 
  // solenoidal
  int index = 0;
  for (int n=0; n<3; n++) {
    cov_temp.setVal(0.0);
    MultiFab::AddProduct(cov_temp,vel_sol_real,n,vel_sol_real,n,0,1,0);
    MultiFab::AddProduct(cov_temp,vel_sol_imag,n,vel_sol_imag,n,0,1,0);

    // copy into a MF with same ba and dm as cov_real/imag/mag
    cov_temp2.ParallelCopy(cov_temp,0,0,1);
        
    if (reset == 1) {
        MultiFab::Copy(cov_real,cov_temp2,0,index,1,0);
    } else {        
        MultiFab::Add(cov_real,cov_temp2,0,index,1,0);
    }

    // Imaginary component of covariance
    cov_temp.setVal(0.0);
    MultiFab::AddProduct(cov_temp,vel_sol_imag,n,vel_sol_real,n,0,1,0);
    cov_temp.mult(-1.0,0);
    MultiFab::AddProduct(cov_temp,vel_sol_real,n,vel_sol_imag,n,0,1,0);

    // copy into a MF with same ba and dm as cov_real/imag/mag
    cov_temp2.ParallelCopy(cov_temp,0,0,1);
    
    if (reset == 1) {
        MultiFab::Copy(cov_imag,cov_temp2,0,index,1,0);
    } else {
        MultiFab::Add(cov_imag,cov_temp2,0,index,1,0);
    }
    index++;
  }
  
  // dilatational
  for (int n=0; n<3; n++) {
    cov_temp.setVal(0.0);
    MultiFab::AddProduct(cov_temp,vel_dil_real,n,vel_dil_real,n,0,1,0);
    MultiFab::AddProduct(cov_temp,vel_dil_imag,n,vel_dil_imag,n,0,1,0);

    // copy into a MF with same ba and dm as cov_real/imag/mag
    cov_temp2.ParallelCopy(cov_temp,0,0,1);
        
    if (reset == 1) {
        MultiFab::Copy(cov_real,cov_temp2,0,index,1,0);
    } else {        
        MultiFab::Add(cov_real,cov_temp2,0,index,1,0);
    }

    // Imaginary component of covariance
    cov_temp.setVal(0.0);
    MultiFab::AddProduct(cov_temp,vel_dil_imag,n,vel_dil_real,n,0,1,0);
    cov_temp.mult(-1.0,0);
    MultiFab::AddProduct(cov_temp,vel_dil_real,n,vel_dil_imag,n,0,1,0);

    // copy into a MF with same ba and dm as cov_real/imag/mag
    cov_temp2.ParallelCopy(cov_temp,0,0,1);
    
    if (reset == 1) {
        MultiFab::Copy(cov_imag,cov_temp2,0,index,1,0);
    } else {
        MultiFab::Add(cov_imag,cov_temp2,0,index,1,0);
    }
    index++;
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
			    const Geometry& geom,
                bool unpack)
{

    BL_PROFILE_VAR("StructFact::ComputeFFT()", ComputeFFT);

#ifdef AMREX_USE_CUDA
    // Print() << "Using cuFFT\n";
#elif AMREX_USE_HIP
    // Print() << "Using rocFFT\n";
#else
    // Print() << "Using FFTW\n";
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
#elif AMREX_USE_HIP
    using FFTplan = rocfft_plan;
    using FFTcomplex = double2;
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

#ifdef AMREX_USE_CUDA // CUDA
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
#elif AMREX_USE_HIP // HIP
                if (is_flattened) {
#if (AMREX_SPACEDIM == 2)
                    const std::size_t lengths[] = {std::size_t(fft_size[0])};
                    rocfft_status result = rocfft_plan_create(&fplan, rocfft_placement_notinplace, 
                                                              rocfft_transform_type_real_forward, rocfft_precision_double,
                                                              1, lengths, 1, nullptr);
                    assert_rocfft_status("rocfft_plan_create", result);
#elif (AMREX_SPACEDIM == 3)
                    const std::size_t lengths[] = {std::size_t(fft_size[0]),std::size_t(fft_size[1])};
                    rocfft_status result = rocfft_plan_create(&fplan, rocfft_placement_notinplace, 
                                                              rocfft_transform_type_real_forward, rocfft_precision_double,
                                                              2, lengths, 1, nullptr);
                    assert_rocfft_status("rocfft_plan_create", result);
#endif
                } else {
#if (AMREX_SPACEDIM == 2)
                    const std::size_t lengths[] = {std::size_t(fft_size[0]),std::size_t(fft_size[1])};
                    rocfft_status result = rocfft_plan_create(&fplan, rocfft_placement_notinplace, 
                                                              rocfft_transform_type_real_forward, rocfft_precision_double,
                                                              2, lengths, 1, nullptr);
                    assert_rocfft_status("rocfft_plan_create", result);
#elif (AMREX_SPACEDIM == 3)
                    const std::size_t lengths[] = {std::size_t(fft_size[0]),std::size_t(fft_size[1]),std::size_t(fft_size[2])};
                    rocfft_status result = rocfft_plan_create(&fplan, rocfft_placement_notinplace, 
                                                              rocfft_transform_type_real_forward, rocfft_precision_double,
                                                              3, lengths, 1, nullptr);
                    assert_rocfft_status("rocfft_plan_create", result);
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
#elif AMREX_USE_HIP
            rocfft_execution_info execinfo = nullptr;
            rocfft_status result = rocfft_execution_info_create(&execinfo);
            assert_rocfft_status("rocfft_execution_info_create", result);

            std::size_t buffersize = 0;
            result = rocfft_plan_get_work_buffer_size(forward_plan[i], &buffersize);
            assert_rocfft_status("rocfft_plan_get_work_buffer_size", result);

            void* buffer = amrex::The_Arena()->alloc(buffersize);
            result = rocfft_execution_info_set_work_buffer(execinfo, buffer, buffersize);
            assert_rocfft_status("rocfft_execution_info_set_work_buffer", result);

            result = rocfft_execution_info_set_stream(execinfo, amrex::Gpu::gpuStream());
            assert_rocfft_status("rocfft_execution_info_set_stream", result);

	        amrex::Real* variables_onegrid_ptr = variables_onegrid[mfi].dataPtr();
	        FFTcomplex* spectral_field_ptr = reinterpret_cast<FFTcomplex*>(spectral_field[i]->dataPtr());
            result = rocfft_execute(forward_plan[i],
                                    (void**) &variables_onegrid_ptr, // in
                                    (void**) &spectral_field_ptr, // out
                                    execinfo);
            assert_rocfft_status("rocfft_execute", result);
            amrex::Gpu::streamSynchronize();
            amrex::The_Arena()->free(buffer);
            result = rocfft_execution_info_destroy(execinfo);
            assert_rocfft_status("rocfft_execution_info_destroy", result);
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
                /*
                  Unpacking rules:

                  For domains from (0,0,0) to (Nx-1,Ny-1,Nz-1)

                  For any cells with i index > Nx/2, these values are complex conjugates of the corresponding
                  entry where (Nx-i,Ny-j,Nz-k) UNLESS that index is zero, in which case you use 0.

                  e.g. for an 8^3 domain, any cell with i index

                  Cell (6,2,3) is complex conjugate of (2,6,5)

                  Cell (4,1,0) is complex conjugate of (4,7,0)  (note that the FFT is computed for 0 <= i <= Nx/2)
                */
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

                    if (unpack) {
                        realpart(i,j,k) =  spectral(iloc,jloc,kloc).real();
                        imagpart(i,j,k) = -spectral(iloc,jloc,kloc).imag();
                    }
                    else {
                        realpart(i,j,k) =  0.0;
                        imagpart(i,j,k) =  0.0;
                    }
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
#elif AMREX_USE_HIP
        rocfft_plan_destroy(forward_plan[i]);
#else
        fftw_destroy_plan(forward_plan[i]);
#endif
    }
//    fftw_mpi_cleanup();
}

void StructFact::InverseFFT(MultiFab& variables,
			    const MultiFab& variables_dft_real, 
			    const MultiFab& variables_dft_imag,
			    const Geometry& geom)
{

    BL_PROFILE_VAR("StructFact::InverseFFT()", InverseFFT);

#ifdef AMREX_USE_CUDA
    // Print() << "Using cuFFT\n";
#elif AMREX_USE_HIP
    // Print() << "Using rocFFT\n";
#else
    // Print() << "Using FFTW\n";
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
    variables_onegrid.define(ba_onegrid, dmap_onegrid, 1, 0);
    
    MultiFab variables_dft_real_onegrid;
    MultiFab variables_dft_imag_onegrid;
    variables_dft_real_onegrid.define(ba_onegrid, dmap_onegrid, 1, 0);
    variables_dft_imag_onegrid.define(ba_onegrid, dmap_onegrid, 1, 0);

//    fftw_mpi_init();

#ifdef AMREX_USE_CUDA
    using FFTplan = cufftHandle;
    using FFTcomplex = cuDoubleComplex;
#elif AMREX_USE_HIP
    using FFTplan = rocfft_plan;
    using FFTcomplex = double2;
#else
    using FFTplan = fftw_plan;
    using FFTcomplex = fftw_complex;
#endif

    // contain to store FFT - note it is shrunk by "half" in x
    Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > > spectral_field;

    Vector<FFTplan> backward_plan;

    // for CUDA builds we only need to build the plan once; track whether we did
    bool built_plan = false;
    
    for (int comp=0; comp<variables_dft_real.nComp(); comp++) {

        // build spectral field from multifabs
        variables_dft_real_onegrid.ParallelCopy(variables_dft_real,comp,0,1);
        variables_dft_imag_onegrid.ParallelCopy(variables_dft_imag,comp,0,1);
        
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

            Array4< GpuComplex<Real> > spectral = (*spectral_field[0]).array();
            Array4<Real> const& realpart = variables_dft_real_onegrid.array(mfi);
            Array4<Real> const& imagpart = variables_dft_imag_onegrid.array(mfi);

            Box bx = mfi.fabbox();

            amrex::ParallelFor(bx, 
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i <= bx.length(0)/2) {
                    GpuComplex<Real> copy(realpart(i,j,k),imagpart(i,j,k));
                    spectral(i,j,k) = copy;
                }
            });
        }

        // build FFTplan if necessary
        for (MFIter mfi(variables_onegrid); mfi.isValid(); ++mfi) {
            
            if (!built_plan) {

                Box realspace_bx = mfi.fabbox();

                IntVect fft_size = realspace_bx.length();

                FFTplan bplan;

#ifdef AMREX_USE_CUDA // CUDA
                if (is_flattened) {
#if (AMREX_SPACEDIM == 2)
                    cufftResult result = cufftPlan1d(&bplan, fft_size[0], CUFFT_Z2D, 1);
                    if (result != CUFFT_SUCCESS) {
                        amrex::AllPrint() << " cufftplan1d forward failed! Error: "
                                          << cufftErrorToString(result) << "\n";
                    }
#elif (AMREX_SPACEDIM == 3)
                    cufftResult result = cufftPlan2d(&bplan, fft_size[1], fft_size[0], CUFFT_Z2D);
                    if (result != CUFFT_SUCCESS) {
                        amrex::AllPrint() << " cufftplan2d forward failed! Error: "
                                          << cufftErrorToString(result) << "\n";
                    }
#endif
                } else {
#if (AMREX_SPACEDIM == 2)
                    cufftResult result = cufftPlan2d(&bplan, fft_size[1], fft_size[0], CUFFT_Z2D);
                    if (result != CUFFT_SUCCESS) {
                        amrex::AllPrint() << " cufftplan2d forward failed! Error: "
                                          << cufftErrorToString(result) << "\n";
                    }
#elif (AMREX_SPACEDIM == 3)
                    cufftResult result = cufftPlan3d(&bplan, fft_size[2], fft_size[1], fft_size[0], CUFFT_Z2D);
                    if (result != CUFFT_SUCCESS) {
                        amrex::AllPrint() << " cufftplan3d forward failed! Error: "
                                          << cufftErrorToString(result) << "\n";
                    }
#endif
                }
#elif AMREX_USE_HIP // HIP
                if (is_flattened) {
#if (AMREX_SPACEDIM == 2)
                    const std::size_t lengths[] = {std::size_t(fft_size[0])};
                    rocfft_status result = rocfft_plan_create(&bplan, rocfft_placement_notinplace, 
                                                              rocfft_transform_type_real_inverse, rocfft_precision_double,
                                                              1, lengths, 1, nullptr);
                    assert_rocfft_status("rocfft_plan_create", result);
#elif (AMREX_SPACEDIM == 3)
                    const std::size_t lengths[] = {std::size_t(fft_size[0]),std::size_t(fft_size[1])};
                    rocfft_status result = rocfft_plan_create(&bplan, rocfft_placement_notinplace, 
                                                              rocfft_transform_type_real_inverse, rocfft_precision_double,
                                                              2, lengths, 1, nullptr);
                    assert_rocfft_status("rocfft_plan_create", result);
#endif
                } else {
#if (AMREX_SPACEDIM == 2)
                    const std::size_t lengths[] = {std::size_t(fft_size[0]),std::size_t(fft_size[1])};
                    rocfft_status result = rocfft_plan_create(&bplan, rocfft_placement_notinplace, 
                                                              rocfft_transform_type_real_inverse, rocfft_precision_double,
                                                              2, lengths, 1, nullptr);
                    assert_rocfft_status("rocfft_plan_create", result);
#elif (AMREX_SPACEDIM == 3)
                    const std::size_t lengths[] = {std::size_t(fft_size[0]),std::size_t(fft_size[1]),std::size_t(fft_size[2])};
                    rocfft_status result = rocfft_plan_create(&bplan, rocfft_placement_notinplace, 
                                                              rocfft_transform_type_real_inverse, rocfft_precision_double,
                                                              3, lengths, 1, nullptr);
                    assert_rocfft_status("rocfft_plan_create", result);
#endif
                }
#else // host

                if (is_flattened) {
#if (AMREX_SPACEDIM == 2)
                    bplan = fftw_plan_dft_c2r_1d(fft_size[0],
                                                 reinterpret_cast<FFTcomplex*>
                                                 (spectral_field.back()->dataPtr()),
                                                 variables_onegrid[mfi].dataPtr(),
                                                 FFTW_ESTIMATE);
#elif (AMREX_SPACEDIM == 3)
                    bplan = fftw_plan_dft_c2r_2d(fft_size[1], fft_size[0],
                                                 reinterpret_cast<FFTcomplex*>
                                                 (spectral_field.back()->dataPtr()),
                                                 variables_onegrid[mfi].dataPtr(),
                                                 FFTW_ESTIMATE);
#endif
                } else {
#if (AMREX_SPACEDIM == 2)
                    bplan = fftw_plan_dft_c2r_2d(fft_size[1], fft_size[0],
                                                 reinterpret_cast<FFTcomplex*>
                                                 (spectral_field.back()->dataPtr()),
                                                 variables_onegrid[mfi].dataPtr(),
                                                 FFTW_ESTIMATE);
#elif (AMREX_SPACEDIM == 3)
                    bplan = fftw_plan_dft_c2r_3d(fft_size[2], fft_size[1], fft_size[0],
                                                 reinterpret_cast<FFTcomplex*>
                                                 (spectral_field.back()->dataPtr()),
                                                 variables_onegrid[mfi].dataPtr(),
                                                 FFTW_ESTIMATE);
#endif
                }
#endif

                backward_plan.push_back(bplan);
            }
	        
            built_plan = true;
        
        } // end MFITer

        ParallelDescriptor::Barrier();

        // InverseTransform
        for (MFIter mfi(variables_onegrid); mfi.isValid(); ++mfi) {
            int i = mfi.LocalIndex();
#ifdef AMREX_USE_CUDA
            cufftSetStream(backward_plan[i], amrex::Gpu::gpuStream());
            cufftResult result = cufftExecZ2D(backward_plan[i],
                                              reinterpret_cast<FFTcomplex*>
                                                  (spectral_field[i]->dataPtr()),
                                              variables_onegrid[mfi].dataPtr());
            if (result != CUFFT_SUCCESS) {
                amrex::AllPrint() << " forward transform using cufftExec failed! Error: "
                                  << cufftErrorToString(result) << "\n";
	        }
#elif AMREX_USE_HIP
            rocfft_execution_info execinfo = nullptr;
            rocfft_status result = rocfft_execution_info_create(&execinfo);
            assert_rocfft_status("rocfft_execution_info_create", result);

            std::size_t buffersize = 0;
            result = rocfft_plan_get_work_buffer_size(backward_plan[i], &buffersize);
            assert_rocfft_status("rocfft_plan_get_work_buffer_size", result);

            void* buffer = amrex::The_Arena()->alloc(buffersize);
            result = rocfft_execution_info_set_work_buffer(execinfo, buffer, buffersize);
            assert_rocfft_status("rocfft_execution_info_set_work_buffer", result);

            result = rocfft_execution_info_set_stream(execinfo, amrex::Gpu::gpuStream());
            assert_rocfft_status("rocfft_execution_info_set_stream", result);

	        amrex::Real* variables_onegrid_ptr = variables_onegrid[mfi].dataPtr();
	        FFTcomplex* spectral_field_ptr = reinterpret_cast<FFTcomplex*>(spectral_field[i]->dataPtr());
            result = rocfft_execute(backward_plan[i],
                                    (void**) &spectral_field_ptr, // in
                                    (void**) &variables_onegrid_ptr, // out
                                    execinfo);
            assert_rocfft_status("rocfft_execute", result);
            amrex::Gpu::streamSynchronize();
            amrex::The_Arena()->free(buffer);
            result = rocfft_execution_info_destroy(execinfo);
            assert_rocfft_status("rocfft_execution_info_destroy", result);
#else
            fftw_execute(backward_plan[i]);
#endif
        }

        variables_onegrid.mult(1.0/sqrtnpts);
        variables.ParallelCopy(variables_onegrid,0,comp,1);
    }

    // destroy fft plan
    for (int i = 0; i < backward_plan.size(); ++i) {
#ifdef AMREX_USE_CUDA
        cufftDestroy(backward_plan[i]);
#elif AMREX_USE_HIP
        rocfft_plan_destroy(backward_plan[i]);
#else
        fftw_destroy_plan(backward_plan[i]);
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

    For each direction d, shift the data in the +d direction by n_cell[d]/2 and then modulo by n_cell[d].

    e.g. for an 8^3 domain

    Cell (7,2,3) is shifted to ( (7+4)%8, (2+4)%8, (3+4)%8 ) =  (3,6,7)

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

      if (zero_avg == 1) {
	  for (MFIter mfi(dft_onegrid); mfi.isValid(); ++mfi) {
	      const Box& bx = mfi.tilebox();
	      const Array4<Real>& dft_temp = dft_onegrid_temp.array(mfi);
	      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept                                                                                                
              {                                                                                                                                                                         
  		  if (i == 0 && j == 0 && k == 0) {
		      dft_temp(i,j,k) = 0.;
		  }
	      });                                                                                                                                                                       
	  }
      }
    
    // Shift DFT by N/2+1 (pi)
    for (MFIter mfi(dft_onegrid); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real>& dft = dft_onegrid.array(mfi);
        const Array4<Real>& dft_temp = dft_onegrid_temp.array(mfi);

        int nx = bx.length(0);
        int nxh = nx/2;
        int ny = bx.length(1);
        int nyh = ny/2;
#if (AMREX_SPACEDIM == 3)
        int nz = bx.length(2);
        int nzh = nz/2;
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
void StructFact::IntegratekShells(const int& step, const Geometry& /*geom*/, const std::string& name) {

    BL_PROFILE_VAR("StructFact::IntegratekShells",IntegratekShells);

    GpuArray<int,AMREX_SPACEDIM> center;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        center[d] = n_cells[d]/2;
    }

    //int npts = n_cells[0]/2-1;
    int npts = n_cells[0]/2;
    //int npts_sq = npts*npts;

    Gpu::DeviceVector<Real> phisum_device(npts);
    Gpu::DeviceVector<int>  phicnt_device(npts);

    Gpu::HostVector<Real> phisum_host(npts);
    
    Real* phisum_ptr = phisum_device.dataPtr();  // pointer to data
    int*  phicnt_ptr = phicnt_device.dataPtr();  // pointer to data

    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
      phisum_ptr[d] = 0.;
      phicnt_ptr[d] = 0;
    });
    

//    Gpu::DeviceVector<Real> phisum_device_large(npts_sq);
//    Gpu::DeviceVector<int>  phicnt_device_large(npts_sq);

//    for (int d=0; d<npts_sq; ++d) {
//        phisum_device_large[d] = 0.;
//        phicnt_device_large[d] = 0;
//    }
    
  //  Real* phisum_large_gpu = phisum_device_large.dataPtr();  // pointer to data
 //   int*  phicnt_large_gpu = phicnt_device_large.dataPtr();  // pointer to data

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
		  amrex::HostDevice::Atomic::Add(&(phisum_ptr[cell]), cov(i,j,k,d));
//		    phisum_large_gpu[idist]  += cov(i,j,k,d);
                }
		amrex::HostDevice::Atomic::Add(&(phicnt_ptr[cell]),1);
 //               ++phicnt_large_gpu[idist];
            }
        });
    }
        
    for (int d=1; d<npts; ++d) {
        ParallelDescriptor::ReduceRealSum(phisum_device[d]);
        ParallelDescriptor::ReduceIntSum(phicnt_device[d]);
    }
        
#if 0
    for (int d=1; d<npts_sq; ++d) {
        ParallelDescriptor::ReduceRealSum(phisum_device_large[d]);
        ParallelDescriptor::ReduceIntSum(phicnt_device_large[d]);
    }

    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream turb_disc;
        std::string turbNamedisc = "turb_disc";
        turbNamedisc += std::to_string(step);
        turbNamedisc += ".txt";
        
        turb_disc.open(turbNamedisc);
        for (int d=1; d<npts_sq; ++d) {
	    if(phicnt_device_large[d]>0) {
		Real dreal = d;
                turb_disc << sqrt(dreal) << " " << 4.*M_PI*d*phisum_device_large[d]/phicnt_device_large[d] << std::endl;
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
            turb_alt << d << " " << phisum_device[d] << std::endl;
        }
    }
#endif

    Real dk = 1.;
    
#if (AMREX_SPACEDIM == 2)
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
      if (d != 0) {
  	  // phisum_ptr[d] *= 2.*M_PI*d*dk*dk/phicnt_ptr[d];
	  phisum_ptr[d] *= 2.*M_PI*(d*dk+.5*dk*dk)/phicnt_ptr[d];
      }
    });
#else
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0) {
	    // phisum_ptr[d] *= 4.*M_PI*(d*d)*dk*dk*dk/phicnt_ptr[d];
	    // phisum_ptr[d] *= 4.*M_PI*(d*d*dk+d*dk*dk+dk*dk*dk/3.)/phicnt_ptr[d];
	    phisum_ptr[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_ptr[d];
	}
    });
#endif

    Gpu::copy(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(), phisum_host.begin());
    
    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream turb;
        std::string turbBaseName = "turb_";
        turbBaseName += name;
        //turbName += std::to_string(step);
        std::string turbName = Concatenate(turbBaseName,step,7);
        turbName += ".txt";
        
        turb.open(turbName);
        for (int d=1; d<npts; ++d) {
            turb << d << " " << phisum_host[d] << std::endl;
        }
    }
}
    
void StructFact::IntegratekShellsDecomp(const int& step, 
                                        const amrex::Geometry& /*geom*/, 
                                        const std::string& name_sol, 
                                        const std::string& name_dil)
{
    BL_PROFILE_VAR("StructFact::IntegratekShellsDecomp",IntegratekShellsDecomp);

    GpuArray<int,AMREX_SPACEDIM> center;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        center[d] = n_cells[d]/2;
    }

    //int npts = n_cells[0]/2-1;
    int npts = n_cells[0]/2;
    //int npts_sq = npts*npts;

    Gpu::DeviceVector<Real> phisum_sol_device(npts);
    Gpu::DeviceVector<Real> phisum_dil_device(npts);
    Gpu::DeviceVector<int>  phicnt_device(npts);

    Gpu::HostVector<Real> phisum_sol_host(npts);
    Gpu::HostVector<Real> phisum_dil_host(npts);
    
    Real* phisum_sol_ptr = phisum_sol_device.dataPtr();  // pointer to data
    Real* phisum_dil_ptr = phisum_dil_device.dataPtr();  // pointer to data
    int*  phicnt_ptr = phicnt_device.dataPtr();  // pointer to data

    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
      phisum_sol_ptr[d] = 0.;
      phisum_dil_ptr[d] = 0.;
      phicnt_ptr[d] = 0;
    });
    
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
            dist = std::sqrt(dist);
            
            if ( dist <= center[0]-0.5) {
	            dist = dist+0.5;
                int cell = int(dist);
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
		            amrex::HostDevice::Atomic::Add(&(phisum_sol_ptr[cell]), cov(i,j,k,d));
		            amrex::HostDevice::Atomic::Add(&(phisum_dil_ptr[cell]), cov(i,j,k,d+3));
                }
		        amrex::HostDevice::Atomic::Add(&(phicnt_ptr[cell]),1);
            }
        });
    }
        
    for (int d=1; d<npts; ++d) {
        ParallelDescriptor::ReduceRealSum(phisum_sol_device[d]);
        ParallelDescriptor::ReduceRealSum(phisum_dil_device[d]);
        ParallelDescriptor::ReduceIntSum(phicnt_device[d]);
    }

    Real dk = 1.;
    
#if (AMREX_SPACEDIM == 2)
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0) {
	        phisum_sol_ptr[d] *= 2.*M_PI*(d*dk+.5*dk*dk)/phicnt_ptr[d];
	        phisum_dil_ptr[d] *= 2.*M_PI*(d*dk+.5*dk*dk)/phicnt_ptr[d];
        }
    });
#else
    amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
    {
        if (d != 0) {
	        phisum_sol_ptr[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_ptr[d];
	        phisum_dil_ptr[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_ptr[d];
	    }
    });
#endif

    Gpu::copy(Gpu::deviceToHost, phisum_sol_device.begin(), 
              phisum_sol_device.end(), phisum_sol_host.begin());

    Gpu::copy(Gpu::deviceToHost, phisum_dil_device.begin(), 
              phisum_dil_device.end(), phisum_dil_host.begin());
    
    if (ParallelDescriptor::IOProcessor()) {
        {
            std::ofstream turb;
            std::string turbBaseName = "turb_";
            turbBaseName += name_sol;
            std::string turbName = Concatenate(turbBaseName,step,7);
            turbName += ".txt";
            
            turb.open(turbName);
            for (int d=1; d<npts; ++d) {
                turb << d << " " << phisum_sol_host[d] << std::endl;
            }
        }
        {
            std::ofstream turb;
            std::string turbBaseName = "turb_";
            turbBaseName += name_dil;
            std::string turbName = Concatenate(turbBaseName,step,7);
            turbName += ".txt";
            
            turb.open(turbName);
            for (int d=1; d<npts; ++d) {
                turb << d << " " << phisum_dil_host[d] << std::endl;
            }
        }
    }
}

// integrate cov_mag over k shells for scalar qtys
void StructFact::IntegratekShellsScalar(const int& step, const Geometry& /*geom*/, const amrex::Vector< std::string >& names) {

    BL_PROFILE_VAR("StructFact::IntegratekShellsMisc",IntegratekShellsMisc);

    int turbvars = NVAR;

    if (names.size() != turbvars) amrex::Error("StructFact::IntegratekShellsMisc() requires names to be of the same length as NVAR");

    GpuArray<int,AMREX_SPACEDIM> center;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        center[d] = n_cells[d]/2;
    }

    //int npts = n_cells[0]/2-1;
    int npts = n_cells[0]/2;

    Gpu::DeviceVector<Real> phisum_device(npts);
    Gpu::DeviceVector<int>  phicnt_device(npts);

    Gpu::HostVector<Real> phisum_host(npts);
    
    Real* phisum_ptr = phisum_device.dataPtr();  // pointer to data
    int*  phicnt_ptr = phicnt_device.dataPtr();  // pointer to data

    for (int var_ind = 0; var_ind < turbvars; var_ind++) {

        amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
        {
          phisum_ptr[d] = 0.;
          phicnt_ptr[d] = 0;
        });
        
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
                dist = std::sqrt(dist);
                
                if ( dist <= center[0]-0.5) {
                    dist = dist+0.5;
                    int cell = int(dist);
                    amrex::HostDevice::Atomic::Add(&(phisum_ptr[cell]), cov(i,j,k,var_ind));
                    amrex::HostDevice::Atomic::Add(&(phicnt_ptr[cell]),1);
                }
            });
        }
            
        for (int d=1; d<npts; ++d) {
            ParallelDescriptor::ReduceRealSum(phisum_device[d]);
            ParallelDescriptor::ReduceIntSum(phicnt_device[d]);
        }
            
        Real dk = 1.;
    
#if (AMREX_SPACEDIM == 2)
#if 0
        for (int d=1; d<npts; ++d) {
          //  phisum_vect[d] *= 2.*M_PI*d*dk*dk/phicnt_vect[d];
            phisum_vect[d] *= 2.*M_PI*(d*dk+.5*dk*dk)/phicnt_vect[d];
        }
#endif
        amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
        {
          if (d != 0) {
          // phisum_ptr[d] *= 2.*M_PI*d*dk*dk/phicnt_ptr[d];
          phisum_ptr[d] *= 2.*M_PI*(d*dk+.5*dk*dk)/phicnt_ptr[d];
          }
        });
#else
#if 0
        for (int d=1; d<npts; ++d) {
          //  phisum_vect[d] *= 4.*M_PI*(d*d)*dk*dk*dk/phicnt_vect[d];
          //  phisum_vect[d] *= 4.*M_PI*(d*d*dk+d*dk*dk+dk*dk*dk/3.)/phicnt_vect[d];
            phisum_vect[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_vect[d];
        }
#endif
        amrex::ParallelFor(npts, [=] AMREX_GPU_DEVICE (int d) noexcept
        {
            if (d != 0) {
            // phisum_ptr[d] *= 4.*M_PI*(d*d)*dk*dk*dk/phicnt_ptr[d];
            // phisum_ptr[d] *= 4.*M_PI*(d*d*dk+d*dk*dk+dk*dk*dk/3.)/phicnt_ptr[d];
            phisum_ptr[d] *= 4.*M_PI*(d*d*dk+dk*dk*dk/12.)/phicnt_ptr[d];
            }
        });
#endif
        Gpu::copy(Gpu::deviceToHost, phisum_device.begin(), phisum_device.end(), phisum_host.begin());
        
        if (ParallelDescriptor::IOProcessor()) {
            std::ofstream turb;
            std::string turbBaseName = "turb_"+names[var_ind];
            std::string turbName = Concatenate(turbBaseName,step,7);
            turbName += ".txt";
            
            turb.open(turbName);
            for (int d=1; d<npts; ++d) {
                turb << d << " " << phisum_host[d] << std::endl;
            }
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


void StructFact::DecomposeVelFourier(const amrex::MultiFab& vel_dft_real, 
                                     const amrex::MultiFab& vel_dft_imag, 
                                     const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("StructFact::DecomposeVelFourier",DecomposeVelFourier);

    const BoxArray& ba = vel_sol_real.boxArray();
    const DistributionMapping& dm = vel_sol_real.DistributionMap();
    MultiFab dft_real, dft_imag;
    dft_real.define(ba, dm, 3, 0);
    dft_imag.define(ba, dm, 3, 0);
    dft_real.ParallelCopy(vel_dft_real,0,0,3);
    dft_imag.ParallelCopy(vel_dft_imag,0,0,3);
    
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for (MFIter mfi(dft_real); mfi.isValid(); ++mfi) {
        
        Box bx = mfi.fabbox();

        Array4<const Real> const& real = dft_real.array(mfi);
        Array4<const Real> const& imag = dft_imag.array(mfi);
        
        Array4<      Real> const& real_sol = vel_sol_real.array(mfi);
        Array4<      Real> const& imag_sol = vel_sol_imag.array(mfi);
        
        Array4<      Real> const& real_dil = vel_dil_real.array(mfi);
        Array4<      Real> const& imag_dil = vel_dil_imag.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int nx = bx.length(0);
            int ny = bx.length(1);
            int nz = bx.length(2);

            Real GxR, GxC, GyR, GyC, GzR, GzC;
            
            if (i <= bx.length(0)/2) { 
                // Gradient Operators
                GxR = (cos(2.0*M_PI*i/nx)-1.0)/dx[0];
                GxC = (sin(2.0*M_PI*i/nx)-0.0)/dx[0];
                GyR = (cos(2.0*M_PI*j/ny)-1.0)/dx[1];
                GyC = (sin(2.0*M_PI*j/ny)-0.0)/dx[1];
                GzR = (cos(2.0*M_PI*k/nz)-1.0)/dx[2];
                GzC = (sin(2.0*M_PI*k/nz)-0.0)/dx[2];
            }
            else { // conjugate
                // Gradient Operators
                GxR = (cos(2.0*M_PI*(nx-i)/nx)-1.0)/dx[0];
                GxC = (sin(2.0*M_PI*(nx-i)/nx)-0.0)/dx[0];
                GyR = (cos(2.0*M_PI*(ny-j)/ny)-1.0)/dx[1];
                GyC = (sin(2.0*M_PI*(ny-j)/ny)-0.0)/dx[1];
                GzR = (cos(2.0*M_PI*(nz-k)/nz)-1.0)/dx[2];
                GzC = (sin(2.0*M_PI*(nz-k)/nz)-0.0)/dx[2];
            }

            // Inverse Laplacian
            Real Lap = GxR*GxR + GxC*GxC + GyR*GyR + GyC*GyC + GzR*GzR + GzC*GzC;

            // Divergence of vel
            Real divR = real(i,j,k,0)*GxR - imag(i,j,k,0)*GxC +
                        real(i,j,k,1)*GyR - imag(i,j,k,1)*GyC +
                        real(i,j,k,2)*GzR - imag(i,j,k,2)*GzC ;
            Real divC = real(i,j,k,0)*GxC + imag(i,j,k,0)*GxR +
                        real(i,j,k,1)*GyC + imag(i,j,k,1)*GyR +
                        real(i,j,k,2)*GzC + imag(i,j,k,2)*GzR ;

            if (Lap < 1.0e-12) { // zero mode for no bulk motion
                real_dil(i,j,k,0) = 0.0;
                real_dil(i,j,k,1) = 0.0;
                real_dil(i,j,k,2) = 0.0;
                imag_dil(i,j,k,0) = 0.0;
                imag_dil(i,j,k,1) = 0.0;
                imag_dil(i,j,k,2) = 0.0;
            }
            else {
                // Dilatational velocity 
                real_dil(i,j,k,0) = (divR*GxR + divC*GxC) / Lap;
                real_dil(i,j,k,1) = (divR*GyR + divC*GyC) / Lap;
                real_dil(i,j,k,2) = (divR*GzR + divC*GzC) / Lap;
                imag_dil(i,j,k,0) = (divC*GxR - divR*GxC) / Lap;
                imag_dil(i,j,k,1) = (divC*GyR - divR*GyC) / Lap;
                imag_dil(i,j,k,2) = (divC*GzR - divR*GzC) / Lap;
                
                // Solenoidal velocity
                real_sol(i,j,k,0) = real(i,j,k,0) - real_dil(i,j,k,0);
                real_sol(i,j,k,1) = real(i,j,k,1) - real_dil(i,j,k,1); 
                real_sol(i,j,k,2) = real(i,j,k,2) - real_dil(i,j,k,2);
                imag_sol(i,j,k,0) = imag(i,j,k,0) - imag_dil(i,j,k,0);
                imag_sol(i,j,k,1) = imag(i,j,k,1) - imag_dil(i,j,k,1);
                imag_sol(i,j,k,2) = imag(i,j,k,2) - imag_dil(i,j,k,2);
            }
        });
    }
}

void StructFact::GetDecompVel(MultiFab& vel_decomp, const Geometry& geom)
{
    BL_PROFILE_VAR("StructFact::GetDecompVel()", GetDecompVel);
    
    if (!decompose) 
        amrex::Error("StructFact::GetDecompVel() is specific for vel decomposition in turbulence");

    const BoxArray& ba_in = vel_decomp.boxArray();
    const DistributionMapping& dmap_in = vel_decomp.DistributionMap();

    MultiFab vel;
    vel.define(ba_in, dmap_in, 3, 0);

    const BoxArray& ba = vel_sol_real.boxArray();
    const DistributionMapping& dm = vel_sol_real.DistributionMap();
    MultiFab dft_real, dft_imag;
    dft_real.define(ba, dm, 3, 0);
    dft_imag.define(ba, dm, 3, 0);
    
    dft_real.ParallelCopy(vel_sol_real,0,0,3);
    dft_imag.ParallelCopy(vel_sol_imag,0,0,3);

    InverseFFT(vel, dft_real, dft_imag, geom);
    vel_decomp.ParallelCopy(vel,0,0,3);

    dft_real.ParallelCopy(vel_dil_real,0,0,3);
    dft_imag.ParallelCopy(vel_dil_imag,0,0,3);

    InverseFFT(vel, dft_real, dft_imag, geom);
    vel_decomp.ParallelCopy(vel,0,3,3);

}
