
#include "common_functions.H"
#include "StructFact.H"

#include <AMReX_MultiFabUtil.H>
#include "AMReX_PlotFileUtil.H"
#include "AMReX_BoxArray.H"

StructFact::StructFact()
{}

StructFact::StructFact(const BoxArray& ba_in, const DistributionMapping& dmap_in,
		       const Vector< std::string >& var_names,
		       const Vector< Real >& var_scaling_in,
		       const Vector< int >& s_pairA_in,
		       const Vector< int >& s_pairB_in,
		       const int& verbosity_in) {

  BL_PROFILE_VAR("StructFact::StructFact()",StructFact);

  verbosity = verbosity_in;
  
  if (s_pairA_in.size() != s_pairA_in.size())
        amrex::Error("StructFact::StructFact() - Must have an equal number of components");

  NVAR = var_names.size();

  //////////////////////////////////////////////////////
  // Reorder pair selecting vectors & error checking

  NCOV = s_pairA_in.size();

  if ( NCOV != var_scaling_in.size() )
      amrex::Error("StructFact::StructFact() Constructor 1 - NCOV != var_scaling_in.size()");

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

  for (int n=0; n<NVARU; n++) {
    Print() << "HACK 1: vector (" << n << ") = " << varu_temp[n] << std::endl;
  }

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
    Print() << "HACK 2: pairs (" << n << ") = " << s_pairA[n] << ", " << s_pairB[n] << std::endl;
  }
  Print() << "HACK: NCOV = " << NCOV << std::endl;

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

StructFact::StructFact(const BoxArray& ba_in, const DistributionMapping& dmap_in,
		       const Vector< std::string >& var_names,
		       const Vector< Real >& var_scaling_in,
		       const int& verbosity_in) {
  
  BL_PROFILE_VAR("StructFact::StructFact()",StructFact);

  verbosity = verbosity_in;
  
  NVAR = var_names.size();
  NCOV = NVAR*(NVAR+1)/2;

  if ( NCOV != var_scaling_in.size() )
      amrex::Error("StructFact::StructFact() Constructor 2 -  NCOV != var_scaling_in.size()");

  scaling.resize(NCOV);
  for (int n=0; n<NCOV; n++) {
      scaling[n] = 1.0/var_scaling_in[n];
  }
  
  s_pairA.resize(NCOV);
  s_pairB.resize(NCOV);
  
  // all variables are selected in this constructor
  NVARU = NVAR;
  var_u.resize(NVARU);
  for (int n=0; n<NVARU; n++) {
    var_u[n] = n;
  }
  
  int index = 0;
  for (int j=0; j<NVAR; j++) {
    for (int i=j; i<NVAR; i++) {
      s_pairA[index] = i;
      s_pairB[index] = j;
      index++;
    }
  }

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

void StructFact::define(const BoxArray& ba_in, const DistributionMapping& dmap_in,
                        const Vector< std::string >& var_names,
                        const Vector< Real >& var_scaling_in,
                        const int& verbosity_in) {
  
  BL_PROFILE_VAR("StructFact::define()",StructFactDefine);

  verbosity = verbosity_in;
  
  NVAR = var_names.size();
  NCOV = NVAR*(NVAR+1)/2;

  if ( NCOV != var_scaling_in.size() )
      amrex::Error("StructFact::define() -  NCOV != var_scaling_in.size()");

  scaling.resize(NCOV);
  for (int n=0; n<NCOV; n++) {
      scaling[n] = 1.0/var_scaling_in[n];
  }
  
  s_pairA.resize(NCOV);
  s_pairB.resize(NCOV);
  
  // all variables are selected in this constructor
  NVARU = NVAR;
  var_u.resize(NVARU);
  for (int n=0; n<NVARU; n++) {
    var_u[n] = n;
  }
  
  int index = 0;
  for (int j=0; j<NVAR; j++) {
    for (int i=j; i<NVAR; i++) {
      s_pairA[index] = i;
      s_pairB[index] = j;
      index++;
    }
  }

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
                               const int& fft_type_in,
                               const int& reset) {

  BL_PROFILE_VAR("StructFact::FortStructure()",FortStructure);

  const BoxArray& ba = variables.boxArray();
  const DistributionMapping& dm = variables.DistributionMap();
  MultiFab variables_dft_real, variables_dft_imag;
  variables_dft_real.define(ba, dm, NVAR, 0);
  variables_dft_imag.define(ba, dm, NVAR, 0);

  if (fft_type_in == 1) {
      Print() << "Using FFTW\n";
      ComputeFFTW(variables, variables_dft_real, variables_dft_imag, geom);
  }
  else if (ba.size() == ParallelDescriptor::NProcs()) {
      Print() << "Using SWFFT\n";
      ComputeSWFFT(variables, variables_dft_real, variables_dft_imag, geom);
  }
  else {
      Print() << "Using SWFFT\n";
      BoxArray ba_temp(geom.Domain());

      ba_temp.maxSize(IntVect(max_grid_size_structfact));

      if (ba_temp.size() != ParallelDescriptor::NProcs()) {
          Print() << "StructFact::FortStructure - number of MPI ranks needs to match the number of grids;\n"
                  << "If this is a full-dimensional dataset, use max_grid_size_structfact;\n"
                  << "If this is a vertically-averaged dataset, use max_grid_projection.\n";
          Abort("");
          exit(0);
      }

      DistributionMapping dm_temp(ba_temp);

      // create variables_temp, variables_dft_real_temp and variables_dft_imag
      // these will have the same number of grids as MPI ranks so they need a different
      // BoxArray and DistributionMapping
      MultiFab variables_temp         (ba_temp, dm_temp, variables.nComp(), 0);
      MultiFab variables_dft_real_temp(ba_temp, dm_temp, NVAR, 0);
      MultiFab variables_dft_imag_temp(ba_temp, dm_temp, NVAR, 0);

      // ParallelCopy variables into variables_temp
      variables_temp.ParallelCopy(variables, 0, 0, variables.nComp());

      ComputeSWFFT(variables_temp, variables_dft_real_temp, variables_dft_imag_temp, geom);

      // ParallelCopy variables_dft_real_temp into variables_dft_real
      // ParallelCopy variables_dft_imag_temp into variables_dft_imag
      variables_dft_real.ParallelCopy(variables_dft_real_temp, 0, 0, NVAR);
      variables_dft_imag.ParallelCopy(variables_dft_imag_temp, 0, 0, NVAR);
  }

  MultiFab cov_temp;
  cov_temp.define(ba, dm, 1, 0);

  int index = 0;
  for (int n=0; n<NCOV; n++) {
    int i = s_pairA[n];
    int j = s_pairB[n];
    
    // Compute temporary real and imaginary components of covariance

    // Real component of covariance
    cov_temp.setVal(0.0);
    MultiFab::AddProduct(cov_temp,variables_dft_real,i,variables_dft_real,j,0,1,0);
    MultiFab::AddProduct(cov_temp,variables_dft_imag,i,variables_dft_imag,j,0,1,0);

    if (reset == 1) {
        MultiFab::Copy(cov_real,cov_temp,0,index,1,0);
    } else {
        MultiFab::Add(cov_real,cov_temp,0,index,1,0);
    }

    // Imaginary component of covariance
    cov_temp.setVal(0.0);
    MultiFab::AddProduct(cov_temp,variables_dft_imag,i,variables_dft_real,j,0,1,0);
    cov_temp.mult(-1.0,0);
    MultiFab::AddProduct(cov_temp,variables_dft_real,i,variables_dft_imag,j,0,1,0);

    if (reset == 1) {
        MultiFab::Copy(cov_imag,cov_temp,0,index,1,0);
    } else {
        MultiFab::Add(cov_imag,cov_temp,0,index,1,0);
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

    cov_real.setVal(0.);
    cov_imag.setVal(0.);
    nsamples = 0;
    
}

void StructFact::ComputeSWFFT(const MultiFab& variables,
                              MultiFab& variables_dft_real, 
                              MultiFab& variables_dft_imag,
                              const Geometry& geom) {

  BL_PROFILE_VAR("StructFact::ComputeSWFFT()", ComputeSWFFT);

  Box domain(geom.Domain());
  const BoxArray& ba = variables.boxArray();
  DistributionMapping dm = variables.DistributionMap();
  
  if (verbosity > 1) {
    amrex::Print() << "BA " << ba << std::endl;
  }

  if (variables_dft_real.nGrow() != 0 || variables.nGrow() != 0) {
    amrex::Error("StructFact::ComputeSWFFT() - Current implementation requires that both variables_temp[0] and variables_dft_real[0] have no ghost cells");
  }

  // We assume that all grids have the same size hence 
  // we have the same nx,ny,nz on all ranks
  int nx = ba[0].size()[0];
  int ny = ba[0].size()[1];
#if (AMREX_SPACEDIM == 2)
  int nz = 1;
#elif (AMREX_SPACEDIM == 3)
  int nz = ba[0].size()[2];
#endif

  int nbx = domain.length(0) / nx;
  int nby = domain.length(1) / ny;
#if (AMREX_SPACEDIM == 2)
  int nbz = 1;
#elif (AMREX_SPACEDIM == 3)
  int nbz = domain.length(2) / nz;
#endif
  int nboxes = nbx * nby * nbz;
  if (verbosity > 1) {
    amrex::Print() << "nx, ny, nz:\t" << nx << ", " << ny << ", " << nz << std::endl;
    amrex::Print() << "Number of boxes:\t" << nboxes << "\tBA size:\t" << ba.size() << std::endl;
  }
  if (nboxes != ba.size())
    amrex::Error("StructFact::ComputeSWFFT() - NBOXES NOT COMPUTED CORRECTLY");

  Vector<int> rank_mapping;
  rank_mapping.resize(nboxes);

  DistributionMapping dmap = variables_dft_real.DistributionMap();

//  Print() << "HACK FFT: " << ba << std::endl << dmap << std::endl << dm << std::endl;

  for (int ib = 0; ib < nboxes; ++ib)
    {
      int i = ba[ib].smallEnd(0) / nx;
      int j = ba[ib].smallEnd(1) / ny;
#if (AMREX_SPACEDIM == 2)
      int k = 0;
#elif (AMREX_SPACEDIM == 3)
      int k = ba[ib].smallEnd(2) / nz;
#endif

      // This would be the "correct" local index if the data wasn't being transformed
      int local_index = k*nbx*nby + j*nbx + i;

      // This is what we [would] pass to dfft to compensate for the Fortran ordering
      //      of amrex data in MultiFabs.
      // int local_index = i*nby*nbz + j*nbz + k;

      rank_mapping[local_index] = dmap[ib];
      if (verbosity > 0)
      	amrex::Print() << "LOADING RANK NUMBER " << dmap[ib] << " FOR GRID NUMBER " << ib 
      		       << " WHICH IS LOCAL NUMBER " << local_index << std::endl;
    }

  // FIXME: Assumes same grid spacing

  // Assume for now that nx = ny = nz
#if (AMREX_SPACEDIM == 2)
  int Ndims[3] = { 1, nby, nbx};
  int     n[3] = { 1, domain.length(1), domain.length(0)};
#elif (AMREX_SPACEDIM == 3)
  int Ndims[3] = { nbz, nby, nbx };
  int     n[3] = { domain.length(2), domain.length(1), domain.length(0)};
#endif
  MPI_Comm comm = ParallelDescriptor::Communicator();
  hacc::Distribution d(comm,n,Ndims,&rank_mapping[0]);

  if (verbosity > 0) {
    Print() << "RANK MAPPING: \n";
    for (int i=0; i<rank_mapping.size(); i++) {
      Print() << "\t" << rank_mapping[i] << std::endl;
    }
  }

  for (int comp=0; comp<NVAR; comp++) {

    bool comp_fft = false;
    for (int i=0; i<NVARU; i++) {
      if (comp == var_u[i]) {
	comp_fft = true;
	break;
      }
    }

    if(comp_fft) {
   
      for (MFIter mfi(variables_dft_real,false); mfi.isValid(); ++mfi) {

        std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;
	std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b;

	a.resize(nx*ny*nz);
	b.resize(nx*ny*nz);

	// Print() << "HACK FFT: got here" << std::endl;

        hacc::Dfft dfft(d);
	dfft.makePlans(&a[0],&b[0],&a[0],&b[0]);

	// *******************************************
	// Copy real data from Rhs into real part of a -- no ghost cells and
	// put into C++ ordering (not Fortran)
	// *******************************************
	complex_t zero(0.0, 0.0);
	size_t local_indx = 0;
	for(size_t k=0; k<(size_t)nz; k++) {
	  for(size_t j=0; j<(size_t)ny; j++) {
	    for(size_t i=0; i<(size_t)nx; i++) {

	      complex_t temp(variables[mfi].dataPtr(comp)[local_indx],0.);
	      a[local_indx] = temp;
	      local_indx++;

	      // Print() << "HACK FFT: a[" << local_indx << "] = \t " << variables[mfi].dataPtr(comp)[local_indx] << std::endl;

	    }
	  }
	}

	//  *******************************************
	//  Compute the forward transform
	//  *******************************************
	dfft.forward(&a[0]);

	// Redistribute data from z-pencils in k-space back to blocks
	d.redistribute_2_to_3(&a[0],&b[0],2);
	
	// Note: Scaling for inverse FFT
	size_t global_size  = dfft.global_size();
	
	// Real pi = 4.0*std::atan(1.0);
	Real fac = sqrt(1.0 / (Real)global_size);

	local_indx = 0;
	for(size_t k=0; k<(size_t)nz; k++) {
	  for(size_t j=0; j<(size_t)ny; j++) {
	    for(size_t i=0; i<(size_t)nx; i++) {

	      // Divide by 2 pi N
	      variables_dft_real[mfi].dataPtr(comp)[local_indx] = fac * std::real(b[local_indx]);
	      variables_dft_imag[mfi].dataPtr(comp)[local_indx] = fac * std::imag(b[local_indx]);
	      local_indx++;
	    }
	  }
	}
      }
    }

    /*
    for (MFIter mfi(variables_dft_real); mfi.isValid(); ++mfi) {

        Box bx = mfi.fabbox();

        Array4<Real> const& realpart = variables_dft_real.array(mfi);
        Array4<Real> const& imagpart = variables_dft_imag.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            std::cout << "HACKFFT " << i << " " << j << " " << k << " "
                      << realpart(i,j,k,comp) << " + " << imagpart(i,j,k,comp) << "i"
                      << std::endl;
        });
    }
    */
  }

  bool write_data = false;
  if (write_data) {
    std::string plotname = "a_DFT_REAL";
    VisMF::Write(variables_dft_real,plotname);
  }
}


void StructFact::ComputeFFTW(const MultiFab& variables,
                             MultiFab& variables_dft_real, 
                             MultiFab& variables_dft_imag,
                             const Geometry& geom) {

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
    
    using FFTplan = fftw_plan;
    using FFTcomplex = fftw_complex;

    // contain to store FFT - note it is shrunk by "half" in x
    Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > > spectral_field;

    Vector<FFTplan> forward_plan;

    for (int comp=0; comp<NVAR; comp++) {

        bool comp_fft = false;
        for (int i=0; i<NVARU; i++) {
            if (comp == var_u[i]) {
                comp_fft = true;
                break;
            }
        }

        variables_onegrid.ParallelCopy(variables,comp,0,1);

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

            forward_plan.push_back(fplan);
        }

        ParallelDescriptor::Barrier();

        // ForwardTransform
        for (MFIter mfi(variables_onegrid); mfi.isValid(); ++mfi) {
            int i = mfi.LocalIndex();
            fftw_execute(forward_plan[i]);
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
        fftw_destroy_plan(forward_plan[i]);
    }

//    fftw_mpi_cleanup();

}

void StructFact::WritePlotFile(const int step, const Real time, const Geometry& geom,
                               std::string plotfile_base,
                               const int& zero_avg) {
  
  BL_PROFILE_VAR("StructFact::WritePlotFile()",WritePlotFile);

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

  if (turbForcing != 1) {
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
  
}

void StructFact::StructOut(MultiFab& struct_out) {

  BL_PROFILE_VAR("StructFact::StructOut()",StructOut);

  if (struct_out.nComp() == cov_mag.nComp()) {
    MultiFab::Copy(struct_out,cov_mag,0,0,cov_mag.nComp(),0);
  } else {
    amrex::Error("StructFact::StructOut() - Must have an equal number of components");
  }
}

void StructFact::Finalize(MultiFab& cov_real_in, MultiFab& cov_imag_in,
                          const Geometry& geom, const int& zero_avg) {

  BL_PROFILE_VAR("StructFact::Finalize()",Finalize);
  
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

void StructFact::ShiftFFT(MultiFab& dft_out, const Geometry& geom, const int& zero_avg) {

  BL_PROFILE_VAR("StructFact::ShiftFFT()",ShiftFFT);

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

        /*
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            std::cout << "HACKFFTSHIFT " << i << " " << j << " " << k << " "
                      << dft(i,j,k) << std::endl;
        });
        */
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
                    phisum_gpu[cell] += cov(i,j,k,d);
//		    phisum_large_gpu[idist]  += cov(i,j,k,d);
                }
                ++phicnt_gpu[cell];
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
