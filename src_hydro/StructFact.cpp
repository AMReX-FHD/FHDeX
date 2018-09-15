
#include "common_functions.H"
#include "hydro_functions_F.H"
#include "StructFact.H"

#include "AMReX_PlotFileUtil.H"
#include "AMReX_BoxArray.H"

StructFact::StructFact(const BoxArray ba_in, const DistributionMapping dmap_in,
		       const Vector< std::string >& var_names) {
  
  BL_PROFILE_VAR("StructFact::StructFact()",StructFact);

  // static_assert(AMREX_SPACEDIM == 3, "3D only");

  N_VAR = var_names.size();
  N_COV = N_VAR*(N_VAR+1)/2;

  // Note that we are defining with NO ghost cells

  cov_real.resize(N_COV);
  cov_imag.resize(N_COV);
  cov_mag.resize(N_COV);

  for (int d=0; d<N_COV; d++) {
    cov_real[d].define(ba_in, dmap_in, 1, 0);
    cov_imag[d].define(ba_in, dmap_in, 1, 0);
    cov_mag[d].define(ba_in, dmap_in, 1, 0);
    cov_mag[d].setVal(0.0);
    cov_real[d].setVal(0.0);
    cov_imag[d].setVal(0.0);
  }

  cov_names.resize(N_COV);
  std::string x;
  int cnt = 0;
  for (int j=0; j<N_VAR; j++) {
    for (int i=j; i<N_VAR; i++) {
      x = "structFact";
      x += '_';
      x += var_names[j];
      x += '_';
      x += var_names[i];
      cov_names[cnt] = x;
      Print() << var_names[j] << "\t" << var_names[i] << "\t" << cov_names[cnt] << std::endl;
      cnt++;
    }
  }
}

void StructFact::FortStructure(const Vector< MultiFab >& umac, const Geometry geom) {

  BL_PROFILE_VAR("StructFact::FortStructure()",FortStructure);

  const BoxArray& ba = cov_real[0].boxArray();
  const DistributionMapping& dm = cov_real[0].DistributionMap();

  // if (ba.size() != 1) {
  //   Abort("Only adapted for single grid");
  //   exit(0);
  // }
  
  Vector< MultiFab > variables_temp;
  variables_temp.resize(N_VAR);
  for (int d=0; d<N_VAR; d++) {
    variables_temp[d].define(ba, dm, 1, 0);
    // Print() << "HACK: variables_temp # ghost cells = " << variables_temp[d].nGrow() << "\n";
  }
  
  for (int d=0; d<N_VAR; d++) { 
    MultiFab::Copy(variables_temp[d],umac[d],0,0,1,0);
  }

  amrex::Vector< MultiFab > variables_dft_real, variables_dft_imag;
  variables_dft_real.resize(N_VAR);
  variables_dft_imag.resize(N_VAR);
  for (int d=0; d<N_VAR; d++) {
    variables_dft_real[d].define(ba, dm, 1, 0);
    variables_dft_imag[d].define(ba, dm, 1, 0);
  }

  ComputeFFT(variables_temp, variables_dft_real, variables_dft_imag, geom);

  MultiFab cov_real_temp;
  MultiFab cov_imag_temp;
  cov_real_temp.define(ba, dm, 1, 0);
  cov_imag_temp.define(ba, dm, 1, 0);

  int index = 0;
  for (int j=0; j<N_VAR; j++) {
    for (int i=j; i<N_VAR; i++) {

      // Compute temporary real and imaginary components of covariance
      cov_real_temp.setVal(0.0);
      MultiFab::AddProduct(cov_real_temp,variables_dft_real[i],0,variables_dft_real[j],0,0,1,0);
      MultiFab::AddProduct(cov_real_temp,variables_dft_imag[i],0,variables_dft_imag[j],0,0,1,0);
      cov_imag_temp.setVal(0.0);
      MultiFab::AddProduct(cov_imag_temp,variables_dft_imag[i],0,variables_dft_real[j],0,0,1,0);
      cov_imag_temp.mult(-1.0,0);
      MultiFab::AddProduct(cov_imag_temp,variables_dft_real[i],0,variables_dft_imag[j],0,0,1,0);
      
      MultiFab::Add(cov_real[index],cov_real_temp,0,0,1,0);
      MultiFab::Add(cov_imag[index],cov_imag_temp,0,0,1,0);
      
      index++;
    }
  }

  bool write_data = false;
  if (write_data) {
    index = 0;
    std::string plotname = "a_COV"; 
    std::string x;
    for (int j=0; j<N_VAR; j++) {
      for (int i=j; i<N_VAR; i++) {
	x = plotname;
	x += '_';
	x += '0' + i;
	x += '0' + j;
	VisMF::Write(cov_real[index],x);
	index++;
      }
    }
  }

  nsamples++;
}

void StructFact::WritePlotFile(const int step, const Real time, const Geometry geom,
			       const Real scale) {
  
  BL_PROFILE_VAR("StructFact::WritePlotFile()",WritePlotFile);

  Real nsamples_inv = 1.0/(double)nsamples;
  
  amrex::Vector< MultiFab > struct_temp;
  struct_temp.resize(N_COV);
  for (int d=0; d<N_COV; d++) {
    struct_temp[d].define(cov_mag[0].boxArray(), cov_mag[0].DistributionMap(), 1, 0);
    struct_temp[d].setVal(0.0);
  }

  StructFinalize();
  for(int d=0; d<N_COV; d++) {   
    MultiFab::Copy(struct_temp[d],cov_mag[d],0,0,1,0);
    struct_temp[d].mult(nsamples_inv,0);
    struct_temp[d].mult(scale,0);
  }

  // Print() << "number of samples = " << nsamples << std::endl;

  ShiftFFT(struct_temp);

  //////////////////////////////////////////////////////////////////////////////////
  // Write out structure factor magnitude to plot file
  //////////////////////////////////////////////////////////////////////////////////
  const std::string plotfilename = "plt_structure_factor";
  int nPlot = N_COV;
  MultiFab plotfile(cov_mag[0].boxArray(), cov_mag[0].DistributionMap(), nPlot, 0);
  Vector<std::string> varNames(nPlot);

  // keep a counter for plotfile variables
  int cnt = 0;

  for (int j=0; j<N_VAR; j++) {
    for (int i=j; i<N_VAR; i++) {
      varNames[cnt++] = cov_names[cnt];
    }
  }

  // reset plotfile variable counter
  cnt = 0;

  // copy structure factor into plotfile
  for (int d=0; d<N_COV; ++d) {
    MultiFab::Copy(plotfile, struct_temp[cnt], 0, cnt, 1, 0);
    cnt++;
  }

  // write a plotfile
  WriteSingleLevelPlotfile(plotfilename,plotfile,varNames,geom,time,step);
  
  //////////////////////////////////////////////////////////////////////////////////
  // Write out real and imaginary components of structure factor to plot file
  //////////////////////////////////////////////////////////////////////////////////
  const std::string plotfilenameReIm = "plt_structure_factor_ReIm";
  nPlot = 2*N_COV;
  MultiFab plotfileReIm(cov_mag[0].boxArray(), cov_mag[0].DistributionMap(), nPlot, 0);
  varNames.resize(nPlot);
  struct_temp.resize(nPlot);
  for (int d=0; d<nPlot; d++) {
    struct_temp[d].define(cov_mag[0].boxArray(), cov_mag[0].DistributionMap(), 1, 0);
    struct_temp[d].setVal(0.0);
  }

  // keep a counter for plotfile variables
  cnt = 0;
  
  int index = 0;
  for (int j=0; j<N_VAR; j++) {
    for (int i=j; i<N_VAR; i++) {
      varNames[cnt] = cov_names[cnt];
      varNames[cnt] += "_real";
      index++;
      cnt++;
    }
  }
  
  index = 0;
  for (int j=0; j<N_VAR; j++) {
    for (int i=j; i<N_VAR; i++) {
      varNames[cnt] = cov_names[index];
      varNames[cnt] += "_imag";
      index++;
      cnt++;
    }
  }

  // reset plotfile variable counter
  cnt = 0;

  for(int d=0; d<N_COV; d++) {   
    MultiFab::Copy(struct_temp[cnt],cov_real[d],0,0,1,0);
    struct_temp[cnt].mult(nsamples_inv,0);
    struct_temp[cnt].mult(scale,0);
    cnt++;
  }

  for(int d=0; d<N_COV; d++) {   
    MultiFab::Copy(struct_temp[cnt],cov_imag[d],0,0,1,0);
    struct_temp[cnt].mult(nsamples_inv,0);
    struct_temp[cnt].mult(scale,0);
    cnt++;
  }

  ShiftFFT(struct_temp);

  // reset plotfile variable counter
  cnt = 0;

  // copy structure factor into plotfile
  for (int d=0; d<nPlot; ++d) {
    MultiFab::Copy(plotfileReIm, struct_temp[cnt], 0, cnt, 1, 0);
    cnt++;
  }

  // write a plotfile
  WriteSingleLevelPlotfile(plotfilenameReIm,plotfileReIm,varNames,geom,time,step);
}

void StructFact::StructOut(amrex::Vector< MultiFab >& struct_out, const amrex::Real scale) {

  BL_PROFILE_VAR("StructFact::StructOut()",StructOut);
  
  struct_out.resize(N_COV);
  for (int d=0; d<N_COV; d++) {
    struct_out[d].define(cov_mag[0].boxArray(), cov_mag[0].DistributionMap(), 1, 0);
    struct_out[d].setVal(0.0);
  }

  StructFinalize();
  for(int d=0; d<N_COV; d++) {   
    MultiFab::Copy(struct_out[d],cov_mag[d],0,0,1,0);
  }

  ShiftFFT(struct_out);
}

void StructFact::StructFinalize(const amrex::Real scale) {

  for(int d=0; d<N_COV; d++) {
    cov_mag[d].setVal(0.0);
    MultiFab::AddProduct(cov_mag[d],cov_real[d],0,cov_real[d],0,0,1,0);
    MultiFab::AddProduct(cov_mag[d],cov_imag[d],0,cov_imag[d],0,0,1,0);
  }

  SqrtMF(cov_mag);

  for(int d=0; d<N_COV; d++) {
    cov_mag[d].mult(scale,0);
  }
}

void StructFact::ComputeFFT(const Vector< MultiFab >& variables_temp,
			    Vector< MultiFab >& variables_dft_real, 
			    Vector< MultiFab >& variables_dft_imag,
			    const Geometry geom) {

  BL_PROFILE_VAR("StructFact::ComputeFFT()",ComputeFFT);

  Box domain(geom.Domain());
  
  // IF
  const BoxArray& ba = variables_dft_real[0].boxArray();
  // const BoxArray& ba_in = variables_dft_real[0].boxArray();
  // BoxArray ba = ba_in;
  // ba.resize(1);
  // ba.refine(2);
  
  // ELSE
  // BoxArray ba;
  // // Initialize the boxarray "ba" from the single box "bx"
  // ba.define(domain);
  // ba.maxSize(IntVect(n_cells));

  DistributionMapping dm(ba);

  Vector< MultiFab > dft_real_temp, dft_imag_temp;
  dft_real_temp.resize(variables_dft_real.size());
  dft_imag_temp.resize(variables_dft_imag.size());
  for (int d=0; d<N_VAR; d++) {
    dft_real_temp[d].define(ba, dm, 1, 0);
    dft_imag_temp[d].define(ba, dm, 1, 0);
  }
  
  amrex::Print() << "BA " << ba << std::endl;

  if (variables_dft_real[0].nGrow() != 0 || variables_temp[0].nGrow() != 0) {
    amrex::Error("Current implementation requires that both variables_temp[0] and variables_dft_real[0] have no ghost cells");
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
  amrex::Print() << "nx, ny, nz:\t" << nx << ", " << ny << ", " << nz << std::endl;
  amrex::Print() << "Number of boxes:\t" << nboxes << "\tBA size:\t" << ba.size() << std::endl;
  if (nboxes != ba.size())
    amrex::Error("NBOXES NOT COMPUTED CORRECTLY");

  Vector<int> rank_mapping;
  rank_mapping.resize(nboxes);

  DistributionMapping dmap = dft_real_temp[0].DistributionMap();

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
      // int local_index = k*nbx*nby + j*nbx + i;

      // This is what we pass to dfft to compensate for the Fortran ordering
      //      of amrex data in MultiFabs.
      int local_index = i*nby*nbz + j*nbz + k;

      rank_mapping[local_index] = dmap[ib];
      // if (verbose)
      // 	amrex::Print() << "LOADING RANK NUMBER " << dmap[ib] << " FOR GRID NUMBER " << ib 
      // 		       << " WHICH IS LOCAL NUMBER " << local_index << std::endl;
    }

  // FIXME: Assumes same grid spacing
  // Real h = geom.CellSize(0);
  // Real hsq = h*h;

  // Assume for now that nx = ny = nz
#if (AMREX_SPACEDIM == 2)
  int Ndims[3] = { nbz, nby};
  int     n[3] = {domain.length(2), domain.length(1)};
#elif (AMREX_SPACEDIM == 3)
  int Ndims[3] = { nbz, nby, nbx };
  int     n[3] = {domain.length(2), domain.length(1), domain.length(0)};
#endif
  hacc::Distribution d(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
  hacc::Dfft dfft(d);

  for (int dim=0; dim<N_VAR; dim++) {
   
    for (MFIter mfi(dft_real_temp[0],false); mfi.isValid(); ++mfi)
      {
	// int gid = mfi.index();
	// size_t local_size  = dfft.local_size();
   
	std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > a;
	std::vector<complex_t, hacc::AlignedAllocator<complex_t, ALIGN> > b;

	a.resize(nx*ny*nz);
	b.resize(nx*ny*nz);

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

	      complex_t temp(variables_temp[dim][mfi].dataPtr()[local_indx],0.);
	      a[local_indx] = temp;
	      local_indx++;

	    }
	  }
	}

	//  *******************************************
	//  Compute the forward transform
	//  *******************************************
	dfft.forward(&a[0]);
	
	// Note: Scaling for inverse FFT
	size_t global_size  = dfft.global_size();
	
	Real pi = 4.0*std::atan(1.0);
	double fac = sqrt(1.0 / global_size / pi);

	local_indx = 0;
	for(size_t k=0; k<(size_t)nz; k++) {
	  for(size_t j=0; j<(size_t)ny; j++) {
	    for(size_t i=0; i<(size_t)nx; i++) {

	      // Divide by 2 pi N
	      dft_real_temp[dim][mfi].dataPtr()[local_indx] = fac * std::real(a[local_indx]);
	      dft_imag_temp[dim][mfi].dataPtr()[local_indx] = fac * std::imag(a[local_indx]);
	      local_indx++;
	    }
	  }
	}
      }

      MultiFab::Copy(variables_dft_real[dim],dft_real_temp[dim],0,0,1,0);
      MultiFab::Copy(variables_dft_imag[dim],dft_imag_temp[dim],0,0,1,0);
  }

  bool write_data = false;
  if (write_data) {
    ShiftFFT(dft_real_temp);
    ShiftFFT(dft_real_temp);
    std::string plotname = "a_DFT"; 
    std::string x;
    for (int i=0; i<N_VAR; i++) {
      x = plotname;
      x += '_';
      x += '0' + i;
      VisMF::Write(variables_dft_real[i],x);
    }
  }
}

void StructFact::ShiftFFT(amrex::Vector< MultiFab >& dft_out) {

  // const BoxArray& ba = dft_out[0].boxArray();
  
  // Box domain(geom.Domain());
  // BoxArray ba;
  // // Initialize the boxarray "ba" from the single box "bx"
  // ba.define(domain);
  // ba.maxSize(IntVect(n_cells));

  const BoxArray& ba_in = dft_out[0].boxArray();
  BoxArray ba = ba_in;
  // ba.resize(1);

  DistributionMapping dm(ba);

  Vector< MultiFab > dft_temp;
  dft_temp.resize(dft_out.size());
  for (int d=0; d<dft_out.size(); d++) {
    dft_temp[d].define(ba, dm, 1, 0);
  }

  // NOT PARALLELIZED
  // Shift DFT by N/2+1 (pi)
  for(int d=0; d<dft_out.size(); d++) {

    MultiFab::Copy(dft_temp[d],dft_out[d],0,0,1,0);

    for (MFIter mfi(dft_out[0]); mfi.isValid(); ++mfi) {
      // Note: Make sure that multifab is cell-centered
      const Box& validBox = mfi.validbox();
      fft_shift(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
  		BL_TO_FORTRAN_ANYD(dft_temp[d][mfi]));
    }

    MultiFab::Copy(dft_out[d],dft_temp[d],0,0,1,0);

  }
}

void StructFact::SqrtMF(amrex::Vector< MultiFab >& struct_out) {
  for(int d=0; d<struct_out.size(); d++) {
    for (MFIter mfi(struct_out[0]); mfi.isValid(); ++mfi) {

      // Note: Make sure that multifab is cell-centered
      const Box& validBox = mfi.validbox();

      sqrt_mf(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
	      BL_TO_FORTRAN_ANYD(struct_out[d][mfi]));

    }
  }
}
