
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

  cov_real.define(ba_in, dmap_in, N_COV, 0);
  cov_imag.define(ba_in, dmap_in, N_COV, 0);
  cov_mag.define( ba_in, dmap_in, N_COV, 0);
  cov_real.setVal(0.0);
  cov_imag.setVal(0.0);
  cov_mag.setVal( 0.0);

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
      cnt++;
    }
  }
}

void StructFact::FortStructure(const MultiFab& variables, const Geometry geom) {

  BL_PROFILE_VAR("StructFact::FortStructure()",FortStructure);

  const BoxArray& ba = variables.boxArray();
  const DistributionMapping& dm = variables.DistributionMap();

  // if (ba.size() != 1) {
  //   Abort("Only adapted for single grid");
  //   exit(0);
  // }

  if (false) {    
    std::string plotname = "a_VAR";
    VisMF::Write(variables,plotname);
  }

  MultiFab variables_dft_real, variables_dft_imag;
  variables_dft_real.define(ba, dm, N_VAR, 0);
  variables_dft_imag.define(ba, dm, N_VAR, 0);

  ComputeFFT(variables, variables_dft_real, variables_dft_imag, geom);

  MultiFab cov_temp;
  cov_temp.define(ba, dm, 1, 0);

  int index = 0;
  for (int j=0; j<N_VAR; j++) {
    for (int i=j; i<N_VAR; i++) {

      // Compute temporary real and imaginary components of covariance

      // Real component of covariance
      cov_temp.setVal(0.0);
      MultiFab::AddProduct(cov_temp,variables_dft_real,i,variables_dft_real,j,0,1,0);
      MultiFab::AddProduct(cov_temp,variables_dft_imag,i,variables_dft_imag,j,0,1,0);

      MultiFab::Add(cov_real,cov_temp,0,index,1,0);

      // Imaginary component of covariance
      cov_temp.setVal(0.0);
      MultiFab::AddProduct(cov_temp,variables_dft_imag,i,variables_dft_real,j,0,1,0);
      cov_temp.mult(-1.0,0);
      MultiFab::AddProduct(cov_temp,variables_dft_real,i,variables_dft_imag,j,0,1,0);
      
      MultiFab::Add(cov_imag,cov_temp,0,index,1,0);
      
      index++;
    }
  }

  bool write_data = false;
  if (write_data) {
    std::string plotname; 
    plotname = "a_COV"; 
    VisMF::Write(cov_real,plotname);
  }

  nsamples++;
}

void StructFact::WritePlotFile(const int step, const Real time, const Geometry geom) {
  
  BL_PROFILE_VAR("StructFact::WritePlotFile()",WritePlotFile);

  MultiFab plotfile;
  Vector<std::string> varNames;
  int nPlot = 1;

  //////////////////////////////////////////////////////////////////////////////////
  // Write out structure factor magnitude to plot file
  //////////////////////////////////////////////////////////////////////////////////
  const std::string plotfilename1 = "plt_structure_factor_mag";
  nPlot = N_COV;
  plotfile.define(cov_mag.boxArray(), cov_mag.DistributionMap(), nPlot, 0);
  varNames.resize(nPlot);

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
    MultiFab::Copy(plotfile, cov_mag, cnt, cnt, 1, 0);
    cnt++;
  }

  // write a plotfile
  WriteSingleLevelPlotfile(plotfilename1,plotfile,varNames,geom,time,step);
  
  //////////////////////////////////////////////////////////////////////////////////
  // Write out real and imaginary components of structure factor to plot file
  //////////////////////////////////////////////////////////////////////////////////
  const std::string plotfilename2 = "plt_structure_factor_real_imag";
  nPlot = 2*N_COV;
  plotfile.define(cov_mag.boxArray(), cov_mag.DistributionMap(), nPlot, 0);
  varNames.resize(nPlot);

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
    MultiFab::Copy(plotfile,cov_real,d,cnt,1,0);
    cnt++;
  }

  for(int d=0; d<N_COV; d++) {   
    MultiFab::Copy(plotfile,cov_imag,d,cnt,1,0);
    cnt++;
  }

  // write a plotfile
  WriteSingleLevelPlotfile(plotfilename2,plotfile,varNames,geom,time,step);
}

void StructFact::StructOut(MultiFab& struct_out) {

  BL_PROFILE_VAR("StructFact::StructOut()",StructOut);
  
  
}

void StructFact::Finalize(const amrex::Real scale) {
  
  Real nsamples_inv = 1.0/(Real)nsamples;
  
  ShiftFFT(cov_real);
  ShiftFFT(cov_imag);

  cov_real.mult(nsamples_inv);
  cov_real.mult(scale);

  cov_imag.mult(nsamples_inv);
  cov_imag.mult(scale);

  cov_mag.setVal(0.0);
  MultiFab::AddProduct(cov_mag,cov_real,0,cov_real,0,0,N_COV,0);
  MultiFab::AddProduct(cov_mag,cov_imag,0,cov_imag,0,0,N_COV,0);

  SqrtMF(cov_mag);

}

void StructFact::ComputeFFT(const MultiFab& variables,
			    MultiFab& variables_dft_real, 
			    MultiFab& variables_dft_imag,
			    const Geometry geom) {

  BL_PROFILE_VAR("StructFact::ComputeFFT()", ComputeFFT);

  Box domain(geom.Domain());
  const BoxArray& ba = variables.boxArray();
  DistributionMapping dm = variables.DistributionMap();

  // BoxArray ba;
  // IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
  // IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
  // Box domain(dom_lo, dom_hi);
  // // Initialize the boxarray "ba" from the single box "bx"
  // ba.define(domain);
  // DistributionMapping dm(ba);

  MultiFab dft_real_temp, dft_imag_temp, variables_temp;
  dft_real_temp.define(ba, dm, 1, 0);
  dft_imag_temp.define(ba, dm, 1, 0);
  variables_temp.define(ba, dm, 1, 0);
  
  amrex::Print() << "BA " << ba << std::endl;

  if (variables_dft_real.nGrow() != 0 || variables.nGrow() != 0) {
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

  DistributionMapping dmap = dft_real_temp.DistributionMap();
  // DistributionMapping dmap = variables.DistributionMap();

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

      // This is what we pass to dfft to compensate for the Fortran ordering
      //      of amrex data in MultiFabs.
      // int local_index = i*nby*nbz + j*nbz + k;

      rank_mapping[local_index] = dmap[ib];
      // if (verbose)
      // 	amrex::Print() << "LOADING RANK NUMBER " << dmap[ib] << " FOR GRID NUMBER " << ib 
      // 		       << " WHICH IS LOCAL NUMBER " << local_index << std::endl;
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
  hacc::Distribution d(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
  hacc::Dfft dfft(d);

  // Print() << "RANK MAPPING: \n";
  // for (int i=0; i<rank_mapping.size(); i++) {
  //   Print() << "\t" << rank_mapping[i] << std::endl;
  // }

  for (int dim=0; dim<N_VAR; dim++) {

    // variables_temp.ParallelCopy(variables, dim, 0, 1);
    MultiFab::Copy(variables_temp,variables,dim,0,1,0);
   
    for (MFIter mfi(dft_real_temp,false); mfi.isValid(); ++mfi)
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

	      complex_t temp(variables_temp[mfi].dataPtr()[local_indx],0.);
	      a[local_indx] = temp;
	      local_indx++;

	    }
	  }
	}

	//  *******************************************
	//  Compute the forward transform
	//  *******************************************
	dfft.forward(&a[0]);

	d.redistribute_2_to_3(&a[0],&b[0],2);
	
	// Note: Scaling for inverse FFT
	size_t global_size  = dfft.global_size();
	
	// Real pi = 4.0*std::atan(1.0);
	double fac = sqrt(1.0 / global_size);

	local_indx = 0;
	for(size_t k=0; k<(size_t)nz; k++) {
	  for(size_t j=0; j<(size_t)ny; j++) {
	    for(size_t i=0; i<(size_t)nx; i++) {

	      // Divide by 2 pi N
	      dft_real_temp[mfi].dataPtr()[local_indx] = fac * std::real(b[local_indx]);
	      dft_imag_temp[mfi].dataPtr()[local_indx] = fac * std::imag(b[local_indx]);
	      local_indx++;
	    }
	  }
	}
      }
    
    // variables_dft_real.ParallelCopy(dft_real_temp, 0, dim, 1);
    // variables_dft_imag.ParallelCopy(dft_imag_temp, 0, dim, 1);
    MultiFab::Copy(variables_dft_real,dft_real_temp,0,dim,1,0);
    MultiFab::Copy(variables_dft_imag,dft_imag_temp,0,dim,1,0);
  }

  bool write_data = false;
  if (write_data) {
    std::string plotname = "a_DFT"; 
    VisMF::Write(variables_dft_real,plotname);
  }
}

void StructFact::ShiftFFT(MultiFab& dft_out) {

  BoxArray ba_onegrid;
  {
      IntVect dom_lo(AMREX_D_DECL(           0,            0,            0));
      IntVect dom_hi(AMREX_D_DECL(n_cells[0]-1, n_cells[1]-1, n_cells[2]-1));
      Box domain(dom_lo, dom_hi);
      
      // Initialize the boxarray "ba" from the single box "bx"
      ba_onegrid.define(domain);
  }

  DistributionMapping dmap_onegrid(ba_onegrid);

  MultiFab dft_onegrid;
  dft_onegrid.define(ba_onegrid, dmap_onegrid, 1, 0);

  for (int d=0; d<N_COV; d++) {
    dft_onegrid.ParallelCopy(dft_out, d, 0, 1);

    // Shift DFT by N/2+1 (pi)
    for (MFIter mfi(dft_onegrid); mfi.isValid(); ++mfi) {
      // Note: Make sure that multifab is cell-centered
      const Box& validBox = mfi.validbox();
      fft_shift(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
		BL_TO_FORTRAN_FAB(dft_onegrid[mfi]));
    }

    dft_out.ParallelCopy(dft_onegrid, 0, d, 1);
  }

}

void StructFact::SqrtMF(MultiFab& struct_out) {
  for (MFIter mfi(struct_out); mfi.isValid(); ++mfi) {

    // Note: Make sure that multifab is cell-centered
    const Box& validBox = mfi.validbox();

    sqrt_mf(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
	    BL_TO_FORTRAN_FAB(struct_out[mfi]));

  }
}
