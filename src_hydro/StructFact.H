#ifndef _StructFact_H
#define _StructFact_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>

// These are for SWFFT
#include <Distribution.H>
#include <AlignedAllocator.h>
#include <Dfft.H>

#include <string>

#include <common_namespace.H>

#define ALIGN 16

using namespace amrex;
using namespace common;

class StructFact {
  
  // Geometry geom;
#if (AMREX_SPACEDIM == 2)
  const int COV_NVAR = 3;
#elif (AMREX_SPACEDIM == 3)
  const int COV_NVAR = 6;
#endif

  // Vectors containing DFT of velocity field components
  amrex::Vector< MultiFab > umac_dft_real;
  amrex::Vector< MultiFab > umac_dft_imag;

  // 2D array of MultiFabs containing inner products (covariances) of velocity DFTs
  amrex::Vector< MultiFab > struct_umac;

public:

  StructFact(BoxArray,DistributionMapping);
  void FortStructure(const std::array< MultiFab, AMREX_SPACEDIM >& umac, Geometry);
  void ComputeFFT(const std::array< MultiFab, AMREX_SPACEDIM >& umac_cc, Geometry);

};

////////////////////////////////////////////////

StructFact::StructFact(BoxArray ba_in, DistributionMapping dmap_in) {

  // static_assert(AMREX_SPACEDIM == 3, "3D only");

  // Note that we are defining with NO ghost cells

  umac_dft_real.resize(AMREX_SPACEDIM);
  umac_dft_imag.resize(AMREX_SPACEDIM);
  struct_umac.resize(COV_NVAR);
  
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    umac_dft_real[d].define(ba_in, dmap_in, 1, 0);
    umac_dft_imag[d].define(ba_in, dmap_in, 1, 0);
  }
  for (int d=0; d<COV_NVAR; d++) {
    struct_umac[d].define(ba_in, dmap_in, 1, 0);
  }

}


void StructFact::FortStructure(const std::array< MultiFab, AMREX_SPACEDIM >& umac, Geometry geom) {

  const BoxArray& ba = umac_dft_real[0].boxArray();
  // const BoxArray& ba = amrex::enclosedCells(ba_nodal);
  const DistributionMapping& dm = umac[0].DistributionMap();
  
  std::array< MultiFab, AMREX_SPACEDIM > umac_cc;
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    umac_cc[d].define(ba, dm, 1, 0);
  }
  
  for (int d=0; d<AMREX_SPACEDIM; d++) { 
    AverageFaceToCC(umac[d], 0, umac_cc[d], 0, 1);
  }

  Real start_time = amrex::second();
  ComputeFFT(umac_cc, geom);
  Real total_time = amrex::second() - start_time;

  bool write_data = false;
  if (write_data) {
    // VisMF::Write(rhs,"RHS");
  }

  // if (verbose) {
  //   // amrex::Print() << "norm is " << diff.norm0() << "\n";
  //   // amrex::Print() << "Time spent in solve: " << total_time << std::endl;    
  // }

}

void StructFact::ComputeFFT(const std::array< MultiFab, AMREX_SPACEDIM >& umac_cc, Geometry geom) {
  
  const BoxArray& ba = umac_dft_real[0].boxArray();
  // amrex::Print() << "BA " << ba << std::endl;
  const DistributionMapping& dm = umac_dft_real[0].DistributionMap();

  if (umac_dft_real[0].nGrow() != 0 || umac_dft_real[0].nGrow() != 0) 
    amrex::Error("Current implementation requires that both umac_cc[0] and umac_dft_real[0] have no ghost cells");

  // Define pi and (two pi) here
  Real  pi = 4 * std::atan(1.0);
  Real tpi = 2 * pi;

  // We assume that all grids have the same size hence 
  // we have the same nx,ny,nz on all ranks
  int nx = ba[0].size()[0];
  int ny = ba[0].size()[1];
  int nz = ba[0].size()[2];

  Box domain(geom.Domain());

  int nbx = domain.length(0) / nx;
  int nby = domain.length(1) / ny;
  int nbz = domain.length(2) / nz;
  int nboxes = nbx * nby * nbz;
  if (nboxes != ba.size()) 
    amrex::Error("NBOXES NOT COMPUTED CORRECTLY");
  amrex::Print() << "Number of boxes:\t" << nboxes << std::endl;

  Vector<int> rank_mapping;
  rank_mapping.resize(nboxes);

  DistributionMapping dmap = umac_cc[0].DistributionMap();

  for (int ib = 0; ib < nboxes; ++ib)
    {
      int i = ba[ib].smallEnd(0) / nx;
      int j = ba[ib].smallEnd(1) / ny;
      int k = ba[ib].smallEnd(2) / nz;

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

  Real h = geom.CellSize(0);
  Real hsq = h*h;

  Real start_time = amrex::second();

  // Assume for now that nx = ny = nz
  int Ndims[3] = { nbz, nby, nbx };
  int     n[3] = {domain.length(2), domain.length(1), domain.length(0)};
  hacc::Distribution d(MPI_COMM_WORLD,n,Ndims,&rank_mapping[0]);
  hacc::Dfft dfft(d);

  for (int dim=0; dim<AMREX_SPACEDIM; dim++) { 
    for (MFIter mfi(umac_cc[0],false); mfi.isValid(); ++mfi)
      {
	int gid = mfi.index();

	size_t local_size  = dfft.local_size();
   
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

	      complex_t temp(umac_cc[0][mfi].dataPtr()[local_indx],0.);
	      a[local_indx] = temp;
	      Print() << temp << "\n";
	      local_indx++;

	    }
	  }
	}

	//  *******************************************
	//  Compute the forward transform
	//  *******************************************
	dfft.forward(&a[0]);

	size_t global_size  = dfft.global_size();
	double fac = hsq / global_size;

	local_indx = 0;
	for(size_t k=0; k<(size_t)nz; k++) {
	  for(size_t j=0; j<(size_t)ny; j++) {
	    for(size_t i=0; i<(size_t)nx; i++) {

	      // Divide by 2 pi N
	      umac_dft_real[dim][mfi].dataPtr()[local_indx] = fac * std::real(a[local_indx]);
	      umac_dft_imag[dim][mfi].dataPtr()[local_indx] = fac * std::imag(a[local_indx]);
	      local_indx++;
	    }
	  }
	}
      }
  }
}

//*/

#endif

