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
  
#if (AMREX_SPACEDIM == 2)
  const int COV_NVAR = 3;
#elif (AMREX_SPACEDIM == 3)
  const int COV_NVAR = 6;
#endif
  
  // Total number of states to average over, updated by FortStructure()
  int nsamples = 0;

  // Vector containing running sums of real and imaginary components 
  // of inner products (covariances) of DFTs
  amrex::Vector< MultiFab > cov_real;
  amrex::Vector< MultiFab > cov_imag;

  // Vector of MultiFabs containing final magnitude of covariances
  amrex::Vector< MultiFab > cov_mag;

public:

  StructFact(const amrex::BoxArray, const amrex::DistributionMapping);
  void FortStructure(const std::array< MultiFab, AMREX_SPACEDIM >&, const amrex::Geometry);
  void WritePlotFile(const int, const amrex::Real, const amrex::Geometry,
		     const amrex::Real scale=1.0);
  void StructOut(amrex::Vector< MultiFab >&, const amrex::Real scale=1.0);
  void StructFinalize(const amrex::Real scale=1.0);
  void ComputeFFT(const std::array< MultiFab, AMREX_SPACEDIM >&, 
		  amrex::Vector< MultiFab >&,
		  amrex::Vector< MultiFab >&, const amrex::Geometry);
  // void ShiftFFT(amrex::Vector< MultiFab >&, const amrex::Geometry);
  void ShiftFFT(amrex::Vector< MultiFab >&);
  void SqrtMF(amrex::Vector< MultiFab >&);

};

#endif

