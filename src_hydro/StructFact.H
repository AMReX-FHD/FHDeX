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

  // Vectors containing DFT of velocity field components
  amrex::Vector< MultiFab > umac_dft_real;
  amrex::Vector< MultiFab > umac_dft_imag;

  // Vector containing running sums of real and imaginary components of covariances
  amrex::Vector< MultiFab > cov_real;
  amrex::Vector< MultiFab > cov_imag;

  // Vector of MultiFabs containing inner products (covariances) of velocity DFTs
  amrex::Vector< MultiFab > struct_umac;

public:

  StructFact(BoxArray,DistributionMapping);
  void FortStructure(const std::array< MultiFab, AMREX_SPACEDIM >&, Geometry);
  void WritePlotFile(const int step, const amrex::Real time, const amrex::Geometry geom);
  void StructOut(amrex::Vector< MultiFab >&);
  void StructFinalize();
  void ComputeFFT(const std::array< MultiFab, AMREX_SPACEDIM >&, Geometry);
  void ShiftFFT(amrex::Vector< MultiFab >&);

};

#endif

