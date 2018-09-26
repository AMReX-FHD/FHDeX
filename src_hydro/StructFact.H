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

  int N_COV = 1;
  int N_VAR = 1;
  
  // Total number of states to average over, updated by FortStructure()
  int nsamples = 0;

  // Vector containing running sums of real and imaginary components 
  // of inner products (covariances) of DFTs
  MultiFab cov_real;
  MultiFab cov_imag;


  // Vector of MultiFabs containing final magnitude of covariances
  MultiFab cov_mag;

  amrex::Vector< std::string > cov_names;

public:

  StructFact(           const amrex::BoxArray, const amrex::DistributionMapping, 
	                const amrex::Vector< std::string >&);
  void FortStructure(   const amrex::MultiFab&, const amrex::Geometry);
  void Finalize(        const amrex::Real scale=1.0);
  void WritePlotFile(   const int, const amrex::Real, const amrex::Geometry);
  void StructOut(       amrex::MultiFab&);
  void ComputeFFT(      const amrex::MultiFab&, 
		        amrex::MultiFab&,
			amrex::MultiFab&, const amrex::Geometry);
  void ShiftFFT(        amrex::MultiFab&);
  void SqrtMF(          amrex::MultiFab&);

};

#endif

