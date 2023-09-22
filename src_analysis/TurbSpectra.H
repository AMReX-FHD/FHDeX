#ifndef _TurbSpectra_H_
#define _TurbSpectra_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>

// These are for heFFTe / FFTW / cuFFT / rocFFT

#if defined(USE_HEFFTE)
#include <heffte.h>
#else
#ifdef AMREX_USE_CUDA
#include <cufft.h>
#elif AMREX_USE_HIP
#  if __has_include(<rocfft/rocfft.h>)  // ROCm 5.3+
#    include <rocfft/rocfft.h>
#  else
#    include <rocfft.h>
#  endif
#else
#include <fftw3.h>
#include <fftw3-mpi.h>
#endif
#endif

#include <AMReX_GpuComplex.H>

#include <string>

#include "common_functions.H"

#define ALIGN 16

using namespace amrex;

#if defined(USE_HEFFTE)
#else
#ifdef AMREX_USE_CUDA
std::string cufftError (const cufftResult& err);
#endif
#ifdef AMREX_USE_HIP
std::string rocfftError (const rocfft_status err);
void rocfft_status (std::string const& name, rocfft_status status);
#endif
#endif

#if defined(USE_HEFFTE)
void IntegrateKScalar(BaseFab<GpuComplex<Real> >& spectral_field,
                      const std::string& name, const Real& scaling,
                      const Box& c_local_box,
                      const Real& sqrtnpts);
void IntegrateKVelocity(BaseFab<GpuComplex<Real> >& spectral_field,
                        const std::string& name, const Real& scaling,
                        const Box& c_local_box);
#else
void IntegrateKScalar(Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > >& spectral_field,
                      const MultiFab& variables_onegrid,
                      const std::string& name,
                      const Real& scaling,
                      const Real& sqrtnpts);
void IntegrateKVelocity(Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > >& spectral_field,
                        const MultiFab& vel_onegrid,
                        const std::string& name,
                        const Real& scaling);
#endif

void TurbSpectrumScalar(const MultiFab& variables, 
                        const amrex::Geometry& geom, 
                        const int& step, 
                        const amrex::Vector<amrex::Real>& var_scaling,
                        const amrex::Vector< std::string >& var_names);

void TurbSpectrumVelDecomp(const MultiFab& vel,
                           MultiFab& vel_decomp,
                           const amrex::Geometry& geom,
                           const int& step,
                           const amrex::Real& var_scaling,
                           const amrex::Vector< std::string >& var_names);


#endif
