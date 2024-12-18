#ifndef _TurbSpectraSingle_H_
#define _TurbSpectraSingle_H_

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

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>

#include <AMReX_GpuComplex.H>

#include <string>

#include "common_functions.H"

#define ALIGN 16

using namespace amrex;

#ifdef AMREX_USE_CUDA
std::string cufftError (const cufftResult& err);
#endif
#ifdef AMREX_USE_HIP
std::string rocfftError (const rocfft_status err);
void Assert_rocfft_status (std::string const& name, rocfft_status status);
#endif

void IntegrateKScalar(const Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > >& spectral_field,
                      const MultiFab& variables_onegrid,
                      const std::string& name,
                      const Real& scaling,
                      const Real& sqrtnpts,
                      const int& step);
void IntegrateKVelocity(const Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > >& spectral_fieldx,
                        const Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > >& spectral_fieldy,
                        const Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > >& spectral_fieldz,
                        const MultiFab& vel_onegrid,
                        const std::string& name,
                        const Real& scaling,
                        const int& step);
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
void InverseFFTVel(Vector<std::unique_ptr<BaseFab<GpuComplex<Real> > > >& spectral_field, 
                   MultiFab& vel_decomp_onegrid, const IntVect& fft_size);


#endif
