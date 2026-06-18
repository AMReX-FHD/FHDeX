#ifndef _spectral_functions_stag_H_
#define _spectral_functions_stag_H_

#include <heffte.h>
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>


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

#include <AMReX_GpuComplex.H>

#include <string>

#define ALIGN 16

using namespace amrex;

#if !defined(HEFFTE_FFTW) && !defined(HEFFTE_CUFFT) && !defined(HEFFTE_ROCFFT)
#ifdef AMREX_USE_CUDA
std::string cufftError (const cufftResult& err);
#endif
#ifdef AMREX_USE_HIP
std::string rocfftError (const rocfft_status err);
void Assert_rocfft_status (std::string const& name, rocfft_status status);
#endif
#endif

void SpectralReadCheckPoint(amrex::Geometry& geom,
                            const amrex::Box& domain,
                            amrex::MultiFab& prim,
                            std::array<MultiFab, 3>& vel,
                            BoxArray& ba, DistributionMapping& dmap,
                            const amrex::Vector<int> n_cells,
                            const int nprimvars,
                            const amrex::Vector<int> max_grid_size,
                            const amrex::IntVect ngc,
                            const int restart);

void SpectralVelDecomp(const MultiFab& vel,
                       MultiFab& vel_decomp_filter,
                       const amrex::Real kmin,
                       const amrex::Real kmax,
                       const amrex::Geometry& geom,
                       const amrex::Vector<int> n_cells);

void SpectralScalarDecomp(const MultiFab& scalar,
                          MultiFab& scalar_filter,
                          const amrex::Real kmin,
                          const amrex::Real kmax,
                          const amrex::Geometry& geom,
                          const amrex::Vector<int> n_cells);

void SpectralWritePlotFile(const int step,
                           const amrex::Real& kmin,
                           const amrex::Real& kmax,
                           const amrex::Geometry& geom,
                           const amrex::MultiFab& vel_decomp_in,
                           const amrex::MultiFab& scalar_in,
                           const amrex::MultiFab& vel_total,
                           const amrex::MultiFab& scalar_total);

void Read_Copy_MF_Checkpoint(amrex::MultiFab& mf, std::string mf_name,
                             const std::string& checkpointname,
                             BoxArray& ba_old, DistributionMapping& dmap_old,
                             int NVARS, int NGC, const amrex::IntVect ngc,
                             int nodal_flag=-1);

void ShiftFaceToCC(const MultiFab& face_in, int face_in_comp,
                   MultiFab& cc_in, int cc_in_comp,
                   int ncomp);

void ComputeGrad(const MultiFab & phi_in, std::array<MultiFab, 3> & gphi,
                 int start_incomp, int start_outcomp, int ncomp, int bccomp, const Geometry & geom,
                 int increment);

void SumStag(const std::array<MultiFab, 3>& m1,
    amrex::Vector<amrex::Real>& sum);

void FCMoments(const std::array<MultiFab, 3>& m1,
    const amrex::Vector<int>& comps,
    std::array<MultiFab, 3>&  mscr,
    const int& power,
    amrex::Vector<amrex::Real>& prod_val);

void SumCC(const amrex::MultiFab& m1,
    const int& comp,
    amrex::Real& sum,
    const bool& divide_by_ncells);

void CCMoments(const amrex::MultiFab& m1,
    const int& comp1,
    amrex::MultiFab& mscr,
    const int& power,
    amrex::Real& prod_val);

#endif