#ifndef _StructFact_H_
#define _StructFact_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>

// These are for FFTW / cuFFT / rocFFT

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
#if AMREX_USE_MPI
#include <fftw3-mpi.h>
#endif
#endif

#include <AMReX_GpuComplex.H>

#include <string>

#include "common_functions.H"

#define ALIGN 16

using namespace amrex;

class StructFact {

    int NVAR = 1;        // Number of variables, as defined by the size of var_names
    int NVARU = NVAR;    // Number of variables we actually need to take the FFT of
    int NCOV = 1;        // Number of covariances

    int verbosity = 0;

    // Total number of states to average over, updated by FortStructure()
    int nsamples = 0;

    // decompose velocity field
    bool decompose = false;

    // Vector containing covariance scaling
    Vector< Real > scaling;

    // Vector containing names of covariances
    amrex::Vector< std::string > cov_names;

    // 2 vectors containing structure factor pairs
    amrex::Vector< int > s_pairA;
    amrex::Vector< int > s_pairB;

    // Define vector of unique selected variables
    amrex::Vector< int > var_u;

public:

    // Vector containing running sums of real and imaginary components
    // of inner products (covariances) of DFTs
    MultiFab cov_real;
    MultiFab cov_imag;

    // Vector of MultiFabs containing final magnitude of covariances
    MultiFab cov_mag;

    // MultiFabs of real/imag for solenoidal/dilatational
    MultiFab vel_sol_real;
    MultiFab vel_sol_imag;
    MultiFab vel_dil_real;
    MultiFab vel_dil_imag;

    StructFact();

    StructFact(const amrex::BoxArray&, const amrex::DistributionMapping&, 
               const amrex::Vector< std::string >&,
               const amrex::Vector< amrex::Real >&,
               const int& verbosity=0);

    StructFact(const amrex::BoxArray&, const amrex::DistributionMapping&, 
               const amrex::Vector< std::string >&,
               const amrex::Vector< amrex::Real >&,
               const amrex::Vector< int >&, const amrex::Vector< int >&,
               const int& verbosity=0);

    void define(const amrex::BoxArray&, const amrex::DistributionMapping&, 
                const amrex::Vector< std::string >&,
                const amrex::Vector< amrex::Real >&,
                const int& verbosity=0);

    void define(const amrex::BoxArray&, const amrex::DistributionMapping&, 
                const amrex::Vector< std::string >&,
                const amrex::Vector< amrex::Real >&,
                const amrex::Vector< int >&, const amrex::Vector< int >&,
                const int& verbosity=0);

    void defineDecomp(const amrex::BoxArray&, const amrex::DistributionMapping&, 
                      const amrex::Vector< std::string >&,
                      const amrex::Vector< amrex::Real >&,
                      const amrex::Vector< int >&,
                      const amrex::Vector< int >&);

    void FortStructure(const amrex::MultiFab&, const amrex::Geometry&,
                       const int& reset=0);

    void FortStructureDecomp(const amrex::MultiFab& vel, const amrex::Geometry& geom,
                             const int& reset=0);

    void DecomposeVelFourier(const amrex::MultiFab& vel_dft_real, 
                             const amrex::MultiFab& vel_dft_imag, 
                             const amrex::Geometry& geom);

    void Reset();
    
    void ComputeFFT(const amrex::MultiFab&, amrex::MultiFab&,
                    amrex::MultiFab&, const amrex::Geometry&,
                    bool unpack=true);
    
    void InverseFFT(amrex::MultiFab&, const amrex::MultiFab&,
                    const amrex::MultiFab&, const amrex::Geometry&);

    void GetDecompVel(amrex::MultiFab&, const amrex::Geometry&); 
    
    void WritePlotFile(const int, const amrex::Real, const amrex::Geometry&, 
                       std::string, const int& zero_avg=1);
    
    void Finalize(amrex::MultiFab&, amrex::MultiFab&, const Geometry& geom,
                  const int& zero_avg=1);

    void CallFinalize(const Geometry& geom, const int& zero_avg=1);
    
    void ShiftFFT(amrex::MultiFab&,  const Geometry& geom,
                  const int& zero_avg=1);

    void IntegratekShells(const int& step, const amrex::Geometry& geom, const std::string& name="");

    void IntegratekShellsScalar(const int& step, const amrex::Geometry& geom, const amrex::Vector< std::string >& names);

    void IntegratekShellsDecomp(const int& step, const amrex::Geometry& geom, 
                                const std::string& name_sol="vel_sol", const std::string& name_dil="vel_dil");

    void AddToExternal(amrex::MultiFab& x_mag, amrex::MultiFab& x_realimag, const amrex::Geometry&, const int& zero_avg=1);

    int get_ncov() const { return NCOV; }

    void WriteCheckPoint(const int& step,
                         const amrex::Real& time,
                         std::string plotfile_base);

    const decltype(cov_names)& get_names() const { return cov_names; }
};

#endif

