#ifndef _StructFact_H
#define _StructFact_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>

// These are for FFTW or cuFFT
#ifdef AMREX_USE_CUDA
#include <cufft.h>
#else
#include <fftw3.h>
#include <fftw3-mpi.h>
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

    void FortStructure(const amrex::MultiFab&, const amrex::Geometry&,
                       const int& reset=0);

    void Reset();
    
    void ComputeFFT(const amrex::MultiFab&, amrex::MultiFab&,
                    amrex::MultiFab&, const amrex::Geometry&);
    
    void WritePlotFile(const int, const amrex::Real, const amrex::Geometry&, 
                       std::string, const int& zero_avg=1);
    
    void Finalize(amrex::MultiFab&, amrex::MultiFab&, const Geometry& geom,
                  const int& zero_avg=1);

    void CallFinalize(const Geometry& geom, const int& zero_avg=1);
    
    void ShiftFFT(amrex::MultiFab&,  const Geometry& geom,
                  const int& zero_avg=1);

    void IntegratekShells(const int& step, const amrex::Geometry& geom);

    void AddToExternal(amrex::MultiFab& x_mag, amrex::MultiFab& x_realimag, const amrex::Geometry&, const int& zero_avg=1);

    int get_ncov() const { return NCOV; }

    const decltype(cov_names)& get_names() const { return cov_names; }
};

#endif

