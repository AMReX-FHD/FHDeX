#ifndef _StructFact_H_
#define _StructFact_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>
#include <AMReX_GpuComplex.H>

#include <string>

#include "common_functions.H"

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

    StructFact(const BoxArray& ba_in,
               const DistributionMapping& dmap_in,
               const Vector< std::string >& var_names,
               const Vector< Real >& var_scaling_in,
               const int& verbosity_in=0);

    StructFact(const BoxArray& ba_in,
               const DistributionMapping& dmap_in,
               const Vector< std::string >& var_names,
               const Vector< Real >& var_scaling_in,
               const Vector< int >& s_pairA_in,
               const Vector< int >& s_pairB_in,
               const int& verbosity_in=0);

    void define(const BoxArray& ba_in,
                const DistributionMapping& dmap_in,
                const Vector< std::string >& var_names,
                const Vector< Real >& var_scaling_in,
                const int& verbosity_in=0);

    void define(const BoxArray& ba_in,
                const DistributionMapping& dmap_in,
                const Vector< std::string >& var_names,
                const Vector< Real >& var_scaling_in,
                const Vector< int >& s_pairA_in,
                const Vector< int >& s_pairB_in,
                const int& verbosity_in=0);

    void FortStructure(const amrex::MultiFab&, const amrex::Geometry&,
                       const int& reset=0);

    void Reset();
    
    void ComputeFFT(const amrex::MultiFab&, amrex::MultiFab&,
                    amrex::MultiFab&, const amrex::Geometry&,
                    bool unpack=true);
    
    void InverseFFT(amrex::MultiFab&, const amrex::MultiFab&,
                    const amrex::MultiFab&, const amrex::Geometry&);

    void WritePlotFile(const int, const amrex::Real, const amrex::Geometry&, 
                       std::string, const int& zero_avg=1);
    
    void Finalize(amrex::MultiFab&, amrex::MultiFab&, const Geometry& geom,
                  const int& zero_avg=1);

    void CallFinalize(const Geometry& geom, const int& zero_avg=1);
    
    void ShiftFFT(amrex::MultiFab&,  const Geometry& geom,
                  const int& zero_avg=1);

    void IntegratekShells(const int& step, const amrex::Geometry& geom, const std::string& name="");

    void IntegratekShellsScalar(const int& step, const amrex::Geometry& geom, const amrex::Vector< std::string >& names);

    void AddToExternal(amrex::MultiFab& x_mag, amrex::MultiFab& x_realimag, const amrex::Geometry&, const int& zero_avg=1);

    int get_ncov() const { return NCOV; }

    void WriteCheckPoint(const int& step,
                         std::string checkfile_base);

    void ReadCheckPoint(std::string checkfile_base,
                        BoxArray& ba_in,
                        DistributionMapping& dmap_in);

    const decltype(cov_names)& get_names() const { return cov_names; }
};

#endif

