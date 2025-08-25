#ifndef _hydro_test_functions_H_
#define _hyrdo_test_functions_H_

#include "common_functions.H"
#include "gmres_functions.H"

#include "TurbForcing.H"

///////////////////////////
// in advance.cpp
void advance(std::array< amrex::MultiFab, AMREX_SPACEDIM >& umac,
	     std::array< amrex::MultiFab, AMREX_SPACEDIM >& umacNew,
             amrex::MultiFab& pres,
	     const std::array< amrex::MultiFab, AMREX_SPACEDIM >& mfluxdiv_stoch,
	     std::array< amrex::MultiFab, AMREX_SPACEDIM >& alpha_fc,
	     amrex::MultiFab& beta, amrex::MultiFab& gamma,
	     std::array< amrex::MultiFab, NUM_EDGE >& beta_ed,
	     const amrex::Geometry& geom, const amrex::Real& dt,
             TurbForcing& turbforce);

///////////////////////////

///////////////////////////
// in Init.cpp

void InitVel(std::array< MultiFab, AMREX_SPACEDIM >& umac,
             const Geometry& geom);

///////////////////////////

///////////////////////////
// in WritePlotFile.cpp

void WritePlotFile(const int step,
                   const amrex::Real time,
                   const amrex::Geometry& geom,
                   std::array< MultiFab, AMREX_SPACEDIM >& umac,
                   const MultiFab& pres);

///////////////////////////

///////////////////////////
// in Checkpoint.cpp

void WriteCheckPoint(const int step,
                     const amrex::Real time,
                     std::array< MultiFab, AMREX_SPACEDIM >& umac,
                     TurbForcing& turbforce);

void ReadCheckPoint(int& step,
                    amrex::Real& time,
                    std::array< MultiFab, AMREX_SPACEDIM >& umac,
                    TurbForcing& turbforce,
                    BoxArray& ba, DistributionMapping& dmap);

void ReadFile(const std::string& filename, Vector<char>& charBuf,
              bool bExitOnError=true);

///////////////////////////

///////////////////////////
// in ExternalForce.cpp

void ExternalForce(std::array< MultiFab, AMREX_SPACEDIM >& gmres_rhs_u,
                   MultiFab& gmres_rhs_p);

#endif
