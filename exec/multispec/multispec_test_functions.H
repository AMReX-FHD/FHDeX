#ifndef _multispec_test_functions_H_
#define _multispec_test_functions_H_

#include "common_functions.H"
#include "gmres_functions.H"

#include "StochMassFlux.H"
#include "StochMomFlux.H"

///////////////////////////
// in AdvanceTimestepInertial.cpp

void AdvanceTimestepInertial(std::array< MultiFab, AMREX_SPACEDIM >& umac,
                             MultiFab& rho_old,
                             MultiFab& rho_new,
                             MultiFab& rhotot_old,
                             MultiFab& rhotot_new,
                             MultiFab& pi,
                             MultiFab& eta,
                             std::array< MultiFab, NUM_EDGE >&  eta_ed,
                             MultiFab& kappa, MultiFab& Temp,
                             std::array< MultiFab, NUM_EDGE >& Temp_ed,
                             MultiFab& diff_mass_fluxdiv,
                             MultiFab& stoch_mass_fluxdiv,
                             std::array< MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                             std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot_old,
                             std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot_new,
                             MultiFab& charge_old,
                             MultiFab& charge_new,
                             MultiFab& Epot,
                             MultiFab& permittivity,
                             StochMassFlux& sMassFlux,
                             StochMomFlux& sMomFlux,
                             const Real& dt,
                             const Real& time,
                             const int& istep,
                             const Geometry& geom);

///////////////////////////
// in AdvanceTimestepBousq.cpp

void AdvanceTimestepBousq(std::array< MultiFab, AMREX_SPACEDIM >& umac,
                          MultiFab& rho_old,
                          MultiFab& rho_new,
                          MultiFab& rhotot_old,
                          MultiFab& rhotot_new,
                          MultiFab& pi,
                          MultiFab& eta,
                          std::array< MultiFab, NUM_EDGE >&  eta_ed,
                          MultiFab& kappa, MultiFab& Temp,
                          std::array< MultiFab, NUM_EDGE >& Temp_ed,
                          MultiFab& diff_mass_fluxdiv,
                          MultiFab& stoch_mass_fluxdiv,
                          std::array< MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                          std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot_old,
                          std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot_new,
                          MultiFab& charge_old,
                          MultiFab& charge_new,
                          MultiFab& Epot,
                          MultiFab& permittivity,
                          StochMassFlux& sMassFlux,
                          StochMomFlux& sMomFlux,
                          const Real& dt,
                          const Real& time,
                          const int& istep,
                          const Geometry& geom);

///////////////////////////
// in Checkpoint.cpp

void WriteCheckPoint(int step,
                     const amrex::Real time,
                     const amrex::Real dt,
                     const MultiFab& rho,
                     const MultiFab& rhotot,
                     const MultiFab& pi,
                     std::array< MultiFab, AMREX_SPACEDIM >& umac,
                     const MultiFab& Epot,
                     std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot);

void ReadCheckPoint(int& step,
                    amrex::Real& time,
                    amrex::Real& dt,
                    MultiFab& rho,
                    MultiFab& rhotot,
                    MultiFab& pi,
                    std::array< MultiFab, AMREX_SPACEDIM >& umac,
                    MultiFab& Epot,
                    std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot,
                    BoxArray& ba,
                    DistributionMapping& dmap);

void ReadFile(const std::string& filename, Vector<char>& charBuf,
              bool bExitOnError=true);

///////////////////////////
// in Init.cpp

void InitRhoUmac(std::array< MultiFab, AMREX_SPACEDIM >& umac,
                 MultiFab& rho_in,
                 const Geometry& geom);

///////////////////////////
// in WritePlotFile.cpp

void WritePlotFile(const int step,
                   const Real time,
                   const Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM >& umac,
                   const MultiFab& rhotot,
                   const MultiFab& rho,
                   const MultiFab& pres,
                   const MultiFab& charge,
                   const MultiFab& Epot);

///////////////////////////

#endif
