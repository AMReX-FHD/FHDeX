#ifndef _INS_functions_H_
#define _INS_functions_H_

#include "common_functions.H"
#include "gmres_functions.H"
#include "FhdParticleContainer.H"

///////////////////////////
// in WritePlotFile.cpp

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   const amrex::Geometry cgeom,
                   const amrex::Geometry egeom,
                   const MultiFab& particleInstant,
                   const MultiFab& particleMeans,
                   FhdParticleContainer& particles,
                   const MultiFab& charge,
                   const MultiFab& chargeM,
                   const MultiFab& potential,
                   const MultiFab& potentialM,
                   const std::array< MultiFab, AMREX_SPACEDIM >& efield);

///////////////////////////

///////////////////////////
// in WritePlotFileHydro.cpp

void WritePlotFileHydro(int step,
                        const amrex::Real time,
                        const amrex::Geometry geom,
                        std::array< MultiFab, AMREX_SPACEDIM >& umac,
                        const MultiFab& pres,
                        std::array< MultiFab, AMREX_SPACEDIM >& umacM);

///////////////////////////

///////////////////////////
// in Checkpoint.cpp

void WriteCheckPoint(int step,
                     amrex::Real time,
                     int statsCount,
                     const std::array< MultiFab, AMREX_SPACEDIM >& umac,
                     const std::array< MultiFab, AMREX_SPACEDIM >& umacM,
                     const MultiFab& pres,
                     const FhdParticleContainer& particles,
                     const MultiFab& particleMeans,
                     const MultiFab& particleVars,
                     const MultiFab& chargeM,
                     const MultiFab& potential,
                     const MultiFab& potentialM);

void ReadCheckPoint(int& step,
                    amrex::Real& time,
                    int& statsCount,
                    std::array< MultiFab, AMREX_SPACEDIM >& umac,
                    std::array< MultiFab, AMREX_SPACEDIM >& umacM,
                    MultiFab& pres,
                    MultiFab& particleMeans,
                    MultiFab& particleVars,
                    MultiFab& chargeM,
                    MultiFab& potential,
                    MultiFab& potentialM);

void ReadCheckPointParticles(FhdParticleContainer& particles, species* particleInfo, const Real* dxp);

void ReadFile(const std::string& filename, Vector<char>& charBuf,
              bool bExitOnError=true);

///////////////////////////

#endif
