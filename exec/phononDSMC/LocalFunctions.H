#ifndef _INS_functions_H_
#define _INS_functions_H_


#include "DsmcParticleContainer.H"

#include "AMReX_PlotFileUtil.H"
#include "AMReX_PlotFileDataImpl.H"

#include <sys/stat.h>

#include "common_functions.H"
#include "common_namespace.H"


///////////////////////////
// in WritePlotFile.cpp


void writePlotFile(
                const MultiFab& mfcuInst,
                const MultiFab& mfcuMeans,
                const MultiFab& mfcuVars,
                FhdParticleContainer& particles,
                const Geometry& geom,
                Real time,
                int step);
///////////////////////////


///////////////////////////
// in Checkpoint.cpp

void WriteCheckPoint(int step,
                     const Real time,
                     const Real dt,
                     int statsCount,
                     const amrex::MultiFab& cuInst,
                     const amrex::MultiFab& cuMeans,
                     const amrex::MultiFab& cuVars,
                     const amrex::MultiFab& primInst,
                     const amrex::MultiFab& primMeans,
                     const amrex::MultiFab& primVars,
                     const amrex::MultiFab& coVars,
                     const FhdParticleContainer& particles,
                     const amrex::MultiFab& spatialCross1D,
                     const int ncross);

void ReadCheckPoint(int& step,
                    Real& time,
                    Real& dt,
                    int& statsCount,
                    amrex::MultiFab& cuInst,
                    amrex::MultiFab& cuMeans,
                    amrex::MultiFab& cuVars,
                    amrex::MultiFab& primInst,
                    amrex::MultiFab& primMeans,
                    amrex::MultiFab& primVars,
                    amrex::MultiFab& coVars,
                    amrex::MultiFab& spatialCross1D,
                    const int ncross);

void ReadCheckPointParticles(FhdParticleContainer& particles);

void ReadFile(const std::string& filename, Vector<char>& charBuf,
              bool bExitOnError);

void globalVars(Real time);

void ReadFile(const std::string& filename, Vector<char>& charBuf,
              bool bExitOnError=true);

///////////////////////////

#endif