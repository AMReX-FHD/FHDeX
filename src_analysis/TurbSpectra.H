#ifndef _TurbSpectraDistributed_H_
#define _TurbSpectraDistributed_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>

#include <AMReX_GpuComplex.H>

#include <string>

#include "common_functions.H"

#define ALIGN 16

using namespace amrex;


void IntegrateKScalar(const MultiFab& cov_mag,
                      const std::string& name,
                      const int& step,
                      const int& comp);

void IntegrateKVelocity(const MultiFab& cov_mag,
                        const std::string& name,
                        const int& step,
                        const int& comp);

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