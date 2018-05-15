#ifndef _common_functions_H_
#define _common_functions_H_

#define MAX_SPECIES 10
#define LOHI 2

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

using namespace amrex;

// copy contents of common_params_module to C++ common namespace
void InitializeCommonNamespace();

///////////////////////////
// in ConvertStag.cpp

void AverageFaceToCC(const MultiFab& face, int face_comp,
                     MultiFab& cc, int cc_comp,
                     int ncomp);

void AverageCCToFace(const MultiFab& cc, int cc_comp,
                     std::array<MultiFab, AMREX_SPACEDIM>& face, int face_comp,
                     int ncomp);

void ShiftFaceToCC(const MultiFab& face, int face_comp,
                     MultiFab& cc, int cc_comp,
                     int ncomp);

///////////////////////////
// in ComputeDivAndGrad.cpp

void ComputeDiv(MultiFab& div, std::array<MultiFab, AMREX_SPACEDIM>& umac, const Geometry geom, int increment);

void ComputeGrad(MultiFab& div, std::array<MultiFab, AMREX_SPACEDIM>& umac, const Geometry geom);

///////////////////////////
// in timeStep.cpp

void eulerStep(const MultiFab& betaCC, const MultiFab& gammaCC, 
                 const MultiFab& betaNodal, const MultiFab& gammaNodal, 
                 std::array<MultiFab, AMREX_SPACEDIM>& umacIn, 
                 std::array<MultiFab, AMREX_SPACEDIM>& umacOut,
                 std::array<MultiFab, AMREX_SPACEDIM>& umacNew,
                 std::array<MultiFab, AMREX_SPACEDIM>& alpha,
                 const Geometry geom,
                 int viscType, Real* dt);

///////////////////////////



#endif
