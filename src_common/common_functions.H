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

void ComputeDiv(MultiFab& div, std::array<MultiFab, AMREX_SPACEDIM>& phi_fc, 
                int start_incomp, int start_outcomp, int ncomp, 
                const Geometry& geom, int increment);

void ComputeGrad(MultiFab& phi, std::array<MultiFab, AMREX_SPACEDIM>& gphi, 
                 int start_incomp, int start_outcomp, int ncomp, const Geometry& geom);


///////////////////////////
// in FindCoords.cpp

void FindFaceCoords(std::array< MultiFab, AMREX_SPACEDIM >& RealFaceCoords, Geometry geom);

void FindNodalValues(const MultiFab& inFab, MultiFab& outFab, const MultiFab& ccFab);

///////////////////////////
// in MultiFABFillRandom.cpp

void MultiFABFillRandom(MultiFab& mf, const int& comp, const amrex::Real& variance, const Geometry& geom);

#endif
