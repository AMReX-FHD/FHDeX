#ifndef _common_functions_H_
#define _common_functions_H_

#define MAX_SPECIES 10
#define LOHI 2

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>

#include "common_functions_F.H"
#include "common_namespace.H"

using namespace common;
using namespace amrex;

// copy contents of common_params_module to C++ common namespace
void InitializeCommonNamespace();

///////////////////////////
// in ConvertStag.cpp

void AverageFaceToCC(const std::array<MultiFab, AMREX_SPACEDIM>& face,
                     MultiFab& cc, int cc_comp);

void AverageCCToFace(const MultiFab& cc, int cc_comp,
                     std::array<MultiFab, AMREX_SPACEDIM>& face, int face_comp,
                     int ncomp);

void ShiftFaceToCC(const MultiFab& face, int face_comp,
                     MultiFab& cc, int cc_comp,
                     int ncomp);

///////////////////////////
// in ComputeDivAndGrad.cpp

void ComputeDiv(MultiFab & div, const std::array<MultiFab, AMREX_SPACEDIM> & phi_fc, 
                int start_incomp, int start_outcomp, int ncomp, 
                const Geometry & geom, int increment);

void ComputeGrad(const MultiFab & phi, std::array<MultiFab, AMREX_SPACEDIM> & gphi, 
                 int start_incomp, int start_outcomp, int ncomp, const Geometry & geom);

void ComputeCentredGrad(const MultiFab& phi, std::array<MultiFab, AMREX_SPACEDIM>& gphi, const Geometry& geom);

///////////////////////////
// in Debug.cpp
void PrintMF(const MultiFab& MF);

void outputMFAscii(const MultiFab& output, std::string filename);

///////////////////////////
// in SqrtMF.cpp
void SqrtMF(MultiFab& MF);

///////////////////////////
// in FindCoords.cpp

void FindFaceCoords(std::array<MultiFab, AMREX_SPACEDIM> & RealFaceCoords,
                    const Geometry & geom);

void FindCenterCoords(MultiFab & RealCenterCoords, const Geometry & geom);

void FindNodalValues(const MultiFab & inFab, MultiFab & outFab, const MultiFab & ccFab);


///////////////////////////
// in ComputeBasicStats.cpp

void ComputeBasicStats(const MultiFab & instant, MultiFab & means, MultiFab & vars,
                const int incomp, const int outcomp, const int steps);

void OutputVolumeMean(const MultiFab & instant, const int comp, const Real domainVol, std::string filename, const Geometry geom);

Real MaskedSum (const MultiFab & inFab,int comp, const Periodicity& period);

Real SumFab(const MultiFab & in, const int ng, const int comp);

void XMeanFab(const MultiFab & in, MultiFab & out, const int ng);

void MaxSpeed(std::array< MultiFab, AMREX_SPACEDIM >& umac);

///////////////////////////
// in ComputeAverages.cpp

void ComputeVerticalAverage(const MultiFab & mf, MultiFab & mf_avg, const Geometry & geom,
                            const int dir, const int incomp, const int ncomp,
                            const int findredist=0);

///////////////////////////
// in MultiFABPhysBC.cpp
// see comments in C++ files for descriptions

void MultiFABPhysBCPres(MultiFab& data, const Geometry& geom);

void MultiFABPhysBCDomainVel(MultiFab& vel, const Geometry& geom, int dim);

void MultiFABPhysBCMacVel(MultiFab& vel, const Geometry& geom, int dim);

void ZeroEdgevalWalls(MultiFab& edge, const Geometry& geom, int scomp, int ncomp);

void MultiFABElectricBC(MultiFab& efieldCC, const Geometry& geom);

void MultiFABPotentialBC(MultiFab& potential, const Geometry& geom);

void MultiFABPotentialBC_solver(MultiFab& potential, const Geometry& geom);

void MultiFABPhysBCCharge(MultiFab& charge, const Geometry& geom);

void MultiFABPhysBCMacStress(MultiFab& stress, const Geometry& geom, int dim);

void MultiFABPhysBCDomainStress(MultiFab& stress, const amrex::Geometry& geom, int dim);

#endif
