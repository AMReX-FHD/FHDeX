#ifndef _common_functions_H_
#define _common_functions_H_

#define MAX_SPECIES 10
#define LOHI 2

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>

#include "common_functions_F.H"
#include "common_namespace.H"

using namespace common;
using namespace amrex;

#ifndef _num_edge_
#define _num_edge_
#if (AMREX_SPACEDIM == 2)
    const int NUM_EDGE = 1;
#elif (AMREX_SPACEDIM == 3)
    const int NUM_EDGE = 3;
#endif
#endif

// copy contents of common_params_module to C++ common namespace
void InitializeCommonNamespace();

///////////////////////////
// in BCPhysToMath.cpp
void BCPhysToMath(int type, Vector<int>& bc_lo, Vector<int>& bc_hi);

///////////////////////////

///////////////////////////
// in ConvertStag.cpp

void AverageFaceToCC(const std::array<MultiFab, AMREX_SPACEDIM>& face_in,
                     MultiFab& cc_in, int cc_in_comp);

void AverageCCToFace(const MultiFab& cc_in,
                     std::array<MultiFab, AMREX_SPACEDIM>& face_in,
                     int scomp, int ncomp, int varType, const Geometry& geom);

void ShiftFaceToCC(const MultiFab& face_in, int face_in_comp,
                   MultiFab& cc_in, int cc_in_comp,
                   int ncomp);

void AverageCCToNode(const MultiFab& cc_in, MultiFab& node_in,
                     int scomp, int ncomp, int varType, const Geometry& geom);

void AverageCCToEdge(const MultiFab& cc_in, std::array<MultiFab, NUM_EDGE>& edge_in,
                     int scomp, int ncomp, int varType, const Geometry& geom);

///////////////////////////
// in ComputeDivAndGrad.cpp

void ComputeDiv(MultiFab & div, const std::array<MultiFab, AMREX_SPACEDIM> & phi_fc, 
                int start_incomp, int start_outcomp, int ncomp, 
                const Geometry & geom, int increment);

void ComputeGrad(const MultiFab & phi_in, std::array<MultiFab, AMREX_SPACEDIM> & gphi, 
                 int start_incomp, int start_outcomp, int ncomp, int varType, const Geometry & geom,
                 int increment=0);

void ComputeCentredGrad(const MultiFab& phi, std::array<MultiFab, AMREX_SPACEDIM>& gphi,
                        const Geometry& geom);

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
// in MultiFabPhysBC.cpp
// see comments in C++ files for descriptions

void MultiFabPhysBC(MultiFab& data, const Geometry& geom, int scomp, int ncomp, int varType);

void MultiFabPhysBCDomainVel(MultiFab& vel, const Geometry& geom, int dim);

void MultiFabPhysBCMacVel(MultiFab& vel, const Geometry& geom, int dim);

void ZeroEdgevalWalls(std::array<MultiFab, AMREX_SPACEDIM>& edge, const Geometry& geom,
                      int scomp, int ncomp);

void ZeroEdgevalPhysical(std::array<MultiFab, AMREX_SPACEDIM>& edge, const Geometry& geom,
                         int scomp, int ncomp);

void MultiFabElectricBC(MultiFab& efieldCC, const Geometry& geom);

void MultiFabPotentialBC(MultiFab& potential, const Geometry& geom);

void MultiFabPotentialBC_solver(MultiFab& potential, const Geometry& geom);

void MultiFabPhysBCCharge(MultiFab& charge, const Geometry& geom);

void MultiFabPhysBCMacStress(MultiFab& stress, const Geometry& geom, int dim);

void MultiFabPhysBCDomainStress(MultiFab& stress, const amrex::Geometry& geom, int dim);

#endif
