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

void ComputeDiv(MultiFab & div, const std::array<MultiFab, AMREX_SPACEDIM> & phi_fc, 
                int start_incomp, int start_outcomp, int ncomp, 
                const Geometry & geom, int increment);

void ComputeGrad(const MultiFab & phi, std::array<MultiFab, AMREX_SPACEDIM> & gphi, 
                 int start_incomp, int start_outcomp, int ncomp, const Geometry & geom);

void ComputeCentredGrad(const MultiFab& phi, std::array<MultiFab, AMREX_SPACEDIM>& gphi, const Geometry& geom);

///////////////////////////
// in Debug.cpp
void PrintMF(const MultiFab& MF);

///////////////////////////
// in FindCoords.cpp

void FindFaceCoords(std::array< MultiFab, AMREX_SPACEDIM >& RealFaceCoords, Geometry geom);

void FindCenterCoords(MultiFab& RealCenterCoords, Geometry geom);

void FindNodalValues(const MultiFab& inFab, MultiFab& outFab, const MultiFab& ccFab);

///////////////////////////
// in MultiFABFillRandom.cpp

void MultiFABFillRandom(MultiFab& mf, const int& comp, const amrex::Real& variance, const Geometry& geom);

void MultiFABFillRandomHack(MultiFab& mf, const int& comp, const amrex::Real& variance, const Geometry& geom);


///////////////////////////
// in MultiFABPhysBC.cpp

/********************************************************************************
 *                                                                              *
 * Functions to fill physical boundary conditions                               *
 *  -> dim_fill_ghost sets for which dimension to fill ghost cells              *
 *                                                                              *
 *******************************************************************************/

void MultiFABPhysBC(MultiFab & pressure, const Geometry & geom);
void MultiFABPhysBC(MultiFab & pressure, int seq_fill_ghost, const Geometry & geom);
void MultiFABPhysBC(MultiFab & pressure, const IntVect & dim_fill_ghost, const Geometry & geom);

void MultiFABPhysBCDomainVel(MultiFab & vel, const Geometry & geom, int dim);
void MultiFABPhysBCDomainVel(MultiFab & vel, int seq_fill_ghost, const Geometry & geom, int dim);
void MultiFABPhysBCDomainVel(MultiFab & vel, const IntVect & dim_fill_ghost, const Geometry & geom, int dim);

void MultiFABPhysBCMacVel(MultiFab & vel, const Geometry & geom, int dim);
void MultiFABPhysBCMacVel(MultiFab & vel, int seq_fill_ghost, const Geometry & geom, int dim);
void MultiFABPhysBCMacVel(MultiFab & vel, const IntVect & dim_fill_ghost, const Geometry & geom, int dim);

//Functions to modify the stress near the boundaries due to peskin particles

void MultiFABPhysBCMacStress(MultiFab & stress, const Geometry & geom, int dim);
void MultiFABPhysBCMacStress(MultiFab & stress, int seq_fill_ghost, const Geometry & geom, int dim);
void MultiFABPhysBCMacStress(MultiFab & stress, const IntVect & dim_fill_ghost, const Geometry & geom, int dim);

void MultiFABPhysBCDomainStress(MultiFab & stress, const amrex::Geometry & geom, int dim);
void MultiFABPhysBCDomainStress(MultiFab & stress, int seq_fill_ghost, const Geometry & geom, int dim);
void MultiFABPhysBCDomainStress(MultiFab & stress, const IntVect & dim_fill_ghost, const Geometry & geom, int dim);

///////////////////////////
// in outputAscii.cpp


void outputMFAscii(const MultiFab& output, std::string filename);


#endif
