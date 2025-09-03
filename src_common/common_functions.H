#ifndef _common_functions_H_
#define _common_functions_H_

#if (AMREX_SPACEDIM == 2)
#define NUM_EDGE 1
#else
#define NUM_EDGE 3
#endif

#define PRES_BC_COMP 0
#define RHO_BC_COMP  1
#define TEMP_BC_COMP 2
#define EPOT_BC_COMP 3
#define TRAN_BC_COMP 4
#define SPEC_BC_COMP 10
#define SPEC_END_BC_COMP 19
#define MOLFRAC_BC_COMP 20
#define MOLFRAC_END_BC_COMP 29

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_PhysBCFunct.H>
#include "common_namespace.H"

using namespace common;
using namespace amrex;

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
                     int scomp, int ncomp, int bccomp, const Geometry& geom);

void ShiftFaceToCC(const MultiFab& face_in, int face_in_comp,
                   MultiFab& cc_in, int cc_in_comp,
                   int ncomp);

void ShiftCCToFace_onegrid(MultiFab& face_in, int face_comp,
                           const MultiFab& cc_in, int cc_comp, int ncomp);

void AverageCCToNode(const MultiFab& cc_in, MultiFab& node_in,
                     int scomp, int ncomp, int bccomp, const Geometry& geom);

void AverageCCToEdge(const MultiFab& cc_in, std::array<MultiFab, NUM_EDGE>& edge_in,
                     int scomp, int ncomp, int bccomp, const Geometry& geom);

///////////////////////////
// in ComputeDivAndGrad.cpp

void ComputeDiv(MultiFab & div, const std::array<MultiFab, AMREX_SPACEDIM> & phi_fc,
                int start_incomp, int start_outcomp, int ncomp,
                const Geometry & geom, int increment=0);

void ComputeGrad(const MultiFab & phi_in, std::array<MultiFab, AMREX_SPACEDIM> & gphi,
                 int start_incomp, int start_outcomp, int ncomp, int bccomp, const Geometry & geom,
                 int increment=0);

void ComputeCentredGrad(const MultiFab& phi, std::array<MultiFab, AMREX_SPACEDIM>& gphi,
                        const Geometry& geom);

void ComputeCentredGradCompDir(const MultiFab & phi,
                               MultiFab& gphi,
                               int dir,
                               int incomp,
                               int outcomp,
                               const Geometry & geom);

void ComputeCentredGradFC(std::array<MultiFab, AMREX_SPACEDIM> & phi,
                          MultiFab & gphi,
                          const Geometry & geom);

void ComputeLap(const MultiFab & phi_in,
                MultiFab& Lphi_in,
                int incomp,
                int outcomp,
                int numcomp,
                const Geometry & geom);

void ComputeStagLap(std::array<MultiFab, AMREX_SPACEDIM> & phi,
                    std::array<MultiFab, AMREX_SPACEDIM> & Lphi,
                    const Geometry & geom);

void ComputeCurlFaceToEdge(std::array<MultiFab, AMREX_SPACEDIM> & umac_in,
                           std::array<MultiFab, NUM_EDGE> & curl,
                           const Geometry & geom);

void ComputeCurlCC(const MultiFab& vel_in,
                   int incomp,
                   MultiFab& curl_in,
                   int outcomp,
                   const Geometry & geom);

void ComputeDivCC(const MultiFab& vel_in,
                   int incomp,
                   MultiFab& div_in,
                   int outcomp,
                   const Geometry & geom);

///////////////////////////
// in Debug.cpp
void PrintMF(const MultiFab& MF, int comp_lo, int comp_hi);

void outputMFAscii(const MultiFab& output, std::string filename);
void outputMFAscii2(const MultiFab& MF, int comp_lo, int comp_hi, std::string filename);
void outputMFAsciiReduce(const MultiFab& MF, int comp, std::string filename, int rat);

///////////////////////////
// in SqrtMF.cpp
void SqrtMF(MultiFab& MF);


///////////////////////////
// in ComputeBasicStats.cpp

Real ComputeSpatialMean(MultiFab& mf, const int& incomp);

Real ComputeSpatialVariance(MultiFab& mf, const int& incomp);

void ComputeBasicStats(MultiFab & instant, MultiFab & means,
                const int incomp, const int outcomp, const int steps);

void OutputVolumeMean(const MultiFab & instant, const int comp, const Real domainVol, std::string filename, const Geometry geom);

Real MaskedSum (const MultiFab & inFab,int comp, const Periodicity& period);

///////////////////////////
// in ComputeAverages.cpp

void WriteHorizontalAverage(const MultiFab& mf_in, const int& dir, const int& incomp,
                            const int& ncomp, const int& step, const Geometry& geom,
                            const std::string& file_prefix = "havg");

void WriteHorizontalAverageToMF(const MultiFab& mf_in, MultiFab& mf_out,
                                const int& dir, const int& incomp,
                                const int& ncomp, int outcomp=-1);

void ComputeVerticalAverage(const MultiFab & mf, MultiFab & mf_flat,
                            const int& dir, const int& incomp, const int& ncomp,
                            const int& slablo=-1, const int& slabhi=99999);

void ExtractSlice(const MultiFab & mf, MultiFab & mf_slice,
                  const int dir, const int slice, const int incomp, const int ncomp);

void ExtractXPencil(const MultiFab& mf, MultiFab& mf_pencil,
                    const int pencily, const int pencilz,
                    const int incomp, const int ncomp);

///////////////////////////
// in MultiFabPhysBC.cpp
// see comments in C++ files for descriptions

void MultiFabPhysBC(MultiFab& data, const Geometry& geom, int scomp, int ncomp, int bccomp, const Real& time=0.);

void MultiFabPhysBCDomainVel(MultiFab& vel, const Geometry& geom, int dim);

void MultiFabPhysBCMacVel(MultiFab& vel, const Geometry& geom, int dim, int is_inhomogeneous=0);

void ZeroEdgevalWalls(std::array<MultiFab, AMREX_SPACEDIM>& edge, const Geometry& geom,
                      int scomp, int ncomp);

void ZeroEdgevalPhysical(std::array<MultiFab, AMREX_SPACEDIM>& edge, const Geometry& geom,
                         int scomp, int ncomp);

void MultiFabElectricBC(MultiFab& efieldCC, const Geometry& geom);

void MultiFabPotentialBC(MultiFab& potential, const Geometry& geom);

void MultiFabPotentialBC_solver(MultiFab& potential, const Geometry& geom);

///////////////////////////
// in NormInnerProduct.cpp

void SumStag(const std::array<MultiFab, AMREX_SPACEDIM> & m1,
    Vector<Real> & sum,
    const bool & divide_by_ncells=false);

void SumEdge(const std::array<MultiFab, NUM_EDGE> & m1,
    Vector<Real> & sum,
    const bool & divide_by_ncells=false);

void SumCC(const MultiFab & m1,
    const int & comp,
    Real & sum,
    const bool & divide_by_ncells=false);

void StagInnerProd(const std::array<MultiFab, AMREX_SPACEDIM> & m1,
    const int & comp1,
    const std::array<MultiFab, AMREX_SPACEDIM> & m2,
    const int & comp2,
    std::array<MultiFab, AMREX_SPACEDIM>& mscr,
    Vector<Real> & prod_val);

void EdgeInnerProd(const std::array<MultiFab, NUM_EDGE>& m1,
    const int& comp1,
    const std::array<MultiFab, NUM_EDGE>& m2,
    const int& comp2,
    std::array<MultiFab, NUM_EDGE>& mscr,
    amrex::Vector<amrex::Real>& prod_val);

void CCInnerProd(const MultiFab & m1,
    const int & comp1,
    const MultiFab & m2,
    const int & comp2,
    amrex::MultiFab& mscr,
    Real & prod_val);

void CCMoments(const MultiFab & m1,
    const int & comp1,
    amrex::MultiFab& mscr,
    const int & power,
    Real & prod_val);

void FCMoments(const std::array<MultiFab, AMREX_SPACEDIM>& m1,
    const amrex::Vector<int>& comps,
    std::array<MultiFab, AMREX_SPACEDIM>&  mscr,
    const int& power,
    amrex::Vector<amrex::Real>& prod_val);

void StagL2Norm(const std::array<MultiFab, AMREX_SPACEDIM> & m1,
    const int & comp,
    std::array<MultiFab, AMREX_SPACEDIM>& mscr,
    Real & norm_l2);

void EdgeL2Norm(const std::array<MultiFab, NUM_EDGE>& m1,
    const int& comp,
    std::array<MultiFab, NUM_EDGE>& mscr,
    Real& norm_l2);

void CCL2Norm(const MultiFab & m1,
    const int & comp,
    amrex::MultiFab& mscr,
    Real & norm_l2);

///////////////////////////
// in InterpCoarsen.cpp
void FaceFillCoarse(Vector<std::array< MultiFab, AMREX_SPACEDIM >>& mf, int map);
void FaceFillCoarse(std::array<MultiFab, AMREX_SPACEDIM>* & mf, int map);
void FaceFillCoarse(std::array<MultiFab*, AMREX_SPACEDIM> & mfCoarse,std::array<MultiFab*, AMREX_SPACEDIM> & mfFine, int map);


void CellFillCoarse(Vector<MultiFab>& mf, Vector<Geometry> geom);
void CellFillCoarse(Vector<MultiFab>& mf, Geometry* geom);
void CellFillCoarse(MultiFab* & mf, Geometry* geom);
void CellFillCoarse(const MultiFab & mfFine, Geometry geomFine, MultiFab & mfCoarse, Geometry geomCoarse);

#endif