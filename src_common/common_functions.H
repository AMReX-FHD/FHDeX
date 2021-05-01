#ifndef _common_functions_H_
#define _common_functions_H_

#define LOHI 2
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
#include <AMReX_ArrayLim.H>
#include <AMReX_BC_TYPES.H>

#include "common_functions_F.H"
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

void AverageCCToNode(const MultiFab& cc_in, MultiFab& node_in,
                     int scomp, int ncomp, int bccomp, const Geometry& geom);

void AverageCCToEdge(const MultiFab& cc_in, std::array<MultiFab, NUM_EDGE>& edge_in,
                     int scomp, int ncomp, int bccomp, const Geometry& geom);

///////////////////////////
// in ComputeDivAndGrad.cpp

void ComputeDiv(MultiFab & div, const std::array<MultiFab, AMREX_SPACEDIM> & phi_fc, 
                int start_incomp, int start_outcomp, int ncomp, 
                const Geometry & geom, Real increment);

void ComputeGrad(const MultiFab & phi_in, std::array<MultiFab, AMREX_SPACEDIM> & gphi, 
                 int start_incomp, int start_outcomp, int ncomp, int bccomp, const Geometry & geom,
                 int increment=0);

void ComputeCentredGrad(const MultiFab& phi, std::array<MultiFab, AMREX_SPACEDIM>& gphi,
                        const Geometry& geom);

void ComputeCentredGradFC(std::array<MultiFab, AMREX_SPACEDIM> & phi,
                          MultiFab & gphi,
                          const Geometry & geom);

///////////////////////////
// in Debug.cpp
void PrintMF(const MultiFab& MF);

void outputMFAscii(const MultiFab& output, std::string filename);

///////////////////////////
// in SqrtMF.cpp
void SqrtMF(MultiFab& MF);


///////////////////////////
// in ComputeBasicStats.cpp

void ComputeBasicStats(MultiFab & instant, MultiFab & means, 
                const int incomp, const int outcomp, const int steps);

void OutputVolumeMean(const MultiFab & instant, const int comp, const Real domainVol, std::string filename, const Geometry geom);

Real MaskedSum (const MultiFab & inFab,int comp, const Periodicity& period);

///////////////////////////
// in ComputeAverages.cpp

void WriteHorizontalAverage(const MultiFab& mf_in, const int& dir, const int& incomp,
                            const int& ncomp, const int& step, const Geometry& geom);

void ComputeVerticalAverage(const MultiFab & mf, MultiFab & mf_avg, const Geometry & geom,
                            const int dir, const int incomp, const int ncomp,
                            const int findredist=0);

///////////////////////////
// in MultiFabPhysBC.cpp
// see comments in C++ files for descriptions

void MultiFabPhysBC(MultiFab& data, const Geometry& geom, int scomp, int ncomp, int bccomp, const Real& time=0.);

void MultiFabPhysBCDomainVel(MultiFab& vel, const Geometry& geom, int dim);

void MultiFabPhysBCMacVel(MultiFab& vel, const Geometry& geom, int dim);

void ZeroEdgevalWalls(std::array<MultiFab, AMREX_SPACEDIM>& edge, const Geometry& geom,
                      int scomp, int ncomp);

void ZeroEdgevalPhysical(std::array<MultiFab, AMREX_SPACEDIM>& edge, const Geometry& geom,
                         int scomp, int ncomp);

void MultiFabElectricBC(MultiFab& efieldCC, const Geometry& geom);

void MultiFabPotentialBC(MultiFab& potential, const Geometry& geom);

void MultiFabPotentialBC_solver(MultiFab& potential, const Geometry& geom);

///////////////////////////
// in NormInnerProduct.cpp

void SumStag(const Geometry& geom,
             const std::array<MultiFab, AMREX_SPACEDIM> & m1,
	     const int & comp,
	     Vector<Real> & sum,
	     const bool & divide_by_ncells=false);

void SumCC(const MultiFab & m1,
	   const int & comp,
	   Real & sum,
	   const bool & divide_by_ncells=false);

void StagInnerProd(const Geometry& geom,
                   const std::array<MultiFab, AMREX_SPACEDIM> & m1,
                   const int & comp1,
                   const std::array<MultiFab, AMREX_SPACEDIM> & m2,
                   const int & comp2,
                   std::array<MultiFab, AMREX_SPACEDIM>& mscr,
                   Vector<Real> & prod_val);

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

void StagL2Norm(const Geometry& geom,
                const std::array<MultiFab, AMREX_SPACEDIM> & m1,
		const int & comp,
                std::array<MultiFab, AMREX_SPACEDIM>& mscr,
		Real & norm_l2);

void CCL2Norm(const MultiFab & m1,
	      const int & comp,
              amrex::MultiFab& mscr,
	      Real & norm_l2);


#endif
