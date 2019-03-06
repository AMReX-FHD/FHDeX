#ifndef _hydro_functions_H
#define _hydro_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_MLABecLaplacian.H>

/////////////////////////////////////////////////////////////////////////////////
// in MkAdvMFluxdiv.cpp

void MkAdvMFluxdiv(const std::array<amrex::MultiFab, AMREX_SPACEDIM>& umac,
		   const std::array<amrex::MultiFab, AMREX_SPACEDIM>& m,
		   std::array<amrex::MultiFab, AMREX_SPACEDIM>& m_update,
		   const amrex::Real* dx,
		   const int& increment=0);

void MkAdvSFluxdiv(const std::array<amrex::MultiFab, AMREX_SPACEDIM>& umac,
		   const amrex::MultiFab& m,
		   amrex::MultiFab& m_update,
		   const amrex::Real* dx,
		   const amrex::Geometry& geom,
		   const int& m_comp=0,
		   const int& increment=0);

/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
// in MacProj.cpp

void MacProj (std::array< amrex::MultiFab, AMREX_SPACEDIM >& umac,
	      const amrex::MultiFab& rho,
	      const amrex::Geometry& geom,
	      const bool& full_solve=true);

void ComputeMACSolverRHS (amrex::MultiFab& solverrhs,
			  const amrex::MultiFab& macrhs,
			  const std::array< amrex::MultiFab, AMREX_SPACEDIM >& umac,
			  const amrex::Geometry& geom);

void AvgFaceBcoeffsInv(std::array< amrex::MultiFab, AMREX_SPACEDIM >& facebcoef,
		       const amrex::MultiFab& rhocc);

void SetMacSolverBCs(amrex::MLABecLaplacian& mlabec);

/////////////////////////////////////////////////////////////////////////////////

// in boundary.cpp
void setBC(amrex::MultiFab & u_mac, amrex::MultiFab & v_mac, amrex::MultiFab & w_mac);


/********************************************************************************
 *                                                                              *
 * Functions to fill physical boundary conditions                               *
 *  -> dim_fill_ghost sets for which dimension to fill ghost cells              *
 *                                                                              *
 *******************************************************************************/

void MultiFABPhysBC(amrex::MultiFab & pressure, const amrex::Geometry & geom);
void MultiFABPhysBC(amrex::MultiFab & pressure, int seq_fill_ghost,
                    const amrex::Geometry & geom);
void MultiFABPhysBC(amrex::MultiFab & pressure, const amrex::IntVect & dim_fill_ghost,
                    const amrex::Geometry & geom);

void MultiFABPhysBCDomainVel(amrex::MultiFab & vel, const amrex::Geometry & geom);
void MultiFABPhysBCDomainVel(amrex::MultiFab & vel, int seq_fill_ghost,
                             const amrex::Geometry & geom);
void MultiFABPhysBCDomainVel(amrex::MultiFab & vel, const amrex::IntVect & dim_fill_ghost,
                             const amrex::Geometry & geom);

void MultiFABPhysBCMacVel(amrex::MultiFab & vel, const amrex::Geometry & geom);
void MultiFABPhysBCMacVel(amrex::MultiFab & vel, int seq_fill_ghost,
                          const amrex::Geometry & geom);
void MultiFABPhysBCMacVel(amrex::MultiFab & vel, const amrex::IntVect & dim_fill_ghost,
                          const amrex::Geometry & geom);

#endif
