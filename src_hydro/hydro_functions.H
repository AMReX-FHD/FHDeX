#ifndef _hydro_functions_H_
#define _hydro_functions_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>
#include <AMReX_MLABecLaplacian.H>

#include "common_functions.H"
#include "gmres_functions.H"

using namespace amrex;

/////////////////////////////////////////////////////////////////////////////////
// in AddMFluctuations.cpp

// used to add noise to an initial momentum field
void addMomFluctuations(std::array< MultiFab, AMREX_SPACEDIM >&,
                        const MultiFab&, const MultiFab&,
                        const Real&,
                        const Geometry& geom);
    
// used to add noise to an initial momentum field
// called by addMomFluctuations
void addMomFluctuations_stag(std::array< MultiFab, AMREX_SPACEDIM >&,
                             const std::array< MultiFab, AMREX_SPACEDIM >&,
                             const std::array< MultiFab, AMREX_SPACEDIM >&,
                             const Real&,
                             const Geometry& geom);

/////////////////////////////////////////////////////////////////////////////////
// in ConvertMToUmac.cpp

void ConvertMToUmac(const std::array< MultiFab, AMREX_SPACEDIM >& s_fc,
                    std::array< MultiFab, AMREX_SPACEDIM >& umac,
                    std::array< MultiFab, AMREX_SPACEDIM >& m,
                    int m_to_umac);

/////////////////////////////////////////////////////////////////////////////////
// in MkAdvMFluxdiv.cpp

void MkAdvMFluxdiv(const std::array<MultiFab, AMREX_SPACEDIM> & umac,
		   const std::array<MultiFab, AMREX_SPACEDIM> & m,
		   std::array<MultiFab, AMREX_SPACEDIM> & m_update,
		   const Real * dx,
		   const int & increment=0);

void MkAdvSFluxdiv_cc(const std::array<MultiFab, AMREX_SPACEDIM> & umac_in,
                      const MultiFab & s_in,
                      MultiFab & s_update_in,
                      const Geometry & geom,
                      const int & scomp,
                      const int & ncomp,
                      const int & increment);

void MkAdvSFluxdiv(const std::array<MultiFab, AMREX_SPACEDIM> & umac_in,
                   const std::array<MultiFab, AMREX_SPACEDIM> & s_fc_in,
                   MultiFab & s_update_in,
                   const Geometry & geom,
                   const int & scomp,
                   const int & ncomp,
                   const int & increment);

/////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
// in MacProj_hydro.cpp

void MacProj_hydro(std::array<MultiFab, AMREX_SPACEDIM> & umac,
                   const MultiFab & rho,
                   const Geometry & geom,
                   const bool & full_solve=true);

void ComputeMACSolverRHS(MultiFab & solverrhs,
			 const MultiFab & macrhs,
			 const std::array<MultiFab, AMREX_SPACEDIM> & umac,
			 const Geometry & geom);

void SetMacSolverBCs(MLABecLaplacian & mlabec);

/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
// in advance.cpp

void advanceStokes(  std::array< MultiFab, AMREX_SPACEDIM >& umac,
	       MultiFab& pres,
	       const std::array< MultiFab, AMREX_SPACEDIM >& stochMfluxdiv,
	       std::array< MultiFab, AMREX_SPACEDIM >& sourceTerms,
	       std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
	       MultiFab& beta,
               MultiFab& gamma,
	       std::array< MultiFab, NUM_EDGE >& beta_ed,
	       const Geometry geom, const Real& dt);

void advanceLowMach(  std::array< MultiFab, AMREX_SPACEDIM >& umac,
	       std::array< MultiFab, AMREX_SPACEDIM >& umacNew,
	       MultiFab& pres, MultiFab& tracer,
	       const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_predict,
	       const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_correct,
	       std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
	       const MultiFab& beta, const MultiFab& gamma,
	       const std::array< MultiFab, NUM_EDGE >& beta_ed,
	       const Geometry geom, const Real& dt);

/////////////////////////////////////////////////////////////////////////////////
// in Vorticity.cpp

void MagVort(const std::array< MultiFab, AMREX_SPACEDIM >& umac,
             MultiFab& magvort,
             const Geometry& geom,
             int outcomp=0);

#endif
