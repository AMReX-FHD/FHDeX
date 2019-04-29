#ifndef _main_driver_H_
#define _main_driver_H_

#include "common_functions.H"
#include "gmres_functions.H"

#include "IBCore.H"


/********************************************************************************
 *                                                                              *
 * Advance Solution                                                             *
 *   - in advance.cpp                                                           *
 *                                                                              *
 *******************************************************************************/

void advance(std::array< MultiFab, AMREX_SPACEDIM > & umac,
             std::array< MultiFab, AMREX_SPACEDIM > & umacNew,
             MultiFab & pres, MultiFab & tracer,
             std::array< MultiFab, AMREX_SPACEDIM > & force_ibm,
             const std::array< MultiFab, AMREX_SPACEDIM > & mfluxdiv_predict,
             const std::array< MultiFab, AMREX_SPACEDIM > & mfluxdiv_correct,
             const std::array< MultiFab, AMREX_SPACEDIM > & alpha_fc,
             const MultiFab & beta, const MultiFab & gamma,
# if (AMREX_SPACEDIM == 2)
             const std::array< MultiFab, 1 > & beta_ed,
# elif (AMREX_SPACEDIM == 3)
             const std::array< MultiFab, 3 > & beta_ed,
# endif
             IBCore & ib_core,
             const Geometry geom, const Real & dt);



/********************************************************************************
 *                                                                              *
 * Output data                                                                  *
 *   - in WritePlotFile.cpp                                                     *
 *                                                                              *
 *******************************************************************************/


void WritePlotFile(const int step,
                   const Real time,
                   const Geometry geom,
                   const std::array< MultiFab, AMREX_SPACEDIM > & umac,
                   const MultiFab & tracer,
                   const MultiFab & pres,
                   const std::array< MultiFab, AMREX_SPACEDIM > & force_ibm,
                   const IBParticleContainer & ib_pc);


#endif
