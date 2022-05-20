#ifndef _main_driver_H_
#define _main_driver_H_

#include "common_functions.H"
#include "gmres_functions.H"


/********************************************************************************
 *                                                                              *
 * Output data                                                                  *
 *   - in WritePlotFile.cpp                                                     *
 *                                                                              *
 *******************************************************************************/

void WritePlotFile(int step,
                   const amrex::Real time,
                   const amrex::Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM >& umac,
                   const MultiFab& rhotot,
                   const MultiFab& rho,
                   const MultiFab& pres,
                   const MultiFab& pressure_jump,
                   const MultiFab& charge,
                   const MultiFab& Epot);


///////////////////////////

#endif
