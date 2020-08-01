#ifndef _main_driver_H_
#define _main_driver_H_

#include "common_functions.H"
#include "gmres_functions.H"
#include "IBMarkerContainer.H"



//_______________________________________________________________________________
// Abbreviations:
using TileIndex       = IBMarkerContainer::TileIndex;
using MarkerListIndex = IBMarkerContainer::MarkerListIndex;
using ParticleType    = IBMarkerContainer::ParticleType;
using AoS             = IBMarkerContainer::AoS;




/********************************************************************************
 *                                                                              *
 * Output data                                                                  *
 *   - in WritePlotFile.cpp                                                     *
 *                                                                              *
 *******************************************************************************/

void WritePlotFile(const int step,
                   const Real time,
                   const Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM >& umac,
                   const MultiFab& pres,
                   const IBMarkerContainer & ib_pc);

///////////////////////////

#endif
