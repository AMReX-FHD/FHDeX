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
 * Molecular Dynamics (MD) helper functions for immersed boundary markers       *
 *   - in ib_md.cpp                                                             *
 *                                                                              *
 *******************************************************************************/

void constrain_ibm_marker(IBMarkerContainer & ib_mc, int ib_lev, int component);

void anchor_coupling_markers(IBMarkerContainer & ib_mc, int ib_lev, int component);

void anchor_first_marker(IBMarkerContainer & ib_mc, int ib_lev, int component);

Real theta(Real amp_ramp, Real time, int i_ib, int index_marker);

void update_bdy_marker(const std::map<std::tuple<int, int>, double> & bond_map,
                       const std::map<int, std::vector<int>> & bond_neighbors,
                       Real time,
                       IBMarkerContainer & ib_mc, int ib_lev,
                       int component, bool pred_pos,
                       const Geometry & geom);

void update_ibm_marker(const RealVect & driv_u, Real driv_amp, Real time,
                       IBMarkerContainer & ib_mc, int ib_lev,
                       int component, bool pred_pos,
                       const Geometry & geom);

void yeax_ibm_marker(Real mot, IBMarkerContainer & ib_mc, int ib_lev,
                     int component_src, int component_dest);

Vector<RealVect> equil_pos(int i_ib, Real time, const Geometry & geom,
    const RealVect & x_anchor, const RealVect & e_anchor);



/********************************************************************************
 *                                                                              *
 * Advance Solution                                                             *
 *   - in advance.cpp and advaance_CN.cpp                                       *
 *                                                                              *
 *******************************************************************************/

void advance(std::array<MultiFab, AMREX_SPACEDIM> & umac,
             std::array<MultiFab, AMREX_SPACEDIM> & umacNew,
             MultiFab & pres,
             IBMarkerContainer & ib_mc,
 //             const std::map<std::tuple<int, int>, double> & bond_map,
 //             const std::map<int, std::vector<int>> & bond_neighbors,
             const std::array<MultiFab, AMREX_SPACEDIM> & mfluxdiv_predict,
             const std::array<MultiFab, AMREX_SPACEDIM> & mfluxdiv_correct,
             std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
             std::array<MultiFab, AMREX_SPACEDIM> & force_ib,
             const MultiFab& beta, const MultiFab & gamma,
# if (AMREX_SPACEDIM == 2)
             const std::array< MultiFab, 1 > & beta_ed,
# endif
# if (AMREX_SPACEDIM == 3)
             const std::array< MultiFab, 3 > & beta_ed,
# endif
             const Geometry geom, const Real & dt, Real time);



void advance_CN(std::array<MultiFab, AMREX_SPACEDIM> & umac,
                std::array<MultiFab, AMREX_SPACEDIM> & umacNew,
                MultiFab & pres,
                IBMarkerContainer & ib_mc,
                const std::array<MultiFab, AMREX_SPACEDIM> & mfluxdiv_predict,
                const std::array<MultiFab, AMREX_SPACEDIM> & mfluxdiv_correct,
                std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                std::array<MultiFab, AMREX_SPACEDIM> & force_ib,
                const MultiFab& beta, const MultiFab & gamma,
# if (AMREX_SPACEDIM == 2)
                const std::array<MultiFab, 1> & beta_ed,
# endif
# if (AMREX_SPACEDIM == 3)
                const std::array<MultiFab, 3> & beta_ed,
# endif
                const Geometry geom, const Real & dt, Real time);



void advance_stokes(std::array<MultiFab, AMREX_SPACEDIM> & umac,
                    std::array<MultiFab, AMREX_SPACEDIM> & umacNew,
                    MultiFab & pres,
                    IBMarkerContainer & ib_mc,
                    const std::map<std::tuple<int, int>, double> & bond_map,
                    const std::map<int, std::vector<int>> & bond_neighbors,
                    const std::array<MultiFab, AMREX_SPACEDIM> & mfluxdiv_predict,
                    const std::array<MultiFab, AMREX_SPACEDIM> & mfluxdiv_correct,
                    std::array<MultiFab, AMREX_SPACEDIM> & alpha_fc,
                    std::array<MultiFab, AMREX_SPACEDIM> & force_ib,
                    const MultiFab& beta, const MultiFab & gamma,
# if (AMREX_SPACEDIM == 2)
                    const std::array<MultiFab, 1> & beta_ed,
# endif
# if (AMREX_SPACEDIM == 3)
                    const std::array<MultiFab, 3> & beta_ed,
# endif
                    const Geometry geom, const Real & dt, Real time);



void constrain_ibm_marker(IBMarkerContainer & ib_mc, int ib_lev, int component);



void update_ibm_marker(const RealVect & driv_u, Real driv_amp, Real time,
                       IBMarkerContainer & ib_mc, int ib_lev,
                       int component, bool pred_pos);



/********************************************************************************
 *                                                                              *
 * Output data                                                                  *
 *   - in WritePlotFile.cpp                                                     *
 *                                                                              *
 *******************************************************************************/

void WritePlotFile(const int step,
                   const Real time,
                   const Geometry geom,
                   std::array< MultiFab, AMREX_SPACEDIM > & umac,
                   std::array< MultiFab, AMREX_SPACEDIM > & umac_avg,
                   std::array< MultiFab, AMREX_SPACEDIM > & force_ib,
                   const MultiFab & pres,
                   const IBMarkerContainer & ib_pc);

#endif
