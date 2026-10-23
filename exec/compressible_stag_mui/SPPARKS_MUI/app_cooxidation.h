/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#ifdef APP_CLASS
AppStyle(cooxidation,AppCooxidation)

#else

#ifndef SPK_APP_COOXIDATION_H
#define SPK_APP_COOXIDATION_H

#include "app_lattice.h"

#if defined(USE_AMREX_MPMD)
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MPMD.H>
#endif

namespace SPPARKS_NS {

class AppCooxidation : public AppLattice {
  friend class DiagCooxidation;

 public:
  AppCooxidation(class SPPARKS *, int, char **);
  ~AppCooxidation();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();
  void reaction_summary_log();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *) {}
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  int firsttime;
  // variables on each lattice site
  int *type,*element,*ac1,*ac2,*ac3,*dc1,*dc2,*dc3;
  double *pressure1,*pressure2,*pressure3,*temp;

  int *esites;
  int *echeck;

  // adding adsorption reactions
  // basically the same structure to first-order reactions
  // however, propensity for adsorption depends on partial pressure of gas phase
  // hence different in each cell and we do not use propensity variable
  // adding desorption reactions is exactly same as the first-order reaction case
  int none,ntwo,nads,ndes,ndissocads,nassocdes,nreaction,rxnsumcount;
  double *srate,*drate,*adsrate,*desrate,*dadsrate,*adesrate,*rxnrate;
  bool *ads_is_rate,*dads_is_rate;
  double *spropensity,*dpropensity,*adespropensity,*rxnpropensity;
  int *stype,**dtype,*adstype,*destype,**dadstype,**adestype,**rxntype;
  int *sinput,**dinput,*adsinput,*desinput,**dadsinput,**adesinput,**rxninput;
  int *soutput,**doutput,*adsoutput,*desoutput,**dadsoutput,**adesoutput,**rxnoutput;
  int *scount,*dcount,*adscount,*descount,*dadscount,*adescount,*rxncount;
  int *adsadsorbate,*desdesorbate,*dadsadsorbate,*adesdesorbate,*reactionsorbate;

  struct Event {           // one event for an owned site
    int style;             // reaction style
    int which;             // which reaction of this type
    int jpartner;          // which J neighbor of I is part of event
    int next;              // index of next event for this site
    double propensity;     // propensity of this event
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  double site_propensity_full(int);
  void clear_events(int);
  void add_event(int, int, int, double, int);
  void grow_reactions(int);

#if defined(USE_AMREX_MPMD)

  void amrex_init_agg ();
  void amrex_push_agg(int,char **);
  void amrex_fetch_agg(int,char **);
  double site_propensity_ads(int);

  double amrex_fhd_lattice_size_x;
  double amrex_fhd_lattice_size_y;
  double amrex_kmc_lattice_offset_x;
  double amrex_kmc_lattice_offset_y;

  int ads_wall_dir;

  int update_full;         // ads-only update or full update
  int *bordersite;         // near a sector border: always full rebuild
  int *cellrep;            // one representative site per FHD cell
  double *adsprop_cell;    // per-cell adsorption propensities [nads x ncell]
  double *dadsprop_cell;   // per-cell dissoc-ads propensities [ndissocads x ncell]

  int dir1, dir2;
  int nlocalFHDcell;               // number of FHD cells overlapping with local domain
  amrex::Vector<double> xFHD;      // x-coord of COM of each overlapping FHD cell region
  amrex::Vector<double> yFHD;      // y-coord of COM of each overlapping FHD cell region
  amrex::Vector<int> intval;       // temp int array for MUI push/fetch
  amrex::Vector<double> dblval;    // temp double array for MUI push/fetch
  amrex::Vector<int> localFHDcell; // map from local KMC site to FHD cell
  // array of nlocalFHDcell for all procs (allocated only for domain->me for
  // debugging purposes)
  amrex::Vector<int> nlocalFHDcell_world;
  amrex::MultiFab mf;
  amrex::iMultiFab imf;
  std::unique_ptr<amrex::MPMD::Copier> mpmd_copier;
  std::unique_ptr<amrex::MultiFab> mf2;
  std::unique_ptr<amrex::iMultiFab> imf2;

    void amrex_send_intval();
    void amrex_recv_dblval();
#endif
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

E: Temperature cannot be 0.0 for app cooxidation

UNDOCUMENTED

*/
