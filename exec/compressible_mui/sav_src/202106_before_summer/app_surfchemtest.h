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
AppStyle(surfchemtest,AppSurfchemtest)

#else

#ifndef SPK_APP_SURFCHEMTEST_H
#define SPK_APP_SURFCHEMTEST_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppSurfchemtest : public AppLattice {
  friend class DiagSurfchemtest;

 public:
  AppSurfchemtest(class SPPARKS *, int, char **);
  ~AppSurfchemtest();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *) {}
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  int engstyle;
  int firsttime;
  // variables on each lattice site
  int *type,*element,*ac1,*ac2,*ac3,*ac4,*ac5,*dc1,*dc2,*dc3,*dc4,*dc5;
  double *density1,*density2,*density3,*density4,*density5,*temp;

  int *esites;
  int *echeck;

  // adding adsorption reactions
  // basically the same structure to first-order reactions
  // however, propensity for adsorption depends on number density of gas phase
  // hence different in each cell and we do not use propensity variable
  // adding desorption reactions is exactly same as the first-order reaction case
  int none,ntwo,nthree,nads,ndes;
  double *srate,*drate,*trate,*adsrate,*desrate;
  double *spropensity,*dpropensity,*tpropensity,*despropensity;
  int *stype,**dtype,**ttype,*adstype,*destype;
  int *sinput,**dinput,**tinput,*adsinput,*desinput;
  int *soutput,**doutput,**toutput,*adsoutput,*desoutput;
  int *scount,*dcount,*tcount,*adscount,*descount;

  struct Event {           // one event for an owned site
    int style;             // reaction style = SINGLE,DOUBLE,TRIPLE
    int which;             // which reaction of this type
    int jpartner,kpartner; // which J,K neighbors of I are part of event
    int next;              // index of next event for this site
    double propensity;     // propensity of this event
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list

  void clear_events(int);
  void add_event(int, int, int, double, int, int);
  void grow_reactions(int);

  void mui_init_agg();
  void mui_push(int,char **);
  void mui_fetch(int,char **);
  void mui_push_agg(int,char **);
  void mui_fetch_agg(int,char **);
  double mui_fhd_lattice_size_x;
  double mui_fhd_lattice_size_y;
  double mui_kmc_lattice_offset_x;
  double mui_kmc_lattice_offset_y;

  int nlocalFHDcell;        // number of FHD cells overlapping with local domain
  double *xFHD;             // x-coord of COM of each overlapping FHD cell region
  double *yFHD;             // y-coord of COM of each overlapping FHD cell region
  int *MUIintval;           // temp int array for MUI push/fetch
  double *MUIdblval;        // temp double array for MUI push/fetch
  int *localFHDcell;        // map from local KMC site to FHD cell
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

E: Temperature cannot be 0.0 for app surfchemtest

UNDOCUMENTED

*/
