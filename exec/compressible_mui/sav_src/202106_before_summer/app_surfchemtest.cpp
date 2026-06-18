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

#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "app_surfchemtest.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "mui.h"
#include <vector>

using namespace mui;

using namespace SPPARKS_NS;

using namespace std;

enum{NOOP,SITEA,SITEB,SITEC};
enum{VACANCY,SPEC1,SPEC2,SPEC3,SPEC4,SPEC5}; // removed ZERO and moved VACANCY to first item // same as DiagSurfchemtest

#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppSurfchemtest::AppSurfchemtest(SPPARKS *spk, int narg, char **arg) :
  AppLattice(spk,narg,arg)
{
  ninteger = 12; // type, element, ac1, ac2, ac3, ac4, ac5, dc1, dc2, dc3, dc4, dc5  (number changes due to ads/des)
  ndouble = 6;  // density1, density2, density3, density4, density5 (number densities in gas phase), temp (local temperature)
  delpropensity = 1;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;

  create_arrays();

  if (narg != 1) error->all(FLERR,"Illegal app_style command");

  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;

  // reaction lists

  none = ntwo = nthree = nads = ndes = 0;
  srate = drate = trate = adsrate = desrate = NULL;
  spropensity = dpropensity = tpropensity = despropensity = NULL;   // no adspropensity
  stype = sinput = soutput = NULL;
  dtype = dinput = doutput = NULL;
  ttype = tinput = toutput = NULL;
  adstype = adsinput = adsoutput = NULL;
  destype = desinput = desoutput = NULL;
  scount = dcount = tcount = adscount = descount = NULL;

  mui_fhd_lattice_size_x = mui_fhd_lattice_size_y = -1.;
  mui_kmc_lattice_offset_x = mui_kmc_lattice_offset_y = 0.;

  xFHD = NULL;
  yFHD = NULL;
  MUIintval = NULL;
  MUIdblval = NULL;
  localFHDcell = NULL;
}

/* ---------------------------------------------------------------------- */

AppSurfchemtest::~AppSurfchemtest()
{
  delete [] esites;
  delete [] echeck;
  memory->sfree(events);
  memory->destroy(firstevent);

  memory->destroy(srate);
  memory->destroy(drate);
  memory->destroy(trate);
  memory->destroy(adsrate);
  memory->destroy(desrate);
  memory->destroy(spropensity);
  memory->destroy(dpropensity);
  memory->destroy(tpropensity);
  memory->destroy(despropensity);
  memory->destroy(stype);
  memory->destroy(sinput);
  memory->destroy(soutput);
  memory->destroy(dtype);
  memory->destroy(dinput);
  memory->destroy(doutput);
  memory->destroy(ttype);
  memory->destroy(tinput);
  memory->destroy(toutput);
  memory->destroy(adstype);
  memory->destroy(adsinput);
  memory->destroy(adsoutput);
  memory->destroy(destype);
  memory->destroy(desinput);
  memory->destroy(desoutput);
  memory->destroy(scount);
  memory->destroy(dcount);
  memory->destroy(tcount);
  memory->destroy(adscount);
  memory->destroy(descount);

  delete [] xFHD;
  delete [] yFHD;
  delete [] MUIintval;
  delete [] MUIdblval;
  delete [] localFHDcell;
}

/* ---------------------------------------------------------------------- */

void AppSurfchemtest::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"event") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal event command");
    int rstyle = atoi(arg[0]);
    grow_reactions(rstyle);

    if (rstyle == 1) {
      if (narg != 5) error->all(FLERR,"Illegal event command");

      if (strcmp(arg[1],"siteA") == 0) stype[none] = SITEA;
      else if (strcmp(arg[1],"siteB") == 0) stype[none] = SITEB;
      else if (strcmp(arg[1],"siteC") == 0) stype[none] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[2],"spec1") == 0) sinput[none] = SPEC1;
      else if (strcmp(arg[2],"spec2") == 0) sinput[none] = SPEC2;
      else if (strcmp(arg[2],"spec3") == 0) sinput[none] = SPEC3;
      else if (strcmp(arg[2],"spec4") == 0) sinput[none] = SPEC4;
      else if (strcmp(arg[2],"spec5") == 0) sinput[none] = SPEC5;
      else if (strcmp(arg[2],"vac") == 0) sinput[none] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      srate[none] = atof(arg[3]);

      if (strcmp(arg[4],"spec1") == 0) soutput[none] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) soutput[none] = SPEC2;
      else if (strcmp(arg[4],"spec3") == 0) soutput[none] = SPEC3;
      else if (strcmp(arg[4],"spec4") == 0) soutput[none] = SPEC4;
      else if (strcmp(arg[4],"spec5") == 0) soutput[none] = SPEC5;
      else if (strcmp(arg[4],"vac") == 0) soutput[none] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      none++;

    } else if (rstyle == 2) {
      if (narg != 8) error->all(FLERR,"Illegal event command");

      if (strcmp(arg[1],"siteA") == 0) dtype[ntwo][0] = SITEA;
      else if (strcmp(arg[1],"siteB") == 0) dtype[ntwo][0] = SITEB;
      else if (strcmp(arg[1],"siteC") == 0) dtype[ntwo][0] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[2],"siteA") == 0) dtype[ntwo][1] = SITEA;
      else if (strcmp(arg[2],"siteB") == 0) dtype[ntwo][1] = SITEB;
      else if (strcmp(arg[2],"siteC") == 0) dtype[ntwo][1] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[3],"spec1") == 0) dinput[ntwo][0] = SPEC1;
      else if (strcmp(arg[3],"spec2") == 0) dinput[ntwo][0] = SPEC2;
      else if (strcmp(arg[3],"spec3") == 0) dinput[ntwo][0] = SPEC3;
      else if (strcmp(arg[3],"spec4") == 0) dinput[ntwo][0] = SPEC4;
      else if (strcmp(arg[3],"spec5") == 0) dinput[ntwo][0] = SPEC5;
      else if (strcmp(arg[3],"vac") == 0) dinput[ntwo][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[4],"spec1") == 0) dinput[ntwo][1] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) dinput[ntwo][1] = SPEC2;
      else if (strcmp(arg[4],"spec3") == 0) dinput[ntwo][1] = SPEC3;
      else if (strcmp(arg[4],"spec4") == 0) dinput[ntwo][1] = SPEC4;
      else if (strcmp(arg[4],"spec5") == 0) dinput[ntwo][1] = SPEC5;
      else if (strcmp(arg[4],"vac") == 0) dinput[ntwo][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      drate[ntwo] = atof(arg[5]);

      if (strcmp(arg[6],"spec1") == 0) doutput[ntwo][0] = SPEC1;
      else if (strcmp(arg[6],"spec2") == 0) doutput[ntwo][0] = SPEC2;
      else if (strcmp(arg[6],"spec3") == 0) doutput[ntwo][0] = SPEC3;
      else if (strcmp(arg[6],"spec4") == 0) doutput[ntwo][0] = SPEC4;
      else if (strcmp(arg[6],"spec5") == 0) doutput[ntwo][0] = SPEC5;
      else if (strcmp(arg[6],"vac") == 0) doutput[ntwo][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[7],"spec1") == 0) doutput[ntwo][1] = SPEC1;
      else if (strcmp(arg[7],"spec2") == 0) doutput[ntwo][1] = SPEC2;
      else if (strcmp(arg[7],"spec3") == 0) doutput[ntwo][1] = SPEC3;
      else if (strcmp(arg[7],"spec4") == 0) doutput[ntwo][1] = SPEC4;
      else if (strcmp(arg[7],"spec5") == 0) doutput[ntwo][1] = SPEC5;
      else if (strcmp(arg[7],"vac") == 0) doutput[ntwo][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      ntwo++;

    } else if (rstyle == 3) {
      if (narg != 11) error->all(FLERR,"Illegal event command");

      if (strcmp(arg[1],"siteA") == 0) ttype[nthree][0] = SITEA;
      else if (strcmp(arg[1],"siteB") == 0) ttype[nthree][0] = SITEB;
      else if (strcmp(arg[1],"siteC") == 0) ttype[nthree][0] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[2],"siteA") == 0) ttype[nthree][1] = SITEA;
      else if (strcmp(arg[2],"siteB") == 0) ttype[nthree][1] = SITEB;
      else if (strcmp(arg[2],"siteC") == 0) ttype[nthree][1] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[3],"siteA") == 0) ttype[nthree][2] = SITEA;
      else if (strcmp(arg[3],"siteB") == 0) ttype[nthree][2] = SITEB;
      else if (strcmp(arg[3],"siteC") == 0) ttype[nthree][2] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[4],"spec1") == 0) tinput[nthree][0] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) tinput[nthree][0] = SPEC2;
      else if (strcmp(arg[4],"spec3") == 0) tinput[nthree][0] = SPEC3;
      else if (strcmp(arg[4],"spec4") == 0) tinput[nthree][0] = SPEC4;
      else if (strcmp(arg[4],"spec5") == 0) tinput[nthree][0] = SPEC5;
      else if (strcmp(arg[4],"vac") == 0) tinput[nthree][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[5],"spec1") == 0) tinput[nthree][1] = SPEC1;
      else if (strcmp(arg[5],"spec2") == 0) tinput[nthree][1] = SPEC2;
      else if (strcmp(arg[5],"spec3") == 0) tinput[nthree][1] = SPEC3;
      else if (strcmp(arg[5],"spec4") == 0) tinput[nthree][1] = SPEC4;
      else if (strcmp(arg[5],"spec5") == 0) tinput[nthree][1] = SPEC5;
      else if (strcmp(arg[5],"vac") == 0) tinput[nthree][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[6],"spec1") == 0) tinput[nthree][2] = SPEC1;
      else if (strcmp(arg[6],"spec2") == 0) tinput[nthree][2] = SPEC2;
      else if (strcmp(arg[6],"spec3") == 0) tinput[nthree][2] = SPEC3;
      else if (strcmp(arg[6],"spec4") == 0) tinput[nthree][2] = SPEC4;
      else if (strcmp(arg[6],"spec5") == 0) tinput[nthree][2] = SPEC5;
      else if (strcmp(arg[6],"vac") == 0) tinput[nthree][2] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      trate[nthree] = atof(arg[7]);

      if (strcmp(arg[8],"spec1") == 0) toutput[nthree][0] = SPEC1;
      else if (strcmp(arg[8],"spec2") == 0) toutput[nthree][0] = SPEC2;
      else if (strcmp(arg[8],"spec3") == 0) toutput[nthree][0] = SPEC3;
      else if (strcmp(arg[8],"spec4") == 0) toutput[nthree][0] = SPEC4;
      else if (strcmp(arg[8],"spec5") == 0) toutput[nthree][0] = SPEC5;
      else if (strcmp(arg[8],"vac") == 0) toutput[nthree][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[9],"spec1") == 0) toutput[nthree][1] = SPEC1;
      else if (strcmp(arg[9],"spec2") == 0) toutput[nthree][1] = SPEC2;
      else if (strcmp(arg[9],"spec3") == 0) toutput[nthree][1] = SPEC3;
      else if (strcmp(arg[9],"spec4") == 0) toutput[nthree][1] = SPEC4;
      else if (strcmp(arg[9],"spec5") == 0) toutput[nthree][1] = SPEC5;
      else if (strcmp(arg[9],"vac") == 0) toutput[nthree][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[10],"spec1") == 0) toutput[nthree][2] = SPEC1;
      else if (strcmp(arg[10],"spec2") == 0) toutput[nthree][2] = SPEC2;
      else if (strcmp(arg[10],"spec3") == 0) toutput[nthree][2] = SPEC3;
      else if (strcmp(arg[10],"spec4") == 0) toutput[nthree][2] = SPEC4;
      else if (strcmp(arg[10],"spec5") == 0) toutput[nthree][2] = SPEC5;
      else if (strcmp(arg[10],"vac") == 0) toutput[nthree][2] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      nthree++;

    } else if (rstyle == 4) {   // adsorption
      if (narg != 5) error->all(FLERR,"Illegal event command");

      if (strcmp(arg[1],"siteA") == 0) adstype[nads] = SITEA;
      else if (strcmp(arg[1],"siteB") == 0) adstype[nads] = SITEB;
      else if (strcmp(arg[1],"siteC") == 0) adstype[nads] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[2],"spec1") == 0)
        error->all(FLERR,"rstyle=4 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[2],"spec2") == 0)
        error->all(FLERR,"rstyle=4 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[2],"spec3") == 0)
        error->all(FLERR,"rstyle=4 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[2],"spec4") == 0)
        error->all(FLERR,"rstyle=4 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[2],"spec5") == 0)
        error->all(FLERR,"rstyle=4 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[2],"vac") == 0) adsinput[nads] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      adsrate[nads] = atof(arg[3]);

      // want to make sure that darray[0], darray[1], darray[2], darray[3], darray[4]
      // correspond to number densities of spec1, spec2, spec3, spec4, spec5
      if (strcmp(arg[4],"spec1") == 0) adsoutput[nads] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) adsoutput[nads] = SPEC2;
      else if (strcmp(arg[4],"spec3") == 0) adsoutput[nads] = SPEC3;
      else if (strcmp(arg[4],"spec4") == 0) adsoutput[nads] = SPEC4;
      else if (strcmp(arg[4],"spec5") == 0) adsoutput[nads] = SPEC5;
      else if (strcmp(arg[4],"vac") == 0)
        error->all(FLERR,"rstyle=4 only allows vac->spec1/2/3/4/5");
      else error->all(FLERR,"Illegal event command");

      nads++;

    } else if (rstyle == 5) {   // desorption
      if (narg != 5) error->all(FLERR,"Illegal event command");

      if (strcmp(arg[1],"siteA") == 0) destype[ndes] = SITEA;
      else if (strcmp(arg[1],"siteB") == 0) destype[ndes] = SITEB;
      else if (strcmp(arg[1],"siteC") == 0) destype[ndes] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[2],"spec1") == 0) desinput[ndes] = SPEC1;
      else if (strcmp(arg[2],"spec2") == 0) desinput[ndes] = SPEC2;
      else if (strcmp(arg[2],"spec3") == 0) desinput[ndes] = SPEC3;
      else if (strcmp(arg[2],"spec4") == 0) desinput[ndes] = SPEC4;
      else if (strcmp(arg[2],"spec5") == 0) desinput[ndes] = SPEC5;
      else if (strcmp(arg[2],"vac") == 0)
        error->all(FLERR,"rstyle=5 only allows spec1/2/3/4/5->vac");
      else error->all(FLERR,"Illegal event command");

      desrate[ndes] = atof(arg[3]);

      // want to make sure that darray[0], darray[1], darray[2], darray[3], darray[4]
      // correspond to number densities of spec1, spec2, spec3, spec4, spec5
      if (strcmp(arg[4],"spec1") == 0)
        error->all(FLERR,"rstyle=5 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[4],"spec2") == 0)
        error->all(FLERR,"rstyle=5 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[4],"spec3") == 0)
        error->all(FLERR,"rstyle=5 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[4],"spec4") == 0)
        error->all(FLERR,"rstyle=5 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[4],"spec5") == 0)
        error->all(FLERR,"rstyle=5 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[4],"vac") == 0) desoutput[ndes] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      ndes++;

    } else error->all(FLERR,"Illegal event command");
  } else if (strcmp(command,"mui_init_agg") == 0) {
    if (narg != 0) error->all(FLERR,"Illegal mui_init_agg command");
    mui_init_agg();
  } else if (strcmp(command,"mui_push") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal mui_push command");
    mui_push(narg,arg);
  } else if (strcmp(command,"mui_fetch") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal mui_fetch command");
    mui_fetch(narg,arg);
  } else if (strcmp(command,"mui_push_agg") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal mui_push_agg command");
    mui_push_agg(narg,arg);
  } else if (strcmp(command,"mui_fetch_agg") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal mui_fetch_agg command");
    mui_fetch_agg(narg,arg);
  } else if (strcmp(command,"mui_fhd_lattice_size") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal mul_fhd_lattice_size command");
    mui_fhd_lattice_size_x = atof(arg[0]);
    mui_fhd_lattice_size_y = atof(arg[1]);
  } else if (strcmp(command,"mui_kmc_lattice_offset") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal mui_kmc_lattice_offset command");
    mui_kmc_lattice_offset_x = atof(arg[0]);
    mui_kmc_lattice_offset_y = atof(arg[1]);
  } else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppSurfchemtest::grow_app()
{
  type = iarray[0];
  element = iarray[1];
  ac1 = iarray[2];
  ac2 = iarray[3];
  ac3 = iarray[4];
  ac4 = iarray[5];
  ac5 = iarray[6];
  dc1 = iarray[7];
  dc2 = iarray[8];
  dc3 = iarray[9];
  dc4 = iarray[10];
  dc5 = iarray[11];
  density1 = darray[0];
  density2 = darray[1];
  density3 = darray[2];
  density4 = darray[3];
  density5 = darray[4];
  temp = darray[5];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppSurfchemtest::init_app()
{
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    memory->create(firstevent,nlocal,"app:firstevent");

    // esites must be large enough for 3 sites and their 1st neighbors

    esites = new int[3 + 3*maxneigh];

    // initializing ac1-ac5 and dc1-dc5

    for (int i = 0; i < nlocal; i++) {
      ac1[i] = 0;
      ac2[i] = 0;
      ac3[i] = 0;
      ac4[i] = 0;
      ac5[i] = 0;
      dc1[i] = 0;
      dc2[i] = 0;
      dc3[i] = 0;
      dc4[i] = 0;
      dc5[i] = 0;
    }

    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: ac1-ac5, dc1-dc5 initialized to zero\n");
      fflush(screen);
    }
  }

  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] < SITEA || type[i] > SITEC) flag = 1;
    if (element[i] < VACANCY || element[i] > SPEC5) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ---------------------------------------------------------------------- */

void AppSurfchemtest::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;

  // clear event list

  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;

  // set propensities from rates
  // propentities for adsorption reactions will be calculated in site_propensity(i)

  if (temperature == 0.0)
    error->all(FLERR,"Temperature cannot be 0.0 for app surfchemtest");

  for (int m = 0; m < none; m++) {
    spropensity[m] = srate[m];
    scount[m] = 0;
  }
  for (int m = 0; m < ntwo; m++) {
    dpropensity[m] = drate[m];
    dcount[m] = 0;
  }
  for (int m = 0; m < nthree; m++) {
    tpropensity[m] = trate[m];
    tcount[m] = 0;
  }
  for (int m = 0; m < nads; m++) {
    adscount[m] = 0;
  }
  for (int m = 0; m < ndes; m++) {
    despropensity[m] = desrate[m];
    descount[m] = 0;
  }
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppSurfchemtest::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppSurfchemtest::site_propensity(int i)
{
  int j,k,m;

  // valid single, double, triple events are in tabulated lists
  // propensity for each event is input by user

  clear_events(i);

  double proball = 0.0;

  // single-site events

  for (m = 0; m < none; m++) {
    if (type[i] != stype[m] || element[i] != sinput[m]) continue;
    add_event(i,1,m,spropensity[m],-1,-1);
    proball += spropensity[m];
  }

  // double-site events

  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (m = 0; m < ntwo; m++) {
      if (type[i] != dtype[m][0] || element[i] != dinput[m][0]) continue;
      if (type[j] != dtype[m][1] || element[j] != dinput[m][1]) continue;
      add_event(i,2,m,dpropensity[m],j,-1);
      proball += dpropensity[m];
    }
  }

  // triple-site events

  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (int kk = 0; kk < numneigh[i]; kk++) {
      if (jj == kk) continue;
      k = neighbor[i][kk];
      for (m = 0; m < nthree; m++) {
    if (type[i] != ttype[m][0] || element[i] != tinput[m][0]) continue;
    if (type[j] != ttype[m][1] || element[j] != tinput[m][1]) continue;
    if (type[k] != ttype[m][2] || element[k] != tinput[m][2]) continue;
    add_event(i,3,m,tpropensity[m],j,k);
    proball += tpropensity[m];
      }
    }
  }

  // adsorption events

  double adspropensity;

  for (m = 0; m < nads; m++) {
    if (type[i] != adstype[m] || element[i] != adsinput[m]) continue;

/*
    // propensity for adsorption = adsrate * num_dens * sqrt(site_temp/sys_temp)
    if (adsoutput[m]==SPEC1) adspropensity = adsrate[m]*density1[i]*sqrt(temp[i]/temperature);
    else if (adsoutput[m]==SPEC2) adspropensity = adsrate[m]*density2[i]*sqrt(temp[i]/temperature);
    else if (adsoutput[m]==SPEC3) adspropensity = adsrate[m]*density3[i]*sqrt(temp[i]/temperature);
    else if (adsoutput[m]==SPEC4) adspropensity = adsrate[m]*density4[i]*sqrt(temp[i]/temperature);
    else if (adsoutput[m]==SPEC5) adspropensity = adsrate[m]*density5[i]*sqrt(temp[i]/temperature);
    add_event(i,4,m,adspropensity,-1,-1);
    proball += adspropensity;
*/

    // propensity for adsorption = adsrate * num_dens
    if (adsoutput[m]==SPEC1) adspropensity = adsrate[m]*density1[i];
    else if (adsoutput[m]==SPEC2) adspropensity = adsrate[m]*density2[i];
    else if (adsoutput[m]==SPEC3) adspropensity = adsrate[m]*density3[i];
    else if (adsoutput[m]==SPEC4) adspropensity = adsrate[m]*density4[i];
    else if (adsoutput[m]==SPEC5) adspropensity = adsrate[m]*density5[i];
    add_event(i,4,m,adspropensity,-1,-1);
    proball += adspropensity;

/*
    // simpler rate
    add_event(i,4,m,adsrate[m],-1,-1);
    proball += adsrate[m];
*/
  }

  // desorption events

  for (m = 0; m < ndes; m++) {
    if (type[i] != destype[m] || element[i] != desinput[m]) continue;
    add_event(i,5,m,despropensity[m],-1,-1);
    proball += despropensity[m];
  }

  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppSurfchemtest::site_event(int i, class RandomPark *random)
{
  int j,k,m,n;

  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event

  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  // perform single, double, triple, or adsoprtion/desorption event

  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;
  k = events[ievent].kpartner;

  if (rstyle == 1) {
    element[i] = soutput[which];
    scount[which]++;
  } else if (rstyle == 2) {
    element[i] = doutput[which][0];
    element[j] = doutput[which][1];
    dcount[which]++;
  } else if (rstyle == 3) {
    element[i] = toutput[which][0];
    element[j] = toutput[which][1];
    element[k] = toutput[which][2];
    tcount[which]++;
  } else if (rstyle == 4) {  // adsorption case
    element[i] = adsoutput[which];
    adscount[which]++;
    if (element[i] == SPEC1) ac1[i]++;
    else if (element[i] == SPEC2) ac2[i]++;
    else if (element[i] == SPEC3) ac3[i]++;
    else if (element[i] == SPEC4) ac4[i]++;
    else if (element[i] == SPEC5) ac5[i]++;
  } else if (rstyle == 5) {  // desorption case
    if (element[i] == SPEC1) dc1[i]++;
    else if (element[i] == SPEC2) dc2[i]++;
    else if (element[i] == SPEC3) dc3[i]++;
    else if (element[i] == SPEC4) dc4[i]++;
    else if (element[i] == SPEC5) dc5[i]++;
    element[i] = desoutput[which];
    descount[which]++;
  }

  // compute propensity changes for participating sites and first neighs
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  int isite = i2site[i];
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  for (n = 0; n < numneigh[i]; n++) {
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
  }

  if (rstyle == 2 || rstyle == 3) {
    for (n = 0; n < numneigh[j]; n++) {
      m = neighbor[j][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
        propensity[isite] = site_propensity(m);
        esites[nsites++] = isite;
        echeck[isite] = 1;
      }
    }
  }

  if (rstyle == 3) {
    for (n = 0; n < numneigh[k]; n++) {
      m = neighbor[k][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
        propensity[isite] = site_propensity(m);
        esites[nsites++] = isite;
        echeck[isite] = 1;
      }
    }
  }

  solve->update(nsites,esites,propensity);

  // clear echeck array

  for (m = 0; m < nsites; m++) echeck[esites[m]] = 0;
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppSurfchemtest::clear_events(int i)
{
  int next;
  int index = firstevent[i];
  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
   add an event to list for site I
   event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppSurfchemtest::add_event(int i, int rstyle, int which, double propensity,
              int jpartner, int kpartner)
{
  // grow event list and setup free list

  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events =
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;

  events[freeevent].style = rstyle;
  events[freeevent].which = which;
  events[freeevent].jpartner = jpartner;
  events[freeevent].kpartner = kpartner;
  events[freeevent].propensity = propensity;

  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   grow list of stored reactions for single, double, or triple
------------------------------------------------------------------------- */

void AppSurfchemtest::grow_reactions(int rstyle)
{
  if (rstyle == 1) {
    int n = none + 1;
    memory->grow(srate,n,"app/surfchemtest:srate");
    memory->grow(spropensity,n,"app/surfchemtest:spropensity");
    memory->grow(stype,n,"app/surfchemtest:stype");
    memory->grow(sinput,n,"app/surfchemtest:sinput");
    memory->grow(soutput,n,"app/surfchemtest:soutput");
    memory->grow(scount,n,"app/surfchemtest:scount");

  } else if (rstyle == 2) {
    int n = ntwo + 1;
    memory->grow(drate,n,"app/surfchemtest:drate");
    memory->grow(dpropensity,n,"app/surfchemtest:dpropensity");
    dtype = memory->grow(dtype,n,2,"app/surfchemtest:dtype");
    dinput = memory->grow(dinput,n,2,"app/surfchemtest:dinput");
    doutput = memory->grow(doutput,n,2,"app/surfchemtest:doutput");
    memory->grow(dcount,n,"app/surfchemtest:dcount");

  } else if (rstyle == 3) {
    int n = nthree + 1;
    memory->grow(trate,n,"app/surfchemtest:trate");
    memory->grow(tpropensity,n,"app/surfchemtest:tpropensity");
    ttype = memory->grow(ttype,n,3,"app/surfchemtest:ttype");
    tinput = memory->grow(tinput,n,3,"app/surfchemtest:tinput");
    toutput = memory->grow(toutput,n,3,"app/surfchemtest:toutput");
    memory->grow(tcount,n,"app/surfchemtest:tcount");
  } else if (rstyle == 4) {
    int n = nads + 1;
    memory->grow(adsrate,n,"app/surfchemtest:adsrate");
    memory->grow(adstype,n,"app/surfchemtest:adstype");
    memory->grow(adsinput,n,"app/surfchemtest:adsinput");
    memory->grow(adsoutput,n,"app/surfchemtest:adsoutput");
    memory->grow(adscount,n,"app/surfchemtest:adscount");
  } else if (rstyle == 5) {
    int n = ndes + 1;
    memory->grow(desrate,n,"app/surfchemtest:desrate");
    memory->grow(despropensity,n,"app/surfchemtest:despropensity");
    memory->grow(destype,n,"app/surfchemtest:destype");
    memory->grow(desinput,n,"app/surfchemtest:desinput");
    memory->grow(desoutput,n,"app/surfchemtest:desoutput");
    memory->grow(descount,n,"app/surfchemtest:descount");
  }
}

/* ----------------------------------------------------------------------
   MUI routines
------------------------------------------------------------------------- */

void AppSurfchemtest::mui_init_agg()
{
  assert(nlocal>0);
  assert(mui_fhd_lattice_size_x>0);
  assert(mui_fhd_lattice_size_y>0);
/*
    // announce span

    double tmp[2];

    tmp[0] = domain->subxlo;
    tmp[1] = domain->subylo;
    point<double,2> send_span_lo(tmp);

    tmp[0] = domain->subxhi;
    tmp[1] = domain->subyhi;
    point<double,2> send_span_hi(tmp);

    mui::geometry::box<config_2d> send_span(send_span_lo,send_span_hi);

    spk->uniface->announce_send_span(0.,1.e10,send_span);

    tmp[0] = domain->subxlo - 0.5*mui_fhd_lattice_size_x;
    tmp[1] = domain->subylo - 0.5*mui_fhd_lattice_size_y;
    point<double,2> recv_span_lo(tmp);

    tmp[0] = domain->subxhi + 0.5*mui_fhd_lattice_size_x;
    tmp[1] = domain->subyhi + 0.5*mui_fhd_lattice_size_y;
    point<double,2> recv_span_hi(tmp);

    mui::geometry::box<config_2d> recv_span(recv_span_lo,recv_span_hi);

    spk->uniface->announce_recv_span(0.,1.e10,recv_span);
*/
    // create maps between KMC sites and FHD cells

    // 1. collect info

    int nFHDcellx = rint((domain->boxxhi-domain->boxxlo)/mui_fhd_lattice_size_x);
    int nFHDcelly = rint((domain->boxyhi-domain->boxylo)/mui_fhd_lattice_size_y);

    vector<vector<int>> cntKMCsite(nFHDcellx,vector<int>(nFHDcelly,0));
    vector<vector<double>> sum1(nFHDcellx,vector<double>(nFHDcelly,0.));
    vector<vector<double>> sum2(nFHDcellx,vector<double>(nFHDcelly,0.));

    for (int i = 0; i < nlocal; i++) {
      int nx = floor((xyz[i][0]+mui_kmc_lattice_offset_x)/mui_fhd_lattice_size_x);
      int ny = floor((xyz[i][1]+mui_kmc_lattice_offset_y)/mui_fhd_lattice_size_y);
      assert(nx>=0 && nx<nFHDcellx && ny>=0 && ny<nFHDcelly);
      cntKMCsite[nx][ny]++;
      sum1[nx][ny] += xyz[i][0]+mui_kmc_lattice_offset_x;
      sum2[nx][ny] += xyz[i][1]+mui_kmc_lattice_offset_y;
    }

    int cntFHDcell = 0;
    for (int nx = 0; nx < nFHDcellx; nx++)
      for (int ny = 0; ny < nFHDcelly; ny++)
        if (cntKMCsite[nx][ny] > 0)
          cntFHDcell++;
    nlocalFHDcell = cntFHDcell;

    // 2. xFHD/yFHD/MUIintval/MUIdblval/localFHDcell

    xFHD = new double[nlocalFHDcell];
    yFHD = new double[nlocalFHDcell];
    MUIintval = new int[nlocalFHDcell];
    MUIdblval = new double[nlocalFHDcell];
    localFHDcell = new int [nlocal];

    vector<vector<int>> FHDcell(nFHDcellx,vector<int>(nFHDcelly,-1));

    cntFHDcell = 0;
    for (int nx = 0; nx < nFHDcellx; nx++) {
      for (int ny = 0; ny < nFHDcelly; ny++) {
        if (cntKMCsite[nx][ny] > 0) {
          FHDcell[nx][ny] = cntFHDcell;
          xFHD[cntFHDcell] = sum1[nx][ny]/cntKMCsite[nx][ny];
          yFHD[cntFHDcell] = sum2[nx][ny]/cntKMCsite[nx][ny];
          cntFHDcell++;
        }
      }
    }
    assert(cntFHDcell==nlocalFHDcell);

    for (int i = 0; i < nlocal; i++) {
      int nx = floor((xyz[i][0]+mui_kmc_lattice_offset_x)/mui_fhd_lattice_size_x);
      int ny = floor((xyz[i][1]+mui_kmc_lattice_offset_y)/mui_fhd_lattice_size_y);
      assert(nx>=0 && nx<nFHDcellx && ny>=0 && ny<nFHDcelly);
      localFHDcell[i] = FHDcell[nx][ny];
    }
}

void AppSurfchemtest::mui_push(int narg, char **arg)
{
  int timestamp = atoi(arg[0]);

  if (domain->me == 0 && screen) {
    fprintf(screen,"** DEBUG: MUI push at timestamp %d\n",timestamp);
    fflush(screen);
  }

  for (int k=1;k<narg;k++)
  {
    if (strcmp(arg[k],"type") == 0) {           // i1
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_type",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},type[i]);
      }
    } else if (strcmp(arg[k],"element") == 0) { // i2
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_element",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},element[i]);
      }
    }
    else if (strcmp(arg[k],"ac1") == 0) {       // i3
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_ac1",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},ac1[i]);
        ac1[i] = 0;
      }
    } else if (strcmp(arg[k],"ac2") == 0) {     // i4
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_ac2",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},ac2[i]);
        ac2[i] = 0;
      }
    } else if (strcmp(arg[k],"ac3") == 0) {     // i5
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_ac3",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},ac3[i]);
        ac3[i] = 0;
      }
    } else if (strcmp(arg[k],"ac4") == 0) {     // i6
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_ac4",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},ac4[i]);
        ac4[i] = 0;
      }
    } else if (strcmp(arg[k],"ac5") == 0) {     // i7
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_ac5",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},ac5[i]);
        ac5[i] = 0;
      }
    } else if (strcmp(arg[k],"dc1") == 0) {     // i8
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_dc1",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},dc1[i]);
        dc1[i] = 0;
      }
    } else if (strcmp(arg[k],"dc2") == 0) {     // i9
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_dc2",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},dc2[i]);
        dc2[i] = 0;
      }
    } else if (strcmp(arg[k],"dc3") == 0) {     // i10
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_dc3",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},dc3[i]);
        dc3[i] = 0;
      }
    } else if (strcmp(arg[k],"dc4") == 0) {     // i11
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_dc4",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},dc4[i]);
        dc4[i] = 0;
      }
    } else if (strcmp(arg[k],"dc5") == 0) {     // i12
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_dc5",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},dc5[i]);
        dc5[i] = 0;
      }
    } else if (strcmp(arg[k],"density1") == 0) {    // d1
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_density1",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},density1[i]);
      }
    } else if (strcmp(arg[k],"density2") == 0) {    // d2
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_density2",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},density2[i]);
      }
    } else if (strcmp(arg[k],"density3") == 0) {    // d3
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_density3",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},density3[i]);
      }
    } else if (strcmp(arg[k],"density4") == 0) {    // d4
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_density4",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},density4[i]);
      }
    } else if (strcmp(arg[k],"density5") == 0) {    // d5
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_density5",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},density5[i]);
      }
    } else if (strcmp(arg[k],"temp") == 0) {        // d6
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_temp",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},temp[i]);
      }
    } else {
      error->all(FLERR,"Illegal mui_push command");
    }

    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: %s pushed\n",arg[k]);
      fflush(screen);
    }
  }

  spk->uniface->commit(timestamp);

  return;
}

void AppSurfchemtest::mui_push_agg(int narg, char **arg)
{
  int timestamp = atoi(arg[0]);

  if (domain->me == 0 && screen) {
    fprintf(screen,"** DEBUG: mui_push_agg at timestamp %d\n",timestamp);
    fflush(screen);
  }

  for (int k=1;k<narg;k++)
  {
    if (strcmp(arg[k],"ac1") == 0) {
      // compute the sum over each FHD domain
      for (int n=0;n<nlocalFHDcell;n++) MUIintval[n] = 0;
      for (int i=0;i<nlocal;i++) {
        MUIintval[localFHDcell[i]] += ac1[i];
        ac1[i] = 0;
      }
      // push for each FHD domain
      for (int n=0;n<nlocalFHDcell;n++)
        spk->uniface->push("CH_ac1",{xFHD[n],yFHD[n]},MUIintval[n]);
    } else if (strcmp(arg[k],"dc1") == 0) {
      // compute the sum over each FHD domain
      for (int n=0;n<nlocalFHDcell;n++) MUIintval[n] = 0;
      for (int i=0;i<nlocal;i++) {
        MUIintval[localFHDcell[i]] += dc1[i];
        dc1[i] = 0;
      }
      // push for each FHD domain
      for (int n=0;n<nlocalFHDcell;n++)
        spk->uniface->push("CH_dc1",{xFHD[n],yFHD[n]},MUIintval[n]);
    } else {
      error->all(FLERR,"Illegal mui_push_agg command");
    }

    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: %s pushed\n",arg[k]);
      fflush(screen);
    }
  }

  spk->uniface->commit(timestamp);

  return;
}

void AppSurfchemtest::mui_fetch(int narg, char **arg)
{
  int timestamp = atoi(arg[0]);

  if (domain->me == 0 && screen) {
    fprintf(screen,"** DEBUG: MUI fetch at timestamp %d\n",timestamp);
    fflush(screen);
  }

  if (mui_fhd_lattice_size_x <= 0. || mui_fhd_lattice_size_y <= 0.)
    error->all(FLERR,"mui_fhd_lattice_size must be set as two positive numbers");

  mui::sampler_kmc_fhd2d<double> s({mui_fhd_lattice_size_x,mui_fhd_lattice_size_y});
  mui::chrono_sampler_exact2d t;

  for (int k=1;k<narg;k++)
  {
    if (strcmp(arg[k],"density1") == 0) {           // d1
      for (int i=0;i<nlocal;i++) {
        density1[i] = spk->uniface->fetch("CH_density1",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},timestamp,s,t);
      }
    } else if (strcmp(arg[k],"density2") == 0) {    // d2
      for (int i=0;i<nlocal;i++) {
        density2[i] = spk->uniface->fetch("CH_density2",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},timestamp,s,t);
      }
    } else if (strcmp(arg[k],"density3") == 0) {    // d3
      for (int i=0;i<nlocal;i++) {
        density3[i] = spk->uniface->fetch("CH_density3",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},timestamp,s,t);
      }
    } else if (strcmp(arg[k],"density4") == 0) {    // d4
      for (int i=0;i<nlocal;i++) {
        density4[i] = spk->uniface->fetch("CH_density4",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},timestamp,s,t);
      }
    } else if (strcmp(arg[k],"density5") == 0) {    // d5
      for (int i=0;i<nlocal;i++) {
        density5[i] = spk->uniface->fetch("CH_density5",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},timestamp,s,t);
      }
    } else if (strcmp(arg[k],"temp") == 0) {        // d6
      for (int i=0;i<nlocal;i++) {
        temp[i] = spk->uniface->fetch("CH_temp",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},timestamp,s,t);
      }
    } else {
      error->all(FLERR,"Illegal mui_fetch command");
    }

    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: %s fetched\n",arg[k]);
      fflush(screen);
    }
  }

  spk->uniface->forget(timestamp);

  return;
}

void AppSurfchemtest::mui_fetch_agg(int narg, char **arg)
{
  int timestamp = atoi(arg[0]);

  if (domain->me == 0 && screen) {
    fprintf(screen,"** DEBUG: mui_fetch_agg at timestamp %d\n",timestamp);
    fflush(screen);
  }

  if (mui_fhd_lattice_size_x <= 0. || mui_fhd_lattice_size_y <= 0.)
    error->all(FLERR,"mui_fhd_lattice_size must be set as two positive numbers");

  mui::sampler_kmc_fhd2d<double> s({mui_fhd_lattice_size_x,mui_fhd_lattice_size_y});
  mui::chrono_sampler_exact2d t;

  for (int k=1;k<narg;k++)
  {
    if (strcmp(arg[k],"density1") == 0) {
      // get info for each FHD cell
      for (int n=0;n<nlocalFHDcell;n++)
        MUIdblval[n] = spk->uniface->fetch("CH_density1",{xFHD[n],yFHD[n]},timestamp,s,t);
      // distribute info to each KMC site
      for (int i=0;i<nlocal;i++)
        density1[i] = MUIdblval[localFHDcell[i]];
    } else {
      error->all(FLERR,"Illegal mui_fetch_agg command");
    }

    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: %s fetched\n",arg[k]);
      fflush(screen);
    }
  }

  spk->uniface->forget(timestamp);

  return;
}
