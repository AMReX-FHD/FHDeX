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
#include "app_lotkavolterra.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "domain.h"

#ifdef MUI
#include "mui.h"
#endif

#if defined(USE_AMREX_MPMD)
#include <AMReX_MPMD.H>
#endif

#include <vector>

#ifdef MUI
using namespace mui;
#endif

using namespace SPPARKS_NS;

using namespace std;

enum{NOOP,SITEA};          // Identical surface site
enum{VACANCY,SPEC1,SPEC2}; // VACANCY: O, SPEC1: A, SPEC2: B

#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppLotkavolterra::AppLotkavolterra(SPPARKS *spk, int narg, char **arg) :
  AppLattice(spk,narg,arg)
{
  ninteger = 6;   // type: site type
                  // element: site element
                  // ac1, ac2 : adsorption count
                  // dc1, dc2 : desorption count
  ndouble = 3;    // pressure1/pressure2: partial of the contacting FHD cell
                  // temp: temperature of the contacting FHD cell
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

  nprey = npredation = npredator = 0;
  preyrate = predationrate = predatorrate = NULL;
  prey_is_rate = NULL;

  preytype = preyinput = preyoutput = NULL;
  predationtype = predationinput = predationoutput = NULL;
  predatortype = predatorinput = predatoroutput = NULL;
  preycount = predationcount = predatorcount = NULL;

  prey = predator = NULL;

#ifdef MUI
  mui_fhd_lattice_size_x = mui_fhd_lattice_size_y = -1.;
  mui_kmc_lattice_offset_x = mui_kmc_lattice_offset_y = 0.;

  xFHD = NULL;
  yFHD = NULL;
  MUIintval = NULL;
  MUIdblval = NULL;
  localFHDcell = NULL;
  nlocalFHDcell_world = NULL;
#endif
}

/* ---------------------------------------------------------------------- */

AppLotkavolterra::~AppLotkavolterra()
{
  delete [] esites;
  delete [] echeck;

  memory->sfree(events);
  memory->destroy(firstevent);
  memory->destroy(preytype);
  memory->destroy(preyinput);
  memory->destroy(preyrate);
  memory->destroy(preyoutput);

  memory->destroy(prey_is_rate);

  memory->destroy(predationtype);
  memory->destroy(predationinput);
  memory->destroy(predationrate);
  memory->destroy(predationoutput);

  memory->destroy(predatortype);
  memory->destroy(predatorinput);
  memory->destroy(predatorrate);
  memory->destroy(predatoroutput);

  memory->destroy(preycount);
  memory->destroy(predationcount);
  memory->destroy(predatorcount);

  memory->destroy(prey);
  memory->destroy(predator);

#ifdef MUI
  delete [] xFHD;
  delete [] yFHD;
  delete [] MUIintval;
  delete [] MUIdblval;
  delete [] localFHDcell;
  if (domain->me == 0) delete [] nlocalFHDcell_world;
#endif
}

/* ---------------------------------------------------------------------- */

void AppLotkavolterra::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"event") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal event command - event style");
    int rstyle = atoi(arg[0]);
    grow_reactions(rstyle);

    if (rstyle == 1) { // event for prey - O + A -> A + A
      if (narg != 8 && narg != 10) error->all(FLERR,"Illegal event command - event keyword");

      if (strcmp(arg[1],"siteA") == 0) preytype[nprey][0] = SITEA;
      else error->all(FLERR,"Illegal event command - site");
      if (strcmp(arg[2],"siteA") == 0) preytype[nprey][1] = SITEA;
      else error->all(FLERR,"Illegal event command - site");

      if (strcmp(arg[3],"vac") == 0) preyinput[nprey][0] = VACANCY;
      else error->all(FLERR,"Illegal event command - vacancy for prey");
      if (strcmp(arg[4],"vac") == 0) preyinput[nprey][1] = VACANCY;
      else if (strcmp(arg[4],"spec1") == 0) preyinput[nprey][1] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) preyinput[nprey][1] = SPEC2;
      else error->all(FLERR,"Illegal event command - input");

      preyrate[nprey] = atof(arg[5]);

      if (strcmp(arg[6],"vac") == 0) preyoutput[nprey][0] = VACANCY;
      else if (strcmp(arg[6],"spec1") == 0) preyoutput[nprey][0] = SPEC1;
      else if (strcmp(arg[6],"spec2") == 0) preyoutput[nprey][0] = SPEC2;
      else error->all(FLERR,"Illegal event command - output");
      if (strcmp(arg[7],"vac") == 0) preyoutput[nprey][1] = VACANCY;
      else if (strcmp(arg[7],"spec1") == 0) preyoutput[nprey][1] = SPEC1;
      else if (strcmp(arg[7],"spec2") == 0) preyoutput[nprey][1] = SPEC2;
      else error->all(FLERR,"Illegal event command - output");

      prey_is_rate[nprey] = true;

      if (narg == 10) {
        if (strcmp(arg[8],"FHD") == 0) {
          prey_is_rate[nprey] = false;
          if (strcmp(arg[9],"spec1") == 0) prey[nprey] = SPEC1;
          else if (strcmp(arg[9],"spec2") == 0) prey[nprey] = SPEC2;
          else error->all(FLERR, "Illegal event command - species");
        }
      }

      nprey++;

    } else if (rstyle == 2) { // event for predation - A + B -> B + B
      if (narg != 8) error->all(FLERR,"Illegal event command - event keyword");

      if (strcmp(arg[1],"siteA") == 0) predationtype[npredation][0] = SITEA;
      else error->all(FLERR,"Illegal event command - site");
      if (strcmp(arg[2],"siteA") == 0) predationtype[npredation][1] = SITEA;
      else error->all(FLERR,"Illegal event command - site");

      if (strcmp(arg[3],"spec1") == 0) predationinput[npredation][0] = SPEC1;
      else error->all(FLERR,"Illegal event command - spec1 for predation");
      if (strcmp(arg[4],"vac") == 0) predationinput[npredation][1] = VACANCY;
      else if (strcmp(arg[4],"spec1") == 0) predationinput[npredation][1] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) predationinput[npredation][1] = SPEC2;
      else error->all(FLERR,"Illegal event command - input");

      predationrate[npredation] = atof(arg[5]);

      if (strcmp(arg[6],"vac") == 0) predationoutput[npredation][0] = VACANCY;
      else if (strcmp(arg[6],"spec1") == 0) predationoutput[npredation][0] = SPEC1;
      else if (strcmp(arg[6],"spec2") == 0) predationoutput[npredation][0] = SPEC2;
      else error->all(FLERR,"Illegal event command - output");
      if (strcmp(arg[7],"vac") == 0) predationoutput[npredation][1] = VACANCY;
      else if (strcmp(arg[7],"spec1") == 0) predationoutput[npredation][1] = SPEC1;
      else if (strcmp(arg[7],"spec2") == 0) predationoutput[npredation][1] = SPEC2;
      else error->all(FLERR,"Illegal event command - output");

      npredation++;

    } else if (rstyle == 3) { // event for predator - B -> O
      if (narg != 5 && narg != 7) error->all(FLERR,"Illegal event command - event keyword");

      if (strcmp(arg[1],"siteA") == 0) predatortype[npredator] = SITEA;
      else error->all(FLERR,"Illegal event command - site");
      if (strcmp(arg[2],"spec2") == 0) predatorinput[npredator] = SPEC2;
      else error->all(FLERR,"Illegal event command - spec2 for predator");

      predatorrate[npredator] = atof(arg[3]);

      if (strcmp(arg[4],"vac") == 0) predatoroutput[npredator] = VACANCY;
      else if (strcmp(arg[4],"spec1") == 0) predatoroutput[npredator] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) predatoroutput[npredator] = SPEC2;
      else error->all(FLERR,"Illegal event command - output");

      if (narg == 7) {
        if (strcmp(arg[5],"FHD") == 0) {
          if (strcmp(arg[6],"spec1") == 0) predator[npredator] = SPEC1;
          else if (strcmp(arg[6],"spec2") == 0) predator[npredator] = SPEC2;
          else error->all(FLERR, "Illegal event command - species");
        }
      }

      npredator++;

    }
  }

#if defined(MUI)
    else if (strcmp(command,"mui_init_agg") == 0) {
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
  } else if (strcmp(command,"mui_commit") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal mui_commit command");
    mui_commit(narg,arg);
  } else if (strcmp(command,"mui_forget") == 0) {
    if (narg < 1) error->all(FLERR,"Illegal mui_forget command");
    mui_forget(narg,arg);
  } else if (strcmp(command,"mui_fhd_lattice_size") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal mul_fhd_lattice_size command");
    mui_fhd_lattice_size_x = atof(arg[0]);
    mui_fhd_lattice_size_y = atof(arg[1]);
  } else if (strcmp(command,"mui_kmc_lattice_offset") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal mui_kmc_lattice_offset command");
    mui_kmc_lattice_offset_x = atof(arg[0]);
    mui_kmc_lattice_offset_y = atof(arg[1]);
  }
#elif defined(USE_AMREX_MPMD)
  else if (strcmp(command,"amrex_init_agg") == 0) {
    if (narg != 0) error->all(FLERR,"Illegal amrex_init_agg command");
    amrex_init_agg();
  } else if (strcmp(command,"amrex_push_agg") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal amrex_push_agg command");
    amrex_push_agg(narg,arg);
  } else if (strcmp(command,"amrex_fetch_agg") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal amrex_fetch_agg command");
    amrex_fetch_agg(narg,arg);
  } else if (strcmp(command,"amrex_fhd_lattice_size") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal amrex_fhd_lattice_size command");
    amrex_fhd_lattice_size_x = atof(arg[0]);
    amrex_fhd_lattice_size_y = atof(arg[1]);
  } else if (strcmp(command,"amrex_kmc_lattice_offset") == 0) {
    if (narg != 2) error->all(FLERR,"Illegal amrex_kmc_lattice_offset command");
    amrex_kmc_lattice_offset_x = atof(arg[0]);
    amrex_kmc_lattice_offset_y = atof(arg[1]);
  }
#endif

  else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void AppLotkavolterra::grow_app()
{
  type = iarray[0];
  element = iarray[1];
  ac1 = iarray[2];
  ac2 = iarray[3];
  dc1 = iarray[4];
  dc2 = iarray[5];
  pressure1 = darray[0];
  pressure2 = darray[1];
  temp = darray[2];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppLotkavolterra::init_app()
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
      dc1[i] = 0;
      dc2[i] = 0;
    }

    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: ac1-2 and dc1-2 initialized to zero\n");
      fflush(screen);
    }
  }

  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] < SITEA || type[i] > SITEA) flag = 1;
    if (element[i] < VACANCY || element[i] > SPEC2) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ---------------------------------------------------------------------- */

void AppLotkavolterra::setup_app()
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
    error->all(FLERR,"Temperature cannot be 0.0 for app lotkavolterra");

  for (int m = 0; m < nprey; m++) {
    preycount[m] = 0;
  }
  for (int m = 0; m < npredation; m++) {
    predationcount[m] = 0;
  }
  for (int m = 0; m < npredator; m++) {
    predatorcount[m] = 0;
  }
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppLotkavolterra::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppLotkavolterra::site_propensity(int i)
{
  clear_events(i);
  double proball = 0.0;

  int j, m;
  double preypropensity;
  double tempratio = temp[i]/temperature;

  if (element[i] == VACANCY) { // prey events
    for (int jj = 0; jj < numneigh[j]; jj++) {
      j = neighbor[i][jj];
      for (m = 0; m < nprey; m++) {
        if (element[j] == preyinput[m][1]) {
          if (prey_is_rate[m]) {
            add_event(i,1,m,preyrate[m],-1,-1);
	    proball += preyrate[m];
	  }
	  else {
            if (prey[m] == SPEC1) preypropensity = preyrate[m]*pressure1[i]*pow(tempratio,-0.5);
            else if (prey[m] == SPEC2) preypropensity = preyrate[m]*pressure2[i]*pow(tempratio,-0.5);
	    add_event(i,1,m,preypropensity,-1,-1);
            proball += preypropensity;
          }
        }
      }
    }
  }

  else if (element[i] == SPEC1) { // predation events
    for (int jj = 0; jj < numneigh[j]; jj++) {
      j = neighbor[i][jj];
      for (m = 0; m < npredation; m++) {
        if (element[j] == predationinput[m][1]) {
          add_event(i,2,m,predationrate[m],-1,-1);
	  proball += predationrate[m];
        }
      }
    }
  }

  else if (element[i] == SPEC2) { // predator events
    for (m = 0; m < npredator; m++) {
      add_event(i,3,m,predatorrate[m],-1,-1);
      proball += predatorrate[m];
    }
  }

  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppLotkavolterra::site_event(int i, class RandomPark *random)
{
  int j, m, n;
  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  // perform prey, predation, predator event

  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;

  if (rstyle == 1) { // prey case
    element[i] = preyoutput[which][0];
    if (prey[which] == SPEC1) ac1[i]++;
    else if (prey[which] == SPEC2) ac2[i]++;
    preycount[which]++;
  }
  else if (rstyle == 2) { // predation case
    element[i] = predationoutput[which][0];
    predationcount[which]++;
  }
  else if (rstyle == 3) { // predator case
    element[i] = predatoroutput[which];
    if (predator[which] == SPEC1) dc1[i]++;
    else if (predator[which] == SPEC2) dc2[i]++;
    predatorcount[which]++;
  }

  // compute propensity changes for participating sites and neighbors
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  int isite = i2site[i];
  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  solve->update(nsites,esites,propensity);

  // clear echeck array

  for (m = 0; m < nsites; m++) echeck[esites[m]] = 0;
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void AppLotkavolterra::clear_events(int i)
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

void AppLotkavolterra::add_event(int i, int rstyle, int which, double propensity,
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
   grow list of stored reactions for prey, predation, predator
------------------------------------------------------------------------- */

void AppLotkavolterra::grow_reactions(int rstyle)
{
  if (rstyle == 1) { // O + A -> A + A
    int n = nprey + 1;
    memory->grow(preyrate,n,"app/lotkavolterra:preyrate");
    preytype = memory->grow(preytype,n,2,"app/lotkavolterra:preytype");
    preyinput = memory->grow(preyinput,n,2,"app/lotkavolterra:preyinput");
    preyoutput = memory->grow(preyoutput,n,2,"app/lotkavolterra:preyoutput");
    memory->grow(preycount,n,"app/lotkavolterra:preycount");
    memory->grow(prey_is_rate,n,"app/lotkavolterra:prey_is_rate");
    memory->grow(prey,n,"app/lotkavolterra:prey");

  } else if (rstyle == 2) { // A + B -> B + B
    int n = npredation + 1;
    memory->grow(predationrate,n,"app/lotkavolterra:predationrate");
    predationtype = memory->grow(predationtype,n,2,"app/lotkavolterra:predationtype");
    predationinput = memory->grow(predationinput,n,2,"app/lotkavolterra:predationinput");
    predationoutput = memory->grow(predationoutput,n,2,"app/lotkavolterra:predationoutput");
    memory->grow(predationcount,n,"app/lotkavolterra:predationcount");

  } else if (rstyle == 3) { // B -> O
    int n = npredator + 1;
    memory->grow(predatorrate,n,"app/lotkavolterra:predatorrate");
    memory->grow(predatortype,n,"app/lotkavolterra:predatortype");
    memory->grow(predatorinput,n,"app/lotkavolterra:predatorinput");
    memory->grow(predatoroutput,n,"app/lotkavolterra:predatoroutput");
    memory->grow(predatorcount,n,"app/lotkavolterra:predatorcount");
    memory->grow(predator,n,"app/lotkavolterra:predator");
  }
}

#if defined(MUI)

/* ----------------------------------------------------------------------
   MUI routines
------------------------------------------------------------------------- */

void AppLotkavolterra::mui_init_agg()
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

    // 1. nlocalFHDcell, nlocalFHDcell_world

    int nFHDcellx = rint((domain->boxxhi-domain->boxxlo)/mui_fhd_lattice_size_x);
    int nFHDcelly = rint((domain->boxyhi-domain->boxylo)/mui_fhd_lattice_size_y);

    if (domain->me == 0) {
        fprintf(logfile,"(boxx, boxy) = %e %e\n",domain->boxxhi-domain->boxxlo,domain->boxyhi-domain->boxylo);
        fprintf(logfile,"(mui_fhd_lattice_size_x, mui_fhd_lattice_size_y) = %e %e\n",mui_fhd_lattice_size_x,mui_fhd_lattice_size_y);
        fprintf(logfile,"(nFHDcellx, nFHDcelly) = %d %d\n",nFHDcellx,nFHDcelly);
    }

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
            if (cntKMCsite[nx][ny] > 0) cntFHDcell++;
    nlocalFHDcell = cntFHDcell;

    if (domain->me == 0) nlocalFHDcell_world = new int[domain->nprocs];

    MPI_Gather(&nlocalFHDcell,1,MPI_INT,nlocalFHDcell_world,1,MPI_INT,0,world);

    if (domain->me == 0) {
        for (int i=0;i<domain->nprocs;i++) {
            fprintf(logfile,"- rank %d: nlocalFHDcell = %d\n",i,nlocalFHDcell_world[i]);
        }
    }

    // 2. xFHD, yFHD

    xFHD = new double[nlocalFHDcell];
    yFHD = new double[nlocalFHDcell];

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

    if (domain->me == 0) {
        // output my info first
        fprintf(logfile,"** rank %d **\n",domain->me);
        for (int j=0;j<nlocalFHDcell;j++) {
            fprintf(logfile,"- (xFHD[%d], yFHD[%d]) = %e %e\n",j,j,xFHD[j],yFHD[j]);
        }
        // output for other procs
        for (int i=1;i<domain->nprocs;i++) {
            double * data1 = new double[nlocalFHDcell_world[i]];
            double * data2 = new double[nlocalFHDcell_world[i]];

            MPI_Recv(data1,nlocalFHDcell_world[i],MPI_DOUBLE,i,0,world,MPI_STATUS_IGNORE);
            MPI_Recv(data2,nlocalFHDcell_world[i],MPI_DOUBLE,i,0,world,MPI_STATUS_IGNORE);

            fprintf(logfile,"** rank %d **\n",i);
            for (int j=0;j<nlocalFHDcell_world[i];j++) {
                fprintf(logfile,"- (xFHD[%d], yFHD[%d]) = %e %e\n",j,j,data1[j],data2[j]);
            }

            delete data1;
            delete data2;
        }
    }
    else {
        MPI_Send(xFHD,nlocalFHDcell,MPI_DOUBLE,0,0,world);
        MPI_Send(yFHD,nlocalFHDcell,MPI_DOUBLE,0,0,world);
    }

    // 3. localFHDcell, MUIintval, MUIdblval

    localFHDcell = new int [nlocal];

    for (int i = 0; i < nlocal; i++) {
        int nx = floor((xyz[i][0]+mui_kmc_lattice_offset_x)/mui_fhd_lattice_size_x);
        int ny = floor((xyz[i][1]+mui_kmc_lattice_offset_y)/mui_fhd_lattice_size_y);
        assert(nx>=0 && nx<nFHDcellx && ny>=0 && ny<nFHDcelly);
        localFHDcell[i] = FHDcell[nx][ny];
    }

    vector<int> cntKMCsite2(nlocalFHDcell,0);
    for (int i = 0; i < nlocal; i++) cntKMCsite2[localFHDcell[i]]++;

    if (domain->me == 0) {
        // output my info first
        fprintf(logfile,"** rank %d **\n",domain->me);
        for (int j=0;j<nlocalFHDcell;j++) {
            fprintf(logfile,"- FHDcell %d has %d KMC sites\n",j,cntKMCsite2[j]);
        }
        // output for other procs
        for (int i=1;i<domain->nprocs;i++) {
            int * data = new int[nlocalFHDcell_world[i]];

            MPI_Recv(data,nlocalFHDcell_world[i],MPI_INT,i,0,world,MPI_STATUS_IGNORE);

            fprintf(logfile,"** rank %d **\n",i);
            for (int j=0;j<nlocalFHDcell_world[i];j++) {
                fprintf(logfile,"- FHDcell %d has %d KMC sites\n",j,data[j]);
            }

            delete data;
        }
    }
    else {
        int * data = &cntKMCsite2[0];
        MPI_Send(data,nlocalFHDcell,MPI_INT,0,0,world);
    }

    MUIintval = new int[nlocalFHDcell];
    MUIdblval = new double[nlocalFHDcell];

    if (domain->me == 0) fflush(logfile);
}

void AppLotkavolterra::mui_print_MUIdblval(int step,const char *str1,const char *str2)
{
    if (domain->me == 0)
    {
        // first print my vals
        fprintf(logfile,"** MUIdblval for %s at step %d **\n",str1,step);
        for (int j=0;j<nlocalFHDcell;j++)
            fprintf(logfile,"%s %d %d %e\n",str2,domain->me,j,MUIdblval[j]);
        // other procs
        for (int i=1;i<domain->nprocs;i++)
        {
            double * data = new double[nlocalFHDcell_world[i]];
            MPI_Recv(data,nlocalFHDcell_world[i],MPI_DOUBLE,i,0,world,MPI_STATUS_IGNORE);
            for (int j=0;j<nlocalFHDcell;j++)
                fprintf(logfile,"%s %d %d %e\n",str2,i,j,data[j]);
            delete [] data;
        }
        fflush(logfile);
    }
    else
    {
        MPI_Send(MUIdblval,nlocalFHDcell,MPI_DOUBLE,0,0,world);
    }
}

void AppLotkavolterra::mui_print_MUIintval(int step,const char *str1,const char *str2)
{
    if (domain->me == 0)
    {
        // first print my vals
        fprintf(logfile,"** MUIintval for %s at step %d **\n",str1,step);
        for (int j=0;j<nlocalFHDcell;j++)
            fprintf(logfile,"%s %d %d %d\n",str2,domain->me,j,MUIintval[j]);
        // other procs
        for (int i=1;i<domain->nprocs;i++)
        {
            int * data = new int[nlocalFHDcell_world[i]];
            MPI_Recv(data,nlocalFHDcell_world[i],MPI_INT,i,0,world,MPI_STATUS_IGNORE);
            for (int j=0;j<nlocalFHDcell;j++)
                fprintf(logfile,"%s %d %d %d\n",str2,i,j,data[j]);
            delete [] data;
        }
        fflush(logfile);
    }
    else
    {
        MPI_Send(MUIintval,nlocalFHDcell,MPI_INT,0,0,world);
    }
}

void AppLotkavolterra::mui_push(int narg, char **arg)
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
    } else if (strcmp(arg[k],"ac1") == 0) {     // i3
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_ac1",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},ac1[i]);
        ac1[i] = 0;
      }
    } else if (strcmp(arg[k],"ac2") == 0) {     // i4
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_ac2",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},ac2[i]);
        ac2[i] = 0;
      }
    } else if (strcmp(arg[k],"dc1") == 0) {     // i5
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_dc1",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},dc1[i]);
        dc1[i] = 0;
      }
    } else if (strcmp(arg[k],"dc2") == 0) {     // i6
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_dc2",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},dc2[i]);
        dc2[i] = 0;
      }
    }

    else if (strcmp(arg[k],"pressure1") == 0) {      // d1
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_pressure1",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},pressure1[i]);
      }
    } else if (strcmp(arg[k],"pressure2") == 0) {    // d2
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_pressure2",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},pressure2[i]);
      }
    } else if (strcmp(arg[k],"temp") == 0) {        // d3
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_temp",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},temp[i]);
      }
    } else if (strcmp(arg[k],"occ1") == 0) {
      for (int i=0;i<nlocal;i++) {
        int is_occ = (element[i]==1) ? 1 : 0;
        spk->uniface->push("CH_occ1",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},is_occ);
      }
    } else if (strcmp(arg[k],"occ2") == 0) {
      for (int i=0;i<nlocal;i++) {
        int is_occ = (element[i]==2) ? 1 : 0;
        spk->uniface->push("CH_occ2",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},is_occ);
      }
    } else {
      error->all(FLERR,"Illegal mui_push command");
    }
    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: %s pushed\n",arg[k]);
      fflush(screen);
    }
  }
  return;
}

void AppLotkavolterra::mui_push_agg(int narg, char **arg)
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
    } else if (strcmp(arg[k],"ac2") == 0) {
      // compute the sum over each FHD domain
      for (int n=0;n<nlocalFHDcell;n++) MUIintval[n] = 0;
      for (int i=0;i<nlocal;i++) {
        MUIintval[localFHDcell[i]] += ac2[i];
        ac2[i] = 0;
      }
      // push for each FHD domain
      for (int n=0;n<nlocalFHDcell;n++)
        spk->uniface->push("CH_ac2",{xFHD[n],yFHD[n]},MUIintval[n]);
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
    } else if (strcmp(arg[k],"dc2") == 0) {
      // compute the sum over each FHD domain
      for (int n=0;n<nlocalFHDcell;n++) MUIintval[n] = 0;
      for (int i=0;i<nlocal;i++) {
        MUIintval[localFHDcell[i]] += dc2[i];
        dc2[i] = 0;
      }
      // push for each FHD domain
      for (int n=0;n<nlocalFHDcell;n++)
        spk->uniface->push("CH_dc2",{xFHD[n],yFHD[n]},MUIintval[n]);
    } else if (strcmp(arg[k],"occ1") == 0) {
      // compute the sum over each FHD domain
      for (int n=0;n<nlocalFHDcell;n++) MUIintval[n] = 0;
      for (int i=0;i<nlocal;i++) {
        int is_occ = (element[i]==1) ? 1 : 0;
        MUIintval[localFHDcell[i]] += is_occ;
      }
      // push for each FHD domain
      for (int n=0;n<nlocalFHDcell;n++)
        spk->uniface->push("CH_occ1",{xFHD[n],yFHD[n]},MUIintval[n]);
    } else if (strcmp(arg[k],"occ2") == 0) {
      // compute the sum over each FHD domain
      for (int n=0;n<nlocalFHDcell;n++) MUIintval[n] = 0;
      for (int i=0;i<nlocal;i++) {
        int is_occ = (element[i]==2) ? 1 : 0;
        MUIintval[localFHDcell[i]] += is_occ;
      }
      // push for each FHD domain
      for (int n=0;n<nlocalFHDcell;n++)
        spk->uniface->push("CH_occ2",{xFHD[n],yFHD[n]},MUIintval[n]);
    } else if (strcmp(arg[k],"one") == 0) {
      // compute the sum over each FHD domain
      for (int n=0;n<nlocalFHDcell;n++) MUIintval[n] = 0;
      for (int i=0;i<nlocal;i++) MUIintval[localFHDcell[i]]++;
      // push for each FHD domain
      for (int n=0;n<nlocalFHDcell;n++)
        spk->uniface->push("CH_one",{xFHD[n],yFHD[n]},MUIintval[n]);
    } else {
      error->all(FLERR,"Illegal mui_push_agg command");
    }

    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: %s pushed\n",arg[k]);
      fflush(screen);
    }
  }

  return;
}

void AppLotkavolterra::mui_commit(int narg, char **arg)
{
  int timestamp = atoi(arg[0]);

  spk->uniface->commit(timestamp);

  if (domain->me == 0 && screen) {
    fprintf(screen,"** DEBUG: mui commit at timestamp %d\n",timestamp);
    fflush(screen);
  }

  return;
}

void AppLotkavolterra::mui_fetch(int narg, char **arg)
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
    if (strcmp(arg[k],"pressure1") == 0) {           // d1
      for (int i=0;i<nlocal;i++) {
        pressure1[i] = spk->uniface->fetch("CH_pressure1",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},timestamp,s,t);
      }
    } else if (strcmp(arg[k],"pressure2") == 0) {    // d2
      for (int i=0;i<nlocal;i++) {
        pressure2[i] = spk->uniface->fetch("CH_pressure2",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},timestamp,s,t);
      }
    } else if (strcmp(arg[k],"temp") == 0) {        // d3
      for (int i=0;i<nlocal;i++) {
        temp[i] = spk->uniface->fetch("CH_temp",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},timestamp,s,t);
      }
    } else error->all(FLERR,"Illegal mui_fetch command");

    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: %s fetched\n",arg[k]);
      fflush(screen);
    }
  }

  return;
}

void AppLotkavolterra::mui_fetch_agg(int narg, char **arg)
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
    if (strcmp(arg[k],"pressure1") == 0) {
      // get info for each FHD cell
      for (int n=0;n<nlocalFHDcell;n++)
        MUIdblval[n] = spk->uniface->fetch("CH_pressure1",{xFHD[n],yFHD[n]},timestamp,s,t);
      // distribute info to each KMC site
      for (int i=0;i<nlocal;i++)
        pressure1[i] = MUIdblval[localFHDcell[i]];
    }
    else if (strcmp(arg[k],"pressure2") == 0) {
      // get info for each FHD cell
      for (int n=0;n<nlocalFHDcell;n++)
        MUIdblval[n] = spk->uniface->fetch("CH_pressure2",{xFHD[n],yFHD[n]},timestamp,s,t);
      // distribute info to each KMC site
      for (int i=0;i<nlocal;i++)
        pressure2[i] = MUIdblval[localFHDcell[i]];
    else if (strcmp(arg[k],"temp") == 0) {
      // get info for each FHD cell
      for (int n=0;n<nlocalFHDcell;n++)
        MUIdblval[n] = spk->uniface->fetch("CH_temp",{xFHD[n],yFHD[n]},timestamp,s,t);
      // distribute info to each KMC site
      for (int i=0;i<nlocal;i++)
        temp[i] = MUIdblval[localFHDcell[i]];
    }
    else error->all(FLERR,"Illegal mui_fetch_agg command");

    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: %s fetched\n",arg[k]);
      fflush(screen);
    }
  }

  return;
}

void AppLotkavolterra::mui_forget(int narg, char **arg)
{
    int timestamp = atoi(arg[0]);

    spk->uniface->forget(timestamp);

    if (domain->me == 0 && screen) {
        fprintf(screen,"** DEBUG: mui forget at timestamp %d\n",timestamp);
        fflush(screen);
    }

    return;
}

#elif defined(USE_AMREX_MPMD)

void AppLotkavolterra::amrex_init_agg ()
{
    AMREX_ASSERT(nlocal>0);
    AMREX_ASSERT(amrex_fhd_lattice_size_x>0);
    AMREX_ASSERT(amrex_fhd_lattice_size_y>0);

    // 1. nlocalFHDcell, nlocalFHDcell_world

    int nFHDcellx = std::rint((domain->boxxhi-domain->boxxlo)/amrex_fhd_lattice_size_x);
    int nFHDcelly = std::rint((domain->boxyhi-domain->boxylo)/amrex_fhd_lattice_size_y);

    if (domain->me == 0) {
        std::fprintf(logfile,"(boxx, boxy) = %e %e\n",domain->boxxhi-domain->boxxlo,domain->boxyhi-domain->boxylo);
        std::fprintf(logfile,"(amrex_fhd_lattice_size_x, amrex_fhd_lattice_size_y) = %e %e\n",amrex_fhd_lattice_size_x,amrex_fhd_lattice_size_y);
        std::fprintf(logfile,"(nFHDcellx, nFHDcelly) = %d %d\n",nFHDcellx,nFHDcelly);
    }

    amrex::Vector<amrex::Vector<int>> cntKMCsite
        (nFHDcellx,amrex::Vector<int>(nFHDcelly,0));
    amrex::Vector<amrex::Vector<double>> sum1
        (nFHDcellx,amrex::Vector<double>(nFHDcelly,0.));
    amrex::Vector<amrex::Vector<double>> sum2
        (nFHDcellx,amrex::Vector<double>(nFHDcelly,0.));

    for (int i = 0; i < nlocal; i++) {
        int nx = std::floor((xyz[i][0]+amrex_kmc_lattice_offset_x)
                            /amrex_fhd_lattice_size_x);
        int ny = std::floor((xyz[i][1]+amrex_kmc_lattice_offset_y)
                            /amrex_fhd_lattice_size_y);
        cntKMCsite[nx][ny]++;
        sum1[nx][ny] += xyz[i][0]+amrex_kmc_lattice_offset_x;
        sum2[nx][ny] += xyz[i][1]+amrex_kmc_lattice_offset_y;
    }

    int cntFHDcell = 0;
    for (int nx = 0; nx < nFHDcellx; nx++) {
        for (int ny = 0; ny < nFHDcelly; ny++) {
            if (cntKMCsite[nx][ny] > 0) cntFHDcell++;
        }
    }
    nlocalFHDcell = cntFHDcell;

    if (domain->me == 0) {
        nlocalFHDcell_world.resize(domain->nprocs);
    }

    MPI_Gather(&nlocalFHDcell,1,MPI_INT,nlocalFHDcell_world.data(),1,MPI_INT,0,world);

    if (domain->me == 0) {
        for (int i=0;i<domain->nprocs;i++) {
            std::fprintf(logfile, "- rank %d: nlocalFHDcell = %d\n",
                         i, nlocalFHDcell_world[i]);
        }
    }

    // 2. xFHD, yFHD

    xFHD.resize(nlocalFHDcell);
    yFHD.resize(nlocalFHDcell);

    amrex::Vector<amrex::Vector<int>> FHDcell(nFHDcellx,amrex::Vector<int>(nFHDcelly,-1));

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
    AMREX_ASSERT(cntFHDcell==nlocalFHDcell);

    if (domain->me == 0) {
        // output my info first
        std::fprintf(logfile,"** rank %d **\n",domain->me);
        for (int j=0;j<nlocalFHDcell;j++) {
            std::fprintf(logfile,"- (xFHD[%d], yFHD[%d]) = %e %e\n",j,j,xFHD[j],yFHD[j]);
        }
        // output for other procs
        for (int i=1;i<domain->nprocs;i++) {
            amrex::Vector<double> data1(nlocalFHDcell_world[i]);
            amrex::Vector<double> data2(nlocalFHDcell_world[i]);

            MPI_Recv(data1.data(),data1.size(),MPI_DOUBLE,i,0,world,MPI_STATUS_IGNORE);
            MPI_Recv(data2.data(),data2.size(),MPI_DOUBLE,i,1,world,MPI_STATUS_IGNORE);

            std::fprintf(logfile,"** rank %d **\n",i);
            for (int j=0;j<nlocalFHDcell_world[i];j++) {
                std::fprintf(logfile,"- (xFHD[%d], yFHD[%d]) = %e %e\n",
                             j,j,data1[j],data2[j]);
            }
        }
    }
    else {
        MPI_Send(xFHD.data(),nlocalFHDcell,MPI_DOUBLE,0,0,world);
        MPI_Send(yFHD.data(),nlocalFHDcell,MPI_DOUBLE,0,1,world);
    }

    // 3. localFHDcell, AMREXintval, AMREXdblval

    localFHDcell.resize(nlocal);

    for (int i = 0; i < nlocal; i++) {
        int nx = std::floor((xyz[i][0]+amrex_kmc_lattice_offset_x)
                            /amrex_fhd_lattice_size_x);
        int ny = std::floor((xyz[i][1]+amrex_kmc_lattice_offset_y)
                            /amrex_fhd_lattice_size_y);
        AMREX_ASSERT(nx>=0 && nx<nFHDcellx && ny>=0 && ny<nFHDcelly);
        localFHDcell[i] = FHDcell[nx][ny];
    }

    amrex::Vector<int> cntKMCsite2(nlocalFHDcell,0);
    for (int i = 0; i < nlocal; i++) cntKMCsite2[localFHDcell[i]]++;

    if (domain->me == 0) {
        // output my info first
        std::fprintf(logfile,"** rank %d **\n",domain->me);
        for (int j=0;j<nlocalFHDcell;j++) {
            std::fprintf(logfile,"- FHDcell %d has %d KMC sites\n",j,cntKMCsite2[j]);
        }
        // output for other procs
        for (int i=1;i<domain->nprocs;i++) {
            amrex::Vector<int> data(nlocalFHDcell_world[i]);

            MPI_Recv(data.data(),data.size(),MPI_INT,i,0,world,MPI_STATUS_IGNORE);

            std::fprintf(logfile,"** rank %d **\n",i);
            for (int j=0;j<nlocalFHDcell_world[i];j++) {
                std::fprintf(logfile,"- FHDcell %d has %d KMC sites\n",j,data[j]);
            }
        }
    }
    else {
        MPI_Send(cntKMCsite2.data(),nlocalFHDcell,MPI_INT,0,0,world);
    }

    intval.resize(nlocalFHDcell);
    dblval.resize(nlocalFHDcell);

    if (domain->me == 0) std::fflush(logfile);

    AMREX_ALWAYS_ASSERT(domain->dimension == 2);

    amrex::Vector<int> proc;
    proc.push_back(amrex::ParallelDescriptor::MyProc());
    amrex::Vector<int> allprocs(amrex::ParallelDescriptor::NProcs());
    amrex::ParallelAllGather::AllGather(proc.data(), 1, allprocs.data(),
                                        amrex::ParallelDescriptor::Communicator());
    amrex::DistributionMapping dmap2(std::move(allprocs));

    double dx = amrex_fhd_lattice_size_x;
    double dy = amrex_fhd_lattice_size_y;
    auto xmm = std::minmax_element(xFHD.begin(), xFHD.end());
    auto ymm = std::minmax_element(yFHD.begin(), yFHD.end());
    int xlo = static_cast<int>(std::floor((*xmm.first -domain->boxxlo)/dx));
    int ylo = static_cast<int>(std::floor((*ymm.first -domain->boxylo)/dy));
    int xhi = static_cast<int>(std::floor((*xmm.second-domain->boxxlo)/dx));
    int yhi = static_cast<int>(std::floor((*ymm.second-domain->boxylo)/dy));
    AMREX_ALWAYS_ASSERT(nlocalFHDcell==(xhi-xlo+1)*(yhi-ylo+1));
    amrex::Vector<amrex::Box> box{amrex::Box(amrex::IntVect(xlo,ylo,0),
                                             amrex::IntVect(xhi,yhi,0))};
    amrex::AllGatherBoxes(box);
    amrex::BoxArray ba2(box.data(), box.size());

    amrex::BoxArray ba;
    amrex::DistributionMapping dmap;
    amrex::Box domainbox = ba2.minimalBox();
    if (domainbox.numPts() == ba2.numPts()) {
        ba = ba2;
        dmap = dmap2;
    } else {
        mf2 = std::make_unique<amrex::MultiFab>(ba2, dmap2, 1, 0,
                                  amrex::MFInfo().SetArena(amrex::The_Cpu_Arena()));
        imf2 = std::make_unique<amrex::iMultiFab>(ba2, dmap2, 1, 0,
                                  amrex::MFInfo().SetArena(amrex::The_Cpu_Arena()));
        ba = amrex::BoxArray(domainbox);
        ba.maxSize(16);
        dmap = amrex::DistributionMapping(ba);
    }

    mf.define(ba, dmap, 1, 0, amrex::MFInfo().SetArena(amrex::The_Cpu_Arena()));
    imf.define(ba, dmap, 1, 0, amrex::MFInfo().SetArena(amrex::The_Cpu_Arena()));

    mpmd_copier = std::make_unique<amrex::MPMD::Copier>(ba, dmap);
}

void AppLotkavolterra::amrex_push_agg(int narg, char **arg)
{
    int timestamp = atoi(arg[0]);

    if (domain->me == 0 && screen) {
        std::fprintf(screen,"** DEBUG: amrex_push_agg at timestamp %d\n",timestamp);
        std::fflush(screen);
    }

    for (int k=1;k<narg;k++)
    {
        if (std::strcmp(arg[k],"ac1") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                intval[localFHDcell[i]] += ac1[i];
                ac1[i] = 0;
            }
            amrex_send_intval();
        } else if (std::strcmp(arg[k],"ac2") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                intval[localFHDcell[i]] += ac2[i];
                ac2[i] = 0;
            }
            amrex_send_intval();
        } else if (std::strcmp(arg[k],"dc1") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                intval[localFHDcell[i]] += dc1[i];
                dc1[i] = 0;
            }
            amrex_send_intval();
        } else if (std::strcmp(arg[k],"dc2") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                intval[localFHDcell[i]] += dc2[i];
                dc2[i] = 0;
            }
            amrex_send_intval();
        } else if (std::strcmp(arg[k],"occ1") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                int is_occ = (element[i]==1) ? 1 : 0;
                intval[localFHDcell[i]] += is_occ;
            }
            amrex_send_intval();
        } else if (std::strcmp(arg[k],"occ2") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                int is_occ = (element[i]==2) ? 1 : 0;
                intval[localFHDcell[i]] += is_occ;
            }
            amrex_send_intval();
        } else if (std::strcmp(arg[k],"one") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) intval[localFHDcell[i]]++;
            amrex_send_intval();
        } else error->all(FLERR,"Illegal amrex_push_agg command");

        if (domain->me == 0 && screen) {
            std::fprintf(screen,"** DEBUG: %s pushed\n",arg[k]);
            std::fflush(screen);
        }
    }
}

void AppLotkavolterra::amrex_fetch_agg(int narg, char **arg)
{
    int timestamp = atoi(arg[0]);

    if (domain->me == 0 && screen) {
        std::fprintf(screen,"** DEBUG: amrex_fetch_agg at timestamp %d\n",timestamp);
        std::fflush(screen);
    }

    if (amrex_fhd_lattice_size_x <= 0. || amrex_fhd_lattice_size_y <= 0.)
        error->all(FLERR,"amrex_fhd_lattice_size must be set as two positive numbers");

    for (int k=1;k<narg;k++) {
        if (std::strcmp(arg[k],"pressure1") == 0) {
            amrex_recv_dblval();
            // distribute info to each KMC site
            for (int i=0;i<nlocal;i++)
                pressure1[i] = dblval[localFHDcell[i]];
      } else if (std::strcmp(arg[k],"pressure2") == 0) {
            amrex_recv_dblval();
            // distribute info to each KMC site
            for (int i=0;i<nlocal;i++)
                pressure2[i] = dblval[localFHDcell[i]];
      } else if (std::strcmp(arg[k],"temp") == 0) {
            amrex_recv_dblval();
            // distribute info to each KMC site
            for (int i=0;i<nlocal;i++)
                temp[i] = dblval[localFHDcell[i]];
      } else error->all(FLERR,"Illegal amrex_fetch_agg command");

      if (domain->me == 0 && screen) {
          std::fprintf(screen,"** DEBUG: %s fetched\n",arg[k]);
          std::fflush(screen);
      }
    }
}

void AppLotkavolterra::amrex_send_intval()
{
    auto& local_imf = imf2 ? *imf2 : imf;

    for (amrex::MFIter mfi(local_imf); mfi.isValid(); ++mfi) {
        amrex::Box const& b = mfi.validbox();
        int const ylen = b.length(1);
        int const offset = b.smallEnd(1) + b.smallEnd(0) * ylen;
        amrex::Array4<int> const& ifab = local_imf.array(mfi);
        int const* p = intval.data();
        amrex::LoopOnCpu(b, [&] (int i, int j, int) noexcept
        {
            ifab(i,j,0) = p[j+i*ylen-offset];
        });
    }

    if (imf2) {
        imf.setVal(0);
        imf.ParallelAdd(*imf2);
    }

    mpmd_copier->send(imf,0,1);
}

void AppLotkavolterra::amrex_recv_dblval()
{
    mpmd_copier->recv(mf,0,1);

    auto& local_mf = mf2 ? *mf2 : mf;
    if (mf2) {
        mf2->ParallelCopy(mf);
    }

    for (amrex::MFIter mfi(local_mf); mfi.isValid(); ++mfi) {
        amrex::Box const& b = mfi.validbox();
        int const ylen = b.length(1);
        int const offset = b.smallEnd(1) + b.smallEnd(0) * ylen;
        amrex::Array4<amrex::Real const> const& fab = local_mf.const_array(mfi);
        double* p = dblval.data();
        amrex::LoopOnCpu(b, [&] (int i, int j, int) noexcept
        {
            p[j+i*ylen-offset] = fab(i,j,0);
        });
    }
}

#endif
