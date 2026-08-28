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

#if defined(USE_AMREX_MPMD)
#include <AMReX_MPMD.H>
#endif

#include <vector>

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
#if defined(USE_AMREX_MPMD)
  ads_wall_dir = 2;
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

#if defined(USE_AMREX_MPMD)
  else if (strcmp(command,"amrex_ads_wall_dir") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal amrex_ads_wall_dir command");
    ads_wall_dir = atoi(arg[0]);
    if (ads_wall_dir < 0 || ads_wall_dir > 2) error->all(FLERR,"Illegal amrex_ads_wall_dir command");
  } else if (strcmp(command,"amrex_init_agg") == 0) {
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
    for (int jj = 0; jj < numneigh[i]; jj++) {
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
    for (int jj = 0; jj < numneigh[i]; jj++) {
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

  for (n = 0; n < numneigh[i]; n++) {
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >=0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
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

#if defined(USE_AMREX_MPMD)

void AppLotkavolterra::amrex_init_agg ()
{
    AMREX_ASSERT(nlocal>0);
    AMREX_ASSERT(amrex_fhd_lattice_size_x>0);
    AMREX_ASSERT(amrex_fhd_lattice_size_y>0);

    // 0. ads_wall_dir
    dir1 = (ads_wall_dir == 0) ? 1 : 0;
    dir2 = (ads_wall_dir == 2) ? 1 : 2;

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
    amrex::IntVect blo(0), bhi(0);
    blo[dir1] = xlo;
    blo[dir2] = ylo;
    bhi[dir1] = xhi;
    bhi[dir2] = yhi;
    amrex::Vector<amrex::Box> box{amrex::Box(blo,bhi)};
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
        int const len2 = b.length(dir2);
        int const offset = b.smallEnd(dir2) + b.smallEnd(dir1) * len2;
        amrex::Array4<int> const& ifab = local_imf.array(mfi);
        int const* p = intval.data();
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) noexcept
        {
            int const idx[3] = {i,j,k};
            ifab(i,j,k) = p[idx[dir2]+idx[dir1]*len2-offset];
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
        int const len2 = b.length(dir2);
        int const offset = b.smallEnd(dir2) + b.smallEnd(dir1) * len2;
        amrex::Array4<amrex::Real const> const& fab = local_mf.const_array(mfi);
        double* p = dblval.data();
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) noexcept
        {
            int const idx[3] = {i,j,k};
            p[idx[dir2]+idx[dir1]*len2-offset] = fab(i,j,k);
        });
    }
}

#endif
