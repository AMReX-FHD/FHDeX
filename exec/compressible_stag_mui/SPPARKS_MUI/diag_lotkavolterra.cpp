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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "diag_lotkavolterra.h"
#include "app.h"
#include "app_lotkavolterra.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

enum{VACANCY,SPEC1,SPEC2};  // removed ZERO and moved VACANCY to first item   // same as AppLotkavolterra
enum{VAC,SP1,SP2,EVENTS,PREY,PREDATION,PREDATOR}; // moved VAC to first item

/* ---------------------------------------------------------------------- */

DiagLotkavolterra::DiagLotkavolterra(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"lotkavolterra") != 0)
    error->all(FLERR,"Diag_style lotkavolterra requires app_style lotkavolterra");

  nlist = 0;

  int iarg = iarg_child;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"list") == 0) {
      nlist = narg - iarg - 1;
      list = new char*[nlist];
      int j = 0;
      for (int i = iarg+1; i < narg; i++) {
        int n = strlen(arg[i]) + 1;
        list[j] = new char[n];
        strcpy(list[j],arg[i]);
        j++;
      }
      iarg = narg;
    } else error->all(FLERR,"Illegal diag_style lotkavolterra command");
  }

  if (nlist == 0) error->all(FLERR,"Illegal diag_style lotkavolterra command");
  which = new int[nlist];
  index = new int[nlist];
  ivector = new int[nlist];
}

/* ---------------------------------------------------------------------- */

DiagLotkavolterra::~DiagLotkavolterra()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] ivector;
}

/* ---------------------------------------------------------------------- */

void DiagLotkavolterra::init()
{
  applotkavolterra = (AppLotkavolterra *) app;

  int nprey = applotkavolterra->nprey;
  int npredation = applotkavolterra->npredation;
  int npredator = applotkavolterra->npredator;

  for (int i = 0; i < nlist; i++) {
    if (strcmp(list[i],"vac") == 0) which[i] = VAC;
    else if (strcmp(list[i],"spec1") == 0) which[i] = SP1;
    else if (strcmp(list[i],"spec2") == 0) which[i] = SP2;
    else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;
    else if (list[i][0] == 'O') {
      which[i] = PREY;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nprey)
        error->all(FLERR,"Invalid value setting in diag_style lotkavolterra");
      index[i] = n - 1;
    } else if (list[i][0] == 'A') {
      which[i] = PREDATION;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > npredation)
        error->all(FLERR,"Invalid value setting in diag_style lotkavolterra");
      index[i] = n - 1;
    } else if (list[i][0] == 'B') {
      which[i] = PREDATOR;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > npredator)
        error->all(FLERR,"Invalid value setting in diag_style lotkavolterra");
      index[i] = n - 1;
    } else error->all(FLERR,"Invalid value setting in diag_style lotkavolterra");
  }

  siteflag = 0;
  for (int i = 0; i < nlist; i++)
    if (which[i] == SP1 || which[i] == SP2 || which[i] == VAC)
      siteflag = 1;

  for (int i = 0; i < nlist; i++) ivector[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DiagLotkavolterra::compute()
{
  int sites[3],ivalue;

  if (siteflag) {
    sites[SPEC1] = sites[SPEC2] = sites[VACANCY] = 0;
    int *element = applotkavolterra->element;
    int nlocal = applotkavolterra->nlocal;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == SP1) ivalue = sites[SPEC1];
    else if (which[i] == SP2) ivalue = sites[SPEC2];
    else if (which[i] == VAC) ivalue = sites[VACANCY];

    else if (which[i] == EVENTS) ivalue = applotkavolterra->nevents;
    else if (which[i] == PREY) ivalue = applotkavolterra->preycount[index[i]];
    else if (which[i] == PREDATION) ivalue = applotkavolterra->predationcount[index[i]];
    else if (which[i] == PREDATOR) ivalue = applotkavolterra->predatorcount[index[i]];
    MPI_Allreduce(&ivalue,&ivector[i],1,MPI_INT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagLotkavolterra::stats(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str,"\t%d",ivector[i]);
    str += strlen(str);
  }
}

/* ---------------------------------------------------------------------- */

void DiagLotkavolterra::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str,"\t%s",list[i]);
    str += strlen(str);
  }
}
