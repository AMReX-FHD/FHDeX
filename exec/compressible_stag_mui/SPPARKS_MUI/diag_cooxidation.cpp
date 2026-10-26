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
#include "diag_cooxidation.h"
#include "app.h"
#include "app_cooxidation.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

enum{VACANCY,SPEC1,SPEC2,SPEC3};  // same as AppCooxidation
enum{VAC,SP1,SP2,SP3,EVENTS,ONE,TWO,ADS,DES,DISSOCADS,ASSOCDES,REACTION,RXNTOTAL};

/* ---------------------------------------------------------------------- */

DiagCooxidation::DiagCooxidation(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"cooxidation") != 0)
    error->all(FLERR,"Diag_style cooxidation requires app_style cooxidation");

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
    } else error->all(FLERR,"Illegal diag_style cooxidation command");
  }

  if (nlist == 0) error->all(FLERR,"Illegal diag_style cooxidation command");
  which = new int[nlist];
  index = new int[nlist];
  ivector = new int[nlist];
}

/* ---------------------------------------------------------------------- */

DiagCooxidation::~DiagCooxidation()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] ivector;
}

/* ---------------------------------------------------------------------- */

void DiagCooxidation::init()
{
  appcooxidation = (AppCooxidation *) app;

  int none = appcooxidation->none;
  int ntwo = appcooxidation->ntwo;
  int nads = appcooxidation->nads;
  int ndes = appcooxidation->ndes;
  int ndissocads = appcooxidation->ndissocads;
  int nassocdes = appcooxidation->nassocdes;
  int nreaction = appcooxidation->nreaction;

  for (int i = 0; i < nlist; i++) {
    if (strcmp(list[i],"spec1") == 0) which[i] = SP1;
    else if (strcmp(list[i],"spec2") == 0) which[i] = SP2;
    else if (strcmp(list[i],"spec3") == 0) which[i] = SP3;
    else if (strcmp(list[i],"vac") == 0) which[i] = VAC;
    else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;
    else if (strcmp(list[i],"reaction") == 0) which[i] = RXNTOTAL;
    else if (list[i][0] == 's') {
      which[i] = ONE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > none)
        error->all(FLERR,"Invalid value setting in diag_style cooxidation");
      index[i] = n - 1;
    } else if (list[i][0] == 'd') {
      which[i] = TWO;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ntwo)
        error->all(FLERR,"Invalid value setting in diag_style cooxidation");
      index[i] = n - 1;
    } else if (list[i][0] == 'A') {
      which[i] = ADS;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nads)
        error->all(FLERR,"Invalid value setting in diag_style cooxidation");
      index[i] = n - 1;
    } else if (list[i][0] == 'D') {
      which[i] = DES;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ndes)
        error->all(FLERR,"Invalid value setting in diag_style cooxidation");
      index[i] = n - 1;
    } else if (list[i][0] == 'B') {
      which[i] = DISSOCADS;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ndissocads)
        error->all(FLERR,"Invalid value setting in diag_style cooxidation");
      index[i] = n - 1;
    } else if (list[i][0] == 'E') {
      which[i] = ASSOCDES;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nassocdes)
        error->all(FLERR,"Invalid value setting in diag_style cooxidation");
      index[i] = n - 1;
    } else if (list[i][0] == 'R') {
      which[i] = REACTION;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nreaction)
        error->all(FLERR,"Invalid value setting in diag_style cooxidation");
      index[i] = n - 1;
    } else error->all(FLERR,"Invalid value setting in diag_style cooxidation");
  }

  siteflag = 0;
  for (int i = 0; i < nlist; i++)
    if (which[i] == SP1 || which[i] == SP2 || which[i] == SP3 || which[i] == VAC)
      siteflag = 1;

  for (int i = 0; i < nlist; i++) ivector[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DiagCooxidation::compute()
{
  int sites[4],ivalue;

  if (siteflag) {
    sites[SPEC1] = sites[SPEC2] = sites[SPEC3] = sites[VACANCY] = 0;
    int *element = appcooxidation->element;
    int nlocal = appcooxidation->nlocal;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == SP1) ivalue = sites[SPEC1];
    else if (which[i] == SP2) ivalue = sites[SPEC2];
    else if (which[i] == SP3) ivalue = sites[SPEC3];
    else if (which[i] == VAC) ivalue = sites[VACANCY];
    else if (which[i] == EVENTS) ivalue = appcooxidation->nevents;
    else if (which[i] == ONE) ivalue = appcooxidation->scount[index[i]];
    else if (which[i] == TWO) ivalue = appcooxidation->dcount[index[i]];
    else if (which[i] == ADS) ivalue = appcooxidation->adscount[index[i]];
    else if (which[i] == DES) ivalue = appcooxidation->descount[index[i]];
    else if (which[i] == DISSOCADS) ivalue = appcooxidation->dadscount[index[i]];
    else if (which[i] == ASSOCDES) ivalue = appcooxidation->adescount[index[i]];
    else if (which[i] == REACTION) ivalue = appcooxidation->rxncount[index[i]];
    else if (which[i] == RXNTOTAL) ivalue = appcooxidation->rxnsumcount;
    MPI_Allreduce(&ivalue,&ivector[i],1,MPI_INT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagCooxidation::stats(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str,"\t%d",ivector[i]);
    str += strlen(str);
  }
}

/* ---------------------------------------------------------------------- */

void DiagCooxidation::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str,"\t%s",list[i]);
    str += strlen(str);
  }
}
