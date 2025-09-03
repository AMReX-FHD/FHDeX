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
#include "diag_surfchemtest.h"
#include "app.h"
#include "app_surfchemtest.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

enum{VACANCY,SPEC1,SPEC2,SPEC3,SPEC4,SPEC5};  // removed ZERO and moved VACANCY to first item   // same as AppSurfchemtest
enum{VAC,SP1,SP2,SP3,SP4,SP5,EVENTS,ONE,TWO,THREE,ADS,DES,DISSOCADS,ASSOCDES,REACTION,RXNTOTAL}; // moved VAC to first item

/* ---------------------------------------------------------------------- */

DiagSurfchemtest::DiagSurfchemtest(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"surfchemtest") != 0)
    error->all(FLERR,"Diag_style surfchemtest requires app_style surfchemtest");

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
    } else error->all(FLERR,"Illegal diag_style surfchemtest command");
  }

  if (nlist == 0) error->all(FLERR,"Illegal diag_style surfchemtest command");
  which = new int[nlist];
  index = new int[nlist];
  ivector = new int[nlist];
}

/* ---------------------------------------------------------------------- */

DiagSurfchemtest::~DiagSurfchemtest()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] ivector;
}

/* ---------------------------------------------------------------------- */

void DiagSurfchemtest::init()
{
  appsurfchemtest = (AppSurfchemtest *) app;

  int none = appsurfchemtest->none;
  int ntwo = appsurfchemtest->ntwo;
  int nthree = appsurfchemtest->nthree;
  int nads = appsurfchemtest->nads;
  int ndes = appsurfchemtest->ndes;
  int ndissocads = appsurfchemtest->ndissocads;
  int nassocdes = appsurfchemtest->nassocdes;
  int nreaction = appsurfchemtest->nreaction;

  for (int i = 0; i < nlist; i++) {
    if (strcmp(list[i],"spec1") == 0) which[i] = SP1;
    else if (strcmp(list[i],"spec2") == 0) which[i] = SP2;
    else if (strcmp(list[i],"spec3") == 0) which[i] = SP3;
    else if (strcmp(list[i],"spec4") == 0) which[i] = SP4;
    else if (strcmp(list[i],"spec5") == 0) which[i] = SP5;
    else if (strcmp(list[i],"vac") == 0) which[i] = VAC;
    else if (strcmp(list[i],"events") == 0) which[i] = EVENTS;
    else if (strcmp(list[i],"reaction") == 0) which[i] = RXNTOTAL;
    else if (list[i][0] == 's') {
      which[i] = ONE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > none)
        error->all(FLERR,"Invalid value setting in diag_style surfchemtest");
      index[i] = n - 1;
    } else if (list[i][0] == 'd') {
      which[i] = TWO;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ntwo)
        error->all(FLERR,"Invalid value setting in diag_style surfchemtest");
      index[i] = n - 1;
    } else if (list[i][0] == 't') {
      which[i] = THREE;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nthree)
        error->all(FLERR,"Invalid value setting in diag_style surfchemtest");
      index[i] = n - 1;
    } else if (list[i][0] == 'A') {
      which[i] = ADS;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nads)
        error->all(FLERR,"Invalid value setting in diag_style surfchemtest");
      index[i] = n - 1;
    } else if (list[i][0] == 'D') {
      which[i] = DES;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ndes)
        error->all(FLERR,"Invalid value setting in diag_style surfchemtest");
      index[i] = n - 1;
    } else if (list[i][0] == 'B') {
      which[i] = DISSOCADS;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > ndissocads)
        error->all(FLERR,"Invalid value setting in diag_style surfchemtest");
      index[i] = n - 1;
    } else if (list[i][0] == 'E') {
      which[i] = ASSOCDES;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nassocdes)
        error->all(FLERR,"Invalid value setting in diag_style surfchemtest");
      index[i] = n - 1;
    } else if (list[i][0] == 'R') {
      which[i] = REACTION;
      int n = atoi(&list[i][1]);
      if (n < 1 || n > nreaction)
        error->all(FLERR,"Invalid value setting in diag_style surfchemtest");
      index[i] = n - 1;
    } else error->all(FLERR,"Invalid value setting in diag_style surfchemtest");
  }

  siteflag = 0;
  for (int i = 0; i < nlist; i++)
    if (which[i] == SP1 || which[i] == SP2 || which[i] == SP3 || which[i] == SP4 || which[i] == SP5 || which[i] == VAC)
      siteflag = 1;

  for (int i = 0; i < nlist; i++) ivector[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DiagSurfchemtest::compute()
{
  int sites[6],ivalue;

  if (siteflag) {
    sites[SPEC1] = sites[SPEC2] = sites[SPEC3] = sites[SPEC4] = sites[SPEC5] = sites[VACANCY] = 0;
    int *element = appsurfchemtest->element;
    int nlocal = appsurfchemtest->nlocal;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == SP1) ivalue = sites[SPEC1];
    else if (which[i] == SP2) ivalue = sites[SPEC2];
    else if (which[i] == SP3) ivalue = sites[SPEC3];
    else if (which[i] == SP4) ivalue = sites[SPEC4];
    else if (which[i] == SP5) ivalue = sites[SPEC5];
    else if (which[i] == VAC) ivalue = sites[VACANCY];
    else if (which[i] == EVENTS) ivalue = appsurfchemtest->nevents;
    else if (which[i] == ONE) ivalue = appsurfchemtest->scount[index[i]];
    else if (which[i] == TWO) ivalue = appsurfchemtest->dcount[index[i]];
    else if (which[i] == THREE) ivalue = appsurfchemtest->tcount[index[i]];
    else if (which[i] == ADS) ivalue = appsurfchemtest->adscount[index[i]];
    else if (which[i] == DES) ivalue = appsurfchemtest->descount[index[i]];
    else if (which[i] == DISSOCADS) ivalue = appsurfchemtest->dadscount[index[i]];
    else if (which[i] == ASSOCDES) ivalue = appsurfchemtest->adescount[index[i]];
    else if (which[i] == REACTION) ivalue = appsurfchemtest->rxncount[index[i]];
    else if (which[i] == RXNTOTAL) ivalue = appsurfchemtest->rxnsumcount;
    MPI_Allreduce(&ivalue,&ivector[i],1,MPI_INT,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagSurfchemtest::stats(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str,"\t%d",ivector[i]);
    str += strlen(str);
  }
}

/* ---------------------------------------------------------------------- */

void DiagSurfchemtest::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str,"\t%s",list[i]);
    str += strlen(str);
  }
}