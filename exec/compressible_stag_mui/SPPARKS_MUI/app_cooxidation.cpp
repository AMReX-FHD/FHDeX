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
#include "app_cooxidation.h"
#include "solve.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "domain.h"


#if defined(USE_AMREX_MPMD)
#include <AMReX_MPMD.H>
#endif


using namespace SPPARKS_NS;

using namespace std;

enum{NOOP,SITEA,SITEB,SITEC};
enum{VACANCY,SPEC1,SPEC2,SPEC3}; // removed ZERO and moved VACANCY to first item

#define DELTAEVENT 100000

static const double BETA = -0.5; // temperature exponent of FHD-fed adsorption rates: k p (T/T0)^BETA

/* ---------------------------------------------------------------------- */

AppCooxidation::AppCooxidation(SPPARKS *spk, int narg, char **arg) :
  AppLattice(spk,narg,arg)
{
  ninteger = 14;  // type: site type
                  // element: site state
                  // ac1/ac2/ac3: adsorption count
                  // dc1/dc2/dc3: desorption count
                  // dac1/dac2/dac3: dissociative adsorption count
                  // adc1/adc2/adc3: associative desorption count
  ndouble = 4;    // pressure1/pressure2/pressure3: partial pressure of the contacting FHD cell
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

  ads_is_rate = dads_is_rate = NULL;

  // reaction lists
  none = ntwo = nads = ndes = ndissocads = nassocdes = nreaction = rxnsumcount = 0;
  srate = drate = adsrate = desrate = dadsrate = adesrate = rxnrate = NULL;
  spropensity = dpropensity = adespropensity = rxnpropensity = NULL;   // no adspropensity/dadspropensity/despropensity
  stype = sinput = soutput = NULL;
  dtype = dinput = doutput = NULL;
  adstype = adsinput = adsoutput = NULL;
  destype = desinput = desoutput = NULL;
  dadstype = dadsinput = dadsoutput = NULL;
  adestype = adesinput = adesoutput = NULL;
  rxntype = rxninput = rxnoutput = NULL;
  scount = dcount = adscount = descount = dadscount = adescount = rxncount = NULL;
  adsadsorbate = dadsadsorbate = adesdesorbate = reactionsorbate = NULL;

#if defined(USE_AMREX_MPMD)
  ads_wall_dir = 2;

  update_full = 0; // update ads events (0) or update all events (1)
  bordersite = NULL;
  cellrep = NULL;
  adsprop_cell = NULL;
  dadsprop_cell = NULL;
#endif
}

/* ---------------------------------------------------------------------- */

AppCooxidation::~AppCooxidation()
{
  delete [] esites;
  delete [] echeck;

  memory->sfree(events);
  memory->destroy(firstevent);
  memory->destroy(srate);
  memory->destroy(drate);
  memory->destroy(adsrate);
  memory->destroy(ads_is_rate);
  memory->destroy(dads_is_rate);
  memory->destroy(desrate);
  memory->destroy(dadsrate);
  memory->destroy(adesrate);
  memory->destroy(rxnrate);
  memory->destroy(spropensity);
  memory->destroy(dpropensity);
  memory->destroy(adespropensity);
  memory->destroy(rxnpropensity);
  memory->destroy(stype);
  memory->destroy(sinput);
  memory->destroy(soutput);
  memory->destroy(dtype);
  memory->destroy(dinput);
  memory->destroy(doutput);
  memory->destroy(adstype);
  memory->destroy(adsinput);
  memory->destroy(adsoutput);
  memory->destroy(destype);
  memory->destroy(desinput);
  memory->destroy(desoutput);
  memory->destroy(dadstype);
  memory->destroy(dadsinput);
  memory->destroy(dadsoutput);
  memory->destroy(adestype);
  memory->destroy(adesinput);
  memory->destroy(adesoutput);
  memory->destroy(rxntype);
  memory->destroy(rxninput);
  memory->destroy(rxnoutput);
  memory->destroy(scount);
  memory->destroy(dcount);
  memory->destroy(adscount);
  memory->destroy(descount);
  memory->destroy(dadscount);
  memory->destroy(adescount);
  memory->destroy(rxncount);
  memory->destroy(adsadsorbate);
  memory->destroy(dadsadsorbate);
  memory->destroy(adesdesorbate);
  memory->destroy(reactionsorbate);

#if defined(USE_AMREX_MPMD)
  delete [] bordersite;
  memory->destroy(cellrep);
  memory->destroy(adsprop_cell);
  memory->destroy(dadsprop_cell);
#endif
}

/* ---------------------------------------------------------------------- */

void AppCooxidation::input_app(char *command, int narg, char **arg)
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
      else if (strcmp(arg[2],"vac") == 0) sinput[none] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      srate[none] = atof(arg[3]);

      if (strcmp(arg[4],"spec1") == 0) soutput[none] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) soutput[none] = SPEC2;
      else if (strcmp(arg[4],"spec3") == 0) soutput[none] = SPEC3;
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
      else if (strcmp(arg[3],"vac") == 0) dinput[ntwo][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[4],"spec1") == 0) dinput[ntwo][1] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) dinput[ntwo][1] = SPEC2;
      else if (strcmp(arg[4],"spec3") == 0) dinput[ntwo][1] = SPEC3;
      else if (strcmp(arg[4],"vac") == 0) dinput[ntwo][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      drate[ntwo] = atof(arg[5]);

      if (strcmp(arg[6],"spec1") == 0) doutput[ntwo][0] = SPEC1;
      else if (strcmp(arg[6],"spec2") == 0) doutput[ntwo][0] = SPEC2;
      else if (strcmp(arg[6],"spec3") == 0) doutput[ntwo][0] = SPEC3;
      else if (strcmp(arg[6],"vac") == 0) doutput[ntwo][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[7],"spec1") == 0) doutput[ntwo][1] = SPEC1;
      else if (strcmp(arg[7],"spec2") == 0) doutput[ntwo][1] = SPEC2;
      else if (strcmp(arg[7],"spec3") == 0) doutput[ntwo][1] = SPEC3;
      else if (strcmp(arg[7],"vac") == 0) doutput[ntwo][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      ntwo++;

    } else if (rstyle == 3) {   // adsorption
      if (narg != 5 && narg != 7) error->all(FLERR,"Illegal event command");
      ads_is_rate[nads] = true;
      adsadsorbate[nads] = -1;
      if (narg == 7) { // FHD-fed: rate constant x pressure x (T/T0)^BETA
        if (strcmp(arg[5],"FHD") != 0) error->all(FLERR,"Illegal event command");
        ads_is_rate[nads] = false;
        if (strcmp(arg[6],"spec1") == 0) adsadsorbate[nads] = SPEC1;
        else if (strcmp(arg[6],"spec2") == 0) adsadsorbate[nads] = SPEC2;
        else if (strcmp(arg[6],"spec3") == 0) adsadsorbate[nads] = SPEC3;
        else error->all(FLERR,"Illegal event command");
      }
      if (strcmp(arg[1],"siteA") == 0) adstype[nads] = SITEA;
      else if (strcmp(arg[1],"siteB") == 0) adstype[nads] = SITEB;
      else if (strcmp(arg[1],"siteC") == 0) adstype[nads] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[2],"spec1") == 0)
        error->all(FLERR,"rstyle=3 only allows vac->spec1/2/3");
      else if (strcmp(arg[2],"spec2") == 0)
        error->all(FLERR,"rstyle=3 only allows vac->spec1/2/3");
      else if (strcmp(arg[2],"spec3") == 0)
        error->all(FLERR,"rstyle=3 only allows vac->spec1/2/3");
      else if (strcmp(arg[2],"vac") == 0) adsinput[nads] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      adsrate[nads] = atof(arg[3]);

      if (strcmp(arg[4],"spec1") == 0) adsoutput[nads] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) adsoutput[nads] = SPEC2;
      else if (strcmp(arg[4],"spec3") == 0) adsoutput[nads] = SPEC3;
      else if (strcmp(arg[4],"vac") == 0)
        error->all(FLERR,"rstyle=3 only allows vac->spec1/2/3");
      else error->all(FLERR,"Illegal event command");

      nads++;

    } else if (rstyle == 4) {   // desorption
      if (narg != 5) error->all(FLERR,"Illegal event command");
      if (strcmp(arg[1],"siteA") == 0) destype[ndes] = SITEA;
      else if (strcmp(arg[1],"siteB") == 0) destype[ndes] = SITEB;
      else if (strcmp(arg[1],"siteC") == 0) destype[ndes] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[2],"spec1") == 0) desinput[ndes] = SPEC1;
      else if (strcmp(arg[2],"spec2") == 0) desinput[ndes] = SPEC2;
      else if (strcmp(arg[2],"spec3") == 0) desinput[ndes] = SPEC3;
      else if (strcmp(arg[2],"vac") == 0)
        error->all(FLERR,"rstyle=4 only allows spec1/2/3->vac");
      else error->all(FLERR,"Illegal event command");

      desrate[ndes] = atof(arg[3]);

      if (strcmp(arg[4],"spec1") == 0)
        error->all(FLERR,"rstyle=4 only allows spec1/2/3->vac");
      else if (strcmp(arg[4],"spec2") == 0)
        error->all(FLERR,"rstyle=4 only allows spec1/2/3->vac");
      else if (strcmp(arg[4],"spec3") == 0)
        error->all(FLERR,"rstyle=4 only allows spec1/2/3->vac");
      else if (strcmp(arg[4],"vac") == 0) desoutput[ndes] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      ndes++;

    } else if (rstyle == 5) { // disassociative adsorption
      if (narg != 8 && narg != 10) error->all(FLERR,"Illegal event command");
      dads_is_rate[ndissocads] = true;
      dadsadsorbate[ndissocads] = -1;
      if (narg == 10) { // FHD-fed: rate constant x pressure x (T/T0)^BETA
        if (strcmp(arg[8],"FHD") != 0) error->all(FLERR,"Illegal event command");
        dads_is_rate[ndissocads] = false;
        if (strcmp(arg[9],"spec1") == 0) dadsadsorbate[ndissocads] = SPEC1;
        else if (strcmp(arg[9],"spec2") == 0) dadsadsorbate[ndissocads] = SPEC2;
        else if (strcmp(arg[9],"spec3") == 0) dadsadsorbate[ndissocads] = SPEC3;
        else error->all(FLERR,"Illegal event command");
      }
      if (strcmp(arg[1],"siteA") == 0) dadstype[ndissocads][0] = SITEA;
      else if (strcmp(arg[1],"siteB") == 0) dadstype[ndissocads][0] = SITEB;
      else if (strcmp(arg[1],"siteC") == 0) dadstype[ndissocads][0] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[2],"siteA") == 0) dadstype[ndissocads][1] = SITEA;
      else if (strcmp(arg[2],"siteB") == 0) dadstype[ndissocads][1] = SITEB;
      else if (strcmp(arg[2],"siteC") == 0) dadstype[ndissocads][1] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[3],"spec1") == 0)
        error->all(FLERR,"rstyle=5 only allows vac->spec1/2/3");
      else if (strcmp(arg[3],"spec2") == 0)
        error->all(FLERR,"rstyle=5 only allows vac->spec1/2/3");
      else if (strcmp(arg[3],"spec3") == 0)
        error->all(FLERR,"rstyle=5 only allows vac->spec1/2/3");
      else if (strcmp(arg[3],"vac") == 0) dadsinput[ndissocads][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[4],"spec1") == 0)
        error->all(FLERR,"rstyle=5 only allows vac->spec1/2/3");
      else if (strcmp(arg[4],"spec2") == 0)
        error->all(FLERR,"rstyle=5 only allows vac->spec1/2/3");
      else if (strcmp(arg[4],"spec3") == 0)
        error->all(FLERR,"rstyle=5 only allows vac->spec1/2/3");
      else if (strcmp(arg[4],"vac") == 0) dadsinput[ndissocads][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      dadsrate[ndissocads] = atof(arg[5]);

      if (strcmp(arg[6],"spec1") == 0) dadsoutput[ndissocads][0] = SPEC1;
      else if (strcmp(arg[6],"spec2") == 0) dadsoutput[ndissocads][0] = SPEC2;
      else if (strcmp(arg[6],"spec3") == 0) dadsoutput[ndissocads][0] = SPEC3;
      else if (strcmp(arg[6],"vac") == 0)
        error->all(FLERR,"rstyle=5 only allows vac->spec1/2/3");
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[7],"spec1") == 0) dadsoutput[ndissocads][1] = SPEC1;
      else if (strcmp(arg[7],"spec2") == 0) dadsoutput[ndissocads][1] = SPEC2;
      else if (strcmp(arg[7],"spec3") == 0) dadsoutput[ndissocads][1] = SPEC3;
      else if (strcmp(arg[7],"vac") == 0)
        error->all(FLERR,"rstyle=5 only allows vac->spec1/2/3");
      ndissocads++;

    } else if (rstyle == 6) { // associative desorption
      if (narg != 9) error->all(FLERR,"Illegal event command");
      if (strcmp(arg[1],"siteA") == 0) adestype[nassocdes][0] = SITEA;
      else if (strcmp(arg[1],"siteB") == 0) adestype[nassocdes][0] = SITEB;
      else if (strcmp(arg[1],"siteC") == 0) adestype[nassocdes][0] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[2],"siteA") == 0) adestype[nassocdes][1] = SITEA;
      else if (strcmp(arg[2],"siteB") == 0) adestype[nassocdes][1] = SITEB;
      else if (strcmp(arg[2],"siteC") == 0) adestype[nassocdes][1] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[3],"spec1") == 0) adesinput[nassocdes][0] = SPEC1;
      else if (strcmp(arg[3],"spec2") == 0) adesinput[nassocdes][0] = SPEC2;
      else if (strcmp(arg[3],"spec3") == 0) adesinput[nassocdes][0] = SPEC3;
      else if (strcmp(arg[3],"vac") == 0)
        error->all(FLERR,"rstyle=6 only allows spec1/2/3->vac");
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[4],"spec1") == 0) adesinput[nassocdes][1] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) adesinput[nassocdes][1] = SPEC2;
      else if (strcmp(arg[4],"spec3") == 0) adesinput[nassocdes][1] = SPEC3;
      else if (strcmp(arg[4],"vac") == 0)
        error->all(FLERR,"rstyle=6 only allows spec1/2/3->vac");
      else error->all(FLERR,"Illegal event command");

      adesrate[nassocdes] = atof(arg[5]);

      if (strcmp(arg[6],"spec1") == 0)
        error->all(FLERR,"rstyle=6 only allows spec1/2/3->vac");
      else if (strcmp(arg[6],"spec2") == 0)
        error->all(FLERR,"rstyle=6 only allows spec1/2/3->vac");
      else if (strcmp(arg[6],"spec3") == 0)
        error->all(FLERR,"rstyle=6 only allows spec1/2/3->vac");
      else if (strcmp(arg[6],"vac") == 0) adesoutput[nassocdes][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[7],"spec1") == 0)
        error->all(FLERR,"rstyle=6 only allows spec1/2/3->vac");
      else if (strcmp(arg[7],"spec2") == 0)
        error->all(FLERR,"rstyle=6 only allows spec1/2/3->vac");
      else if (strcmp(arg[7],"spec3") == 0)
        error->all(FLERR,"rstyle=6 only allows spec1/2/3->vac");
      else if (strcmp(arg[7],"vac") == 0) adesoutput[nassocdes][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[8],"spec1") == 0) adesdesorbate[nassocdes] = SPEC1;
      else if (strcmp(arg[8],"spec2") == 0) adesdesorbate[nassocdes] = SPEC2;
      else if (strcmp(arg[8],"spec3") == 0) adesdesorbate[nassocdes] = SPEC3;
      else if (strcmp(arg[8],"vac") == 0)
        error->all(FLERR,"rstyle=6 only allows desorption of spec1/2/3");
      else error->all(FLERR,"Illegal event command");

      nassocdes++;

    } else if (rstyle == 7) { // (reaction) associative desorption
      if (narg != 9) error->all(FLERR,"Illegal event command");
      if (strcmp(arg[1],"siteA") == 0) rxntype[nreaction][0] = SITEA;
      else if (strcmp(arg[1],"siteB") == 0) rxntype[nreaction][0] = SITEB;
      else if (strcmp(arg[1],"siteC") == 0) rxntype[nreaction][0] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[2],"siteA") == 0) rxntype[nreaction][1] = SITEA;
      else if (strcmp(arg[2],"siteB") == 0) rxntype[nreaction][1] = SITEB;
      else if (strcmp(arg[2],"siteC") == 0) rxntype[nreaction][1] = SITEC;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[3],"spec1") == 0) rxninput[nreaction][0] = SPEC1;
      else if (strcmp(arg[3],"spec2") == 0) rxninput[nreaction][0] = SPEC2;
      else if (strcmp(arg[3],"spec3") == 0) rxninput[nreaction][0] = SPEC3;
      else if (strcmp(arg[3],"vac") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3->vac");
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[4],"spec1") == 0) rxninput[nreaction][1] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) rxninput[nreaction][1] = SPEC2;
      else if (strcmp(arg[4],"spec3") == 0) rxninput[nreaction][1] = SPEC3;
      else if (strcmp(arg[4],"vac") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3->vac");
      else error->all(FLERR,"Illegal event command");

      rxnrate[nreaction] = atof(arg[5]);

      if (strcmp(arg[6],"spec1") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3->vac");
      else if (strcmp(arg[6],"spec2") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3->vac");
      else if (strcmp(arg[6],"spec3") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3->vac");
      else if (strcmp(arg[6],"vac") == 0) rxnoutput[nreaction][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[7],"spec1") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3->vac");
      else if (strcmp(arg[7],"spec2") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3->vac");
      else if (strcmp(arg[7],"spec3") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3->vac");
      else if (strcmp(arg[7],"vac") == 0) rxnoutput[nreaction][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[8],"spec1") == 0) reactionsorbate[nreaction] = SPEC1;
      else if (strcmp(arg[8],"spec2") == 0) reactionsorbate[nreaction] = SPEC2;
      else if (strcmp(arg[8],"spec3") == 0) reactionsorbate[nreaction] = SPEC3;
      else if (strcmp(arg[8],"vac") == 0)
        error->all(FLERR,"rstyle=7 only allows desorption of spec1/2/3");
      else error->all(FLERR,"Illegal event command");

      nreaction++;

    }  else error->all(FLERR,"Illegal event command");
  }

#if defined(USE_AMREX_MPMD)
  else if (strcmp(command,"amrex_update_full") == 0) {
    if (narg != 1) error->all(FLERR,"Illegal amrex_update_full command");
    update_full = atoi(arg[0]);
    if (update_full < 0 || update_full > 1) error->all(FLERR,"Illegal amrex_update_full command");
  } else if (strcmp(command,"amrex_ads_wall_dir") == 0) {
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

void AppCooxidation::grow_app()
{
  type = iarray[0];
  element = iarray[1];
  ac1 = iarray[2];
  ac2 = iarray[3];
  ac3 = iarray[4];
  dc1 = iarray[5];
  dc2 = iarray[6];
  dc3 = iarray[7];
  dac1 = iarray[8];
  dac2 = iarray[9];
  dac3 = iarray[10];
  adc1 = iarray[11];
  adc2 = iarray[12];
  adc3 = iarray[13];
  pressure1 = darray[0];
  pressure2 = darray[1];
  pressure3 = darray[2];
  temp = darray[3];
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppCooxidation::init_app()
{
  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    memory->create(firstevent,nlocal,"app:firstevent");

#if defined(USE_AMREX_MPMD)
    bordersite = new int[nlocal];
    for (int i = 0; i < nlocal; i++) bordersite[i] = 0;
#endif

    // esites must be large enough for 2 sites and their 1st neighbors

    esites = new int[2 + 2*maxneigh];

    // initializing ac1-ac3 and dc1-dc3

    for (int i = 0; i < nlocal; i++) {
      ac1[i] = 0;
      ac2[i] = 0;
      ac3[i] = 0;
      dc1[i] = 0;
      dc2[i] = 0;
      dc3[i] = 0;
      dac1[i] = 0;
      dac2[i] = 0;
      dac3[i] = 0;
      adc1[i] = 0;
      adc2[i] = 0;
      adc3[i] = 0;
    }

    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: ac1-ac3, dc1-dc3, dac1-dac3, adc1-adc3 initialized to zero\n");
      fflush(screen);
    }
  }

  // site validity

  int flag = 0;
  for (int i = 0; i < nlocal; i++) {
    if (type[i] < SITEA || type[i] > SITEC) flag = 1;
    if (element[i] < VACANCY || element[i] > SPEC3) flag = 1;
  }
  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");
}

/* ---------------------------------------------------------------------- */

void AppCooxidation::setup_app()
{
  // in ads-only update mode (amrex_update_full 0, default), event lists
  // persist across coupling cycles (the per-cycle refresh updates only
  // the gas-fed adsorption entries), so they are built from scratch only
  // on the first setup

#if defined(USE_AMREX_MPMD)
  // sites whose event lists can be invalidated by communication (sector
  // borders and their first neighbors) always take the full-rebuild path

  if (update_full == 0)
    for (int s = 0; s < nset; s++)
      for (int m = 0; m < set[s].nborder; m++) {
        int b = set[s].border[m];
        bordersite[b] = 1;
        for (int k = 0; k < numneigh[b]; k++)
          if (neighbor[b][k] < nlocal) bordersite[neighbor[b][k]] = 1;
      }

  if (update_full || nevents == 0)
#endif
  {
    for (int i = 0; i < nlocal; i++) echeck[i] = 0;

    // clear event list

    nevents = 0;
    for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
    for (int i = 0; i < maxevent; i++) events[i].next = i+1;
    freeevent = 0;
  }

  // set propensities from rates
  // propentities for adsorption reactions will be calculated in site_propensity(i)

  if (temperature == 0.0)
    error->all(FLERR,"Temperature cannot be 0.0 for app cooxidation");

  for (int m = 0; m < none; m++) {
    spropensity[m] = srate[m];
    scount[m] = 0;
  }
  for (int m = 0; m < ntwo; m++) {
    dpropensity[m] = drate[m];
    dcount[m] = 0;
  }
  for (int m = 0; m < nads; m++) {
    adscount[m] = 0;
  }
  for (int m = 0; m < ndes; m++) {
    descount[m] = 0;
  }
  for (int m = 0; m < ndissocads; m++) {
    dadscount[m] = 0;
  }
  for (int m = 0; m < nassocdes; m++) {
    adespropensity[m] = adesrate[m];
    adescount[m] = 0;
  }
  for (int m = 0; m < nreaction; m++) {
    rxnpropensity[m] = rxnrate[m];
    rxncount[m] = 0;
  }
}

/* ----------------------------------------------------------------------
   compute energy of site
------------------------------------------------------------------------- */

double AppCooxidation::site_energy(int i)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppCooxidation::site_propensity(int i)
{
  // ads-only fast path for the per-cycle refresh: the configuration
  // around i is unchanged since its event list was built (site_event
  // rebuilds affected sites via site_propensity_full, and border sites
  // always rebuild), so only the gas-fed adsorption entries need new
  // values, inserted from the per-cell table filled at fetch.
  // firstevent[i] >= 0 keeps never-built sites (first setup) on the
  // full path; cellrep != NULL means the coupling was initialized

#if defined(USE_AMREX_MPMD)
  if (update_full == 0 && cellrep != NULL && firstevent[i] >= 0 && !bordersite[i])
    return site_propensity_ads(i);
#endif

  return site_propensity_full(i);
}

/* ---------------------------------------------------------------------- */

double AppCooxidation::site_propensity_full(int i)
{
  int j,m;

  // valid single and double events are in tabulated lists
  // propensity for each event is input by user

  clear_events(i);

  double proball = 0.0;

  // single-site events

  for (m = 0; m < none; m++) {
    if (type[i] != stype[m] || element[i] != sinput[m]) continue;
    add_event(i,1,m,spropensity[m],-1);
    proball += spropensity[m];
  }

  // double-site events

  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (m = 0; m < ntwo; m++) {
      if (type[i] != dtype[m][0] || element[i] != dinput[m][0]) continue;
      if (type[j] != dtype[m][1] || element[j] != dinput[m][1]) continue;

      add_event(i,2,m,dpropensity[m],j);
      proball += dpropensity[m];
    }
  }

  // adsorption events

  double adspropensity;
  double tempratio = temp[i]/temperature;

  for (m = 0; m < nads; m++) {
    if (type[i] != adstype[m] || element[i] != adsinput[m]) continue;

    if (ads_is_rate[m]) adspropensity = adsrate[m];
    else
    {
      if (adsadsorbate[m]==SPEC1) adspropensity = adsrate[m]*pressure1[i]*pow(tempratio,BETA);
      else if (adsadsorbate[m]==SPEC2) adspropensity = adsrate[m]*pressure2[i]*pow(tempratio,BETA);
      else if (adsadsorbate[m]==SPEC3) adspropensity = adsrate[m]*pressure3[i]*pow(tempratio,BETA);
    }

    add_event(i,3,m,adspropensity,-1);
    proball += adspropensity;
  }

  // desorption events

  double despropensity;

  for (m = 0; m < ndes; m++) {
    if (type[i] != destype[m] || element[i] != desinput[m]) continue;

    despropensity = desrate[m];
    add_event(i,4,m,despropensity,-1);
    proball += despropensity;
  }

  // dissociative adsorption events

  double dadspropensity;

  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (m = 0; m < ndissocads; m++) {
      if (type[i] != dadstype[m][0] || element[i] != dadsinput[m][0]) continue;
      if (type[j] != dadstype[m][1] || element[j] != dadsinput[m][1]) continue;
      if (dads_is_rate[m]) dadspropensity = dadsrate[m];
      else {
        if (dadsadsorbate[m]==SPEC1) dadspropensity = dadsrate[m]*pressure1[i]*pow(tempratio,BETA);
        else if (dadsadsorbate[m]==SPEC2) dadspropensity = dadsrate[m]*pressure2[i]*pow(tempratio,BETA);
        else if (dadsadsorbate[m]==SPEC3) dadspropensity = dadsrate[m]*pressure3[i]*pow(tempratio,BETA);
      }
      add_event(i,5,m,dadspropensity,j);
      proball += dadspropensity;
    }
  }

  // associative desorption events

  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (m = 0; m < nassocdes; m++) {
      if (type[i] != adestype[m][0] || element[i] != adesinput[m][0]) continue;
      if (type[j] != adestype[m][1] || element[j] != adesinput[m][1]) continue;
      add_event(i,6,m,adespropensity[m],j);
      proball += adespropensity[m];
    }
  }

  // reaction events

  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (m = 0; m < nreaction; m++) {
      if (type[i] != rxntype[m][0] || element[i] != rxninput[m][0]) continue;
      if (type[j] != rxntype[m][1] || element[j] != rxninput[m][1]) continue;
      add_event(i,7,m,rxnpropensity[m],j);
      proball += rxnpropensity[m];
    }
  }

  return proball;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void AppCooxidation::site_event(int i, class RandomPark *random)
{
  int j,m,n;

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

  // perform single, double, or adsorption/desorption event

  int rstyle = events[ievent].style;
  int which = events[ievent].which;
  j = events[ievent].jpartner;

  if (rstyle == 1) {
    element[i] = soutput[which];
    scount[which]++;
  } else if (rstyle == 2) {
    element[i] = doutput[which][0];
    element[j] = doutput[which][1];
    dcount[which]++;
  } else if (rstyle == 3) {  // adsorption case
    element[i] = adsoutput[which];
    adscount[which]++;
    if (element[i] == SPEC1) ac1[i]++;
    else if (element[i] == SPEC2) ac2[i]++;
    else if (element[i] == SPEC3) ac3[i]++;
  } else if (rstyle == 4) {  // desorption case
    if (element[i] == SPEC1) dc1[i]++;
    else if (element[i] == SPEC2) dc2[i]++;
    else if (element[i] == SPEC3) dc3[i]++;
    element[i] = desoutput[which];
    descount[which]++;
  } else if (rstyle == 5) { // dissociative adsorption case
    element[i] = dadsoutput[which][0];
    element[j] = dadsoutput[which][1];
    int adsorbate = dadsadsorbate[which];
    if (adsorbate == SPEC1) dac1[i]++;
    else if (adsorbate == SPEC2) dac2[i]++;
    else if (adsorbate == SPEC3) dac3[i]++;
    dadscount[which]++;
  } else if (rstyle == 6) { // associative desorption case
    element[i] = adesoutput[which][0];
    element[j] = adesoutput[which][1];
    int desorbate = adesdesorbate[which];
    if (desorbate == SPEC1) adc1[i]++;
    else if (desorbate == SPEC2) adc2[i]++;
    else if (desorbate == SPEC3) adc3[i]++;
    adescount[which]++;
  } else if (rstyle == 7) { // (reaction) associative desorption case
    element[i] = rxnoutput[which][0];
    element[j] = rxnoutput[which][1];
    int rxnsorbate = reactionsorbate[which];
    if (rxnsorbate == SPEC1) adc1[i]++;
    else if (rxnsorbate == SPEC2) adc2[i]++;
    else if (rxnsorbate == SPEC3) adc3[i]++;
    rxncount[which]++;
    rxnsumcount++;
  }

  // compute propensity changes for participating sites and first neighs
  // ignore update of sites with isite < 0
  // use echeck[] to avoid resetting propensity of same site

  int nsites = 0;

  int isite = i2site[i];
  propensity[isite] = site_propensity_full(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  for (n = 0; n < numneigh[i]; n++) {
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity_full(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
  }

  if (rstyle == 2 || rstyle == 5 || rstyle == 6 || rstyle == 7) {
    for (n = 0; n < numneigh[j]; n++) {
      m = neighbor[j][n];
      isite = i2site[m];
      if (isite >= 0 && echeck[isite] == 0) {
        propensity[isite] = site_propensity_full(m);
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

void AppCooxidation::clear_events(int i)
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

void AppCooxidation::add_event(int i, int rstyle, int which, double propensity,
              int jpartner)
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
  events[freeevent].propensity = propensity;

  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   grow list of stored reactions for single or double
------------------------------------------------------------------------- */

void AppCooxidation::grow_reactions(int rstyle)
{
  if (rstyle == 1) {
    int n = none + 1;
    memory->grow(srate,n,"app/cooxidation:srate");
    memory->grow(spropensity,n,"app/cooxidation:spropensity");
    memory->grow(stype,n,"app/cooxidation:stype");
    memory->grow(sinput,n,"app/cooxidation:sinput");
    memory->grow(soutput,n,"app/cooxidation:soutput");
    memory->grow(scount,n,"app/cooxidation:scount");

  } else if (rstyle == 2) {
    int n = ntwo + 1;
    memory->grow(drate,n,"app/cooxidation:drate");
    memory->grow(dpropensity,n,"app/cooxidation:dpropensity");
    dtype = memory->grow(dtype,n,2,"app/cooxidation:dtype");
    dinput = memory->grow(dinput,n,2,"app/cooxidation:dinput");
    doutput = memory->grow(doutput,n,2,"app/cooxidation:doutput");
    memory->grow(dcount,n,"app/cooxidation:dcount");

  } else if (rstyle == 3) {
    int n = nads + 1;
    memory->grow(adsrate,n,"app/cooxidation:adsrate");
    memory->grow(ads_is_rate,n,"app/cooxidation:ads_is_rate");
    memory->grow(adsadsorbate,n,"app/cooxidation:adsadsorbate");
    memory->grow(adstype,n,"app/cooxidation:adstype");
    memory->grow(adsinput,n,"app/cooxidation:adsinput");
    memory->grow(adsoutput,n,"app/cooxidation:adsoutput");
    memory->grow(adscount,n,"app/cooxidation:adscount");
  } else if (rstyle == 4) {
    int n = ndes + 1;
    memory->grow(desrate,n,"app/cooxidation:desrate");
    memory->grow(destype,n,"app/cooxidation:destype");
    memory->grow(desinput,n,"app/cooxidation:desinput");
    memory->grow(desoutput,n,"app/cooxidation:desoutput");
    memory->grow(descount,n,"app/cooxidation:descount");
  } else if (rstyle == 5) {
    int n = ndissocads + 1;
    memory->grow(dadsrate,n,"app/cooxidation:dadsrate");
    memory->grow(dads_is_rate,n,"app/cooxidation:dads_is_rate");
    dadstype = memory->grow(dadstype,n,2,"app/cooxidation:dadstype");
    dadsinput = memory->grow(dadsinput,n,2,"app/cooxidation:dadsinput");
    dadsoutput = memory->grow(dadsoutput,n,2,"app/cooxidation:dadsoutput");
    memory->grow(dadscount,n,"app/cooxidation:dadscount");
    memory->grow(dadsadsorbate,n,"app/cooxidation:dadsadsorbate");
  } else if (rstyle == 6) {
    int n = nassocdes + 1;
    memory->grow(adesrate,n,"app/cooxidation:adesrate");
    memory->grow(adespropensity,n,"app/cooxidation:adespropensity");
    adestype = memory->grow(adestype,n,2,"app/cooxidation:adestype");
    adesinput = memory->grow(adesinput,n,2,"app/cooxidation:adesinput");
    adesoutput = memory->grow(adesoutput,n,2,"app/cooxidation:adesoutput");
    memory->grow(adescount,n,"app/cooxidation:adescount");
    memory->grow(adesdesorbate,n,"app/cooxidation:adesdesorbate");
  } else if (rstyle == 7) {
    int n = nreaction + 1;
    memory->grow(rxnrate,n,"app/cooxidation:rxnrate");
    memory->grow(rxnpropensity,n,"app/cooxidation:rxnpropensity");
    rxntype = memory->grow(rxntype,n,2,"app/cooxidation:rxntype");
    rxninput = memory->grow(rxninput,n,2,"app/cooxidation:rxninput");
    rxnoutput = memory->grow(rxnoutput,n,2,"app/cooxidation:rxnoutput");
    memory->grow(rxncount,n,"app/cooxidation:rxncount");
    memory->grow(reactionsorbate,n,"app/cooxidation:reactionsorbate");
  }
}

#if defined(USE_AMREX_MPMD)

void AppCooxidation::amrex_init_agg ()
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

    // 4. cellrep, adsprop_cell, dadsprop_cell
    //    (per-cell adsorption rate table; filled at each amrex_fetch_agg,
    //     zero until the first fetch)

    memory->create(cellrep,nlocalFHDcell,"app/cooxidation:cellrep");
    for (int c = 0; c < nlocalFHDcell; c++) cellrep[c] = -1;
    for (int i = 0; i < nlocal; i++)
        if (cellrep[localFHDcell[i]] < 0) cellrep[localFHDcell[i]] = i;

    memory->create(adsprop_cell,nads*nlocalFHDcell,
                   "app/cooxidation:adsprop_cell");
    for (int k = 0; k < nads*nlocalFHDcell; k++) adsprop_cell[k] = 0.;
    memory->create(dadsprop_cell,ndissocads*nlocalFHDcell,
                   "app/cooxidation:dadsprop_cell");
    for (int k = 0; k < ndissocads*nlocalFHDcell; k++) dadsprop_cell[k] = 0.;
}

void AppCooxidation::amrex_push_agg(int narg, char **arg)
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
        } else if (std::strcmp(arg[k],"dac1") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                intval[localFHDcell[i]] += dac1[i];
                dac1[i] = 0;
            }
            amrex_send_intval();
        } else if (std::strcmp(arg[k],"dac2") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                intval[localFHDcell[i]] += dac2[i];
                dac2[i] = 0;
            }
            amrex_send_intval();
        } else if (std::strcmp(arg[k],"dac3") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                intval[localFHDcell[i]] += dac3[i];
                dac3[i] = 0;
            }
            amrex_send_intval();
        } else if (std::strcmp(arg[k],"adc1") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                intval[localFHDcell[i]] += adc1[i];
                adc1[i] = 0;
            }
            amrex_send_intval();
        } else if (std::strcmp(arg[k],"adc2") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                intval[localFHDcell[i]] += adc2[i];
                adc2[i] = 0;
            }
            amrex_send_intval();
        } else if (std::strcmp(arg[k],"adc3") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                intval[localFHDcell[i]] += adc3[i];
                adc3[i] = 0;
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
        } else if (std::strcmp(arg[k],"occ3") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                int is_occ = (element[i]==3) ? 1 : 0;
                intval[localFHDcell[i]] += is_occ;
            }
            amrex_send_intval();
        } else if (std::strcmp(arg[k],"one") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) intval[localFHDcell[i]]++;
            amrex_send_intval();
        } else {
            error->all(FLERR,"Illegal amrex_push_agg command");
        }

        if (domain->me == 0 && screen) {
            std::fprintf(screen,"** DEBUG: %s pushed\n",arg[k]);
            std::fflush(screen);
        }
    }
}

void AppCooxidation::amrex_fetch_agg(int narg, char **arg)
{
    int timestamp = atoi(arg[0]);

    if (domain->me == 0 && screen) {
        std::fprintf(screen,"** DEBUG: amrex_fetch_agg at timestamp %d\n",timestamp);
        std::fflush(screen);
    }

    if (amrex_fhd_lattice_size_x <= 0. || amrex_fhd_lattice_size_y <= 0.)
        error->all(FLERR,"amrex_fhd_lattice_size must be set as two positive numbers");

    for (int k=1;k<narg;k++)
    {
        if (std::strcmp(arg[k],"pressure1") == 0) {
            amrex_recv_dblval();
            // distribute info to each KMC site
            for (int i=0;i<nlocal;i++)
                pressure1[i] = dblval[localFHDcell[i]];
        }
        else if (std::strcmp(arg[k],"pressure2") == 0) {
            amrex_recv_dblval();
            // distribute info to each KMC site
            for (int i=0;i<nlocal;i++)
                pressure2[i] = dblval[localFHDcell[i]];
        }
        else if (std::strcmp(arg[k],"pressure3") == 0) {
            amrex_recv_dblval();
            // distribute info to each KMC site
            for (int i=0;i<nlocal;i++)
                pressure3[i] = dblval[localFHDcell[i]];
        }
        else if (std::strcmp(arg[k],"temp") == 0) {
            amrex_recv_dblval();
            // distribute info to each KMC site
            for (int i=0;i<nlocal;i++)
                temp[i] = dblval[localFHDcell[i]];
        }
        else {
            error->all(FLERR,"Illegal amrex_fetch_agg command");
        }

        if (domain->me == 0 && screen) {
            std::fprintf(screen,"** DEBUG: %s fetched\n",arg[k]);
            std::fflush(screen);
        }
    }

    // the fetched fields are new: recompute the per-cell adsorption rates
    // (pressure and temperature are constant within an FHD cell, so each
    // (event, cell) value is computed once per cycle and the per-site
    // refresh only inserts it; the table is allocated in amrex_init_agg,
    // section 4)

    if (update_full == 0) {

        for (int c = 0; c < nlocalFHDcell; c++) {
            int rep = cellrep[c];
            double tempratio = temp[rep]/temperature;
            for (int m = 0; m < nads; m++) {
                double p;
                if (ads_is_rate[m]) adsprop_cell[m*nlocalFHDcell+c] = adsrate[m];
                else {
                    double p;
                    if (adsadsorbate[m] == SPEC1) p = pressure1[rep];
                    else if (adsadsorbate[m] == SPEC2) p = pressure2[rep];
                    else p = pressure3[rep];
                    adsprop_cell[m*nlocalFHDcell+c] = adsrate[m]*p*pow(tempratio,BETA);
                }
            }
            for (int m = 0; m < ndissocads; m++) {
                double p;
                if (dads_is_rate[m]) dadsprop_cell[m*nlocalFHDcell+c] = dadsrate[m];
                else {
                    double p;
                    if (dadsadsorbate[m] == SPEC1) p = pressure1[rep];
                    else if (dadsadsorbate[m] == SPEC2) p = pressure2[rep];
                    else p = pressure3[rep];
                    dadsprop_cell[m*nlocalFHDcell+c] = dadsrate[m]*p*pow(tempratio,BETA);
                }
            }
        }

    }
}

void AppCooxidation::amrex_send_intval()
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

void AppCooxidation::amrex_recv_dblval()
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

/* ----------------------------------------------------------------------
   ads-only update: walk the existing event list of site i and update
   only the gas-fed entries (adsorption, style 3; dissociative
   adsorption, style 5) by lookup in the per-cell rate table; all other
   entries depend only on the configuration, which site_event keeps
   current
------------------------------------------------------------------------- */

double AppCooxidation::site_propensity_ads(int i)
{
  int c = localFHDcell[i];
  double proball = 0.0;

  for (int ev = firstevent[i]; ev >= 0; ev = events[ev].next) {
    if (events[ev].style == 3)
      events[ev].propensity = adsprop_cell[events[ev].which*nlocalFHDcell+c];
    else if (events[ev].style == 5)
      events[ev].propensity = dadsprop_cell[events[ev].which*nlocalFHDcell+c];
    proball += events[ev].propensity;
  }

  return proball;
}

#endif
