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

enum{NOOP,SITEA,SITEB,SITEC};
enum{VACANCY,SPEC1,SPEC2,SPEC3,SPEC4,SPEC5}; // removed ZERO and moved VACANCY to first item // same as DiagSurfchemtest

#define DELTAEVENT 100000

/* ---------------------------------------------------------------------- */

AppSurfchemtest::AppSurfchemtest(SPPARKS *spk, int narg, char **arg) :
  AppLattice(spk,narg,arg)
{
  ninteger = 22;  // type: site type
                  // element: site state
                  // ac1/ac2/ac3/ac4/ac5: adsorption count
                  // dc1/dc2/dc3/dc4/dc5: desorption count
                  // dac1/dac2/dac3/dac4/dac5: dissociative adsorption count
                  // adc1/adc2/adc3/adc4/adc5: associative desorption count
  ndouble = 7;    // density1/density2/density3/density4/density5: number density of the contacting FHD cell
                  // temp: temperature of the contacting FHD cell
                  // Vz: normal velocity of the contacting FHD cell
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

  ads_is_rate = false;
  dads_is_rate = false;

  // reaction lists
  none = ntwo = nthree = nads = ndes = ndissocads = nassocdes = nreaction = rxnsumcount = 0;
  srate = drate = trate = adsrate = ads_beta = dads_beta = desrate = dadsrate = adesrate = rxnrate = NULL; // beta implementation
  spropensity = dpropensity = tpropensity = adespropensity = rxnpropensity = NULL;   // no adspropensity/dadspropensity/despropensity
  stype = sinput = soutput = NULL;
  dtype = dinput = doutput = NULL;
  ttype = tinput = toutput = NULL;
  adstype = adsinput = adsoutput = NULL;
  destype = desinput = desoutput = NULL;
  dadstype = dadsinput = dadsoutput = NULL;
  adestype = adesinput = adesoutput = NULL;
  rxntype = rxninput = rxnoutput = NULL;
  scount = dcount = tcount = adscount = descount = dadscount = adescount = rxncount = NULL;
  dadsadsorbate = adesdesorbate = reactionsorbate = NULL;

  neighboring_diff = neighboring_des = neighboring_ades = neighboring_rxn = NULL;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 6; k++) {
        for (int l = 0; l < 6; l++) {
          V_neighbor[i][j][k][l] = 0.0;
        }
      }
    }
  }

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
  memory->destroy(ads_beta);
  memory->destroy(dads_beta);
  memory->destroy(desrate);
  memory->destroy(dadsrate);
  memory->destroy(adesrate);
  memory->destroy(rxnrate);
  memory->destroy(spropensity);
  memory->destroy(dpropensity);
  memory->destroy(tpropensity);
  memory->destroy(adespropensity);
  memory->destroy(rxnpropensity);
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
  memory->destroy(tcount);
  memory->destroy(adscount);
  memory->destroy(descount);
  memory->destroy(dadscount);
  memory->destroy(adescount);
  memory->destroy(rxncount);
  memory->destroy(dadsadsorbate);
  memory->destroy(adesdesorbate);
  memory->destroy(reactionsorbate);

  memory->destroy(neighboring_diff);
  memory->destroy(neighboring_des);
  memory->destroy(neighboring_ades);
  memory->destroy(neighboring_rxn);

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

void AppSurfchemtest::input_app(char *command, int narg, char **arg)
{
  if (strcmp(command,"neighbor") == 0) {
    if (narg != 5) error->all(FLERR,"Illegal event command");
    int i,j,k,l;
    if (strcmp(arg[0],"siteA") == 0) i = SITEA;
    else if (strcmp(arg[0],"siteB") == 0) i = SITEB;
    else if (strcmp(arg[0],"siteC") == 0) i = SITEC;
    else error->all(FLERR,"Illegal event command");
    if (strcmp(arg[1],"siteA") == 0) j = SITEA;
    else if (strcmp(arg[1],"siteB") == 0) j = SITEB;
    else if (strcmp(arg[1],"siteC") == 0) j = SITEC;
    else error->all(FLERR,"Illegal event command");
    if (strcmp(arg[2],"spec1") == 0) k = SPEC1;
    else if (strcmp(arg[2],"spec2") == 0) k = SPEC2;
    else if (strcmp(arg[2],"sepc3") == 0) k = SPEC3;
    else if (strcmp(arg[2],"sepc4") == 0) k = SPEC4;
    else if (strcmp(arg[2],"sepc5") == 0) k = SPEC5;
    else if (strcmp(arg[2],"vac") == 0) k = VACANCY;
    else error->all(FLERR,"Illegal event command");
    if (strcmp(arg[3],"spec1") == 0) l = SPEC1;
    else if (strcmp(arg[3],"spec2") == 0) l = SPEC2;
    else if (strcmp(arg[3],"sepc3") == 0) l = SPEC3;
    else if (strcmp(arg[3],"sepc4") == 0) l = SPEC4;
    else if (strcmp(arg[3],"sepc5") == 0) l = SPEC5;
    else if (strcmp(arg[3],"vac") == 0) l = VACANCY;
    else error->all(FLERR,"Illegal event command");
    V_neighbor[i][j][k][l] = atof(arg[4]);
  }

  else if (strcmp(command,"event") == 0) {
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
      if (narg < 8 || narg > 9) error->all(FLERR,"Illegal event command");
      if (narg == 9) {
        if (strcmp(arg[8],"neighbor") == 0) neighboring_diff[ntwo] = true;
        else error->all(FLERR,"Illegal event command");
      }
      else {
        neighboring_diff[ntwo] = false;
      }
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
      if (narg < 5 || narg > 7) error->all(FLERR,"Illegal event command");
      if (narg == 6) {
        if (strcmp(arg[5],"rate") == 0) ads_is_rate = true;
        else error->all(FLERR,"Illegal event command");
      }
      ads_beta[nads] = 0.;
      if (narg == 7) { // beta implementation
        if (strcmp(arg[5],"beta") == 0) ads_beta[nads] = atof(arg[6]);
        else error->all(FLERR,"Illegal event command");
      }
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
      if (narg < 5 || narg > 6) error->all(FLERR,"Illegal event command");
      if (narg == 6) {
        if (strcmp(arg[5],"neighbor") == 0) neighboring_des[ndes] = true;
        else error->all(FLERR,"Illegal event command");
      }
      else {
        neighboring_des[ndes] = false;
      }
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

    } else if (rstyle == 6) { // disassociative adsorption
      if (narg < 9 || narg > 11) error->all(FLERR,"Illegal event command");
      if (narg == 10) {
        if (strcmp(arg[9],"rate") == 0) dads_is_rate = true;
        else error->all(FLERR,"Illegal event command");
      }
      dads_beta[ndissocads] = 0.;
      if (narg == 11) { // beta implementation
        if (strcmp(arg[9],"beta") == 0) dads_beta[ndissocads] = atof(arg[10]);
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
        error->all(FLERR,"rstyle=6 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[3],"spec2") == 0)
        error->all(FLERR,"rstyle=6 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[3],"spec3") == 0)
        error->all(FLERR,"rstyle=6 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[3],"spec4") == 0)
        error->all(FLERR,"rstyle=6 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[3],"spec5") == 0)
        error->all(FLERR,"rstyle=6 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[3],"vac") == 0) dadsinput[ndissocads][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[4],"spec1") == 0)
        error->all(FLERR,"rstyle=6 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[4],"spec2") == 0)
        error->all(FLERR,"rstyle=6 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[4],"spec3") == 0)
        error->all(FLERR,"rstyle=6 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[4],"spec4") == 0)
        error->all(FLERR,"rstyle=6 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[4],"spec5") == 0)
        error->all(FLERR,"rstyle=6 only allows vac->spec1/2/3/4/5");
      else if (strcmp(arg[4],"vac") == 0) dadsinput[ndissocads][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");

      dadsrate[ndissocads] = atof(arg[5]);

      if (strcmp(arg[6],"spec1") == 0) dadsoutput[ndissocads][0] = SPEC1;
      else if (strcmp(arg[6],"spec2") == 0) dadsoutput[ndissocads][0] = SPEC2;
      else if (strcmp(arg[6],"spec3") == 0) dadsoutput[ndissocads][0] = SPEC3;
      else if (strcmp(arg[6],"spec4") == 0) dadsoutput[ndissocads][0] = SPEC4;
      else if (strcmp(arg[6],"spec5") == 0) dadsoutput[ndissocads][0] = SPEC5;
      else if (strcmp(arg[6],"vac") == 0)
        error->all(FLERR,"rstyle=6 only allows vac->spec1/2/3/4/5");
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[7],"spec1") == 0) dadsoutput[ndissocads][1] = SPEC1;
      else if (strcmp(arg[7],"spec2") == 0) dadsoutput[ndissocads][1] = SPEC2;
      else if (strcmp(arg[7],"spec3") == 0) dadsoutput[ndissocads][1] = SPEC3;
      else if (strcmp(arg[7],"spec4") == 0) dadsoutput[ndissocads][1] = SPEC4;
      else if (strcmp(arg[7],"spec5") == 0) dadsoutput[ndissocads][1] = SPEC5;
      else if (strcmp(arg[6],"vac") == 0)
        error->all(FLERR,"rstyle=6 only allows vac->spec1/2/3/4/5");
      if (strcmp(arg[8],"spec1") == 0) dadsadsorbate[ndissocads] = SPEC1;
      else if (strcmp(arg[8],"spec2") == 0) dadsadsorbate[ndissocads] = SPEC2;
      else if (strcmp(arg[8],"spec3") == 0) dadsadsorbate[ndissocads] = SPEC3;
      else if (strcmp(arg[8],"spec4") == 0) dadsadsorbate[ndissocads] = SPEC4;
      else if (strcmp(arg[8],"spec5") == 0) dadsadsorbate[ndissocads] = SPEC5;
      else if (strcmp(arg[8],"vac") == 0)
        error->all(FLERR,"rstyle=6 only allows adsorption of spec1/2/3/4/5");

      else error->all(FLERR,"Illegal event command");

      ndissocads++;

    } else if (rstyle == 7) { // associative desorption
      if (narg < 9 || narg > 10) error->all(FLERR,"Illegal event command");
      if (narg == 10) {
        if (strcmp(arg[9],"neighbor") == 0) neighboring_ades[nassocdes] = true;
        else error->all(FLERR,"Illegal event command");
      }
      else {
        neighboring_ades[nassocdes] = false;
      }
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
      else if (strcmp(arg[3],"spec4") == 0) adesinput[nassocdes][0] = SPEC4;
      else if (strcmp(arg[3],"spec5") == 0) adesinput[nassocdes][0] = SPEC5;
      else if (strcmp(arg[3],"vac") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3/4/5->vac");
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[4],"spec1") == 0) adesinput[nassocdes][1] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) adesinput[nassocdes][1] = SPEC2;
      else if (strcmp(arg[4],"spec3") == 0) adesinput[nassocdes][1] = SPEC3;
      else if (strcmp(arg[4],"spec4") == 0) adesinput[nassocdes][1] = SPEC4;
      else if (strcmp(arg[4],"spec5") == 0) adesinput[nassocdes][1] = SPEC5;
      else if (strcmp(arg[4],"vac") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3/4/5->vac");
      else error->all(FLERR,"Illegal event command");

      adesrate[nassocdes] = atof(arg[5]);

      if (strcmp(arg[6],"spec1") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[6],"spec2") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[6],"spec3") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[6],"spec4") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[6],"spec5") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[6],"vac") == 0) adesoutput[nassocdes][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[7],"spec1") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[7],"spec2") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[7],"spec3") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[7],"spec4") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[7],"spec5") == 0)
        error->all(FLERR,"rstyle=7 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[7],"vac") == 0) adesoutput[nassocdes][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[8],"spec1") == 0) adesdesorbate[nassocdes] = SPEC1;
      else if (strcmp(arg[8],"spec2") == 0) adesdesorbate[nassocdes] = SPEC2;
      else if (strcmp(arg[8],"spec3") == 0) adesdesorbate[nassocdes] = SPEC3;
      else if (strcmp(arg[8],"spec4") == 0) adesdesorbate[nassocdes] = SPEC4;
      else if (strcmp(arg[8],"spec5") == 0) adesdesorbate[nassocdes] = SPEC5;
      else if (strcmp(arg[8],"vac") == 0)
        error->all(FLERR,"rstyle=7 only allows desorption of spec1/2/3/4/5");
      else error->all(FLERR,"Illegal event command");

      nassocdes++;

    } else if (rstyle == 8) { // (reaction) associative desorption
      if (narg < 9 || narg > 10) error->all(FLERR,"Illegal event command");
      if (narg == 10) {
        if (strcmp(arg[9],"neighbor") == 0) neighboring_rxn[nreaction] = true;
        else error->all(FLERR,"Illegal event command");
      }
      else {
        neighboring_rxn[nreaction] = false;
      }
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
      else if (strcmp(arg[3],"spec4") == 0) rxninput[nreaction][0] = SPEC4;
      else if (strcmp(arg[3],"spec5") == 0) rxninput[nreaction][0] = SPEC5;
      else if (strcmp(arg[3],"vac") == 0)
        error->all(FLERR,"rstyle=8 only allows spec1/2/3/4/5->vac");
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[4],"spec1") == 0) rxninput[nreaction][1] = SPEC1;
      else if (strcmp(arg[4],"spec2") == 0) rxninput[nreaction][1] = SPEC2;
      else if (strcmp(arg[4],"spec3") == 0) rxninput[nreaction][1] = SPEC3;
      else if (strcmp(arg[4],"spec4") == 0) rxninput[nreaction][1] = SPEC4;
      else if (strcmp(arg[4],"spec5") == 0) rxninput[nreaction][1] = SPEC5;
      else if (strcmp(arg[4],"vac") == 0)
        error->all(FLERR,"rstyle=8 only allows spec1/2/3/4/5->vac");
      else error->all(FLERR,"Illegal event command");

      rxnrate[nreaction] = atof(arg[5]);

      if (strcmp(arg[6],"spec1") == 0)
        error->all(FLERR,"rstyle=8 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[6],"spec2") == 0)
        error->all(FLERR,"rstyle=8 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[6],"spec3") == 0)
        error->all(FLERR,"rstyle=8 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[6],"spec4") == 0)
        error->all(FLERR,"rstyle=8 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[6],"spec5") == 0)
        error->all(FLERR,"rstyle=8 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[6],"vac") == 0) rxnoutput[nreaction][0] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[7],"spec1") == 0)
        error->all(FLERR,"rstyle=8 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[7],"spec2") == 0)
        error->all(FLERR,"rstyle=8 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[7],"spec3") == 0)
        error->all(FLERR,"rstyle=8 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[7],"spec4") == 0)
        error->all(FLERR,"rstyle=8 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[7],"spec5") == 0)
        error->all(FLERR,"rstyle=8 only allows spec1/2/3/4/5->vac");
      else if (strcmp(arg[7],"vac") == 0) rxnoutput[nreaction][1] = VACANCY;
      else error->all(FLERR,"Illegal event command");
      if (strcmp(arg[8],"spec1") == 0) reactionsorbate[nreaction] = SPEC1;
      else if (strcmp(arg[8],"spec2") == 0) reactionsorbate[nreaction] = SPEC2;
      else if (strcmp(arg[8],"spec3") == 0) reactionsorbate[nreaction] = SPEC3;
      else if (strcmp(arg[8],"spec4") == 0) reactionsorbate[nreaction] = SPEC4;
      else if (strcmp(arg[8],"spec5") == 0) reactionsorbate[nreaction] = SPEC5;
      else if (strcmp(arg[8],"vac") == 0)
        error->all(FLERR,"rstyle=8 only allows desorption of spec1/2/3/4/5");
      else error->all(FLERR,"Illegal event command");

      nreaction++;

    }  else error->all(FLERR,"Illegal event command");
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
  dac1 = iarray[12];
  dac2 = iarray[13];
  dac3 = iarray[14];
  dac4 = iarray[15];
  dac5 = iarray[16];
  adc1 = iarray[17];
  adc2 = iarray[18];
  adc3 = iarray[19];
  adc4 = iarray[20];
  adc5 = iarray[21];
  density1 = darray[0];
  density2 = darray[1];
  density3 = darray[2];
  density4 = darray[3];
  density5 = darray[4];
  temp = darray[5];
  Vz = darray[6];
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
      dac1[i] = 0;
      dac2[i] = 0;
      dac3[i] = 0;
      dac4[i] = 0;
      dac5[i] = 0;
      adc1[i] = 0;
      adc2[i] = 0;
      adc3[i] = 0;
      adc4[i] = 0;
      adc5[i] = 0;
    }

    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: ac1-ac5, dc1-dc5, dac1-dac5, adc1-adc5 initialized to zero\n");
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

  reaction_summary_log();
}

void AppSurfchemtest::reaction_summary_log()
{
  if (domain->me == 0) {
    fprintf(logfile,"** reaction summary **\n");

    fprintf(logfile,"- nads = %d\n",nads);
    for (int i=0;i<nads;i++) {
      char sitechar = 'A'+adstype[i]-1;
      fprintf(logfile,"ads%d: [%c] %d -> %d, %e",i,sitechar,adsinput[i],adsoutput[i],adsrate[i]);
      if (ads_is_rate) fprintf(logfile," (rate), ");
      else fprintf(logfile," (rate const), ");
      fprintf(logfile,"beta= %f\n",ads_beta[i]);
    }

    fprintf(logfile,"- ndes = %d\n",ndes);
    for (int i=0;i<ndes;i++) {
      char sitechar = 'A'+destype[i]-1;
      fprintf(logfile,"des%d: [%c] %d -> %d, %e\n",i,sitechar,desinput[i],desoutput[i],desrate[i]);
    }

    fprintf(logfile,"- none = %d\n",none);
    fprintf(logfile,"- ntwo = %d\n",ntwo);
    fprintf(logfile,"- nthree = %d\n",nthree);
    fprintf(logfile,"- ndissocads = %d\n",ndissocads);
    fprintf(logfile,"- nassocdes = %d\n",nassocdes);
    fprintf(logfile,"- nreaction = %d\n",nreaction);
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
  int j,k,l,m;

  // valid single, double, triple events are in tabulated lists
  // propensity for each event is input by user

  clear_events(i);

  double proball = 0.0;
  double V_int;
  double kB = 0.00008618; // eV/K

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

      if (neighboring_diff[m]) { // diffusion
        V_int = 0.0;
        for (int x = 0; x < numneigh[i]; x++) {
          k = neighbor[i][x];
          V_int += -2*V_neighbor[type[i]][type[k]][element[i]][element[k]];
        }
        for (int y = 0; y < numneigh[j]; y++) {
          l = neighbor[j][y];
          V_int += 2*V_neighbor[type[j]][type[l]][element[i]][element[l]];
        }
        if (V_int > 0) {
          add_event(i,2,m,dpropensity[m]*exp(-V_int/(kB*temperature)),j,-1);
          proball += dpropensity[m]*exp(-V_int/(kB*temperature));
        }
        else {
          add_event(i,2,m,dpropensity[m],j,-1);
          proball += dpropensity[m];
        }
      }
      else {
        add_event(i,2,m,dpropensity[m],j,-1);
        proball += dpropensity[m];
      }
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
  double tempratio = temp[i]/temperature;

  for (m = 0; m < nads; m++) {
    if (type[i] != adstype[m] || element[i] != adsinput[m]) continue;

    if (ads_is_rate) adspropensity = adsrate[m];
    else
    {
      if (adsoutput[m]==SPEC1) adspropensity = adsrate[m]*density1[i]*pow(tempratio,ads_beta[m]); // beta
      else if (adsoutput[m]==SPEC2) adspropensity = adsrate[m]*density2[i]*pow(tempratio,ads_beta[m]);
      else if (adsoutput[m]==SPEC3) adspropensity = adsrate[m]*density3[i]*pow(tempratio,ads_beta[m]);
      else if (adsoutput[m]==SPEC4) adspropensity = adsrate[m]*density4[i]*pow(tempratio,ads_beta[m]);
      else if (adsoutput[m]==SPEC5) adspropensity = adsrate[m]*density5[i]*pow(tempratio,ads_beta[m]);
    }

    add_event(i,4,m,adspropensity,-1,-1);
    proball += adspropensity;
  }

  // desorption events

  double despropensity;

  for (m = 0; m < ndes; m++) {
    if (type[i] != destype[m] || element[i] != desinput[m]) continue;

    despropensity = desrate[m];

    if (neighboring_des[m]) {
      V_int = 0.0;
      for (int k = 0; k < numneigh[i]; k++) {
        j = neighbor[i][k];
        V_int += 2*V_neighbor[type[i]][type[j]][element[i]][element[j]];
      }
      despropensity = desrate[m]*exp(V_int/(kB*temperature));
    }
    add_event(i,5,m,despropensity,-1,-1);
    proball += despropensity;
  }

  // dissociative adsorption events

  double dadspropensity;

  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (m = 0; m < ndissocads; m++) {
      if (type[i] != dadstype[m][0] || element[i] != dadsinput[m][0]) continue;
      if (type[j] != dadstype[m][1] || element[j] != dadsinput[m][1]) continue;
      if (dads_is_rate) dadspropensity = dadsrate[m];
      else {
        if (dadsadsorbate[m]==SPEC1) dadspropensity = dadsrate[m]*density1[i]*pow(tempratio,dads_beta[m]);
        else if (dadsadsorbate[m]==SPEC2) dadspropensity = dadsrate[m]*density2[i]*pow(tempratio,dads_beta[m]);
        else if (dadsadsorbate[m]==SPEC3) dadspropensity = dadsrate[m]*density3[i]*pow(tempratio,dads_beta[m]);
        else if (dadsadsorbate[m]==SPEC4) dadspropensity = dadsrate[m]*density4[i]*pow(tempratio,dads_beta[m]);
        else if (dadsadsorbate[m]==SPEC5) dadspropensity = dadsrate[m]*density5[i]*pow(tempratio,dads_beta[m]);
      }
      add_event(i,6,m,dadspropensity,j,-1);
      proball += dadspropensity;
    }
  }

  // associative desorption events

  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (m = 0; m < nassocdes; m++) {
      if (type[i] != adestype[m][0] || element[i] != adesinput[m][0]) continue;
      if (type[j] != adestype[m][1] || element[j] != adesinput[m][1]) continue;
      if (neighboring_ades[m]) {
        V_int = 0.0;
        for (int x = 0; x < numneigh[i]; x++) {
          k = neighbor[i][x];
          V_int += V_neighbor[type[i]][type[k]][element[i]][element[k]];
          V_int += V_neighbor[type[k]][type[i]][element[k]][element[i]];
        }
        for (int y = 0; y < numneigh[j]; y++) {
          l = neighbor[j][y];
          V_int += V_neighbor[type[j]][type[l]][element[j]][element[l]];
          V_int += V_neighbor[type[l]][type[j]][element[l]][element[j]];
        }
        V_int -= 2*V_neighbor[type[i]][type[j]][element[i]][element[j]];
        V_int -= 2*V_neighbor[type[j]][type[i]][element[j]][element[i]];
        add_event(i,7,m,adespropensity[m]*exp(V_int/(kB*temperature)),j,-1);
        proball += adespropensity[m]*exp(V_int/(kB*temperature));
      }
      else {
        add_event(i,7,m,adespropensity[m],j,-1);
        proball += adespropensity[m];
      }
    }
  }

  // reaction events

  for (int jj = 0; jj < numneigh[i]; jj++) {
    j = neighbor[i][jj];
    for (m = 0; m < nreaction; m++) {
      if (type[i] != rxntype[m][0] || element[i] != rxninput[m][0]) continue;
      if (type[j] != rxntype[m][1] || element[j] != rxninput[m][1]) continue;
      if (neighboring_rxn[m]) {
        V_int = 0.0;
        for (int x = 0; x < numneigh[i]; x++) {
          k = neighbor[i][x];
          V_int += V_neighbor[type[i]][type[k]][element[i]][element[k]];
          V_int += V_neighbor[type[k]][type[i]][element[k]][element[i]];
        }
        for (int y = 0; y < numneigh[j]; y++) {
          l = neighbor[j][y];
          V_int += V_neighbor[type[j]][type[l]][element[j]][element[l]];
          V_int += V_neighbor[type[l]][type[j]][element[l]][element[j]];
        }
        V_int -= 2*V_neighbor[type[i]][type[j]][element[i]][element[j]];
        V_int -= 2*V_neighbor[type[j]][type[i]][element[j]][element[i]];
        add_event(i,8,m,rxnpropensity[m]*exp(V_int/(kB*temperature)),j,-1);
        proball += rxnpropensity[m]*exp(V_int/(kB*temperature));
      }
      else {
        add_event(i,8,m,rxnpropensity[m],j,-1);
        proball += rxnpropensity[m];
      }
    }
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
  } else if (rstyle == 6) { // dissociative adsorption case
    element[i] = dadsoutput[which][0];
    element[j] = dadsoutput[which][1];
    int adsorbate = dadsadsorbate[which];
    if (adsorbate == SPEC1) dac1[i]++;
    else if (adsorbate == SPEC2) dac2[i]++;
    else if (adsorbate == SPEC3) dac3[i]++;
    else if (adsorbate == SPEC4) dac4[i]++;
    else if (adsorbate == SPEC5) dac5[i]++;
    dadscount[which]++;
  } else if (rstyle == 7) { // associative desorption case
    element[i] = adesoutput[which][0];
    element[j] = adesoutput[which][1];
    int desorbate = adesdesorbate[which];
    if (desorbate == SPEC1) adc1[i]++;
    else if (desorbate == SPEC2) adc2[i]++;
    else if (desorbate == SPEC3) adc3[i]++;
    else if (desorbate == SPEC4) adc4[i]++;
    else if (desorbate == SPEC5) adc5[i]++;
    adescount[which]++;
  } else if (rstyle == 8) { // (reaction) associative desorption case
    element[i] = rxnoutput[which][0];
    element[j] = rxnoutput[which][1];
    int rxnsorbate = reactionsorbate[which];
    if (rxnsorbate == SPEC1) adc1[i]++;
    else if (rxnsorbate == SPEC2) adc2[i]++;
    else if (rxnsorbate == SPEC3) adc3[i]++;
    else if (rxnsorbate == SPEC4) adc4[i]++;
    else if (rxnsorbate == SPEC5) adc5[i]++;
    rxncount[which]++;
    rxnsumcount++;
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

  if (rstyle == 2 || rstyle == 3 || rstyle == 6 || rstyle == 7 || rstyle == 8) {
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
    memory->grow(neighboring_diff,n,"app/surfchemtest:neighboring_diff");

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
    memory->grow(ads_beta,n,"app/surfchemtest:ads_beta");
    memory->grow(adstype,n,"app/surfchemtest:adstype");
    memory->grow(adsinput,n,"app/surfchemtest:adsinput");
    memory->grow(adsoutput,n,"app/surfchemtest:adsoutput");
    memory->grow(adscount,n,"app/surfchemtest:adscount");
  } else if (rstyle == 5) {
    int n = ndes + 1;
    memory->grow(desrate,n,"app/surfchemtest:desrate");
    memory->grow(destype,n,"app/surfchemtest:destype");
    memory->grow(desinput,n,"app/surfchemtest:desinput");
    memory->grow(desoutput,n,"app/surfchemtest:desoutput");
    memory->grow(descount,n,"app/surfchemtest:descount");
    memory->grow(neighboring_des,n,"app/surfchemtest:neighboring_des");
  } else if (rstyle == 6) {
    int n = ndissocads + 1;
    memory->grow(dadsrate,n,"app/surfchemtest:dadsrate");
    memory->grow(dads_beta,n,"app/surfchemtest:dads_beta");
    dadstype = memory->grow(dadstype,n,2,"app/surfchemtest:dadstype");
    dadsinput = memory->grow(dadsinput,n,2,"app/surfchemtest:dadsinput");
    dadsoutput = memory->grow(dadsoutput,n,2,"app/surfchemtest:dadsoutput");
    memory->grow(dadscount,n,"app/surfchemtest:dadscount");
    memory->grow(dadsadsorbate,n,"app/surfchemtest:dadsadsorbate");
  } else if (rstyle == 7) {
    int n = nassocdes + 1;
    memory->grow(adesrate,n,"app/surfchemtest:adesrate");
    memory->grow(adespropensity,n,"app/surfchemtest:adespropensity");
    adestype = memory->grow(adestype,n,2,"app/surfchemtest:adestype");
    adesinput = memory->grow(adesinput,n,2,"app/surfchemtest:adesinput");
    adesoutput = memory->grow(adesoutput,n,2,"app/surfchemtest:adesoutput");
    memory->grow(adescount,n,"app/surfchemtest:adescount");
    memory->grow(adesdesorbate,n,"app/surfchemtest:adesdesorbate");
    memory->grow(neighboring_ades,n,"app/surfchemtest:neighboring_ades");
  } else if (rstyle == 8) {
    int n = nreaction + 1;
    memory->grow(rxnrate,n,"app/surfchemtest:adesrate");
    memory->grow(rxnpropensity,n,"app/surfchemtest:adespropensity");
    rxntype = memory->grow(rxntype,n,2,"app/surfchemtest:adestype");
    rxninput = memory->grow(rxninput,n,2,"app/surfchemtest:adesinput");
    rxnoutput = memory->grow(rxnoutput,n,2,"app/surfchemtest:adesoutput");
    memory->grow(rxncount,n,"app/surfchemtest:adescount");
    memory->grow(reactionsorbate,n,"app/surfchemtest:adesdesorbate");
    memory->grow(neighboring_rxn,n,"app/surfchemtest:neighboring_ades");
  }
}

#if defined(MUI)

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

void AppSurfchemtest::mui_print_MUIdblval(int step,const char *str1,const char *str2)
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

void AppSurfchemtest::mui_print_MUIintval(int step,const char *str1,const char *str2)
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
    } else if (strcmp(arg[k],"dac1") == 0) {     // i13
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_dac1",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},dac1[i]);
        dac1[i] = 0;
      }
    } else if (strcmp(arg[k],"dac2") == 0) {     // i14
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_dac2",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},dac2[i]);
        dac2[i] = 0;
      }
    } else if (strcmp(arg[k],"dac3") == 0) {     // i15
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_dac3",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},dac3[i]);
        dac3[i] = 0;
      }
    } else if (strcmp(arg[k],"dac4") == 0) {     // i16
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_dac4",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},dac4[i]);
        dac4[i] = 0;
      }
    } else if (strcmp(arg[k],"dac5") == 0) {     // i17
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_dac5",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},dac5[i]);
        dac5[i] = 0;
      }
    } else if (strcmp(arg[k],"adc1") == 0) {     // i18
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_adc1",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},adc1[i]);
        adc1[i] = 0;
      }
    } else if (strcmp(arg[k],"adc2") == 0) {     // i19
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_adc2",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},adc2[i]);
        adc2[i] = 0;
      }
    } else if (strcmp(arg[k],"adc3") == 0) {     // i20
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_adc3",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},adc3[i]);
        adc3[i] = 0;
      }
    } else if (strcmp(arg[k],"adc4") == 0) {     // i21
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_adc4",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},adc4[i]);
        adc4[i] = 0;
      }
    } else if (strcmp(arg[k],"adc5") == 0) {     // i22
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_adc5",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},adc5[i]);
        adc5[i] = 0;
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
    } else if (strcmp(arg[k],"Vz") == 0) {          // d7
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_Vz",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},Vz[i]);
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
    } else if (strcmp(arg[k],"occ3") == 0) {
      for (int i=0;i<nlocal;i++) {
        int is_occ = (element[i]==3) ? 1 : 0;
        spk->uniface->push("CH_occ3",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},is_occ);
      }
    } else if (strcmp(arg[k],"occ4") == 0) {
      for (int i=0;i<nlocal;i++) {
        int is_occ = (element[i]==4) ? 1 : 0;
        spk->uniface->push("CH_occ4",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},is_occ);
      }
    } else if (strcmp(arg[k],"occ5") == 0) {
      for (int i=0;i<nlocal;i++) {
        int is_occ = (element[i]==5) ? 1 : 0;
        spk->uniface->push("CH_occ5",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},is_occ);
      }
    } else if (strcmp(arg[k],"one") == 0) {
      for (int i=0;i<nlocal;i++) {
        spk->uniface->push("CH_one",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},1);
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

void AppSurfchemtest::mui_commit(int narg, char **arg)
{
  int timestamp = atoi(arg[0]);

  spk->uniface->commit(timestamp);

  if (domain->me == 0 && screen) {
    fprintf(screen,"** DEBUG: mui commit at timestamp %d\n",timestamp);
    fflush(screen);
  }

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
    } else if (strcmp(arg[k],"Vz") == 0) {          // d7
      for (int i=0;i<nlocal;i++) {
        Vz[i] = spk->uniface->fetch("CH_Vz",{xyz[i][0]+mui_kmc_lattice_offset_x,xyz[i][1]+mui_kmc_lattice_offset_y},timestamp,s,t);
      }
    } else {
      error->all(FLERR,"Illegal mui_fetch command");
    }

    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: %s fetched\n",arg[k]);
      fflush(screen);
    }
  }

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
    }
    else if (strcmp(arg[k],"temp") == 0) {
      // get info for each FHD cell
      for (int n=0;n<nlocalFHDcell;n++)
        MUIdblval[n] = spk->uniface->fetch("CH_temp",{xFHD[n],yFHD[n]},timestamp,s,t);
      // distribute info to each KMC site
      for (int i=0;i<nlocal;i++)
        temp[i] = MUIdblval[localFHDcell[i]];
    }
    else if (strcmp(arg[k],"Vz") == 0) {
      // get info for each FHD cell
      for (int n=0;n<nlocalFHDcell;n++)
        MUIdblval[n] = spk->uniface->fetch("CH_Vz",{xFHD[n],yFHD[n]},timestamp,s,t);
      // distribute info to each KMC site
      for (int i=0;i<nlocal;i++)
        Vz[i] = MUIdblval[localFHDcell[i]];
    }
    else {
      error->all(FLERR,"Illegal mui_fetch_agg command");
    }

    if (domain->me == 0 && screen) {
      fprintf(screen,"** DEBUG: %s fetched\n",arg[k]);
      fflush(screen);
    }
  }

  return;
}

void AppSurfchemtest::mui_forget(int narg, char **arg)
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

void AppSurfchemtest::amrex_init_agg ()
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

void AppSurfchemtest::amrex_push_agg(int narg, char **arg)
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
        } else if (std::strcmp(arg[k],"adc4") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                intval[localFHDcell[i]] += adc4[i];
                adc4[i] = 0;
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
        } else if (std::strcmp(arg[k],"occ4") == 0) {
            // compute the sum over each FHD domain
            for (int n=0;n<nlocalFHDcell;n++) intval[n] = 0;
            for (int i=0;i<nlocal;i++) {
                int is_occ = (element[i]==4) ? 1 : 0;
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

void AppSurfchemtest::amrex_fetch_agg(int narg, char **arg)
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
        if (std::strcmp(arg[k],"density1") == 0) {
            amrex_recv_dblval();
            // distribute info to each KMC site
            for (int i=0;i<nlocal;i++)
                density1[i] = dblval[localFHDcell[i]];
        }
        else if (std::strcmp(arg[k],"density2") == 0) {
            amrex_recv_dblval();
            // distribute info to each KMC site
            for (int i=0;i<nlocal;i++)
                density2[i] = dblval[localFHDcell[i]];
        }
        else if (std::strcmp(arg[k],"density3") == 0) {
            amrex_recv_dblval();
            // distribute info to each KMC site
            for (int i=0;i<nlocal;i++)
                density3[i] = dblval[localFHDcell[i]];
        }
        else if (std::strcmp(arg[k],"temp") == 0) {
            amrex_recv_dblval();
            // distribute info to each KMC site
            for (int i=0;i<nlocal;i++)
                temp[i] = dblval[localFHDcell[i]];
        }
        else if (std::strcmp(arg[k],"Vz") == 0) {
            amrex_recv_dblval();
            // distribute info to each KMC site
            for (int i=0;i<nlocal;i++)
                Vz[i] = dblval[localFHDcell[i]];
        }
        else {
            error->all(FLERR,"Illegal amrex_fetch_agg command");
        }

        if (domain->me == 0 && screen) {
            std::fprintf(screen,"** DEBUG: %s fetched\n",arg[k]);
            std::fflush(screen);
        }
    }
}

void AppSurfchemtest::amrex_send_intval()
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

void AppSurfchemtest::amrex_recv_dblval()
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
