#ifndef _species_H_
#define _species_H_

typedef struct {

  int type;
  int total;
  int ppb;
  double m;
  double d;
  double T;
  double R;
  double q;
  double gamma1;
  double mu;
  double n0;
  double mfp;
  double P;
  double Neff;
  double cp;
  double propulsion;
  double totalDiff;
  double dryDiff;
  double wetDiff;
  double sigma;
  double eepsilon;

} species;

typedef struct {

  int P;
  double mfp;
  double nu;
  double n0;
  double T;
  double R;

} kineticProperties;

#endif
