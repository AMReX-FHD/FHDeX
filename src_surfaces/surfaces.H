#ifndef _surfaces_H_
#define _surfaces_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

//Parametric plane

typedef struct {

#if (BL_SPACEDIM == 3)
    double x0;
    double y0;
    double z0;

    double ux;
    double uy;
    double uz;

    double vx;
    double vy;
    double vz;

    double lnx;
    double lny;
    double lnz;

    double rnx;
    double rny;
    double rnz;

    double uTop;
    double vTop;

    double cosThetaLeft;
    double sinThetaLeft;
    double cosPhiLeft;
    double sinPhiLeft;

    double cosThetaRight;
    double sinThetaRight;
    double cosPhiRight;
    double sinPhiRight;


#endif

#if (BL_SPACEDIM == 2)
    double x0;
    double y0;

    double ux;
    double uy;

    double lnx;
    double lny;

    double rnx;
    double rny;

    double uTop;

    double cosThetaLeft;
    double sinThetaLeft;

    double cosThetaRight;
    double sinThetaRight;
#endif

    double porosityLeft;
    double specularityLeft;
    double temperatureLeft;

    double porosityRight;
    double specularityRight;
    double temperatureRight;

} surface;

void BuildSurfaces(surface* surfaceList, const int surfaces, const Real* domainLo, const Real* domainHi);

#endif
