#ifndef _paramPlane_H_
#define _paramPlane_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <common_namespace.H>

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

    double velx;
    double vely;
    double velz;

    double c0;
    double resomg;
    double agraph;
    double bgraph;
    double a0graph;
    double b0graph;
    double coltime;

    double besslist[100];
    double dbesslist[100];

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
    double momentumConsLeft;

    double porosityRight;
    double specularityRight;
    double temperatureRight;
    double momentumConsRight;
    double densityRight[MAX_SPECIES];

    double massFluxRight[MAX_SPECIES];
    double massFluxLeft[MAX_SPECIES];

    double sourceLeft;
    double sourceRight;

    double sinkLeft;
    double sinkRight;

    double periodicity;
    double area;

    double fxLeft;
    double fyLeft;
    double fzLeft;

    double fxRight;
    double fyRight;
    double fzRight;

    double fxLeftAv;
    double fyLeftAv;
    double fzLeftAv;

    double fxRightAv;
    double fyRightAv;
    double fzRightAv;

    int boundary;

    double species[5];

    double tension;
    double stiffness;
    double temperature;

} paramPlane;

void BuildParamplanes(paramPlane* paramPlaneList, const int paramplanes, const Real* domainLo, const Real* domainHi);

double getTheta(double nx, double ny, double nz);
double getPhi(double nx, double ny, double nz);

#endif
