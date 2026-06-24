#ifndef _paramPlane_H_
#define _paramPlane_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <common_namespace.H>

using namespace amrex;

//Parametric plane

typedef struct {

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

    double porosityLeft;
    double specularityLeft;
    double temperatureLeft;
    double momentumConsLeft;
    double densityLeft[MAX_SPECIES];

    double porosityRight;
    double specularityRight;
    double temperatureRight;
    double momentumConsRight;
    double densityRight[MAX_SPECIES];

    double massFluxRight[MAX_SPECIES];
    double massFluxLeft[MAX_SPECIES];

    double xMomFluxRight[MAX_SPECIES];
    double xMomFluxLeft[MAX_SPECIES];
    double yMomFluxRight[MAX_SPECIES];
    double yMomFluxLeft[MAX_SPECIES];
    double zMomFluxRight[MAX_SPECIES];
    double zMomFluxLeft[MAX_SPECIES];

    double energyFluxRight[MAX_SPECIES];
    double energyFluxLeft[MAX_SPECIES];

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

    double* xVelRecLeft;
    double* yVelRecLeft;
    double* zVelRecLeft;

    double* xVelRecRight;
    double* yVelRecRight;
    double* zVelRecRight;

    double* xPosRecLeft;
    double* yPosRecLeft;
    double* zPosRecLeft;

    double* xPosRecRight;
    double* yPosRecRight;
    double* zPosRecRight;

    double* freqRecRight;
    double* freqRecLeft;

    double* timeRecRight;
    double* timeRecLeft;

    int recCountRight;
    int recCountLeft;

} paramPlane;

void BuildParamplanes(paramPlane* paramPlaneList, const int paramplanes, const Real* domainLo, const Real* domainHi);
void BuildParamplanesPhonon(paramPlane* paramPlaneList, const int paramplanes, const Real* domainLo, const Real* domainHi);

void SetBoundaryCells(paramPlane* paramPlaneList, const int paramplanes, const Real* domainLo, const Real* domainHi, const Geometry & Geom, iMultiFab& bCell, int paramPlaneCount);

double getTheta(double nx, double ny, double nz);
double getPhi(double nx, double ny, double nz);

#endif
