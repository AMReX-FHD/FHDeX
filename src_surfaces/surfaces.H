#ifndef _surfaces_H_
#define _surfaces_H_

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

    double nx;
    double ny;
    double nz;

    double utop;
    double vtop;

    double porosity;
    double specularity;

    double temperature;

} surface;

void BuildSurfaces(surface* surfaceList, int surfaces);

#endif
