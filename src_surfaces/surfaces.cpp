#include "surfaces.H"

using namespace amrex;

void BuildSurfaces(surface* surfaceList, const int surfaces, const Real* domainLo, const Real* domainHi)
{
#if (BL_SPACEDIM == 3)

//Domain boundaries

//domainLo x plane
    surfaceList[0].x0 = domainLo[0];
    surfaceList[0].y0 = domainLo[1];
    surfaceList[0].z0 = domainLo[2];

    surfaceList[0].ux = 0;
    surfaceList[0].uy = 1;
    surfaceList[0].uz = 0;

    surfaceList[0].vx = 0;
    surfaceList[0].vy = 0;
    surfaceList[0].vz = 1;

    surfaceList[0].uTop = domainHi[1];
    surfaceList[0].vTop = domainHi[2];

    surfaceList[0].lnx = 1;
    surfaceList[0].lny = 0;
    surfaceList[0].lnz = 0;
    
    surfaceList[0].rnx = -1;
    surfaceList[0].rny = 0;
    surfaceList[0].rnz = 0;

    surfaceList[0].porosityLeft = 1;
    surfaceList[0].specularityLeft = 1;
    surfaceList[0].temperatureLeft = 273;

    surfaceList[0].porosityRight = 1;
    surfaceList[0].specularityRight = 1;
    surfaceList[0].temperatureRight = 273;


//domianHi x plane
    surfaceList[1].x0 = domainHi[0];
    surfaceList[1].y0 = domainLo[1];
    surfaceList[1].z0 = domainLo[2];

    surfaceList[1].ux = 0;
    surfaceList[1].uy = 1;
    surfaceList[1].uz = 0;

    surfaceList[1].vx = 0;
    surfaceList[1].vy = 0;
    surfaceList[1].vz = 1;

    surfaceList[1].uTop = domainHi[1];
    surfaceList[1].vTop = domainHi[2];

    surfaceList[1].lnx = 1;
    surfaceList[1].lny = 0;
    surfaceList[1].lnz = 0;
    
    surfaceList[1].rnx = -1;
    surfaceList[1].rny = 0;
    surfaceList[1].rnz = 0;

    surfaceList[1].porosityLeft = 1;
    surfaceList[1].specularityLeft = 1;
    surfaceList[1].temperatureLeft = 273;

    surfaceList[1].porosityRight = 1;
    surfaceList[1].specularityRight = 1;
    surfaceList[1].temperatureRight = 273;

//domainLo y plane
    surfaceList[2].x0 = domainLo[0];
    surfaceList[2].y0 = domainLo[1];
    surfaceList[2].z0 = domainLo[2];

    surfaceList[2].ux = 1;
    surfaceList[2].uy = 0;
    surfaceList[2].uz = 0;

    surfaceList[2].vx = 0;
    surfaceList[2].vy = 0;
    surfaceList[2].vz = 1;

    surfaceList[2].uTop = domainHi[0];
    surfaceList[2].vTop = domainHi[2];

    surfaceList[2].lnx = 0;
    surfaceList[2].lny = 1;
    surfaceList[2].lnz = 0;
    
    surfaceList[2].rnx = 0;
    surfaceList[2].rny = -1;
    surfaceList[2].rnz = 0;

    surfaceList[2].porosityLeft = 1;
    surfaceList[2].specularityLeft = 1;
    surfaceList[2].temperatureLeft = 273;

    surfaceList[2].porosityRight = 1;
    surfaceList[2].specularityRight = 1;
    surfaceList[2].temperatureRight = 273;


//domianHi y plane
    surfaceList[3].x0 = domainLo[0];
    surfaceList[3].y0 = domainHi[1];
    surfaceList[3].z0 = domainLo[2];

    surfaceList[3].ux = 1;
    surfaceList[3].uy = 0;
    surfaceList[3].uz = 0;

    surfaceList[3].vx = 0;
    surfaceList[3].vy = 0;
    surfaceList[3].vz = 1;

    surfaceList[3].uTop = domainHi[0];
    surfaceList[3].vTop = domainHi[2];

    surfaceList[3].lnx = 0;
    surfaceList[3].lny = 1;
    surfaceList[3].lnz = 0;
    
    surfaceList[3].rnx = 0;
    surfaceList[3].rny = -1;
    surfaceList[3].rnz = 0;

    surfaceList[3].porosityLeft = 1;
    surfaceList[3].specularityLeft = 1;
    surfaceList[3].temperatureLeft = 273;

    surfaceList[3].porosityRight = 1;
    surfaceList[3].specularityRight = 1;
    surfaceList[3].temperatureRight = 273;

//domainLo z plane
    surfaceList[4].x0 = domainLo[0];
    surfaceList[4].y0 = domainLo[1];
    surfaceList[4].z0 = domainLo[2];

    surfaceList[4].ux = 1;
    surfaceList[4].uy = 0;
    surfaceList[4].uz = 0;

    surfaceList[4].vx = 0;
    surfaceList[4].vy = 1;
    surfaceList[4].vz = 0;

    surfaceList[4].uTop = domainHi[0];
    surfaceList[4].vTop = domainHi[1];

    surfaceList[4].lnx = 0;
    surfaceList[4].lny = 0;
    surfaceList[4].lnz = 1;
    
    surfaceList[4].rnx = 0;
    surfaceList[4].rny = 0;
    surfaceList[4].rnz = -1;

    surfaceList[4].porosityLeft = 1;
    surfaceList[4].specularityLeft = 1;
    surfaceList[4].temperatureLeft = 273;

    surfaceList[4].porosityRight = 1;
    surfaceList[4].specularityRight = 1;
    surfaceList[4].temperatureRight = 273;


//domianHi z plane
    surfaceList[5].x0 = domainLo[0];
    surfaceList[5].y0 = domainLo[1];
    surfaceList[5].z0 = domainHi[3];

    surfaceList[5].ux = 1;
    surfaceList[5].uy = 0;
    surfaceList[5].uz = 0;

    surfaceList[5].vx = 0;
    surfaceList[5].vy = 1;
    surfaceList[5].vz = 0;

    surfaceList[5].uTop = domainHi[0];
    surfaceList[5].vTop = domainHi[1];

    surfaceList[5].lnx = 0;
    surfaceList[5].lny = 0;
    surfaceList[5].lnz = 1;
    
    surfaceList[5].rnx = 0;
    surfaceList[5].rny = 0;
    surfaceList[5].rnz = -1;

    surfaceList[5].porosityLeft = 1;
    surfaceList[5].specularityLeft = 1;
    surfaceList[5].temperatureLeft = 273;

    surfaceList[5].porosityRight = 1;
    surfaceList[5].specularityRight = 1;
    surfaceList[5].temperatureRight = 273;

#endif

#if (BL_SPACEDIM == 2)

//Domain boundaries

//domainLo x
    surfaceList[0].x0 = domainLo[0];
    surfaceList[0].y0 = domainLo[1];

    surfaceList[0].ux = 0;
    surfaceList[0].uy = 1;

    surfaceList[0].uTop = domainHi[1];

    surfaceList[0].lnx = 1;
    surfaceList[0].lny = 0;
    
    surfaceList[0].rnx = -1;
    surfaceList[0].rny = 0;

    surfaceList[0].porosityLeft = 1;
    surfaceList[0].specularityLeft = 1;
    surfaceList[0].temperatureLeft = 273;

    surfaceList[0].porosityRight = 1;
    surfaceList[0].specularityRight = 1;
    surfaceList[0].temperatureRight = 273;


//domianHi x plane
    surfaceList[0].x0 = domainHi[0];
    surfaceList[0].y0 = domainLo[1];

    surfaceList[0].ux = 0;
    surfaceList[0].uy = 1;

    surfaceList[0].uTop = domainHi[1];

    surfaceList[0].lnx = 1;
    surfaceList[0].lny = 0;
    
    surfaceList[0].rnx = -1;
    surfaceList[0].rny = 0;

    surfaceList[0].porosityLeft = 1;
    surfaceList[0].specularityLeft = 1;
    surfaceList[0].temperatureLeft = 273;

    surfaceList[0].porosityRight = 1;
    surfaceList[0].specularityRight = 1;
    surfaceList[0].temperatureRight = 273;

//domainLo y plane
    surfaceList[0].x0 = domainLo[0];
    surfaceList[0].y0 = domainLo[1];

    surfaceList[0].ux = 1;
    surfaceList[0].uy = 0;

    surfaceList[0].uTop = domainHi[0];

    surfaceList[0].lnx = 0;
    surfaceList[0].lny = 1;
    
    surfaceList[0].rnx = 0;
    surfaceList[0].rny = -1;

    surfaceList[0].porosityLeft = 1;
    surfaceList[0].specularityLeft = 1;
    surfaceList[0].temperatureLeft = 273;

    surfaceList[0].porosityRight = 1;
    surfaceList[0].specularityRight = 1;
    surfaceList[0].temperatureRight = 273;

//domianHi y plane
    surfaceList[0].x0 = domainLo[0];
    surfaceList[0].y0 = domainHi[1];

    surfaceList[0].ux = 1;
    surfaceList[0].uy = 0;

    surfaceList[0].uTop = domainHi[0];

    surfaceList[0].lnx = 0;
    surfaceList[0].lny = 1;
    
    surfaceList[0].rnx = 0;
    surfaceList[0].rny = -1;

    surfaceList[0].porosityLeft = 1;
    surfaceList[0].specularityLeft = 1;
    surfaceList[0].temperatureLeft = 273;

    surfaceList[0].porosityRight = 1;
    surfaceList[0].specularityRight = 1;
    surfaceList[0].temperatureRight = 273;

#endif


}
