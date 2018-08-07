#include "surfaces.H"
#include <math.h>

using namespace amrex;

double getTheta(double nx, double ny, double nz)
{

#if (BL_SPACEDIM == 3)
    double r = sqrt(nx*nx+ny*ny+nz*nz);
    return acos(nz/r);
#endif
#if (BL_SPACEDIM == 2)
    double r = sqrt(nx*nx+ny*ny+nz*nz);
    return acos(nx/r);
#endif
}

double getPhi(double nx, double ny, double nz)
{
    return atan2(ny,nx);
}

void BuildSurfaces(surface* surfaceList, const int surfaces, const Real* domainLo, const Real* domainHi)
{
#if (BL_SPACEDIM == 3)

    double theta, phi;
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

    surfaceList[0].lnx = -1;
    surfaceList[0].lny = 0;
    surfaceList[0].lnz = 0;
    
    surfaceList[0].rnx = 1;
    surfaceList[0].rny = 0;
    surfaceList[0].rnz = 0;

    surfaceList[0].porosityLeft = 0;
    surfaceList[0].specularityLeft = 0;
    surfaceList[0].temperatureLeft = 273;
    surfaceList[0].momentumConsLeft = 1;

    surfaceList[0].porosityRight = 0;
    surfaceList[0].specularityRight = 0;
    surfaceList[0].temperatureRight = 273;
    surfaceList[0].momentumConsRight = 1;

    surfaceList[0].periodicity = 0;
    surfaceList[0].boundary = 1;
        

    theta = getTheta(surfaceList[0].lnx, surfaceList[0].lny, surfaceList[0].lnz);
    phi   = getPhi(surfaceList[0].lnx, surfaceList[0].lny, surfaceList[0].lnz);

    surfaceList[0].cosThetaLeft = cos(theta);
    surfaceList[0].sinThetaLeft = sin(theta);
    surfaceList[0].cosPhiLeft = cos(phi);
    surfaceList[0].sinPhiLeft = sin(phi);

    theta = getTheta(surfaceList[0].rnx, surfaceList[0].rny, surfaceList[0].rnz);
    phi   = getPhi(surfaceList[0].rnx, surfaceList[0].rny, surfaceList[0].rnz);

    surfaceList[0].cosThetaRight = cos(theta);
    surfaceList[0].sinThetaRight = sin(theta);
    surfaceList[0].cosPhiRight = cos(phi);
    surfaceList[0].sinPhiRight = sin(phi);

    surfaceList[0].fxLeftAv = 0;
    surfaceList[0].fyLeftAv = 0;
    surfaceList[0].fzLeftAv = 0;

    surfaceList[0].fxRightAv = 0;
    surfaceList[0].fyRightAv = 0;
    surfaceList[0].fzRightAv = 0;

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

    surfaceList[1].lnx = -1;
    surfaceList[1].lny = 0;
    surfaceList[1].lnz = 0;
    
    surfaceList[1].rnx = 1;
    surfaceList[1].rny = 0;
    surfaceList[1].rnz = 0;

    surfaceList[1].porosityLeft = 0;
    surfaceList[1].specularityLeft = 0;
    surfaceList[1].temperatureLeft = 819;
    surfaceList[1].momentumConsLeft = 1;

    surfaceList[1].porosityRight = 0;
    surfaceList[1].specularityRight = 0;
    surfaceList[1].temperatureRight = 819;
    surfaceList[1].momentumConsRight = 1;

    surfaceList[1].periodicity = 0;
    surfaceList[1].boundary = 2;

    theta = getTheta(surfaceList[1].lnx, surfaceList[1].lny, surfaceList[1].lnz);
    phi   = getPhi(surfaceList[1].lnx, surfaceList[1].lny, surfaceList[1].lnz);

    surfaceList[1].cosThetaLeft = cos(theta);
    surfaceList[1].sinThetaLeft = sin(theta);
    surfaceList[1].cosPhiLeft = cos(phi);
    surfaceList[1].sinPhiLeft = sin(phi);

    theta = getTheta(surfaceList[1].rnx, surfaceList[1].rny, surfaceList[1].rnz);
    phi   = getPhi(surfaceList[1].rnx, surfaceList[1].rny, surfaceList[1].rnz);

    surfaceList[1].cosThetaRight = cos(theta);
    surfaceList[1].sinThetaRight = sin(theta);
    surfaceList[1].cosPhiRight = cos(phi);
    surfaceList[1].sinPhiRight = sin(phi);

    surfaceList[1].fxLeftAv = 0;
    surfaceList[1].fyLeftAv = 0;
    surfaceList[1].fzLeftAv = 0;

    surfaceList[1].fxRightAv = 0;
    surfaceList[1].fyRightAv = 0;
    surfaceList[1].fzRightAv = 0;

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
    surfaceList[2].temperatureLeft = 300;
    surfaceList[2].momentumConsLeft = 1;

    surfaceList[2].porosityRight = 1;
    surfaceList[2].specularityRight = 1;
    surfaceList[2].temperatureRight = 300;
    surfaceList[2].momentumConsRight = 1;

    surfaceList[2].periodicity = 1;
    surfaceList[2].boundary = 3;

    theta = getTheta(surfaceList[2].lnx, surfaceList[2].lny, surfaceList[2].lnz);
    phi   = getPhi(surfaceList[2].lnx, surfaceList[2].lny, surfaceList[2].lnz);

    surfaceList[2].cosThetaLeft = cos(theta);
    surfaceList[2].sinThetaLeft = sin(theta);
    surfaceList[2].cosPhiLeft = cos(phi);
    surfaceList[2].sinPhiLeft = sin(phi);

    theta = getTheta(surfaceList[2].rnx, surfaceList[2].rny, surfaceList[2].rnz);
    phi   = getPhi(surfaceList[2].rnx, surfaceList[2].rny, surfaceList[2].rnz);

    surfaceList[2].cosThetaRight = cos(theta);
    surfaceList[2].sinThetaRight = sin(theta);
    surfaceList[2].cosPhiRight = cos(phi);
    surfaceList[2].sinPhiRight = sin(phi);

    surfaceList[2].fxLeftAv = 0;
    surfaceList[2].fyLeftAv = 0;
    surfaceList[2].fzLeftAv = 0;

    surfaceList[2].fxRightAv = 0;
    surfaceList[2].fyRightAv = 0;
    surfaceList[2].fzRightAv = 0;

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
    surfaceList[3].temperatureLeft = 700;
    surfaceList[3].momentumConsLeft = 1;

    surfaceList[3].porosityRight = 1;
    surfaceList[3].specularityRight = 0;
    surfaceList[3].temperatureRight = 700;
    surfaceList[3].momentumConsRight = 1;

    surfaceList[3].periodicity = 1;
    surfaceList[3].boundary = 4;

    theta = getTheta(surfaceList[3].lnx, surfaceList[3].lny, surfaceList[3].lnz);
    phi   = getPhi(surfaceList[3].lnx, surfaceList[3].lny, surfaceList[3].lnz);

    surfaceList[3].cosThetaLeft = cos(theta);
    surfaceList[3].sinThetaLeft = sin(theta);
    surfaceList[3].cosPhiLeft = cos(phi);
    surfaceList[3].sinPhiLeft = sin(phi);

    theta = getTheta(surfaceList[3].rnx, surfaceList[3].rny, surfaceList[3].rnz);
    phi   = getPhi(surfaceList[3].rnx, surfaceList[3].rny, surfaceList[3].rnz);

    surfaceList[3].cosThetaRight = cos(theta);
    surfaceList[3].sinThetaRight = sin(theta);
    surfaceList[3].cosPhiRight = cos(phi);
    surfaceList[3].sinPhiRight = sin(phi);

    surfaceList[3].fxLeftAv = 0;
    surfaceList[3].fyLeftAv = 0;
    surfaceList[3].fzLeftAv = 0;

    surfaceList[3].fxRightAv = 0;
    surfaceList[3].fyRightAv = 0;
    surfaceList[3].fzRightAv = 0;

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
    surfaceList[4].momentumConsLeft = 1;

    surfaceList[4].porosityRight = 1;
    surfaceList[4].specularityRight = 0;
    surfaceList[4].temperatureRight = 273;
    surfaceList[4].momentumConsRight = 1;

    surfaceList[4].periodicity = 1;
    surfaceList[4].boundary = 5;

    theta = getTheta(surfaceList[4].lnx, surfaceList[4].lny, surfaceList[4].lnz);
    phi   = getPhi(surfaceList[4].lnx, surfaceList[4].lny, surfaceList[4].lnz);

    surfaceList[4].cosThetaLeft = cos(theta);
    surfaceList[4].sinThetaLeft = sin(theta);
    surfaceList[4].cosPhiLeft = cos(phi);
    surfaceList[4].sinPhiLeft = sin(phi);

    theta = getTheta(surfaceList[4].rnx, surfaceList[4].rny, surfaceList[4].rnz);
    phi   = getPhi(surfaceList[4].rnx, surfaceList[4].rny, surfaceList[4].rnz);

    surfaceList[4].cosThetaRight = cos(theta);
    surfaceList[4].sinThetaRight = sin(theta);
    surfaceList[4].cosPhiRight = cos(phi);
    surfaceList[4].sinPhiRight = sin(phi);

    surfaceList[4].fxLeftAv = 0;
    surfaceList[4].fyLeftAv = 0;
    surfaceList[4].fzLeftAv = 0;

    surfaceList[4].fxRightAv = 0;
    surfaceList[4].fyRightAv = 0;
    surfaceList[4].fzRightAv = 0;

//domianHi z plane
    surfaceList[5].x0 = domainLo[0];
    surfaceList[5].y0 = domainLo[1];
    surfaceList[5].z0 = domainHi[2];

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
    surfaceList[5].momentumConsLeft = 1;

    surfaceList[5].porosityRight = 1;
    surfaceList[5].specularityRight = 0;
    surfaceList[5].temperatureRight = 273;
    surfaceList[5].momentumConsRight = 1;

    surfaceList[5].periodicity = 1;
    surfaceList[5].boundary = 6;

    theta = getTheta(surfaceList[5].lnx, surfaceList[5].lny, surfaceList[5].lnz);
    phi   = getPhi(surfaceList[5].lnx, surfaceList[5].lny, surfaceList[5].lnz);

    surfaceList[5].cosThetaLeft = cos(theta);
    surfaceList[5].sinThetaLeft = sin(theta);
    surfaceList[5].cosPhiLeft = cos(phi);
    surfaceList[5].sinPhiLeft = sin(phi);

    theta = getTheta(surfaceList[5].rnx, surfaceList[5].rny, surfaceList[5].rnz);
    phi   = getPhi(surfaceList[5].rnx, surfaceList[5].rny, surfaceList[5].rnz);

    surfaceList[5].cosThetaRight = cos(theta);
    surfaceList[5].sinThetaRight = sin(theta);
    surfaceList[5].cosPhiRight = cos(phi);
    surfaceList[5].sinPhiRight = sin(phi);

    surfaceList[5].fxLeftAv = 0;
    surfaceList[5].fyLeftAv = 0;
    surfaceList[5].fzLeftAv = 0;

    surfaceList[5].fxRightAv = 0;
    surfaceList[5].fyRightAv = 0;
    surfaceList[5].fzRightAv = 0;

//Membrane
    surfaceList[6].x0 = 3.753e-7;
    surfaceList[6].y0 = domainLo[1];
    surfaceList[6].z0 = domainLo[2];

    surfaceList[6].ux = 0;
    surfaceList[6].uy = 1;
    surfaceList[6].uz = 0;

    surfaceList[6].vx = 0;
    surfaceList[6].vy = 0;
    surfaceList[6].vz = 1;

    surfaceList[6].uTop = domainHi[1];
    surfaceList[6].vTop = domainHi[2];

    surfaceList[6].lnx = -1;
    surfaceList[6].lny = 0;
    surfaceList[6].lnz = 0;
    
    surfaceList[6].rnx = 1;
    surfaceList[6].rny = 0;
    surfaceList[6].rnz = 0;

    surfaceList[6].porosityLeft = 0.2;
    surfaceList[6].specularityLeft = 1;
    surfaceList[6].temperatureLeft = 273;
    surfaceList[6].momentumConsLeft = 1;

    surfaceList[6].porosityRight = 0.2;
    surfaceList[6].specularityRight = 1;
    surfaceList[6].temperatureRight = 273;
    surfaceList[6].momentumConsRight = 1;

    surfaceList[6].periodicity = 0;
    surfaceList[6].boundary = 0;
        
    theta = getTheta(surfaceList[6].lnx, surfaceList[6].lny, surfaceList[6].lnz);
    phi   = getPhi(surfaceList[6].lnx, surfaceList[6].lny, surfaceList[6].lnz);

    surfaceList[6].cosThetaLeft = cos(theta);
    surfaceList[6].sinThetaLeft = sin(theta);
    surfaceList[6].cosPhiLeft = cos(phi);
    surfaceList[6].sinPhiLeft = sin(phi);

    theta = getTheta(surfaceList[6].rnx, surfaceList[6].rny, surfaceList[6].rnz);
    phi   = getPhi(surfaceList[6].rnx, surfaceList[6].rny, surfaceList[6].rnz);

    surfaceList[6].cosThetaRight = cos(theta);
    surfaceList[6].sinThetaRight = sin(theta);
    surfaceList[6].cosPhiRight = cos(phi);
    surfaceList[6].sinPhiRight = sin(phi);

    surfaceList[6].fxLeftAv = 0;
    surfaceList[6].fyLeftAv = 0;
    surfaceList[6].fzLeftAv = 0;

    surfaceList[6].fxRightAv = 0;
    surfaceList[6].fyRightAv = 0;
    surfaceList[6].fzRightAv = 0;

#endif

#if (BL_SPACEDIM == 2)

    double theta;

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

    surfaceList[0].porosityLeft = 0;
    surfaceList[0].specularityLeft = 0;
    surfaceList[0].temperatureLeft = 273;
    surfaceList[0].momentumConsLeft = 1;

    surfaceList[0].porosityRight = 0;
    surfaceList[0].specularityRight = 0;
    surfaceList[0].temperatureRight = 273;
    surfaceList[0].momentumConsRight = 1;

    surfaceList[0].periodicity = 0;
    surfaceList[0].boundary = 1;

    theta = getTheta(surfaceList[0].lnx, surfaceList[0].lny, 0);

    surfaceList[0].cosThetaLeft = cos(theta);
    surfaceList[0].sinThetaLeft = sin(theta);

    theta = getTheta(surfaceList[0].rnx, surfaceList[0].rny, 0);

    surfaceList[0].cosThetaRight = cos(theta);
    surfaceList[0].sinThetaRight = sin(theta);

    surfaceList[0].fxLeftAv = 0;
    surfaceList[0].fyLeftAv = 0;
    surfaceList[0].fzLeftAv = 0;

    surfaceList[0].fxRightAv = 0;
    surfaceList[0].fyRightAv = 0;
    surfaceList[0].fzRightAv = 0;


//domianHi x plane
    surfaceList[1].x0 = domainHi[0];
    surfaceList[1].y0 = domainLo[1];

    surfaceList[1].ux = 0;
    surfaceList[1].uy = 1;

    surfaceList[1].uTop = domainHi[1];

    surfaceList[1].lnx = 1;
    surfaceList[1].lny = 0;
    
    surfaceList[1].rnx = -1;
    surfaceList[1].rny = 0;

    surfaceList[1].porosityLeft = 0;
    surfaceList[1].specularityLeft = 0;
    surfaceList[1].temperatureLeft = 819;
    surfaceList[1].momentumConsLeft = 1;

    surfaceList[1].porosityRight = 0;
    surfaceList[1].specularityRight = 0;
    surfaceList[1].temperatureRight = 819;
    surfaceList[1].momentumConsRight = 1;

    surfaceList[1].periodicity = 0;
    surfaceList[1].boundary = 2;

    theta = getTheta(surfaceList[1].lnx, surfaceList[1].lny, 0);

    surfaceList[1].cosThetaLeft = cos(theta);
    surfaceList[1].sinThetaLeft = sin(theta);

    theta = getTheta(surfaceList[1].rnx, surfaceList[1].rny, 0);

    surfaceList[1].cosThetaRight = cos(theta);
    surfaceList[1].sinThetaRight = sin(theta);

    surfaceList[1].fxLeftAv = 0;
    surfaceList[1].fyLeftAv = 0;
    surfaceList[1].fzLeftAv = 0;

    surfaceList[1].fxRightAv = 0;
    surfaceList[1].fyRightAv = 0;
    surfaceList[1].fzRightAv = 0;

//domainLo y plane
    surfaceList[2].x0 = domainLo[0];
    surfaceList[2].y0 = domainLo[1];

    surfaceList[2].ux = 1;
    surfaceList[2].uy = 0;

    surfaceList[2].uTop = domainHi[0];

    surfaceList[2].lnx = 0;
    surfaceList[2].lny = 1;
    
    surfaceList[2].rnx = 0;
    surfaceList[2].rny = -1;

    surfaceList[2].porosityLeft = 1;
    surfaceList[2].specularityLeft = 1;
    surfaceList[2].temperatureLeft = 273;
    surfaceList[2].momentumConsLeft = 1;

    surfaceList[2].porosityRight = 1;
    surfaceList[2].specularityRight = 1;
    surfaceList[2].temperatureRight = 273;
    surfaceList[2].momentumConsRight = 1;

    surfaceList[2].periodicity = 1;
    surfaceList[2].boundary = 3;

    theta = getTheta(surfaceList[2].lnx, surfaceList[2].lny, 0);

    surfaceList[2].cosThetaLeft = cos(theta);
    surfaceList[2].sinThetaLeft = sin(theta);

    theta = getTheta(surfaceList[2].rnx, surfaceList[2].rny, 0);

    surfaceList[2].cosThetaRight = cos(theta);
    surfaceList[2].sinThetaRight = sin(theta);

    surfaceList[2].fxLeftAv = 0;
    surfaceList[2].fyLeftAv = 0;
    surfaceList[2].fzLeftAv = 0;

    surfaceList[2].fxRightAv = 0;
    surfaceList[2].fyRightAv = 0;
    surfaceList[2].fzRightAv = 0;

//domianHi y plane
    surfaceList[3].x0 = domainLo[0];
    surfaceList[3].y0 = domainHi[1];

    surfaceList[3].ux = 1;
    surfaceList[3].uy = 0;

    surfaceList[3].uTop = domainHi[0];

    surfaceList[3].lnx = 0;
    surfaceList[3].lny = 1;
    
    surfaceList[3].rnx = 0;
    surfaceList[3].rny = -1;

    surfaceList[3].porosityLeft = 1;
    surfaceList[3].specularityLeft = 1;
    surfaceList[3].temperatureLeft = 273;
    surfaceList[3].momentumConsLeft = 1;

    surfaceList[3].porosityRight = 1;
    surfaceList[3].specularityRight = 1;
    surfaceList[3].temperatureRight = 273;
    surfaceList[3].momentumConsRight = 1;

    surfaceList[3].periodicity = 1;
    surfaceList[3].boundary = 4;

    theta = getTheta(surfaceList[3].lnx, surfaceList[3].lny, 0);

    surfaceList[3].cosThetaLeft = cos(theta);
    surfaceList[3].sinThetaLeft = sin(theta);

    theta = getTheta(surfaceList[3].rnx, surfaceList[3].rny, 0);

    surfaceList[3].cosThetaRight = cos(theta);
    surfaceList[3].sinThetaRight = sin(theta);

    surfaceList[3].fxLeftAv = 0;
    surfaceList[3].fyLeftAv = 0;
    surfaceList[3].fzLeftAv = 0;

    surfaceList[3].fxRightAv = 0;
    surfaceList[3].fyRightAv = 0;
    surfaceList[3].fzRightAv = 0;

//Membrane
    surfaceList[4].x0 = 3.753e-7;
    surfaceList[4].y0 = domainLo[1];

    surfaceList[4].ux = 0;
    surfaceList[4].uy = 1;

    surfaceList[4].uTop = domainHi[1];

    surfaceList[4].lnx = -1;
    surfaceList[4].lny = 0;
    
    surfaceList[4].rnx = 1;
    surfaceList[4].rny = 0;

    surfaceList[4].porosityLeft = 0.05;
    surfaceList[4].specularityLeft = 1;
    surfaceList[4].temperatureLeft = 273;
    surfaceList[4].momentumConsLeft = 1;

    surfaceList[4].porosityRight = 0.05;
    surfaceList[4].specularityRight = 1;
    surfaceList[4].temperatureRight = 273;
    surfaceList[4].momentumConsRight = 1;

    surfaceList[4].periodicity = 0;
    surfaceList[4].boundary = 0;
        
    theta = getTheta(surfaceList[4].lnx, surfaceList[4].lny, 0);

    surfaceList[4].cosThetaLeft = cos(theta);
    surfaceList[4].sinThetaLeft = sin(theta);

    theta = getTheta(surfaceList[4].rnx, surfaceList[6].rny, 0);

    surfaceList[4].cosThetaRight = cos(theta);
    surfaceList[4].sinThetaRight = sin(theta);

    surfaceList[4].fxLeftAv = 0;
    surfaceList[4].fyLeftAv = 0;
    surfaceList[4].fzLeftAv = 0;

    surfaceList[4].fxRightAv = 0;
    surfaceList[4].fyRightAv = 0;
    surfaceList[4].fzRightAv = 0;

#endif

}
