#include "paramPlane.H"
#include <math.h>

#include "common_functions.H"


using namespace amrex;

double getTheta(double nx, double ny, double nz)
{

    double r = sqrt(nx*nx+ny*ny+nz*nz);
    return acos(nz/r);


}

double getPhi(double nx, double ny, double nz)
{
    return atan2(ny,nx);
}

void BuildParamplanes(paramPlane* paramPlaneList, const int paramplanes, const Real* domainLo, const Real* domainHi)
{

    double theta, phi;
//Domain boundaries

//domainLo x plane
    paramPlaneList[0].x0 = domainLo[0];
    paramPlaneList[0].y0 = domainLo[1];
    paramPlaneList[0].z0 = domainLo[2];

    paramPlaneList[0].ux = 0;
    paramPlaneList[0].uy = 1;
    paramPlaneList[0].uz = 0;

    paramPlaneList[0].vx = 0;
    paramPlaneList[0].vy = 0;
    paramPlaneList[0].vz = 1;

    paramPlaneList[0].uTop = domainHi[1] - domainLo[1];
    paramPlaneList[0].vTop = domainHi[2] - domainLo[2];

    paramPlaneList[0].lnx = -1;
    paramPlaneList[0].lny = 0;
    paramPlaneList[0].lnz = 0;
    
    paramPlaneList[0].rnx = 1;
    paramPlaneList[0].rny = 0;
    paramPlaneList[0].rnz = 0;

    paramPlaneList[0].specularityLeft = 0;
    paramPlaneList[0].temperatureLeft = T_init[0];
    paramPlaneList[0].momentumConsLeft = 1;
    paramPlaneList[0].densityLeft[0] = 0;
    paramPlaneList[0].sourceLeft = 0;
    paramPlaneList[0].sinkLeft = 0;

    paramPlaneList[0].specularityRight = 0;
    paramPlaneList[0].temperatureRight = T_init[0];
    paramPlaneList[0].momentumConsRight = 1;
    paramPlaneList[0].densityRight[0] = 0;
    paramPlaneList[0].sourceRight = 0;
    paramPlaneList[0].sinkRight = 0;

    if(bc_vel_lo[0] == -1)
    {   
        paramPlaneList[0].periodicity = 1;
        paramPlaneList[0].porosityLeft = 1;
        paramPlaneList[0].porosityRight = 1;
    }
    else if(bc_vel_lo[0] == 2)
    {
        paramPlaneList[0].periodicity = 0;
        paramPlaneList[0].porosityLeft = 0;
        paramPlaneList[0].porosityRight = 0;
        paramPlaneList[0].specularityLeft = 1;
        paramPlaneList[0].specularityRight = 1;

       // paramPlaneList[0].x0 = domainLo[0] + sigma[0]/2.0;
    }
    else{

        paramPlaneList[0].periodicity = 0;
        paramPlaneList[0].porosityLeft = 0;
        paramPlaneList[0].porosityRight = 0;
        paramPlaneList[0].specularityLeft = 1;
        paramPlaneList[0].specularityRight = 1;
    }

    paramPlaneList[0].boundary = 1;
        
    theta = getTheta(paramPlaneList[0].lnx, paramPlaneList[0].lny, paramPlaneList[0].lnz);
    phi   = getPhi(paramPlaneList[0].lnx, paramPlaneList[0].lny, paramPlaneList[0].lnz);

    paramPlaneList[0].cosThetaLeft = cos(theta);
    paramPlaneList[0].sinThetaLeft = sin(theta);
    paramPlaneList[0].cosPhiLeft = cos(phi);
    paramPlaneList[0].sinPhiLeft = sin(phi);

    theta = getTheta(paramPlaneList[0].rnx, paramPlaneList[0].rny, paramPlaneList[0].rnz);
    phi   = getPhi(paramPlaneList[0].rnx, paramPlaneList[0].rny, paramPlaneList[0].rnz);

    paramPlaneList[0].cosThetaRight = cos(theta);
    paramPlaneList[0].sinThetaRight = sin(theta);
    paramPlaneList[0].cosPhiRight = cos(phi);
    paramPlaneList[0].sinPhiRight = sin(phi);

    paramPlaneList[0].fxLeftAv = 0;
    paramPlaneList[0].fyLeftAv = 0;
    paramPlaneList[0].fzLeftAv = 0;

    paramPlaneList[0].fxRightAv = 0;
    paramPlaneList[0].fyRightAv = 0;
    paramPlaneList[0].fzRightAv = 0;

//domainHi x plane
    paramPlaneList[1].x0 = domainHi[0];
    paramPlaneList[1].y0 = domainLo[1];
    paramPlaneList[1].z0 = domainLo[2];

    paramPlaneList[1].ux = 0;
    paramPlaneList[1].uy = 1;
    paramPlaneList[1].uz = 0;

    paramPlaneList[1].vx = 0;
    paramPlaneList[1].vy = 0;
    paramPlaneList[1].vz = 1;

    paramPlaneList[1].uTop = domainHi[1] - domainLo[1];
    paramPlaneList[1].vTop = domainHi[2] - domainLo[2];

    paramPlaneList[1].lnx = -1;
    paramPlaneList[1].lny = 0;
    paramPlaneList[1].lnz = 0;
    
    paramPlaneList[1].rnx = 1;
    paramPlaneList[1].rny = 0;
    paramPlaneList[1].rnz = 0;

    paramPlaneList[1].specularityLeft = 0;
    paramPlaneList[1].temperatureLeft = T_init[0];
    paramPlaneList[1].densityLeft[0] = 0;
    paramPlaneList[1].momentumConsLeft = 1;
    paramPlaneList[1].sourceLeft = 0;
    paramPlaneList[1].sinkLeft = 0;

    paramPlaneList[1].specularityRight = 0;
    paramPlaneList[1].temperatureRight = T_init[0];
    paramPlaneList[1].momentumConsRight = 1;
    paramPlaneList[1].densityRight[0] = 0;
    paramPlaneList[1].sourceRight = 0;
    paramPlaneList[1].sinkRight = 0;

    if(bc_vel_hi[0] == -1)
    {   
        paramPlaneList[1].periodicity = 1;
        paramPlaneList[1].porosityLeft = 1;
        paramPlaneList[1].porosityRight = 1;
    }
    else if (bc_vel_hi[0] == 2) 
    {
        paramPlaneList[1].periodicity = 0;
        paramPlaneList[1].porosityLeft = 0;
        paramPlaneList[1].porosityRight = 0;
        paramPlaneList[1].specularityLeft = 1;
        paramPlaneList[1].specularityRight = 1;

      //  paramPlaneList[1].x0 = domainHi[0] - sigma[0]/2.0;
    }
    else{
        paramPlaneList[1].periodicity = 0;
        paramPlaneList[1].porosityLeft = 0;
        paramPlaneList[1].porosityRight = 0;
        paramPlaneList[1].specularityLeft = 1;
        paramPlaneList[1].specularityRight = 1;

    }

    paramPlaneList[1].boundary = 2;

    theta = getTheta(paramPlaneList[1].lnx, paramPlaneList[1].lny, paramPlaneList[1].lnz);
    phi   = getPhi(paramPlaneList[1].lnx, paramPlaneList[1].lny, paramPlaneList[1].lnz);

    paramPlaneList[1].cosThetaLeft = cos(theta);
    paramPlaneList[1].sinThetaLeft = sin(theta);
    paramPlaneList[1].cosPhiLeft = cos(phi);
    paramPlaneList[1].sinPhiLeft = sin(phi);

    theta = getTheta(paramPlaneList[1].rnx, paramPlaneList[1].rny, paramPlaneList[1].rnz);
    phi   = getPhi(paramPlaneList[1].rnx, paramPlaneList[1].rny, paramPlaneList[1].rnz);

    paramPlaneList[1].cosThetaRight = cos(theta);
    paramPlaneList[1].sinThetaRight = sin(theta);
    paramPlaneList[1].cosPhiRight = cos(phi);
    paramPlaneList[1].sinPhiRight = sin(phi);

    paramPlaneList[1].fxLeftAv = 0;
    paramPlaneList[1].fyLeftAv = 0;
    paramPlaneList[1].fzLeftAv = 0;

    paramPlaneList[1].fxRightAv = 0;
    paramPlaneList[1].fyRightAv = 0;
    paramPlaneList[1].fzRightAv = 0;

//domainLo y plane
    paramPlaneList[2].x0 = domainLo[0];
    paramPlaneList[2].y0 = domainLo[1];
    paramPlaneList[2].z0 = domainLo[2];

    paramPlaneList[2].ux = 1;
    paramPlaneList[2].uy = 0;
    paramPlaneList[2].uz = 0;

    paramPlaneList[2].vx = 0;
    paramPlaneList[2].vy = 0;
    paramPlaneList[2].vz = 1;

    paramPlaneList[2].uTop = domainHi[0] - domainLo[0];
    paramPlaneList[2].vTop = domainHi[2] - domainLo[2];

    paramPlaneList[2].lnx = 0;
    paramPlaneList[2].lny = 1;
    paramPlaneList[2].lnz = 0;
    
    paramPlaneList[2].rnx = 0;
    paramPlaneList[2].rny = -1;
    paramPlaneList[2].rnz = 0;

    paramPlaneList[2].specularityLeft = 0;
    paramPlaneList[2].temperatureLeft = T_init[0];
    paramPlaneList[2].momentumConsLeft = 1;
    paramPlaneList[2].densityLeft[0] = 0;
    paramPlaneList[2].sourceLeft = 0;
    paramPlaneList[2].sinkLeft = 0;

    paramPlaneList[2].specularityRight = 0;
    paramPlaneList[2].temperatureRight = T_init[0];
    paramPlaneList[2].momentumConsRight = 1;
    paramPlaneList[2].densityRight[0] = 0;
    paramPlaneList[2].sourceRight = 0;
    paramPlaneList[2].sinkRight = 0;

    if(bc_vel_lo[1] == -1)
    {   
        paramPlaneList[2].periodicity = 1;
        paramPlaneList[2].porosityLeft = 1;
        paramPlaneList[2].porosityRight = 1;
    }
    else if (bc_vel_lo[1] == 2) 
    {
        paramPlaneList[2].periodicity = 0;
        paramPlaneList[2].porosityLeft = 0;
        paramPlaneList[2].porosityRight = 0;
        paramPlaneList[2].specularityLeft = 1;
        paramPlaneList[2].specularityRight = 1;

       // paramPlaneList[2].y0 = domainLo[1] + sigma[0]/2.0;
    }
    else{
        paramPlaneList[2].periodicity = 0;
        paramPlaneList[2].porosityLeft = 0;
        paramPlaneList[2].porosityRight = 0;
        paramPlaneList[2].specularityLeft = 1;
        paramPlaneList[2].specularityRight = 1;
    }

    paramPlaneList[2].boundary = 3;

    theta = getTheta(paramPlaneList[2].lnx, paramPlaneList[2].lny, paramPlaneList[2].lnz);
    phi   = getPhi(paramPlaneList[2].lnx, paramPlaneList[2].lny, paramPlaneList[2].lnz);

    paramPlaneList[2].cosThetaLeft = cos(theta);
    paramPlaneList[2].sinThetaLeft = sin(theta);
    paramPlaneList[2].cosPhiLeft = cos(phi);
    paramPlaneList[2].sinPhiLeft = sin(phi);

    theta = getTheta(paramPlaneList[2].rnx, paramPlaneList[2].rny, paramPlaneList[2].rnz);
    phi   = getPhi(paramPlaneList[2].rnx, paramPlaneList[2].rny, paramPlaneList[2].rnz);

    paramPlaneList[2].cosThetaRight = cos(theta);
    paramPlaneList[2].sinThetaRight = sin(theta);
    paramPlaneList[2].cosPhiRight = cos(phi);
    paramPlaneList[2].sinPhiRight = sin(phi);

    paramPlaneList[2].fxLeftAv = 0;
    paramPlaneList[2].fyLeftAv = 0;
    paramPlaneList[2].fzLeftAv = 0;

    paramPlaneList[2].fxRightAv = 0;
    paramPlaneList[2].fyRightAv = 0;
    paramPlaneList[2].fzRightAv = 0;

//domianHi y plane
    paramPlaneList[3].x0 = domainLo[0];
    paramPlaneList[3].y0 = domainHi[1];
    paramPlaneList[3].z0 = domainLo[2];

    paramPlaneList[3].ux = 1;
    paramPlaneList[3].uy = 0;
    paramPlaneList[3].uz = 0;

    paramPlaneList[3].vx = 0;
    paramPlaneList[3].vy = 0;
    paramPlaneList[3].vz = 1;

    paramPlaneList[3].uTop = domainHi[0] - domainLo[0];
    paramPlaneList[3].vTop = domainHi[2] - domainLo[2];

    paramPlaneList[3].lnx = 0;
    paramPlaneList[3].lny = 1;
    paramPlaneList[3].lnz = 0;
    
    paramPlaneList[3].rnx = 0;
    paramPlaneList[3].rny = -1;
    paramPlaneList[3].rnz = 0;

    paramPlaneList[3].specularityLeft = 0;
    paramPlaneList[3].temperatureLeft = T_init[0];
    paramPlaneList[3].momentumConsLeft = 1;
    paramPlaneList[3].densityLeft[0] = 0;
    paramPlaneList[3].sourceLeft = 0;
    paramPlaneList[3].sinkLeft = 0;

    paramPlaneList[3].specularityRight = 0;
    paramPlaneList[3].temperatureRight = T_init[0];
    paramPlaneList[3].momentumConsRight = 1;
    paramPlaneList[3].densityRight[0] = 0;
    paramPlaneList[3].sourceRight = 0;
    paramPlaneList[3].sinkRight = 0;


    if(bc_vel_hi[1] == -1)
    {   
        paramPlaneList[3].periodicity = 1;
        paramPlaneList[3].porosityLeft = 1;
        paramPlaneList[3].porosityRight = 1;
    }
    else if (bc_vel_hi[1] == 2) 
    {
        paramPlaneList[3].periodicity = 0;
        paramPlaneList[3].porosityLeft = 0;
        paramPlaneList[3].porosityRight = 0;
        paramPlaneList[3].specularityLeft = 1;
        paramPlaneList[3].specularityRight = 1;

       // paramPlaneList[3].y0 = domainHi[1] - sigma[0]/2.0;
    }
    else{
        paramPlaneList[3].periodicity = 0;
        paramPlaneList[3].porosityLeft = 0;
        paramPlaneList[3].porosityRight = 0;
        paramPlaneList[3].specularityLeft = 1;
        paramPlaneList[3].specularityRight = 1;

    }

    paramPlaneList[3].boundary = 4;

    theta = getTheta(paramPlaneList[3].lnx, paramPlaneList[3].lny, paramPlaneList[3].lnz);
    phi   = getPhi(paramPlaneList[3].lnx, paramPlaneList[3].lny, paramPlaneList[3].lnz);

    paramPlaneList[3].cosThetaLeft = cos(theta);
    paramPlaneList[3].sinThetaLeft = sin(theta);
    paramPlaneList[3].cosPhiLeft = cos(phi);
    paramPlaneList[3].sinPhiLeft = sin(phi);

    theta = getTheta(paramPlaneList[3].rnx, paramPlaneList[3].rny, paramPlaneList[3].rnz);
    phi   = getPhi(paramPlaneList[3].rnx, paramPlaneList[3].rny, paramPlaneList[3].rnz);

    paramPlaneList[3].cosThetaRight = cos(theta);
    paramPlaneList[3].sinThetaRight = sin(theta);
    paramPlaneList[3].cosPhiRight = cos(phi);
    paramPlaneList[3].sinPhiRight = sin(phi);

    paramPlaneList[3].fxLeftAv = 0;
    paramPlaneList[3].fyLeftAv = 0;
    paramPlaneList[3].fzLeftAv = 0;

    paramPlaneList[3].fxRightAv = 0;
    paramPlaneList[3].fyRightAv = 0;
    paramPlaneList[3].fzRightAv = 0;

//domainLo z plane
    paramPlaneList[4].x0 = domainLo[0];
    paramPlaneList[4].y0 = domainLo[1];
    paramPlaneList[4].z0 = domainLo[2];

    paramPlaneList[4].ux = 1;
    paramPlaneList[4].uy = 0;
    paramPlaneList[4].uz = 0;

    paramPlaneList[4].vx = 0;
    paramPlaneList[4].vy = 1;
    paramPlaneList[4].vz = 0;

    paramPlaneList[4].uTop = domainHi[0] - domainLo[0];
    paramPlaneList[4].vTop = domainHi[1] - domainLo[1];

    paramPlaneList[4].lnx = 0;
    paramPlaneList[4].lny = 0;
    paramPlaneList[4].lnz = 1;
    
    paramPlaneList[4].rnx = 0;
    paramPlaneList[4].rny = 0;
    paramPlaneList[4].rnz = -1;

    paramPlaneList[4].specularityLeft = 0;
    paramPlaneList[4].temperatureLeft = T_init[0];
    paramPlaneList[4].momentumConsLeft = 1;
    paramPlaneList[4].sourceLeft = 0;
    paramPlaneList[4].densityLeft[0] = 0;
    paramPlaneList[4].sinkLeft = 0;

    paramPlaneList[4].specularityRight = 0;
    paramPlaneList[4].temperatureRight = T_init[0];
    paramPlaneList[4].momentumConsRight = 1;
    paramPlaneList[4].densityRight[0] = 0;
    paramPlaneList[4].sourceRight = 0;
    paramPlaneList[4].sinkRight = 0;

    if(bc_vel_lo[2] == -1)
    {   
        paramPlaneList[4].periodicity = 1;
        paramPlaneList[4].porosityLeft = 1;
        paramPlaneList[4].porosityRight = 1;
    }
    else if (bc_vel_lo[2] == 2) 
    {
        paramPlaneList[4].periodicity = 0;
        paramPlaneList[4].porosityLeft = 0;
        paramPlaneList[4].porosityRight = 0;
        paramPlaneList[4].specularityLeft = 1;
        paramPlaneList[4].specularityRight = 1;

       // paramPlaneList[4].z0 = domainLo[2] + sigma[0]/2.0;
    }
    else{
        paramPlaneList[4].periodicity = 0;
        paramPlaneList[4].porosityLeft = 0;
        paramPlaneList[4].porosityRight = 0;
        paramPlaneList[4].specularityLeft = 1;
        paramPlaneList[4].specularityRight = 1;


    }

    paramPlaneList[4].boundary = 5;

    theta = getTheta(paramPlaneList[4].lnx, paramPlaneList[4].lny, paramPlaneList[4].lnz);
    phi   = getPhi(paramPlaneList[4].lnx, paramPlaneList[4].lny, paramPlaneList[4].lnz);

    paramPlaneList[4].cosThetaLeft = cos(theta);
    paramPlaneList[4].sinThetaLeft = sin(theta);
    paramPlaneList[4].cosPhiLeft = cos(phi);
    paramPlaneList[4].sinPhiLeft = sin(phi);

    theta = getTheta(paramPlaneList[4].rnx, paramPlaneList[4].rny, paramPlaneList[4].rnz);
    phi   = getPhi(paramPlaneList[4].rnx, paramPlaneList[4].rny, paramPlaneList[4].rnz);

    paramPlaneList[4].cosThetaRight = cos(theta);
    paramPlaneList[4].sinThetaRight = sin(theta);
    paramPlaneList[4].cosPhiRight = cos(phi);
    paramPlaneList[4].sinPhiRight = sin(phi);

    paramPlaneList[4].fxLeftAv = 0;
    paramPlaneList[4].fyLeftAv = 0;
    paramPlaneList[4].fzLeftAv = 0;

    paramPlaneList[4].fxRightAv = 0;
    paramPlaneList[4].fyRightAv = 0;
    paramPlaneList[4].fzRightAv = 0;

//domianHi z plane
    paramPlaneList[5].x0 = domainLo[0];
    paramPlaneList[5].y0 = domainLo[1];
    paramPlaneList[5].z0 = domainHi[2];

    paramPlaneList[5].ux = 1;
    paramPlaneList[5].uy = 0;
    paramPlaneList[5].uz = 0;

    paramPlaneList[5].vx = 0;
    paramPlaneList[5].vy = 1;
    paramPlaneList[5].vz = 0;

    paramPlaneList[5].uTop = domainHi[0] - domainLo[0];
    paramPlaneList[5].vTop = domainHi[1] - domainLo[1];

    paramPlaneList[5].lnx = 0;
    paramPlaneList[5].lny = 0;
    paramPlaneList[5].lnz = 1;
    
    paramPlaneList[5].rnx = 0;
    paramPlaneList[5].rny = 0;
    paramPlaneList[5].rnz = -1;

    paramPlaneList[5].specularityLeft = 0;
    paramPlaneList[5].temperatureLeft = T_init[0];
    paramPlaneList[5].momentumConsLeft = 1;
    paramPlaneList[5].densityLeft[0] = 0;
    paramPlaneList[5].sourceLeft = 0;
    paramPlaneList[5].sinkLeft = 0;

    paramPlaneList[5].specularityRight = 0;
    paramPlaneList[5].temperatureRight = T_init[0];
    paramPlaneList[5].momentumConsRight = 1;
    paramPlaneList[5].densityRight[0] = 0;
    paramPlaneList[5].sourceRight = 0;
    paramPlaneList[5].sinkRight = 0;


    paramPlaneList[5].velx=0;
    paramPlaneList[5].vely=0;
    paramPlaneList[5].velz=0;

    paramPlaneList[5].c0=0;
    paramPlaneList[5].resomg=0;
    paramPlaneList[5].agraph=0;
    paramPlaneList[5].bgraph=0;
    paramPlaneList[5].a0graph=0;
    paramPlaneList[5].b0graph=0;
    paramPlaneList[5].coltime=0;

    paramPlaneList[5].besslist;
    paramPlaneList[5].dbesslist;

    if(bc_vel_hi[2] == -1)
    {   
        paramPlaneList[5].periodicity = 1;
        paramPlaneList[5].porosityLeft = 1;
        paramPlaneList[5].porosityRight = 1;
    }
    else if (bc_vel_hi[2] == 2) 
    {
        paramPlaneList[5].periodicity = 0;
        paramPlaneList[5].porosityLeft = 0;
        paramPlaneList[5].porosityRight = 0;
        paramPlaneList[5].specularityLeft = 1;
        paramPlaneList[5].specularityRight = 1;

       // paramPlaneList[5].z0 = domainHi[2] - sigma[0]/2.0;
    }
    else{
        paramPlaneList[5].periodicity = 0;
        paramPlaneList[5].porosityLeft = 0;
        paramPlaneList[5].porosityRight = 0;
        paramPlaneList[5].specularityLeft = 1;
        paramPlaneList[5].specularityRight = 1;

    }

    paramPlaneList[5].boundary = 6;

    theta = getTheta(paramPlaneList[5].lnx, paramPlaneList[5].lny, paramPlaneList[5].lnz);
    phi   = getPhi(paramPlaneList[5].lnx, paramPlaneList[5].lny, paramPlaneList[5].lnz);

    paramPlaneList[5].cosThetaLeft = cos(theta);
    paramPlaneList[5].sinThetaLeft = sin(theta);
    paramPlaneList[5].cosPhiLeft = cos(phi);
    paramPlaneList[5].sinPhiLeft = sin(phi);

    theta = getTheta(paramPlaneList[5].rnx, paramPlaneList[5].rny, paramPlaneList[5].rnz);
    phi   = getPhi(paramPlaneList[5].rnx, paramPlaneList[5].rny, paramPlaneList[5].rnz);

    paramPlaneList[5].cosThetaRight = cos(theta);
    paramPlaneList[5].sinThetaRight = sin(theta);
    paramPlaneList[5].cosPhiRight = cos(phi);
    paramPlaneList[5].sinPhiRight = sin(phi);

    paramPlaneList[5].fxLeftAv = 0;
    paramPlaneList[5].fyLeftAv = 0;
    paramPlaneList[5].fzLeftAv = 0;

    paramPlaneList[5].fxRightAv = 0;
    paramPlaneList[5].fyRightAv = 0;
    paramPlaneList[5].fzRightAv = 0;

    std::ifstream planeFile("paramplanes.dat");
    int fileCount = 0;
    if(planeFile.good())
    {
        planeFile >> fileCount;
    }
    
    int totalCount = 6+fileCount;

    for(int i=6; i<totalCount; i++)
    {

        planeFile >> paramPlaneList[i].x0;

        Print() << "surface " << i << " xo " << paramPlaneList[i].x0 << "\n";

        planeFile >> paramPlaneList[i].y0;
        planeFile >> paramPlaneList[i].z0;

        planeFile >> paramPlaneList[i].ux;
        planeFile >> paramPlaneList[i].uy;
        planeFile >> paramPlaneList[i].uz;

        planeFile >> paramPlaneList[i].vx;
        planeFile >> paramPlaneList[i].vy;
        planeFile >> paramPlaneList[i].vz;

        planeFile >> paramPlaneList[i].uTop;
        planeFile >> paramPlaneList[i].vTop;

        planeFile >> paramPlaneList[i].rnx;
        planeFile >> paramPlaneList[i].rny;
        planeFile >> paramPlaneList[i].rnz;

        planeFile >> paramPlaneList[i].lnx;
        planeFile >> paramPlaneList[i].lny;
        planeFile >> paramPlaneList[i].lnz;
        
        planeFile >> paramPlaneList[i].porosityRight;
        planeFile >> paramPlaneList[i].specularityRight;
        planeFile >> paramPlaneList[i].temperatureRight;
        for(int j=0; j < nspecies; j++)
        {
            planeFile >> paramPlaneList[i].densityRight[j];
        }
        planeFile >> paramPlaneList[i].sourceRight;
        planeFile >> paramPlaneList[i].sinkRight;
        planeFile >> paramPlaneList[i].momentumConsRight;

        planeFile >> paramPlaneList[i].porosityLeft;
        planeFile >> paramPlaneList[i].specularityLeft;
        planeFile >> paramPlaneList[i].temperatureLeft;
        for(int j=0; j < nspecies; j++)
        {
           planeFile >> paramPlaneList[i].densityLeft[j];
        }
        planeFile >> paramPlaneList[i].sourceLeft;
        planeFile >> paramPlaneList[i].sinkLeft;
        planeFile >> paramPlaneList[i].momentumConsLeft;

        planeFile >> paramPlaneList[i].periodicity;

        planeFile >> paramPlaneList[i].area;

        paramPlaneList[i].boundary = i+1;
            
        theta = getTheta(paramPlaneList[i].lnx, paramPlaneList[i].lny, paramPlaneList[i].lnz);
        phi   = getPhi(paramPlaneList[i].lnx, paramPlaneList[i].lny, paramPlaneList[i].lnz);

        paramPlaneList[i].cosThetaLeft = cos(theta);
        paramPlaneList[i].sinThetaLeft = sin(theta);
        paramPlaneList[i].cosPhiLeft = cos(phi);
        paramPlaneList[i].sinPhiLeft = sin(phi);

        theta = getTheta(paramPlaneList[i].rnx, paramPlaneList[i].rny, paramPlaneList[i].rnz);
        phi   = getPhi(paramPlaneList[i].rnx, paramPlaneList[i].rny, paramPlaneList[i].rnz);

        paramPlaneList[i].cosThetaRight = cos(theta);
        paramPlaneList[i].sinThetaRight = sin(theta);
        paramPlaneList[i].cosPhiRight = cos(phi);
        paramPlaneList[i].sinPhiRight = sin(phi);

        paramPlaneList[i].fxLeftAv = 0;
        paramPlaneList[i].fyLeftAv = 0;
        paramPlaneList[i].fzLeftAv = 0;

        paramPlaneList[i].fxRightAv = 0;
        paramPlaneList[i].fyRightAv = 0;
        paramPlaneList[i].fzRightAv = 0;        

    }

    planeFile.close();
}
