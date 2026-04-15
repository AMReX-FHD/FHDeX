#ifndef INIT_METRIC_H
#define INIT_METRIC_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Print.H>

void init_dk_metric(int i, int j, int k,
               amrex::Array4<amrex::Real> const& gmet,
               amrex::Array4<amrex::Real> const& sqrmet,
               amrex::Array4<amrex::Real> const& detg,
               amrex::Array4<amrex::Real> const& dhdx,
               amrex::Array4<amrex::Real> const& dhdy,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo)
{

    amrex::Print() << " init metric i j " << i << " " <<  j << " " << dhdx(i,j,k,0) << " " << dhdy(i,j,k,0) << std::endl; 
    amrex::Real p = dhdx(i,j,k,0);
    amrex::Real q = dhdy(i,j,k,0);

    amrex::Real det = 1. + p*p + q*q;

    detg(i,j,k,0) = det;
    detg(i,j,k,1) = std::sqrt(det);
    detg(i,j,k,2) = std::sqrt(detg(i,j,k,1));

    gmet(i,j,k,0) = (1.+ q*q)/det;
    gmet(i,j,k,1) = -p*q/det;
    gmet(i,j,k,2) = -p*q/det;
    gmet(i,j,k,3) = (1.+ p*p)/det;

    amrex::Real factor = 1./(detg(i,j,k,1)*(1.+detg(i,j,k,1)));

    sqrmet(i,j,k,0) = (1.+q*q+detg(i,j,k,1))*factor;
    sqrmet(i,j,k,1) =  - p*q*factor;
    sqrmet(i,j,k,2) =  - p*q*factor;
    sqrmet(i,j,k,3) = (1.+p*p+detg(i,j,k,1))*factor;


#if 0
    amrex::Real amp = 4.0;

    amrex::Real caps, capsy, capr, caprx;

    caps  = 0.;
    capsy = 0.;

    capr  = std::pow(std::cos(0.5*(xloc-PI)),4);
    caprx = -2.*std::sin(0.5*(xloc-PI))*std::pow(std::cos(0.5*(xloc-PI)),3);

    for (int iterm  = 0; iterm < 3; iterm++)
    {
       amrex::Real shift = PI/3. + 2.*iterm*PI/3.;
       caps += std::pow(std::cos(0.5*(yloc-shift)),16);
       capsy += -8*std::sin(0.5*(yloc-shift))*std::pow(std::cos(0.5*(yloc-shift)),15);
    }


    amrex::Real p,q;
    p = amp*caps*caprx;
    q = amp*capsy*capr;

    amrex::Real det = 1. + p*p + q*q;

    detg(i,j,k,0) = det;
    detg(i,j,k,1) = std::sqrt(det);
    detg(i,j,k,2) = std::sqrt(detg(i,j,k,1));

    gmet(i,j,k,0) = (1.+ q*q)/det;
    gmet(i,j,k,1) = -p*q/det;
    gmet(i,j,k,2) = -p*q/det;
    gmet(i,j,k,3) = (1.+ p*p)/det;

    amrex::Real factor = 1./(detg(i,j,k,1)*(1.+detg(i,j,k,1)));

    sqrmet(i,j,k,0) = (1.+q*q+detg(i,j,k,1))*factor;
    sqrmet(i,j,k,1) =  - p*q*factor;
    sqrmet(i,j,k,2) =  - p*q*factor;
    sqrmet(i,j,k,3) = (1.+p*p+detg(i,j,k,1))*factor;

#endif

#if 0
    amrex::Print() << " in metric " << i << " " << j << std::endl;
    amrex::Print() << "det " << detg(i,j,k,0) << " " << detg(i,j,k,1) << " " << detg(i,j,k,2) << std::endl; 
    for(int n=0; n<4; n++){
    amrex::Print() << "gmet " << n  << " " << gmet(i,j,k,n) << std::endl;
    }
    for(int n=0; n<4; n++){
    amrex::Print() << "sqrgmet " << n  << " " << sqrmet(i,j,k,n) << std::endl;
    }
#endif

}


#endif
