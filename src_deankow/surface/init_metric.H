#ifndef INIT_METRIC_H
#define INIT_METRIC_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Print.H>

void init_dk_metric(int i, int j, int k,
               amrex::Array4<amrex::Real> const& gmet,
               amrex::Array4<amrex::Real> const& sqrmet,
               amrex::Array4<amrex::Real> const& detg,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo)
{

    constexpr double  TWOPI = 2.*3.14159265358979323846264338327950288;
    constexpr double  PI = 3.14159265358979323846264338327950288;

    amrex::Real xloc = prob_lo[0] + (i+0.5) * dx[0] ;
    amrex::Real yloc = prob_lo[1] + (j+0.5) * dx[1] ;

#if 1

    amrex::Real amp = 3.;
    amrex::Real sinx,siny,cosx,cosy,det,sig11,sig12,sig21,sig22;

    sinx = std::sin(xloc);
    siny = std::sin(yloc);
    cosx = std::cos(xloc);
    cosy = std::cos(yloc);

    amrex::Real p,q;
    p = amp*cosx*siny;
    q = amp*cosy*sinx;

    det = 1. + p*p + q*q;

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

/*
    det = 1. + amp*amp*(cosx*cosx*siny*siny+sinx*sinx*cosy*cosy);

    detg(i,j,k,0) = det;
    detg(i,j,k,1) = std::sqrt(det);
    detg(i,j,k,2) = std::sqrt(detg(i,j,k,1));

    gmet(i,j,k,0) = (1.+ amp*amp*cosy*cosy*sinx*sinx)/det;
    gmet(i,j,k,1) = -(amp*amp*cosy*cosx*sinx*siny)/det;
    gmet(i,j,k,2) = -(amp*amp*cosy*cosx*sinx*siny)/det;
    gmet(i,j,k,3) = (1.+ amp*amp*cosx*cosx*siny*siny)/det;


//    amrex::Real detm1 = det - 1.;
//    amrex::Real cfac = (1. - 1./std::sqrt(1+detm1))/detm1;
    amrex::Real cfac = 1./(detg(i,j,k,0)+detg(i,j,k,1));

    sqrmet(i,j,k,0) = 1. - cfac * amp*amp * cosx*cosx*siny*siny;
    sqrmet(i,j,k,1) =  - cfac * amp*amp * cosy*cosx*sinx*siny;
    sqrmet(i,j,k,2) =  - cfac * amp*amp * cosy*cosx*sinx*siny;
    sqrmet(i,j,k,3) = 1. - cfac * amp*amp * cosy*cosy*sinx*sinx;
*/
#endif
#if 0
    amrex::Real amp = 4.;
    amrex::Real sinx,siny,cosx,cosy,det,sig11,sig12,sig21,sig22;

    sinx = std::sin(xloc);
    siny = std::sin(yloc);
    cosx = std::cos(xloc);
    cosy = std::cos(yloc);

    det = 1. + 4.* amp*amp*(cosx*cosx*siny*siny+sinx*sinx*cosy*cosy)*sinx*sinx*siny*siny;

    detg(i,j,k,0) = det;
    detg(i,j,k,1) = std::sqrt(det);
    detg(i,j,k,2) = std::sqrt(detg(i,j,k,1));

    gmet(i,j,k,0) = (1.+ 4.*amp*amp*cosy*cosy*siny*siny*sinx*sinx*sinx*sinx)/det;
    gmet(i,j,k,1) = -4.*(amp*amp*cosy*cosx*sinx*sinx*sinx*siny*siny*siny)/det;
    gmet(i,j,k,2) = -4.*(amp*amp*cosy*cosx*sinx*sinx*sinx*siny*siny*siny)/det;
    gmet(i,j,k,3) = (1.+ 4.*amp*amp*cosx*cosx*sinx*sinx*siny*siny*siny*siny)/det;

    sqrmet(i,j,k,0) = std::sqrt(gmet(i,j,k,0));
    sqrmet(i,j,k,1) =  0.;
    sqrmet(i,j,k,2) =  gmet(i,j,k,2)/sqrmet(i,j,k,0);
    sqrmet(i,j,k,3) = std::sqrt(gmet(i,j,k,3)-gmet(i,j,k,2)*gmet(i,j,k,2)/gmet(i,j,k,0));

#endif
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
