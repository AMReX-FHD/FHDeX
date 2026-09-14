#ifndef UPDATE_METRIC_H
#define UPDATE_METRIC_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Print.H>

void new_dk_metric(int i, int j, int k,
               amrex::Array4<amrex::Real> const& gmet,
               amrex::Array4<amrex::Real> const& sqrmet,
               amrex::Array4<amrex::Real> const& detg,
               amrex::Array4<amrex::Real> const& newgmet,
               amrex::Array4<amrex::Real> const& newsqrmet,
               amrex::Array4<amrex::Real> const& newdetg,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
               amrex::Real time)
{

    constexpr double  TWOPI = 2.*3.14159265358979323846264338327950288;
    constexpr double  PI = 3.14159265358979323846264338327950288;

    amrex::Real xloc = prob_lo[0] + (i+0.5) * dx[0] ;
    amrex::Real yloc = prob_lo[1] + (j+0.5) * dx[1] ;

#if 0

    amrex::Real amp = 3.;
    amrex::Real sinx,siny,cosx,cosy,det,sig11,sig12,sig21,sig22;

    sinx = std::sin(xloc);
    siny = std::sin(yloc);
    cosx = std::cos(xloc);
    cosy = std::cos(yloc);

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

#endif

    amrex::Real amp = 4.;
    amrex::Real tscale = 0.3;
    amrex::Real sinx,siny,cosx,cosy,det,sig11,sig12,sig21,sig22;

    sinx = std::sin(xloc);
    siny = std::sin(yloc);
    cosx = std::cos(xloc);
    cosy = std::cos(yloc);

    amrex::Real p,q,p0,q0, p1, q1, aoft;

    detg(i,j,k,0) = newdetg(i,j,k,0);
    detg(i,j,k,1) = newdetg(i,j,k,1);
    detg(i,j,k,2) = newdetg(i,j,k,2);

    gmet(i,j,k,0) = newgmet(i,j,k,0);
    gmet(i,j,k,1) = newgmet(i,j,k,1);
    gmet(i,j,k,2) = newgmet(i,j,k,2);
    gmet(i,j,k,3) = newgmet(i,j,k,3);

    sqrmet(i,j,k,0) = newsqrmet(i,j,k,0);
    sqrmet(i,j,k,1) = newsqrmet(i,j,k,1);
    sqrmet(i,j,k,2) = newsqrmet(i,j,k,2);
    sqrmet(i,j,k,3) = newsqrmet(i,j,k,3);

    p0 = 2.*amp*cosx*sinx*siny*siny;
    q0 = 2.*amp*cosy*sinx*sinx*siny;

    p1 = -2.*amp*cosx*sinx*cosy*cosy;
    q1 = -2.*amp*cosy*cosx*cosx*siny;

    aoft = std::cos(0.5*PI * std::min(time,tscale) / tscale)*std::cos(0.5*PI * std::min(time,tscale) / tscale);

    p = aoft*p0 + (1.-aoft)*p1;
    q = aoft*q0 + (1.-aoft)*q1;


    det = 1. + p*p + q*q;

    newdetg(i,j,k,0) = det;
    newdetg(i,j,k,1) = std::sqrt(det);
    newdetg(i,j,k,2) = std::sqrt(detg(i,j,k,1));

    newgmet(i,j,k,0) = (1.+ q*q)/det;
    newgmet(i,j,k,1) = -p*q/det;
    newgmet(i,j,k,2) = -p*q/det;
    newgmet(i,j,k,3) = (1.+ p*p)/det;

    amrex::Real factor = 1./(detg(i,j,k,1)*(1.+detg(i,j,k,1)));

    newsqrmet(i,j,k,0) = (1.+q*q+detg(i,j,k,1))*factor;
    newsqrmet(i,j,k,1) =  - p*q*factor;
    newsqrmet(i,j,k,2) =  - p*q*factor;
    newsqrmet(i,j,k,3) = (1.+p*p+detg(i,j,k,1))*factor;

}

#endif
