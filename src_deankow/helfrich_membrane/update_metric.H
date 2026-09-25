#ifndef UPDATE_METRIC_H
#define UPDATE_METRIC_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Print.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void new_dk_metric(int i, int j, int k,
               amrex::Array4<amrex::Real> const& gmet,
               amrex::Array4<amrex::Real> const& sqrmet,
               amrex::Array4<amrex::Real> const& detg,
               amrex::Array4<amrex::Real> const& newgmet,
               amrex::Array4<amrex::Real> const& newsqrmet,
               amrex::Array4<amrex::Real> const& newdetg,
               amrex::Array4<amrex::Real> const& dhdx,
               amrex::Array4<amrex::Real> const& dhdy,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
               amrex::Real time)
{
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

//    amrex::Print() << " new metric i j " << i << " " <<  j << " " << dhdx(i,j,k,0) << " " << dhdy(i,j,k,0) << std::endl;

    amrex::Real p = dhdx(i,j,k,0);
    amrex::Real q = dhdy(i,j,k,0);

    amrex::Real det = 1. + p*p + q*q;

    newdetg(i,j,k,0) = det;
    newdetg(i,j,k,1) = std::sqrt(det);
    newdetg(i,j,k,2) = std::sqrt(newdetg(i,j,k,1));

    newgmet(i,j,k,0) = (1.+ q*q)/det;
    newgmet(i,j,k,1) = -p*q/det;
    newgmet(i,j,k,2) = -p*q/det;
    newgmet(i,j,k,3) = (1.+ p*p)/det;

    // The inverse square root of g = I + v v^T (v = (p,q)) is
    // I - v v^T/(s(s+1)) with s = sqrt(det).  s must be the one that belongs
    // to the p,q used in the same expressions, i.e. the NEW metric -- detg
    // holds the previous time level at this point.
    amrex::Real factor = 1./(newdetg(i,j,k,1)*(1.+newdetg(i,j,k,1)));

    newsqrmet(i,j,k,0) = (1.+q*q+newdetg(i,j,k,1))*factor;
    newsqrmet(i,j,k,1) =  - p*q*factor;
    newsqrmet(i,j,k,2) =  - p*q*factor;
    newsqrmet(i,j,k,3) = (1.+p*p+newdetg(i,j,k,1))*factor;

}

#endif
