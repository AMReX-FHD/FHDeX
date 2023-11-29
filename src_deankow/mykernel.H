#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void init_phi (int i, int j, int k,
               amrex::Array4<amrex::Real> const& phi,
               GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
               amrex::Real npts_scale,int Ncomp)
{
    using amrex::Real;;

    Real x = prob_lo[0] + (i+Real(0.5)) * dx[0];
    Real y = prob_lo[1] + (j+Real(0.5)) * dx[1];
#if (AMREX_SPACEDIM > 2)
    Real z = prob_lo[2] + (k+Real(0.5)) * dx[2];
    Real r2 = ((x-Real(0.25))*(x-Real(0.25))+(y-Real(0.25))*(y-Real(0.25))+(z-Real(0.25))*(z-Real(0.25)))/Real(0.01);
#else
    Real z = Real(0.);
    Real r2 = ((x-Real(0.25))*(x-Real(0.25))+(y-Real(0.25))*(y-Real(0.25)))/Real(0.01);
#endif
    phi(i,j,k,0) = Real(1.) + std::exp(-r2);
    if( x > 0.5){
    phi(i,j,k,0) = 1.5*npts_scale;
    } else {
    phi(i,j,k,0) = 0.5*npts_scale;
    }
    if(i==4 && j ==4)
    {
       phi(i,j,k,0) = 1./(dx[0]*dx[1]);
    }
    else
    {
       phi(i,j,k,0) = 0;;
    }
    if(Ncomp == 2){
        phi(i,j,k,1) = phi(i,j,k,0);
    }
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_x (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxx,
                     amrex::Array4<amrex::Real> const& stochfluxx,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dxinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi,int Ncomp)
{
    if (lo == dom_lo &&
        (bc_lo == BCType::foextrap ||
         bc_lo == BCType::ext_dir))
    {
        if(i == lo)
        {
            fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv / 0.5;
        }
        else
        {
            fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv;
        }
    }
    else if (hi == dom_hi &&
             (bc_hi == BCType::foextrap ||
              bc_hi == BCType::ext_dir))
    {
        if(i == hi+1)
        {
            fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv / 0.5;
        }
        else
        {
            fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv;
        }
    }
    else
    {
//	amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i-1,j,k))/(phi(i,j,k)+phi(i-1,j,k)+1.e-16);
        Real phip = std::max(phi(i,j,k,1),0.);
        Real phim = std::max(phi(i-1,j,k,1),0.);
	amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
        fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv
	  + stochfluxx(i,j,k) * std::sqrt(phiavg);
        fluxx(i,j,k,1) = 0.5*(phi(i,j,k,1)-phi(i-1,j,k,1)) * dxinv;
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_y (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxy,
                     amrex::Array4<amrex::Real> const& stochfluxy,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dyinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi, int Ncomp)
{
    if (lo == dom_lo &&
        (bc_lo == BCType::foextrap ||
         bc_lo == BCType::ext_dir))
    {
        if(j == lo)
        {
            fluxy(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j-1,k,0)) * dyinv / 0.5;
        }
        else
        {
            fluxy(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j-1,k,0)) * dyinv;
        }
    }
    else if (hi == dom_hi &&
             (bc_hi == BCType::foextrap ||
              bc_hi == BCType::ext_dir))
    {
        if(j == hi+1)
        {
            fluxy(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j-1,k,0)) * dyinv / 0.5;
        }
        else
        {
            fluxy(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j-1,k,0)) * dyinv;
        }
    }
    else
    {
//	amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i,j-1,k))/(phi(i,j,k)+phi(i,j-1,k)+1.e-16);
        Real phip = std::max(phi(i,j,k,1),0.);
        Real phim = std::max(phi(i,j-1,k,1),0.);
	amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
        fluxy(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j-1,k,0)) * dyinv
	+ stochfluxy(i,j,k)* std::sqrt(phiavg);
        fluxy(i,j,k,1) = 0.5*(phi(i,j,k,1)-phi(i,j-1,k,1)) * dyinv;
    }
}


#if (AMREX_SPACEDIM > 2)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_z (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxz,
                     amrex::Array4<amrex::Real> const& stochfluxz,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dzinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi)
{
    if (lo == dom_lo &&
        (bc_lo == BCType::foextrap ||
         bc_lo == BCType::ext_dir))
    {
        if(k == lo)
        {
            fluxz(i,j,k) = 0.5*(phi(i,j,k)-phi(i,j,k-1)) * dzinv / 0.5;
        }
        else
        {
            fluxz(i,j,k) = 0.5*(phi(i,j,k)-phi(i,j,k-1)) * dzinv;
        }
    }
    else if (hi == dom_hi &&
             (bc_hi == BCType::foextrap ||
              bc_hi == BCType::ext_dir))
    {
        if(k == hi+1)
        {
            fluxz(i,j,k) = 0.5*(phi(i,j,k)-phi(i,j,k-1)) * dzinv / 0.5;
        }
        else
        {
            fluxz(i,j,k) = 0.5*(phi(i,j,k)-phi(i,j,k-1)) * dzinv;
        }
    }
    else
    {
//	amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i,j,k-1))/(phi(i,j,k)+phi(i,j,k-1)+1.e-16);
        Real phip = std::max(phi(i,j,k),0.);
        Real phim = std::max(phi(i,j,k-1),0.);
	amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
        fluxz(i,j,k) = 0.5*(phi(i,j,k)-phi(i,j,k-1)) * dzinv
	 +stochfluxz(i,j,k) * std::sqrt(phiavg);
    }
}
#endif

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void update_phi (int i, int j, int k,
                 amrex::Array4<amrex::Real const> const& phiold,
                 amrex::Array4<amrex::Real      > const& phinew,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fluxx,
                              amrex::Array4<amrex::Real const> const& fluxy,
                              amrex::Array4<amrex::Real const> const& fluxz),
                 amrex::Real dt,
                 AMREX_D_DECL(amrex::Real dxinv,
                              amrex::Real dyinv,
                              amrex::Real dzinv),
                 int Ncomp)
{
  for (int n = 0; n < Ncomp; n++){
    phinew(i,j,k,n) = phiold(i,j,k,n)
        +  dt * dxinv * (fluxx(i+1,j  ,k,n  ) - fluxx(i,j,k,n))
        +  dt * dyinv * (fluxy(i  ,j+1,k,n  ) - fluxy(i,j,k,n))
#if (AMREX_SPACEDIM > 2)
        +  dt * dzinv * (fluxz(i  ,j  ,k+1,n) - fluxz(i,j,k,n));
#else
        ;
#endif
  }
}

#endif
