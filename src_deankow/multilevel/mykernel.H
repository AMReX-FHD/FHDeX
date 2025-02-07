#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void init_phi (int i, int j, int k,
               amrex::Array4<amrex::Real> const& phi,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
               amrex::Real /*npts_scale*/, int Ncomp)
{
    using amrex::Real;;
    Real cellvol = dx[0]*dx[1];

#if(AMREX_SPACEDIM > 2)
    cellvol *=dx[2];
#endif


#if 0
    Real rad = Real(0.1);
    Real xcent = Real(0.5);
    Real ycent = Real(0.5);
    Real x = prob_lo[0] + (i+Real(0.5)) * dx[0] - xcent;
    Real y = prob_lo[1] + (j+Real(0.5)) * dx[1] - ycent;

#if (AMREX_SPACEDIM > 2)
    Real zcent = Real(0.5);
    Real z = prob_lo[2] + (k+Real(0.5)) * dx[2] - zcent;
    Real r2 = (x*x + y*y + z*z) / (rad*rad);
#else
    Real z = Real(0.);
    Real r2 = (x*x + y*y) / (rad*rad);
#endif

    // Gaussian 
    phi(i,j,k,0) = Real(1.) + std::exp(-r2);
#endif

#if 0
    // Jump in x-direction
    Real x = prob_lo[0] + (i+Real(0.5)) * dx[0];
    Real y = prob_lo[1] + (j+Real(0.5)) * dx[1];
    if ( x > 0.5) {
        phi(i,j,k,0) = 1.5*npts_scale;
    } else {
        phi(i,j,k,0) = 0.5*npts_scale;
    }
#endif

#if 0
    // Delta function 
    if (i==4 && j ==4)
    {
       phi(i,j,k,0) = 1./(dx[0]*dx[1]);
    }
    else
    {
       phi(i,j,k,0) = 0;;
    }
#endif

#if 0
    // 2D Test for multilevel
    Real rad = Real(0.2);
    Real xcent = Real(0.5);
    Real ycent = Real(0.5);
    Real x = prob_lo[0] + (i+Real(0.5)) * dx[0] - xcent;
    Real y = prob_lo[1] + (j+Real(0.5)) * dx[1] - ycent;
    Real r2 = (x*x + y*y);
    r2 = 4.*(x-y)*(x-y)+ .25*(x+y)*(x+y);
#if (AMREX_SPACEDIM == 3)
    Real zcent = Real(0.75);
    Real z = prob_lo[2] + (k+Real(0.5)) * dx[2] - zcent;
    r2 = (x*x + y*y+z*z);
#endif

    if ( r2 > rad*rad) {
        phi(i,j,k,0) = 15./cellvol;
    } else {
        phi(i,j,k,0) = 0./cellvol;
    }
#endif

#if 1
    // 2D Test for multilevel
    Real rad = Real(0.4);
    Real radi = Real(0.1);
    Real xcent = Real(0.5);
    Real ycent = Real(0.5);
    Real x = prob_lo[0] + (i+Real(0.5)) * dx[0] - xcent;
    Real y = prob_lo[1] + (j+Real(0.5)) * dx[1] - ycent;
    Real r2 = (x*x + y*y);
    r2 = 2.*(x-y)*(x-y)+ .5*(x+y)*(x+y);
    Real r3 = 0.5*(x-y)*(x-y)+ 2.*(x+y)*(x+y);
#if (AMREX_SPACEDIM == 3)
    Real zcent = Real(0.5);
    Real z = prob_lo[2] + (k+Real(0.5)) * dx[2] - zcent;
    r3 = (x*x + y*y+z*z);
    r2 = (x-.1)*(x-.1)+(y-.2)*(y-.2)+z*z;
#endif

    if ( r3 > rad*rad) {
        //phi(i,j,k,0) = 50./cellvol;
        phi(i,j,k,0) = 30./cellvol;
    } else if (r2 < radi*radi) {
        //phi(i,j,k,0) = 20./cellvol;
        phi(i,j,k,0) = 15./cellvol;
    } else {
        phi(i,j,k,0) = 0./cellvol;
    }
/*

    //if(std::abs(x) <= .1 && y > 0.){
    if(x >= .3 && y > .3){
       phi(i,j,k,0) = 0.;
    } else {
       phi(i,j,k,0) = 30./cellvol;
    }
    int ioff = 16;
    int joff = 16;
    ioff = 0;
    joff = 0;

    if( i-ioff >= 0 && i - ioff < 32 && j-joff >= 0 && j - joff < 32 ) {
       phi(i,j,k,0) = 0./cellvol;
    } else {
       phi(i,j,k,0) = 30./cellvol;
    }

*/

#endif

    if (Ncomp == 2){
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
    int coeff_index = Ncomp-1;
    if (lo == dom_lo &&
        (bc_lo == amrex::BCType::foextrap ||
         bc_lo == amrex::BCType::ext_dir))
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
             (bc_hi == amrex::BCType::foextrap ||
              bc_hi == amrex::BCType::ext_dir))
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
        amrex::Real phip = std::max(phi(i,j,k,coeff_index),0.);
        amrex::Real phim = std::max(phi(i-1,j,k,coeff_index),0.);
	//amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
	//amrex::Real phiavg = 0.5*(phip+phim);

	// HACK TO COMMENT OUT STOCHASTIC TERMS
	// phiavg = 0.;

//        fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv
//    	                    + stochfluxx(i,j,k) * std::sqrt(phiavg);

	  amrex::Real sqphiavg = 0.5*(std::sqrt(phip)+std::sqrt(phim));
          fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv
      	                    + stochfluxx(i,j,k) * sqphiavg;
        if(Ncomp == 2){
           fluxx(i,j,k,1) = 0.5*(phi(i,j,k,1)-phi(i-1,j,k,1)) * dxinv;
        }
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_y (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxy,
                     amrex::Array4<amrex::Real> const& stochfluxy,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dyinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi, int Ncomp)
{
    int coeff_index = Ncomp-1;
    if (lo == dom_lo &&
        (bc_lo == amrex::BCType::foextrap ||
         bc_lo == amrex::BCType::ext_dir))
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
             (bc_hi == amrex::BCType::foextrap ||
              bc_hi == amrex::BCType::ext_dir))
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
        amrex::Real phip = std::max(phi(i,j,k,coeff_index),0.);
        amrex::Real phim = std::max(phi(i,j-1,k,coeff_index),0.);
	//amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
	//amrex::Real phiavg = 0.5*(phip+phim);

	// HACK TO COMMENT OUT STOCHASTIC TERMS
	// phiavg = 0.;

        //fluxy(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j-1,k,0)) * dyinv
	//+ stochfluxy(i,j,k)* std::sqrt(phiavg);

	amrex::Real sqphiavg = 0.5*(std::sqrt(phip)+std::sqrt(phim));
        fluxy(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j-1,k,0)) * dyinv
	+ stochfluxy(i,j,k)* sqphiavg;
        if(Ncomp == 2){
           fluxy(i,j,k,1) = 0.5*(phi(i,j,k,1)-phi(i,j-1,k,1)) * dyinv;
        }
    }
}


#if (AMREX_SPACEDIM > 2)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_z (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxz,
                     amrex::Array4<amrex::Real> const& stochfluxz,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dzinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi, int Ncomp)
{
    int coeff_index = Ncomp-1;
    if (lo == dom_lo &&
        (bc_lo == amrex::BCType::foextrap ||
         bc_lo == amrex::BCType::ext_dir))
    {
        if(k == lo)
        {
            fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv / 0.5;
        }
        else
        {
            fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv;
        }
    }
    else if (hi == dom_hi &&
             (bc_hi == amrex::BCType::foextrap ||
              bc_hi == amrex::BCType::ext_dir))
    {
        if(k == hi+1)
        {
            fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv / 0.5;
        }
        else
        {
            fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv;
        }
    }
    else
    {
//	amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i,j,k-1))/(phi(i,j,k)+phi(i,j,k-1)+1.e-16);
        amrex::Real phip = std::max(phi(i,j,k,coeff_index),0.);
        amrex::Real phim = std::max(phi(i,j,k-1,coeff_index),0.);
	//amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
        //fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv
	// +stochfluxz(i,j,k) * std::sqrt(phiavg);
	amrex::Real sqphiavg = 0.5*(std::sqrt(phip)+std::sqrt(phim));
        fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv
	 +stochfluxz(i,j,k) * sqphiavg;
        if(Ncomp == 2){
           fluxz(i,j,k,1) = 0.5*(phi(i,j,k,1)-phi(i,j,k-1,1)) * dzinv;
        }
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
