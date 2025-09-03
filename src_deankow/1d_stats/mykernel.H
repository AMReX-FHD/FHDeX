#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real reg_sqrt (amrex::Real phi, amrex::Real deltaN) noexcept
{
    amrex::Real phiabs = std::abs(phi);
    amrex::Real smsqrt;
    if(phiabs <= deltaN/2.) {
        smsqrt = phiabs/std::sqrt(deltaN);
    } else if (phiabs <= deltaN) {
        smsqrt = std::sqrt(deltaN)/2. - 1.5*phiabs/std::sqrt(deltaN) + 4.*phiabs*phiabs*std::pow(deltaN,-1.5)-2.*std::pow(phiabs,3.)*std::pow(deltaN,-2.5);
    } else {
        smsqrt=std::sqrt(phiabs);
    }
    return (  (phi > 0.) ? smsqrt : -smsqrt);
    //return (  phi/std::sqrt(deltaN));
}



AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void init_phi (int i, int j, int k,
               amrex::Array4<amrex::Real> const& phi,
               GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
               amrex::Real npts_scale,int Ncomp, amrex::Real phileft, amrex::Real phiright)
{
    using amrex::Real;;

    amrex::Real slope = phiright-phileft;


    Real x = prob_lo[0] + (i+Real(0.5)) * dx[0];
    Real y = prob_lo[1] + (j+Real(0.5)) * dx[1];
#if (AMREX_SPACEDIM > 2)
    Real z = prob_lo[2] + (k+Real(0.5)) * dx[2];
    Real r2 = ((x-Real(0.25))*(x-Real(0.25))+(y-Real(0.25))*(y-Real(0.25))+(z-Real(0.25))*(z-Real(0.25)))/Real(0.01);
#else
    Real z = Real(0.);
    Real r2 = ((x-Real(0.25))*(x-Real(0.25))+(y-Real(0.25))*(y-Real(0.25)))/Real(0.01);
#endif
//    phi(i,j,k,0) = 50./cellvol;
//    phi(i,j,k,0) = .1;

    phi(i,j,k,0) = phileft+x*slope;
    if(Ncomp == 2){
        phi(i,j,k,1) = phi(i,j,k,0);
    }
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_x (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxx,
                     amrex::Array4<amrex::Real> const& stochfluxx,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dxinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi,int Ncomp,
                     amrex::Real phileft, amrex::Real phiright, int avg_type)
{
    //    amrex::Real cellvol = 1./dxinv;
    //    phileft = .1;
    //    phiright = .1;
    //amrex::Real deltaN = std::pow(0.5*(phileft+phiright),2./3.);
    amrex::Real deltaN = std::pow(0.5*(phileft+phiright),1./3.);
    //amrex::Print() << "deltaN " << deltaN << std::endl;
    //amrex::Real val,flag;
    //for (int n=0; n<201; n++){
    //   val=(n-100.)*deltaN/50.;
    //   flag = (val>0) ? 1.:-1.;
    //   amrex::Print() << "curve " << val << " " << reg_sqrt(val,deltaN) << " " << flag*std::sqrt(std::abs(val)) << std::endl;
    //}

    int coeff_index = Ncomp-1;
    if (i == dom_lo && bc_lo == BCType::foextrap )
    {
        fluxx(i,j,k,0) = 0.;
        if(Ncomp == 2){
            fluxx(i,j,k,1) = 0.;
        }
    } else if ( i == dom_lo && bc_lo == BCType::ext_dir )
    {
        fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phileft) * dxinv / 0.5 + std::sqrt(2.*phileft)*stochfluxx(i,j,k);
        if(Ncomp == 2){
            fluxx(i,j,k,1) = 0.5*(phi(i,j,k,1)-phileft) * dxinv / 0.5 ;
        }
    } else if ( i == dom_hi+1 && bc_lo == BCType::foextrap )
    {
        fluxx(i,j,k,0) = 0.;
        if(Ncomp == 2){
            fluxx(i,j,k,1) = 0.;
        }
    } else if ( i == dom_hi+1 && bc_lo == BCType::ext_dir )
    {
        fluxx(i,j,k,0) = 0.5*(phiright-phi(i-1,j,k,0)) * dxinv / 0.5 + std::sqrt(2.*phiright)*stochfluxx(i,j,k);
        if(Ncomp == 2){
            fluxx(i,j,k,1) = 0.5*(phiright-phi(i-1,j,k,1)) * dxinv / 0.5 ;
        }
    } else
    {
        //amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i-1,j,k))/(phi(i,j,k)+phi(i-1,j,k)+1.e-16);
        amrex::Real phiavg, phip, phim, noise_coeff;
        if(avg_type == 0 || avg_type == 2 ) {

            phip = std::max(phi(i,j,k,coeff_index),0.);
            phim = std::max(phi(i-1,j,k,coeff_index),0.);
            phiavg = 0.5*(phip+phim);
            noise_coeff = std::sqrt(phiavg);

        } else if (avg_type == 1) {

            phip = std::max(phi(i,j,k,coeff_index),0.);
            phim = std::max(phi(i-1,j,k,coeff_index),0.);
            phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
            noise_coeff = std::sqrt(phiavg);

        } else if (avg_type == 3) {

            phip = phi(i,j,k,coeff_index);
            phim = phi(i-1,j,k,coeff_index);
            phiavg = 0.5*(phip+phim);
            noise_coeff = reg_sqrt(phiavg,deltaN);
        }
        //phiavg = .1;
        //amrex::Print() << "test " << std::sqrt(phiavg) << " " << phiavg/std::sqrt(deltaN) << " " << deltaN << std::endl;
        //phiavg = 1600.;
        fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv
            + stochfluxx(i,j,k) * noise_coeff;
        //  + stochfluxx(i,j,k) * reg_sqrt(phiavg,deltaN);
        //  + stochfluxx(i,j,k) * phiavg/std::sqrt(deltaN);
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
        //amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i,j-1,k))/(phi(i,j,k)+phi(i,j-1,k)+1.e-16);
        Real phip = std::max(phi(i,j,k,coeff_index),0.);
        Real phim = std::max(phi(i,j-1,k,coeff_index),0.);
        amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
        fluxy(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j-1,k,0)) * dyinv
            + stochfluxy(i,j,k)* std::sqrt(phiavg);
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
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi)
{
    int coeff_index = Ncomp-1;
    if (lo == dom_lo &&
        (bc_lo == BCType::foextrap ||
         bc_lo == BCType::ext_dir))
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
             (bc_hi == BCType::foextrap ||
              bc_hi == BCType::ext_dir))
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
        //amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i,j,k-1))/(phi(i,j,k)+phi(i,j,k-1)+1.e-16);
        Real phip = std::max(phi(i,j,k,coeff_index),0.);
        Real phim = std::max(phi(i,j,k-1,coeff_index),0.);
        amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
        fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv
            + stochfluxz(i,j,k) * std::sqrt(phiavg);
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