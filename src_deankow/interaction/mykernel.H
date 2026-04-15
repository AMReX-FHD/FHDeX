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
    constexpr double  PI = 3.14159265358979323846264338327950288;
    using amrex::Real;;
    Real cellvol = dx[0]*dx[1];

#if(AMREX_SPACEDIM > 2)
    cellvol *=dx[2];
#endif

//     phi(i,j,k,0) = 20./cellvol;
//     phi(i,j,k,0) += std::cos(2.*PI*(i)*dx[0])*std::cos(2.*PI*(j)*dx[1])/cellvol;
     //phi(i,j,k,0) += std::cos(16.*2.*PI*(i)*dx[0])*std::cos(16*2.*PI*j*dx[1])/cellvol;
     phi(i,j,k,0) = 20480.;
     phi(i,j,k,0) += std::sin(2.*PI*(i)*dx[0])*1024;
     //phi(i,j,k,0) += std::sin(2.*PI*(j)*dx[1])*1024;
//     phi(i,j,k,0) += std::cos(2.*PI*(i)*dx[0])*std::cos(2.*PI*(j)*dx[1])*1024;
//     phi(i,j,k,0) += std::cos(8.*PI*(i)*dx[0]+.7)*std::cos(10.*PI*(j)*dx[1]-.3)*1024;
     // phi(i,j,k,0) += std::cos(2.*3.14156265*(i+0.5)*dx[0])*std::cos(2.*3.14156265*(j+0.5)*dx[1])/cellvol;

    if (Ncomp == 2){
        phi(i,j,k,1) = phi(i,j,k,0);
    }
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void init_int_pot (int i, int j, int k,
               amrex::Array4<amrex::Real> const& U,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_hi)
{
    constexpr double  PI = 3.14159265358979323846264338327950288;
    amrex::Real eps, R, alpha;
    eps = 0.0333*20000.;
    R = 0.1;
    alpha = 3.;

//    shiftinging definit of U by 1/2 for fft
    amrex::Real x,y,z,xlen,ylen,zlen;
    
    x = (i) * dx[0] ;
    xlen = prob_hi[0]-prob_lo[0];
    if(x + .5*dx[0] > 0.5*xlen){
       x = x-xlen;
    }

    y = (j) * dx[1] ;
    ylen = prob_hi[1]-prob_lo[1];
    if(y  + .5*dx[1]> 0.5*ylen){
       y = y-ylen;
    }

#if (AMREX_SPACEDIM > 2)
    z = (k) * dx[2] ;
    zlen = prob_hi[2]-prob_lo[2];
    if(z + .5*dx[2] > 0.5*zlen){
       z = z-zlen;
    }
    amrex::Real r = std::sqrt(x*x+y*y+z*z);
#else
    amrex::Real r = std::sqrt(x*x+y*y);
#endif

    U(i,j,k) = eps * std::exp( - std::pow(r/R,alpha));
//    U(i,j,k) = std::cos(2*PI*i*dx[0]);
/*
    U(i,j,k) = 0.;

    if (i == 0 and j == 0 ){
       U(i,j,k) = 1.;
    } 
    if (i == 1 and j == 0 ){
       U(i,j,k) = .5;
    } 
    if (i == 31 and j == 0 ){
       U(i,j,k) = .5;
    } 
    if (i == 0 and j == 1 ){
       U(i,j,k) = .5;
    } 
    if (i == 1 and j == 1 ){
       U(i,j,k) = .25;
    } 
    if (i == 31 and j == 1 ){
       U(i,j,k) = .25;
    } 
    if (i == 0 and j == 31 ){
       U(i,j,k) = .5;
    } 
    if (i == 1 and j == 31 ){
       U(i,j,k) = .25;
    } 
    if (i == 31 and j == 31 ){
       U(i,j,k) = .25;
    } 
*/


//        amrex::Print() << "i j U " << i << " " << j << " " << U(i,j,k) << std::endl;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_x (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxx,
                     amrex::Array4<amrex::Real> const& stochfluxx,
                     amrex::Array4<amrex::Real> const& C,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dxinv, amrex::Real diff_coeff,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi,int Ncomp)
{
    int coeff_index = Ncomp-1;
    amrex::Real phibdry = 0.;
    amrex::Real stofact = 1.;

    int ext_pot = 0;
    int inter_pot = 1;
    amrex::Real alpha = .25;
    amrex::Real beta = .75;
    amrex::Real gamma = 5.e-4;
    gamma = 1.1e-3;

    amrex::Real phip,phim, phipn, phimn;

/*
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
*/
    //amrex::Print() << "BC info " << i << " " << dom_lo << " " << bc_lo << " " << dom_hi << " " << bc_hi << std::endl;
    if (i == dom_lo && (bc_lo == amrex::BCType::foextrap || bc_lo == amrex::BCType::ext_dir)) {

        stofact = std::sqrt(2.);

        phim = 2.*phibdry - phi(i,j,k,0);
        phip = phi(i,j,k,0);

        //phimn = 2.*phibdry - phi(i,j,k,coeff_index);
        //phipn = phi(i,j,k,coeff_index);
        phimn = phibdry;
        phipn = phibdry;
        //amrex::Print() << " A " << std::endl;

    } else if (i == dom_hi+1 && (bc_hi == amrex::BCType::foextrap || bc_hi == amrex::BCType::ext_dir)) {

        stofact = std::sqrt(2.);
        phip = 2.*phibdry - phi(i-1,j,k,0);
        phim = phi(i-1,j,k,0);
        //phipn = 2.*phibdry - phi(i-1,j,k,coeff_index);
        //phimn = phi(i-1,j,k,coeff_index);
        phimn = phibdry;
        phipn = phibdry;
        //amrex::Print() << " B " << std::endl;
    } else {

        phip = phi(i,j,k,0);
        phim = phi(i-1,j,k,0);
        phipn = std::abs(phi(i,j,k,coeff_index));
        phimn = std::abs(phi(i-1,j,k,coeff_index));
        stofact = 1.;
        //amrex::Print() << " C " << std::endl;

    }

//        amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i-1,j,k))/(phi(i,j,k)+phi(i-1,j,k)+1.e-16);
        //amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
        //amrex::Real phiavg = 0.5*(phip+phim);

        // HACK TO COMMENT OUT STOCHASTIC TERMS
        // phiavg = 0.;

        //fluxx(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i-1,j,k,0)) * dxinv
            //                    + stochfluxx(i,j,k) * std::sqrt(phiavg);
          //amrex::Print() << " point " << i << " " << j << " " <<  dxinv << std::endl;
          //amrex::Print() << "phip " << phip << std::endl;
          //amrex::Print() << "phib " << phibdry << std::endl;
          //amrex::Print() << "phi " << phi(i,j,k,0) << std::endl;
          //amrex::Print() << "phim " << phim << std::endl;

          fluxx(i,j,k,0) = diff_coeff*(phip-phim) * dxinv;

          if(ext_pot == 1) {

             amrex::Real xloc;
             xloc = i / dxinv;
             amrex::Real phiavg = 0.5*(phim+phip);

             amrex::Real Vsubx;

             Vsubx = 2.*(xloc - beta) * (xloc - alpha)* (2.*xloc - alpha - beta) / gamma;

             fluxx(i,j,k,0) += phiavg * Vsubx;

          }

          if(inter_pot == 1) {

             amrex::Real phiavg = 0.5*(phim+phip);

             fluxx(i,j,k,0) += phiavg * (C(i,j,k) - C(i-1,j,k)) * dxinv;

          }

        if(Ncomp == 2){
           fluxx(i,j,k,1) = fluxx(i,j,k,0);
        }

          amrex::Real sqphiavg = 0.5*(std::sqrt(phipn)+std::sqrt(phimn));
          fluxx(i,j,k,0) += stofact * stochfluxx(i,j,k) * sqphiavg;
//          amrex::Print() << " x flux " << i << " " << j << " " << fluxx(i,j,k,0) << std::endl;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_y (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxy,
                     amrex::Array4<amrex::Real> const& stochfluxy,
                     amrex::Array4<amrex::Real> const& C,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dyinv, amrex::Real diff_coeff,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi, int Ncomp)
{
    int coeff_index = Ncomp-1;
    amrex::Real phibdry = 0.;
    amrex::Real stofact = 1.;

    int ext_pot = 0;
    int inter_pot = 1;
    amrex::Real alpha = .25;
    amrex::Real beta = .75;
    amrex::Real gamma = 5.e-4;
    gamma = 1.1e-3;

    amrex::Real phip,phim, phipn, phimn;

    if (j == dom_lo && (bc_lo == amrex::BCType::foextrap || bc_lo == amrex::BCType::ext_dir))
    {

        stofact = std::sqrt(2.);
        phim = 2.*phibdry - phi(i,j,k,0);
        phip = phi(i,j,k,0);
        //phimn = 2.*phibdry - std::abs(phi(i,j,k,coeff_index));
        //phipn = std::abs( phi(i,j,k,coeff_index));
        phipn = phibdry;
        phimn = phibdry;

    } else if (j == dom_hi+1 && (bc_hi == amrex::BCType::foextrap || bc_hi == amrex::BCType::ext_dir)) {

        stofact = std::sqrt(2.);
        phip = 2.*phibdry - phi(i,j-1,k,0);
        phim = phi(i,j-1,k,0);
        //phipn = 2.*phibdry - std::abs(phi(i,j-1,k,coeff_index));
        //phimn = std::abs(phi(i,j-1,k,coeff_index));
        phipn = phibdry;
        phimn = phibdry;

     } else {

        phip = phi(i,j,k,0);
        phim = phi(i,j-1,k,0);
        phipn = std::abs(phi(i,j,k,coeff_index));
        phimn = std::abs(phi(i,j-1,k,coeff_index));
        stofact = 1.;

    }
//        amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i,j-1,k))/(phi(i,j,k)+phi(i,j-1,k)+1.e-16);
        //amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
        //amrex::Real phiavg = 0.5*(phip+phim);

        // HACK TO COMMENT OUT STOCHASTIC TERMS
        // phiavg = 0.;

        //fluxy(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j-1,k,0)) * dyinv
        //+ stochfluxy(i,j,k)* std::sqrt(phiavg);

        fluxy(i,j,k,0) = diff_coeff*(phip-phim) * dyinv;

          if(ext_pot == 1) {

             amrex::Real yloc;
             yloc = j / dyinv;
             amrex::Real phiavg = 0.5*(phip+phim);

             amrex::Real Vsuby;

             Vsuby =0.5* 4.*(yloc - .5)*(yloc - .5)*(yloc - .5) / gamma;

             fluxy(i,j,k,0) += phiavg * Vsuby;

          }

          if(inter_pot == 1) {

             amrex::Real phiavg = 0.5*(phim+phip);

             fluxy(i,j,k,0) += phiavg * (C(i,j,k) - C(i,j-1,k)) * dyinv;

          }

        if(Ncomp == 2){
           fluxy(i,j,k,1) = fluxy(i,j,k,0);
        }

        amrex::Real sqphiavg = 0.5*(std::sqrt(phipn)+std::sqrt(phimn));

        fluxy(i,j,k,0) += stofact * stochfluxy(i,j,k)* sqphiavg ;
//        amrex::Print() << " y flux " << i << " " << j << " " << fluxy(i,j,k,0) << std::endl;
}


#if (AMREX_SPACEDIM > 2)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_z (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxz,
                     amrex::Array4<amrex::Real> const& stochfluxz,
                     amrex::Array4<amrex::Real> const& C,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dzinv, amrex::Real diff_coeff,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi, int Ncomp)
{
    int coeff_index = Ncomp-1;
    amrex::Real phibdry = 0.;
    amrex::Real stofact = 1.;

    int ext_pot = 0;
    int inter_pot = 1;
    amrex::Real alpha = .25;
    amrex::Real beta = .75;
    amrex::Real gamma = 5.e-4;
    gamma = 1.1e-3;

    amrex::Real phip,phim, phipn, phimn;

    if (k == dom_lo && (bc_lo == amrex::BCType::foextrap || bc_lo == amrex::BCType::ext_dir))
    {

        stofact = std::sqrt(2.);
        phim = 2.*phibdry - phi(i,j,k,0);
        phip = phi(i,j,k,0);
        //phimn = 2.*phibdry - std::abs(phi(i,j,k,coeff_index));
        //phipn = std::abs( phi(i,j,k,coeff_index));
        phipn = phibdry;
        phimn = phibdry;

    } else if (k == dom_hi+1 && (bc_hi == amrex::BCType::foextrap || bc_hi == amrex::BCType::ext_dir)) {

        stofact = std::sqrt(2.);
        phip = 2.*phibdry - phi(i,j,k-1,0);
        phim = phi(i,j,k-1,0);
        //phipn = 2.*phibdry - std::abs(phi(i,j-1,k,coeff_index));
        //phimn = std::abs(phi(i,j-1,k,coeff_index));
        phipn = phibdry;
        phimn = phibdry;

     } else {

        phip = phi(i,j,k,0);
        phim = phi(i,j,k-1,0);
        phipn = std::abs(phi(i,j,k,coeff_index));
        phimn = std::abs(phi(i,j,k-1,coeff_index));
        stofact = 1.;

    }
//        amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i,j-1,k))/(phi(i,j,k)+phi(i,j-1,k)+1.e-16);
        //amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
        //amrex::Real phiavg = 0.5*(phip+phim);

        // HACK TO COMMENT OUT STOCHASTIC TERMS
        // phiavg = 0.;

        //fluxy(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j-1,k,0)) * dyinv
        //+ stochfluxy(i,j,k)* std::sqrt(phiavg);

        fluxz(i,j,k,0) = diff_coeff*(phip-phim) * dzinv;

          if(ext_pot == 1) {

             amrex::Real zloc;
             zloc = k / dzinv;
             amrex::Real phiavg = 0.5*(phip+phim);

             amrex::Real Vsubz;

             Vsubz =0.5* 4.*(zloc - .5)*(zloc - .5)*(zloc - .5) / gamma;

             fluxz(i,j,k,0) += phiavg * Vsuby;

          }

          if(inter_pot == 1) {

             amrex::Real phiavg = 0.5*(phim+phip);

             fluxz(i,j,k,0) += phiavg * (C(,j,k) - C(i,j,k-1)) * dzinv;

          }

        if(Ncomp == 2){
           fluxz(i,j,k,1) = fluxz(i,j,k,0);
        }

        amrex::Real sqphiavg = 0.5*(std::sqrt(phipn)+std::sqrt(phimn));

        fluxz(i,j,k,0) += stofact * stochfluxz(i,j,k)* sqphiavg;

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
