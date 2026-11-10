#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>
#include "PotParams.H"

// Initial condition, in particles per cell (see ICParams in PotParams.H)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void init_phi (int i, int j, int k,
               amrex::Array4<amrex::Real> const& phi,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
               ICParams const& ic, PotParams const& pot, int Ncomp)
{
    using amrex::Real;
    Real cellvol = dx[0]*dx[1];
#if(AMREX_SPACEDIM > 2)
    cellvol *= dx[2];
#endif

    Real npc = ic.npc;
    if (ic.type == 1) {
        // Boltzmann distribution of the external potential, mean ic.npc per cell
        const Real x = prob_lo[0] + (i + Real(0.5))*dx[0];
        const Real y = prob_lo[1] + (j + Real(0.5))*dx[1];
        npc = ic.npc*std::exp(-2.*pot.V(x,y)/pot.gamma)/ic.norm;
    } else if (ic.type == 2) {
        // step: ic.npc_in inside the index box, ic.npc outside
        const bool inside = (i >= ic.box_lo[0] && i <= ic.box_hi[0] &&
                             j >= ic.box_lo[1] && j <= ic.box_hi[1]
#if (AMREX_SPACEDIM > 2)
                             && k >= ic.box_lo[2] && k <= ic.box_hi[2]
#endif
                             );
        if (inside) { npc = ic.npc_in; }
    }
    phi(i,j,k,0) = npc/cellvol;

    if (Ncomp == 2){
        phi(i,j,k,1) = phi(i,j,k,0);
    }
}


AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_x (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxx,
                     amrex::Array4<amrex::Real> const& stochfluxx,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dxinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi,int Ncomp,
                     amrex::Real xlo, PotParams const& pot)
{
    int coeff_index = Ncomp-1;
    amrex::Real phibdry = 0.;
    amrex::Real stofact = 1.;

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
        //phipn = std::abs(phi(i,j,k,coeff_index));
        //phimn = std::abs(phi(i-1,j,k,coeff_index));
        phipn = std::max(phi(i,j,k,coeff_index),0.);
        phimn = std::max(phi(i-1,j,k,coeff_index),0.);
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

          fluxx(i,j,k,0) = 0.5*(phip-phim) * dxinv;

          if (pot.ext_pot == 1) {
             // drift flux phi dV/dx / gamma at the face x_i
             const amrex::Real xloc = xlo + i / dxinv;
             const amrex::Real phiavg = 0.5*(phim+phip);
             fluxx(i,j,k,0) += phiavg * pot.dVx(xloc) / pot.gamma;
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
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dyinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi, int Ncomp,
                     amrex::Real ylo, PotParams const& pot)
{
    int coeff_index = Ncomp-1;
    amrex::Real phibdry = 0.;
    amrex::Real stofact = 1.;

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
        // phipn = std::abs(phi(i,j,k,coeff_index));
        // phimn = std::abs(phi(i,j-1,k,coeff_index));
        phipn = std::max(phi(i,j,k,coeff_index),0.);
        phimn = std::max(phi(i,j-1,k,coeff_index),0.);
        stofact = 1.;

    }
//        amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i,j-1,k))/(phi(i,j,k)+phi(i,j-1,k)+1.e-16);
        //amrex::Real phiavg = 2.*(phip*phim)/(phip+phim+1.e-16);
        //amrex::Real phiavg = 0.5*(phip+phim);

        // HACK TO COMMENT OUT STOCHASTIC TERMS
        // phiavg = 0.;

        //fluxy(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j-1,k,0)) * dyinv
        //+ stochfluxy(i,j,k)* std::sqrt(phiavg);

        fluxy(i,j,k,0) = 0.5*(phip-phim) * dyinv;

          if (pot.ext_pot == 1) {
             // drift flux phi dV/dy / gamma at the face y_j
             const amrex::Real yloc = ylo + j / dyinv;
             const amrex::Real phiavg = 0.5*(phip+phim);
             fluxy(i,j,k,0) += phiavg * pot.dVy(yloc) / pot.gamma;
          }

        if(Ncomp == 2){
           fluxy(i,j,k,1) = fluxy(i,j,k,0);
        }

        amrex::Real sqphiavg = 0.5*(std::sqrt(phipn)+std::sqrt(phimn));

        fluxy(i,j,k,0) += stofact * stochfluxy(i,j,k)* sqphiavg;
//        amrex::Print() << " y flux " << i << " " << j << " " << fluxy(i,j,k,0) << std::endl;
}


#if (AMREX_SPACEDIM > 2)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void compute_flux_z (int i, int j, int k,
                     amrex::Array4<amrex::Real> const& fluxz,
                     amrex::Array4<amrex::Real> const& stochfluxz,
                     amrex::Array4<amrex::Real const> const& phi, amrex::Real dzinv,
                     int lo, int hi, int dom_lo, int dom_hi, int bc_lo, int bc_hi, int Ncomp)
{
    amrex::ignore_unused(lo, hi);
    int coeff_index = Ncomp-1;
    // Note: the test has to be on the face index k, not on the box bounds lo/hi,
    // otherwise every face of a box that touches the domain boundary is treated
    // as a boundary face.
    if (k == dom_lo &&
        (bc_lo == amrex::BCType::foextrap ||
         bc_lo == amrex::BCType::ext_dir))
    {
        fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv / 0.5;
        if(Ncomp == 2){
            fluxz(i,j,k,1) = 0.5*(phi(i,j,k,1)-phi(i,j,k-1,1)) * dzinv / 0.5;
        }
    }
    else if (k == dom_hi+1 &&
             (bc_hi == amrex::BCType::foextrap ||
              bc_hi == amrex::BCType::ext_dir))
    {
        fluxz(i,j,k,0) = 0.5*(phi(i,j,k,0)-phi(i,j,k-1,0)) * dzinv / 0.5;
        if(Ncomp == 2){
            fluxz(i,j,k,1) = 0.5*(phi(i,j,k,1)-phi(i,j,k-1,1)) * dzinv / 0.5;
        }
    }
    else
    {
//        amrex::Real phiavg = 2.*(phi(i,j,k)*phi(i,j,k-1))/(phi(i,j,k)+phi(i,j,k-1)+1.e-16);
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
