#include "hydro_functions.H"
#include "hydro_functions_F.H"

#include "common_functions.H"

using namespace amrex;

void MkAdvSFluxdiv(const std::array<MultiFab, AMREX_SPACEDIM>& umac,
		   const MultiFab& s,
		   MultiFab& s_update,
		   const amrex::Real* dx,
		   const Geometry& geom,
		   const int& comp,
		   const int& increment)
{

     BL_PROFILE_VAR("MkAdvSFluxdiv()",MkAdvSFluxdiv);

     Real dxinv = 1./dx[0];
     int ncomp = 1;

     // Loop over boxes
     for (MFIter mfi(s); mfi.isValid(); ++mfi) {

         // Create cell-centered box from semi-nodal box
         const Box& bx = mfi.validbox();

         Array4<Real const> const& s_fab = s.array(mfi);

         Array4<Real> const& s_update_fab = s_update.array(mfi);

         AMREX_D_TERM(Array4<Real const> const& umac_fab = umac[0].array(mfi);,
                      Array4<Real const> const& vmac_fab = umac[1].array(mfi);,
                      Array4<Real const> const& wmac_fab = umac[2].array(mfi););

         if (increment == 1) {

         }
         else {
             AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
             {
                 s_update_fab(i,j,k,n+comp) = 
                     - dxinv*( 0.5*(s_fab(i+1,j,k,n+comp)+s_fab(i,j,k,n+comp))*umac_fab(i+1,j,k) - 0.5*(s_fab(i,j,k,n+comp)+s_fab(i-1,j,k,n+comp))*umac_fab(i,j,k) )
                     - dxinv*( 0.5*(s_fab(i,j+1,k,n+comp)+s_fab(i,j,k,n+comp))*vmac_fab(i,j+1,k) - 0.5*(s_fab(i,j,k,n+comp)+s_fab(i,j-1,k,n+comp))*vmac_fab(i,j,k) )
#if (AMREX_SPACEDIM == 3)
                     - dxinv*( 0.5*(s_fab(i,j,k+1,n+comp)+s_fab(i,j,k,n+comp))*wmac_fab(i,j,k+1) - 0.5*(s_fab(i,j,k,n+comp)+s_fab(i,j,k-1,n+comp))*wmac_fab(i,j,k) )
#endif
                     ;
             });
         }
     }
}
