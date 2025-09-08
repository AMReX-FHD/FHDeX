#include "compressible_functions.H"
#include "common_functions.H"

using namespace common;
using namespace compressible;

#if defined(PELEPHYSICS)
void calculateTransportCoeffs(const MultiFab& prim_in, 
			                  MultiFab& eta_in, MultiFab& zeta_in, MultiFab& kappa_in,
			                  MultiFab& chi_in, MultiFab& Dij_in,
                              pele::physics::PeleParams<pele::physics::transport::TransParm<
                              pele::physics::PhysicsType::eos_type,
                              pele::physics::PhysicsType::transport_type>>& trans_parms)
{
    BL_PROFILE_VAR("calculateTransportCoeffs()",calculateTransportCoeffs);

    // see comments in conservedPrimitiveConversions.cpp regarding alternate ways of declaring
    // thread shared and thread private arrays on GPUs
    // if the size is not known at compile time, alternate approaches are required
    // here we know the size at compile time

    amrex::IntVect ng_temp;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        ng_temp[d] = (do_1D == 1) ? 1 : ngc[d];
    }

    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        // grow the box by ngc
        const Box& bx = amrex::grow(mfi.tilebox(), ng_temp);

        const Array4<const Real>& prim = prim_in.array(mfi);

        const Array4<Real>& eta   =   eta_in.array(mfi);
        const Array4<Real>& zeta  =  zeta_in.array(mfi);
        const Array4<Real>& kappa = kappa_in.array(mfi);
        const Array4<Real>& chi   =   chi_in.array(mfi);
        const Array4<Real>& Dij   =   Dij_in.array(mfi);

        auto const* ltransparm = trans_parms.device_parm();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
        
            GpuArray<Real,MAX_SPECIES> Yk_fixed;
            GpuArray<Real,MAX_SPECIES> Xk_fixed;
            
            Real sumYk = 0.;
            for (int n=0; n<nspecies; ++n) {
	        if (prim(i,j,k,6+n) <= 0.0) amrex::Abort("Negative mass fraction encountered");
   	        if (prim(i,j,k,6+n) >= 1.0) amrex::Abort("Greater than unity mass fraction encountered");
                Yk_fixed[n] = amrex::max(0.,amrex::min(1.,prim(i,j,k,6+n)));
                sumYk += Yk_fixed[n];
            }

            for (int n=0; n<nspecies; ++n) {
                Yk_fixed[n] /= sumYk;
            }

            // compute mole fractions from mass fractions
            GetMolfrac(Yk_fixed, Xk_fixed);

            amrex::Real Tloc = prim(i,j,k,4);
            amrex::Real rholoc = prim(i,j,k,0);
            amrex::Real Yloc[MAX_SPECIES] = {0.0};
            amrex::Real Xloc[MAX_SPECIES] = {0.0};
            for (int ns = 0; ns < nspecies; ++ns) {
                Yloc[ns] = Yk_fixed[ns];
                Xloc[ns] = Xk_fixed[ns];
            }

            const bool get_xi = true, get_mu = true, get_lam = true,
                  get_Ddiag = true, get_chi = true;
            amrex::Real muloc, xiloc, lamloc;
            amrex::Real Dbinloc[MAX_SPECIES*MAX_SPECIES], chi_mixloc[MAX_SPECIES] = {0.0};

            auto trans = pele::physics::PhysicsType::transport();
            trans.transport_full(get_xi, get_mu, get_lam, get_Ddiag, get_chi, 
                                 Tloc, rholoc, Yloc, 
                                 Dbinloc, chi_mixloc, muloc, xiloc, lamloc, ltransparm);
            amrex::GpuArray<amrex::Real,MAX_SPECIES*MAX_SPECIES> Dloc;
            D_GIO1(Dbinloc,Yk_fixed,Xk_fixed,Dloc,nspecies);
            eta(i,j,k) = muloc;
            kappa(i,j,k) = lamloc;
            zeta(i,j,k) = xiloc;
            for (int kk=0; kk<nspecies; ++kk) {
                chi(i,j,k,kk) = chi_mixloc[kk]/Xloc[kk]; // re-scale thermal diffusion coefficients (pg. 24, Eq. 2.5.24, Giovangigli)
            }
            
            // want this multiplied by rho for all times (rho*D_tilde = rho*Y*D)
            for (int kk=0; kk<nspecies; ++kk) {
                for (int ll=0; ll<nspecies; ++ll) {
                    int n = kk*nspecies + ll;
                    Dij(i,j,k,n) = Dloc[n]*prim(i,j,k,0);
                }
            }

            // set bulk viscosity
            if (amrex::Math::abs(visc_type) == 3) {
                zeta(i,j,k) = zeta_ratio * eta(i,j,k);
            }


        });
    }
}
#else
void calculateTransportCoeffs(const MultiFab& prim_in,
                              MultiFab& eta_in, MultiFab& zeta_in, MultiFab& kappa_in,
                              MultiFab& chi_in, MultiFab& Dij_in)
{
    BL_PROFILE_VAR("calculateTransportCoeffs()",calculateTransportCoeffs);

    // see comments in conservedPrimitiveConversions.cpp regarding alternate ways of declaring
    // thread shared and thread private arrays on GPUs
    // if the size is not known at compile time, alternate approaches are required
    // here we know the size at compile time

    amrex::IntVect ng_temp;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        ng_temp[d] = (do_1D == 1) ? 1 : ngc[d];
    }

    // Loop over boxes
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        // grow the box by ngc
        const Box& bx = amrex::grow(mfi.tilebox(), ng_temp);

        const Array4<const Real>& prim = prim_in.array(mfi);

        const Array4<Real>& eta   =   eta_in.array(mfi);
        const Array4<Real>& zeta  =  zeta_in.array(mfi);
        const Array4<Real>& kappa = kappa_in.array(mfi);
        const Array4<Real>& chi   =   chi_in.array(mfi);
        const Array4<Real>& Dij   =   Dij_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            GpuArray<Real,MAX_SPECIES> Yk_fixed;
            GpuArray<Real,MAX_SPECIES> Xk_fixed;

            Real sumYk = 0.;
            for (int n=0; n<nspecies; ++n) {
                if (prim(i,j,k,6+n) <= 0.0) amrex::Abort("Negative mass fraction encountered");
                if (prim(i,j,k,6+n) >= 1.0) amrex::Abort("Greater than unity mass fraction encountered");
                Yk_fixed[n] = amrex::max(0.,amrex::min(1.,prim(i,j,k,6+n)));
                sumYk += Yk_fixed[n];
            }

            for (int n=0; n<nspecies; ++n) {
                Yk_fixed[n] /= sumYk;
            }

            // compute mole fractions from mass fractions
            GetMolfrac(Yk_fixed, Xk_fixed);

            amrex::GpuArray<amrex::Real,MAX_SPECIES*MAX_SPECIES> Dloc;
            amrex::GpuArray<amrex::Real,MAX_SPECIES> chiloc;

            if (transport_type == 1) { // Giovangigli
                IdealMixtureTransportGIO(prim(i,j,k,0), prim(i,j,k,4), prim(i,j,k,5),
                                         Yk_fixed, eta(i,j,k), kappa(i,j,k), zeta(i,j,k),
                                         Dloc, chiloc);
            }

            else if (transport_type == 2) { // Waldmann-Valk
                IdealMixtureTransportVW(prim(i,j,k,0), prim(i,j,k,4), prim(i,j,k,5),
                                      Yk_fixed, Xk_fixed, eta(i,j,k), kappa(i,j,k), zeta(i,j,k),
                                      Dloc, chiloc);
            }
            else if (transport_type == 3) { // Hirschfelder-Curtiss-Bird for binary mixtures
                IdealMixtureTransportHCBBin(prim(i,j,k,0), prim(i,j,k,4), prim(i,j,k,5),
                                            Yk_fixed, Xk_fixed, eta(i,j,k), kappa(i,j,k), zeta(i,j,k),
                                            Dloc, chiloc);
            }
            for (int kk=0; kk<nspecies; ++kk) {
                chi(i,j,k,kk) = chiloc[kk];
            }
            // want this multiplied by rho for all times (rho*D_tilde = rho*Y*D)
            for (int kk=0; kk<nspecies; ++kk) {
                for (int ll=0; ll<nspecies; ++ll) {
                    int n = kk*nspecies + ll;
                    Dij(i,j,k,n) = Dloc[n]*prim(i,j,k,0);
                }
            }

            // set bulk viscosity
            if (amrex::Math::abs(visc_type) == 3) {
                zeta(i,j,k) = zeta_ratio * eta(i,j,k);
            }


        });
    }
}
#endif
