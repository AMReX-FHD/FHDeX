#include "common_functions.H"
#include "chemistry_functions.H"
#include "rng_functions.H"

void RK3step_chem_only(MultiFab& rho_old, MultiFab& rho_new,
                       const amrex::Geometry geom, const amrex::Real dt)
{
    if (reaction_type!=1) amrex::Abort("RK3step_chem_only assumes reaction_type=1");

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    BoxArray ba = rho_old.boxArray();
    DistributionMapping dm = rho_old.DistributionMap();

    MultiFab rhop(ba,dm,nspecies,0);
    MultiFab rhop2(ba,dm,nspecies,0);
    MultiFab rhop3(ba,dm,nspecies,0);;

    MultiFab source(ba,dm,nspecies,0);;

    MultiFab ranchem(ba,dm,nreaction,0);
    MultiFab ranchem_A(ba,dm,nreaction,0);
    MultiFab ranchem_B(ba,dm,nreaction,0);
    
    // weights for stochastic fluxes; swgt2 changes each stage
    amrex::Real swgt1, swgt2;
    swgt1 = 1.;

    // initialize white noise fields
    for (int m=0;m<nreaction;m++) {
        MultiFabFillRandom(ranchem_A,m,1.,geom);
        MultiFabFillRandom(ranchem_B,m,1.,geom);
    }

    // stage1
    swgt2 = (2.*std::sqrt(2.)+std::sqrt(3.))/5.;

    MultiFab::LinComb(ranchem,
                swgt1, ranchem_A, 0,
                swgt2, ranchem_B, 0,
                0, nreaction, 0);

    compute_chemistry_source_CLE_2(dt,dx[0]*dx[1]*dx[2],rho_old,0,source,0,ranchem);

    for ( MFIter mfi(rho_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & rho_old_fab = rho_old.array(mfi);
        const Array4<Real> & rhop_fab = rhop.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            rhop_fab(i,j,k,n) = rho_old_fab(i,j,k,n) + dt*source_fab(i,j,k,n);
        });
    }

    // stage2
    swgt2 = (-4.*std::sqrt(2.)+3.*std::sqrt(3.))/5.;

    MultiFab::LinComb(ranchem,
                swgt1, ranchem_A, 0,
                swgt2, ranchem_B, 0,
                0, nreaction, 0);

    compute_chemistry_source_CLE_2(dt,dx[0]*dx[1]*dx[2],rhop,0,source,0,ranchem);

    for ( MFIter mfi(rho_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & rho_old_fab = rho_old.array(mfi);
        const Array4<Real> & rhop_fab = rhop.array(mfi);
        const Array4<Real> & rhop2_fab = rhop2.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            rhop2_fab(i,j,k,n) = 0.25*( 3.*rho_old_fab(i,j,k,n) + rhop_fab(i,j,k,n) + dt*source_fab(i,j,k,n) );
        });
    }

    // stage3
    swgt2 = (std::sqrt(2.)-2.*std::sqrt(3.))/10.;

    MultiFab::LinComb(ranchem,
                swgt1, ranchem_A, 0,
                swgt2, ranchem_B, 0,
                0, nreaction, 0);

    compute_chemistry_source_CLE_2(dt,dx[0]*dx[1]*dx[2],rhop2,0,source,0,ranchem);

    for ( MFIter mfi(rho_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & rho_old_fab = rho_old.array(mfi);
        const Array4<Real> & rho_new_fab = rho_new.array(mfi);
        const Array4<Real> & rhop2_fab = rhop2.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            rho_new_fab(i,j,k,n) = (2./3.)*( 0.5*rho_old_fab(i,j,k,n) + rhop2_fab(i,j,k,n) + dt*source_fab(i,j,k,n) );
        });
    }
}
