#include "common_functions.H"
#include "chemistry_functions.H"
#include "rng_functions.H"

void EMstep_chem_only(MultiFab& rho_old, MultiFab& rho_new,
                       const amrex::Geometry geom, const amrex::Real dt)
{
    if (reaction_type!=1) amrex::Abort("EMstep_chem_only assumes reaction_type=1");

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    BoxArray ba = rho_old.boxArray();
    DistributionMapping dm = rho_old.DistributionMap();

    MultiFab source(ba,dm,nspecies,0);;
    
    MultiFab ranchem(ba,dm,nreaction,0);

    // initialize white noise field
    for (int m=0;m<nreaction;m++) {
        MultiFabFillRandom(ranchem,m,1.,geom);
    }

    compute_chemistry_source_CLE_2(dt,dx[0]*dx[1]*dx[2],rho_old,0,source,0,ranchem);

    for ( MFIter mfi(rho_old,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & rho_old_fab = rho_old.array(mfi);
        const Array4<Real> & rho_new_fab = rho_new.array(mfi);
        const Array4<Real> & source_fab = source.array(mfi);

        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            rho_new_fab(i,j,k,n) = rho_old_fab(i,j,k,n) + dt*source_fab(i,j,k,n);
        });
    }

}
