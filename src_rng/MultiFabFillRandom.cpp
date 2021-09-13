#include "common_functions.H"

#include "rng_functions.H"

void MultiFabFillRandom(MultiFab& mf, const int& comp, const amrex::Real& variance,
                        const Geometry& geom, const int& ng)
{
    BL_PROFILE_VAR("MultiFabFillRandom()",MultiFabFillRandom);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = (ng==0) ? mfi.validbox() : mfi.growntilebox(ng);
        const Array4<Real>& mf_fab = mf.array(mfi);
        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
        {
	  mf_fab(i,j,k,comp) = amrex::RandomNormal(0.,1.,engine);
        });
    }

//----------------------------------------

    // Scale standard gaussian samples by standard deviation
    mf.mult(sqrt(variance), comp, 1, 0);

    // sync up random numbers of faces/nodes that are at the same physical location
    mf.OverrideSync(geom.periodicity());

    // fill interior and periodic ghost cells
    mf.FillBoundary(geom.periodicity());

//----------------------------------------
}
