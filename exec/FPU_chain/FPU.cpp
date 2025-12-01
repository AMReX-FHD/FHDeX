#include "FPU.H"

using namespace amrex;

void init_p(MultiFab& state,
            const Real& beta) {

    Real sigma = std::sqrt(1./beta);

    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<Real>& state_fab = state.array(mfi);

        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
        {
            state_fab(i,j,k,0) = amrex::RandomNormal(0.,sigma,engine);
        });
    }

    
}
