#include "reactDiff_functions.H"

#include "AMReX_MLMG.H"
#include <AMReX_MLABecLaplacian.H>

#include <random>

void MultinomialDiffusion(MultiFab& n_in,
                          const std::array< MultiFab, AMREX_SPACEDIM >& diff_coef_face,
                          const Geometry& geom,
                          const Real& dt) {

    const GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    Real dv = (AMREX_SPACEDIM==2) ? dx[0]*dx[1]*cell_depth : dx[0]*dx[1]*dx[2]*cell_depth;

    for (MFIter mfi(n_in); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real> & n_arr = n_in.array(mfi);
        
        AMREX_D_TERM(const Array4<const Real> & diffx = diff_coef_face[0].array(mfi);,
                     const Array4<const Real> & diffy = diff_coef_face[1].array(mfi);,
                     const Array4<const Real> & diffz = diff_coef_face[2].array(mfi););

        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {

            GpuArray<Real,2*AMREX_SPACEDIM> p;
            GpuArray<Real,2*AMREX_SPACEDIM> fluxes;

            p[0] = diffx(i,j,k,n)*dt/(dx[0]*dx[0]);
            p[1] = diffx(i+1,j,k,n)*dt/(dx[0]*dx[0]);
            p[2] = diffy(i,j,k,n)*dt/(dx[0]*dx[1]);
            p[3] = diffy(i,j+1,k,n)*dt/(dx[0]*dx[1]);
#if (AMREX_SPACEDIM == 3)
            p[4] = diffz(i,j,k,n)*dt/(dx[0]*dx[2]);
            p[5] = diffz(i,j,k+1,n)*dt/(dx[0]*dx[2]);
#endif

            int N = std::max(0., n_arr(i,j,k,n)*dv);
            
            multinomial_rng(fluxes, N, p);

        });
    }
}

AMREX_GPU_HOST_DEVICE void multinomial_rng(GpuArray<Real,2*AMREX_SPACEDIM>& samples,
                                           const int& N,
                                           GpuArray<Real,2*AMREX_SPACEDIM>& p)
{
#if (AMREX_USE_CUDA)
    Abort("MultinomialRNG not supported for CUDA");
#else

    Real sum_p = 0;
    for (int sample=0; sample<2*AMREX_SPACEDIM; ++sample) {
        sum_p += p[sample];
    }
    if (sum_p > 1.) {
        Abort("multinomial_rng: probabilities must sum to 1 or less");
    }

    std::default_random_engine generator;
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());

    sum_p = 0.;
    int sum_n = 0;
    
    for (int sample=0; sample<2*AMREX_SPACEDIM; ++sample) {
        std::binomial_distribution<int> distribution(N-sum_n, p[sample]/(1.-sum_p));
        samples[sample] = distribution(generator);

        sum_n += samples[sample];
        sum_p += p[sample];
    }

#endif
}
