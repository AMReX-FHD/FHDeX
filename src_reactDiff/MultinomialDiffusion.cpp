#include "reactDiff_functions.H"

#include "AMReX_MLMG.H"
#include <AMReX_MLABecLaplacian.H>

#include <random>

void MultinomialDiffusion(MultiFab& n_old,
                          MultiFab& n_new,
                          const std::array< MultiFab, AMREX_SPACEDIM >& diff_coef_face,
                          const Geometry& geom,
                          const Real& dt,
                          const Real& time)
{
#if (AMREX_USE_CUDA)
    Abort("std::MultinomailDiffusion not supported for CUDA (need sum reductions)");
#endif

    BoxArray ba = n_old.boxArray();
    DistributionMapping dmap = n_old.DistributionMap();

    MultiFab cell_update(ba, dmap, nspecies, 1);
    cell_update.setVal(0.);

    // set new state to zero everywhere, including ghost cells
    n_new.setVal(0.);

    // copy old state into new in valid region only
    MultiFab::Copy(n_new,n_old,0,0,nspecies,0);

    const GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    Real dv = (AMREX_SPACEDIM==2) ? dx[0]*dx[1]*cell_depth : dx[0]*dx[1]*dx[2]*cell_depth;

    for (MFIter mfi(n_new); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real> & n_arr = n_new.array(mfi);

        const Array4<Real> & update = cell_update.array(mfi);

        AMREX_D_TERM(const Array4<const Real> & diffx = diff_coef_face[0].array(mfi);,
                     const Array4<const Real> & diffy = diff_coef_face[1].array(mfi);,
                     const Array4<const Real> & diffz = diff_coef_face[2].array(mfi););

        amrex::ParallelForRNG(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n, amrex::RandomEngine const& engine) noexcept
        {

            GpuArray<Real,2*AMREX_SPACEDIM> p;
            GpuArray<Real,2*AMREX_SPACEDIM> fluxes;

            p[0] = diffx(i  ,j  ,k,n)*dt/(dx[0]*dx[0]);
            p[1] = diffx(i+1,j  ,k,n)*dt/(dx[0]*dx[0]);
            p[2] = diffy(i  ,j  ,k,n)*dt/(dx[1]*dx[1]);
            p[3] = diffy(i  ,j+1,k,n)*dt/(dx[1]*dx[1]);
#if (AMREX_SPACEDIM == 3)
            p[4] = diffz(i  ,j  ,k  ,n)*dt/(dx[2]*dx[2]);
            p[5] = diffz(i  ,j  ,k+1,n)*dt/(dx[2]*dx[2]);
#endif

            int N = std::max(0., std::round(n_arr(i,j,k,n)*dv));

            multinomial_rng(fluxes, N, p, engine);

            // lo-x face
            update(i  ,j,k,n) -= fluxes[0];
            update(i-1,j,k,n) += fluxes[0];

            // hi-x face
            update(i  ,j,k,n) -= fluxes[1];
            update(i+1,j,k,n) += fluxes[1];

            // lo-y face
            update(i,j,  k,n) -= fluxes[2];
            update(i,j-1,k,n) += fluxes[2];

            // hi-y face
            update(i,j  ,k,n) -= fluxes[3];
            update(i,j+1,k,n) += fluxes[3];

#if (AMREX_SPACEDIM == 3)
            // lo-z face
            update(i,j,k,  n) -= fluxes[4];
            update(i,j,k-1,n) += fluxes[4];

            // hi-z face
            update(i,j,k,  n) -= fluxes[5];
            update(i,j,k+1,n) += fluxes[5];
#endif
        });
    }

    for (MFIter mfi(n_new); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(1);

        const Array4<Real> & n_arr = n_new.array(mfi);

        const Array4<Real> & update = cell_update.array(mfi);

        amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            n_arr(i,j,k,n) += update(i,j,k,n) / dv;
        });
    }

    n_new.SumBoundary(geom.periodicity());
    n_new.FillBoundary(geom.periodicity());
    MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);
}

AMREX_GPU_HOST_DEVICE void multinomial_rng(GpuArray<Real,2*AMREX_SPACEDIM>& samples,
                                           const int& N,
                                           GpuArray<Real,2*AMREX_SPACEDIM>& p,
                                           const amrex::RandomEngine& engine)
{
    Real sum_p = 0;
    for (int sample=0; sample<2*AMREX_SPACEDIM; ++sample) {
        sum_p += p[sample];
    }
    if (sum_p > 1.) {
        printf("sum_p = %f",sum_p);
        Abort("multinomial_rng: probabilities must sum to 1 or less");
    }

    // brute force multinomial
    for (int sample=0; sample<2*AMREX_SPACEDIM; ++sample) {
        samples[sample] = 0.;
    }
    for (int n=0; n<N; ++n) {
        Real x = amrex::Random(engine); // uniform over [0,1)
        Real sum_p = 0.;
        // find the multinomial bin the RNG lands in
        for (int sample=0; sample<2*AMREX_SPACEDIM; ++sample) {
            sum_p += p[sample];
            if (x <= sum_p) {
                samples[sample] += 1.;
                break;
            }
        }
    }

#if 0
    // not sure why std:: binomial_distribition gives grid artifacts
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
