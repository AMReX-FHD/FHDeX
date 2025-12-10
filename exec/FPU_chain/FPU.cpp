#include "FPU.H"

#include <AMReX_MultiFabUtil.H>

using namespace amrex;

void init(MultiFab& state,
          const Real& beta,
          const Real& pressure,
          const Real& a,
          const Real& b,
          const Real& c,
          const Real& r0,
          const int& burn,
          const Real& step_size,
          const int& n_particles,
          const int& n_ensembles,
          const Geometry& geom) {

/*
  Initialize r
  
  Overdamped Langevin sampler for π(r) ∝ exp[-β (V(r) + P r)],
  where V(r) = (a/2) r^2 + (b/3) r^3 + (c/4) r^4.

  Parameters
  ----------
  beta : real
    Inverse temperature β.
  P : real
    Linear coefficient in the exponent (acts like a constant force term).
  a, b, c : real
    Coefficients in the quartic potential V(r).
  r0 : real
    Initial position of the chain.
  step_size : real
    Time step Δt for the Langevin discretization.

  Returns
  -------
  samples : np.ndarray, shape (n_steps,)
  Samples approximately distributed according to π(r).
*/

    Real sqrt_2dt = std::sqrt(2.*step_size);

    for ( MFIter mfi(state); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<Real> state_fab = state.array(mfi);

        for (auto k = lo.z; k <= hi.z; ++k) {
        for (auto j = lo.y; j <= hi.y; ++j) {

            Real r = r0;
            Real random;

            // discard burn random numbers
            for (auto aaa=0; aaa<burn; ++aaa) {
                random = amrex::RandomNormal(0.,1.);
            }

            for (auto i = lo.x; i <= hi.x; ++i) {
                Real grad_U = beta * (a*r + b*r*r + c*r*r*r + pressure);
                r += -step_size * grad_U + sqrt_2dt * amrex::RandomNormal(0.,1.);
                state_fab(i,j,k,0) = r;
            }
        }
        }

    }

    // Initialize p

    Real sigma = std::sqrt(1./beta);

    for (MFIter mfi(state); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real>& state_fab = state.array(mfi);

        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
        {
            state_fab(i,j,k,1) = amrex::RandomNormal(0.,sigma,engine);
        });
    }

    state.FillBoundary(geom.periodicity());
    

    // compute means
    BoxArray ba = state.boxArray();
    Box domain = ba.minimalBox().enclosedCells();
    
    Gpu::HostVector<Real> sum_r(n_ensembles);
    sum_r = sumToLine(state, 0, 1, domain, 1);

    for (int i=0; i<n_ensembles; ++i) {
        Print() << "For ensemble " << i << " The mean stretch r is: " << sum_r[i]/n_particles << std::endl;
    }

    Gpu::HostVector<Real> sum_p(n_ensembles);
    sum_p = sumToLine(state, 1, 1, domain, 1);

    for (int i=0; i<n_ensembles; ++i) {
        Print() << "For ensemble " << i << " The mean momentum p is: " << sum_p[i]/n_particles << std::endl;
    }
}

void init_p(MultiFab& state_p,
            const Real& beta,
            const int& n_particles,
            const int& n_ensembles,
            const Geometry& geom) {

    
}
