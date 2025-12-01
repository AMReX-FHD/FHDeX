#include "FPU.H"

#include <AMReX_MultiFabUtil.H>

using namespace amrex;

void init_r(MultiFab& state_r,
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

      """
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
    """

    # Set up RNG
    if isinstance(random_state, np.random.Generator):
        rng = random_state
    else:
        rng = np.random.default_rng(random_state)

    def grad_U(r):
        """
        Gradient of U(r) = β (V(r) + P r)
        with V'(r) = a r + b r^2 + c r^3.
        """
        dV = a * r + b * r**2 + c * r**3
        return beta * (dV + P)

    r = float(r0)
    total_steps = burn_in + n_steps * thin
    samples = np.empty(n_steps, dtype=float)

    k_sample = 0
    sqrt_2dt = np.sqrt(2.0 * step_size)

    for k in range(total_steps):
        # Langevin update
        r += -step_size * grad_U(r) + sqrt_2dt * rng.normal()

        # Store samples after burn-in, with thinning
        if k >= burn_in and ((k - burn_in) % thin == 0):
            samples[k_sample] = r
            k_sample += 1
            if k_sample >= n_steps:
                break

    return samples
*/
    
    Real sqrt_2dt = std::sqrt(2.*step_size);

    for ( MFIter mfi(state_r); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        const Array4<Real> r_fab = state_r.array(mfi);

        for (auto k = lo.z; k <= hi.z; ++k) {
        for (auto j = lo.y; j <= hi.y; ++j) {

            Real r = r0;
            Real random;
            for (auto aaa=0; aaa<burn; ++aaa) {
                random = amrex::RandomNormal(0.,1.);
            }
            
            for (auto i = lo.x; i <= hi.x; ++i) {
                Real grad_U = beta * (a*r + b*r*r + c*r*r*r + pressure);
                r += -step_size * grad_U + sqrt_2dt * amrex::RandomNormal(0.,1.);
                r_fab(i,j,k) = r;
            }
        }
        }
        
    }
        
    state_r.FillBoundary(geom.periodicity());

    // compute mean
    BoxArray ba = state_r.boxArray();
    Box domain = ba.minimalBox().enclosedCells();
    
    Gpu::HostVector<Real> sum_r(n_ensembles);
    sum_r = sumToLine(state_r, 0, 1, domain, 1);

    for (int i=0; i<n_ensembles; ++i) {
        Print() << "For ensemble " << i << " The mean stretch r is: " << sum_r[i]/n_particles << std::endl;
    }
}

void init_p(MultiFab& state_p,
            const Real& beta,
            const int& n_particles,
            const int& n_ensembles,
            const Geometry& geom) {

    Real sigma = std::sqrt(1./beta);

    for (MFIter mfi(state_p); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real>& p_fab = state_p.array(mfi);

        amrex::ParallelForRNG(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
        {
            p_fab(i,j,k) = amrex::RandomNormal(0.,sigma,engine);
        });
    }

    state_p.FillBoundary(geom.periodicity());
    
    // compute mean
    BoxArray ba = state_p.boxArray();
    Box domain = ba.minimalBox().enclosedCells();

    Gpu::HostVector<Real> sum_p(n_ensembles);
    sum_p = sumToLine(state_p, 0, 1, domain, 1);

    for (int i=0; i<n_ensembles; ++i) {
        Print() << "For ensemble " << i << " The mean momentum p is: " << sum_p[i]/n_particles << std::endl;
    }

    
}
