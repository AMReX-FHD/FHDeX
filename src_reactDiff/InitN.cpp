#include "reactDiff_functions.H"

void InitN(MultiFab& n_in,
           const Geometry& geom,
           const Real& time) {

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    for ( MFIter mfi(n_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.tilebox();

        const Array4<Real> & n_init = n_in.array(mfi);

        if (prob_type == 0) {
            //============================================================
            // Thermodynamic equilibrium
            //============================================================

            amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                n_init(i,j,k,n) = n_init_in(0,n);
            });

        } else if (prob_type == 5) {
            //=================================================================
            // bubble having radius = 0.5*perturb_width*dx(1)
            // n_init = n_init_in(1,:) inside, n_init = n_init_in (2,:) outside
            // can be discontinous or smooth depending on smoothing_width
            //=================================================================

            Real rad = 0.5*perturb_width*dx[0];

            amrex::ParallelFor(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                Real x = prob_lo[0] + (i+0.5)*dx[0] - 0.5*(prob_lo[0]+prob_hi[0]);
                Real y = prob_lo[1] + (j+0.5)*dx[1] - 0.5*(prob_lo[1]+prob_hi[1]);
                Real r = std::sqrt(x*x + y*y);
#if (AMREX_SPACEDIM == 3)
                Real z = prob_lo[2] + (k+0.5)*dx[2] - 0.5*(prob_lo[2]+prob_hi[2]);
                r = std::sqrt(x*x + y*y + z*z);
#endif

                if (smoothing_width == 0.) {
                    // discontinuous interface
                    if (r < rad) {
                        n_init(i,j,k,n) = n_init_in(0,n);
                    } else {
                        n_init(i,j,k,n) = n_init_in(1,n);
                    }
                } else {
                    // smooth interface
                    n_init(i,j,k,n) = n_init_in(0,n) + (n_init_in(1,n) - n_init_in(0,n))* 0.5*(1. + std::tanh((r-rad)/(smoothing_width*dx[0])));
                }

            });

        } else {
            Abort("prob_type not implemented yet");
        }

    }

    if (integer_populations == 1) { // Ensure that the initial number of molecules are integers

        Real dv = (AMREX_SPACEDIM == 3) ? dx[0]*dx[1]*dx[2]*cell_depth : dx[0]*dx[1]*cell_depth;

        if (initial_variance_mass < 0.) { // Distribute the particles on the box using a multinomial sampler

            Abort("integer_populations=1 with initial_variance_mass < 0. not supported yet");

        } else if (initial_variance_mass > 0.) { // Make the number of molecules in each cell Poisson distributed with desired mean

            for ( MFIter mfi(n_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

                const Box& bx = mfi.tilebox();

                const Array4<Real> & n_init = n_in.array(mfi);

                amrex::ParallelForRNG(bx, nspecies, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n, amrex::RandomEngine const& engine) noexcept
                {
                    // Generate the initial fluctuations using a Poisson random number generator
                    // This assumes that the distribution of initial conditions is a product Poisson measure
                    int nparticles = RandomPoisson(n_init(i,j,k,n)*dv, engine);
                    n_init(i,j,k,n) = nparticles / dv;
                });
            }

        }
    }

    n_in.FillBoundary(geom.periodicity());
    MultiFabPhysBC(n_in, geom, 0, nspecies, SPEC_BC_COMP, time);
}
