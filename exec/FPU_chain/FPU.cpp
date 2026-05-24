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

        const Array4<Real>& state_fab = state.array(mfi);

        for (auto k = lo.z; k <= hi.z; ++k) {
        for (auto j = lo.y; j <= hi.y; ++j) {

            Real r = r0;
            Real random;

            // discard burn random numbers
            for (auto aaa=0; aaa<burn; ++aaa) {
                random = amrex::RandomNormal(0.,1.);
            }

            for (auto i = lo.x; i <= hi.x; ++i) {
                for (int kkk=0; kkk<1000; kkk++) {
                    Real grad_U = (a*r + b*r*r + c*r*r*r + pressure);
                    r += -step_size * beta * grad_U + sqrt_2dt * amrex::RandomNormal(0.,1.);
                }
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
}



void FPU_RK4(MultiFab& state,
             const Real& a,
             const Real& b,
             const Real& c,
             const Real& dt,
             const int& n_particles,
             const int& n_ensembles,
             const Geometry& geom) {

    BoxArray ba = state.boxArray();
    DistributionMapping dm = state.DistributionMap();

    // we only need a ghost cell for the x-direction
    IntVect ng_vect(1,0);

    // storage for intermediate states
    MultiFab state2(ba,dm,2,ng_vect);
    MultiFab state3(ba,dm,2,ng_vect);
    MultiFab state4(ba,dm,2,ng_vect);

    // storage for rhs at each stage
    MultiFab rhs1(ba,dm,2,0);
    MultiFab rhs2(ba,dm,2,0);
    MultiFab rhs3(ba,dm,2,0);
    MultiFab rhs4(ba,dm,2,0);

    // compute rhs1 = f(rhs)
    rhs(rhs1,state,a,b,c,geom);

    // state2 = state + 0.5*dt*rhs1
    MultiFab::LinComb(state2,1.0,state,0,0.5*dt,rhs1,0,0,2,0);
    state2.FillBoundary(geom.periodicity());

    // compute rhs2 = f(rhs2)
    rhs(rhs2,state2,a,b,c,geom);

    // state3 = state + 0.5*dt*rhs2
    MultiFab::LinComb(state3,1.0,state,0,0.5*dt,rhs2,0,0,2,0);
    state3.FillBoundary(geom.periodicity());

    // compute rhs3 = f(rhs3)
    rhs(rhs3,state3,a,b,c,geom);

    // state4 = state + dt*rhs3
    MultiFab::LinComb(state4,1.0,state,0,dt,rhs3,0,0,2,0);
    state4.FillBoundary(geom.periodicity());

    // compute rhs4 = f(rhs4)
    rhs(rhs4,state4,a,b,c,geom);

    // RK4 update: state += (dt/6.0) * (rhs1 + 2*rhs2 + 2*rhs3 + rhs4)
    MultiFab::Saxpy(state,dt/6.0,rhs1,0,0,2,0);
    MultiFab::Saxpy(state,dt/3.0,rhs2,0,0,2,0);
    MultiFab::Saxpy(state,dt/3.0,rhs3,0,0,2,0);
    MultiFab::Saxpy(state,dt/6.0,rhs4,0,0,2,0);

    state.FillBoundary(geom.periodicity());

}

void rhs(MultiFab& rhs,
         MultiFab& state,
         const Real& a,
         const Real& b,
         const Real& c,
         const Geometry& geom) {

    BoxArray ba = state.boxArray();
    DistributionMapping dm = state.DistributionMap();

    MultiFab V_prime(ba,dm,1,1);

    // compute V'
    for (MFIter mfi(state); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<const Real>&   state_fab = state.array(mfi);
        const Array4<      Real>& V_prime_fab = V_prime.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real r = state_fab(i,j,k,0);
            V_prime_fab(i,j,k,0) = a*r + b*r*r + c*r*r*r;
        });
    }

    V_prime.FillBoundary(geom.periodicity());

    // compute rhs
    for (MFIter mfi(state); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<const Real>&   state_fab = state.array(mfi);
        const Array4<const Real>& V_prime_fab = V_prime.array(mfi);
        const Array4<      Real>&     rhs_fab = rhs.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // dr_dt = p(i+1) - p(i)
            rhs_fab(i,j,k,0) = state_fab(i+1,j,k,1) - state_fab(i,j,k,1);

            // dp_dt = V'(i) - V'(i-1)
            rhs_fab(i,j,k,1) = V_prime_fab(i,j,k,0) - V_prime_fab(i-1,j,k,0);
        });
    }
}

void compute_means(MultiFab& state,
                   const int& n_particles,
                   const int& n_ensembles,
                   const int& step) {

    BoxArray ba = state.boxArray();
    Box domain = ba.minimalBox().enclosedCells();

    // compute mean r
    Gpu::HostVector<Real> sum_r(n_ensembles);
    sum_r = sumToLine(state, 0, 1, domain, 1);

    // compute mean p
    Gpu::HostVector<Real> sum_p(n_ensembles);
    sum_p = sumToLine(state, 1, 1, domain, 1);

    // compute mean energy
    Gpu::HostVector<Real> sum_e(n_ensembles);
    sum_e = sumToLine(state, 2, 1, domain, 1);

    const std::string meanfile = amrex::Concatenate("mean",step,7);
    amrex::Print() << "Writing means " << meanfile << std::endl;

    std::ofstream meanout;
    if (ParallelDescriptor::IOProcessor()) {
        meanout.open(meanfile, std::ios::out);
        for (int i=0; i<n_ensembles; ++i) {
            meanout << "For ensemble " << i << " the mean r/p/e is: "
                    << sum_r[i]/n_particles << " "
                    << sum_p[i]/n_particles << " "
                    << sum_e[i]/n_particles << "\n";
        }
    }
}

void compute_energy(MultiFab& state,
                    const Real& a,
                    const Real& b,
                    const Real& c) {

    // compute energy
    for (MFIter mfi(state); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real>&  state_fab = state.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real r = state_fab(i,j,k,0);
            state_fab(i,j,k,2) = 0.5*state_fab(i,j,k,1)*state_fab(i,j,k,1) +
                (1./2.)*a*r*r + (1./3.)*b*r*r*r + (1./4.)*c*r*r*r*r;
        });
    }

}

//void compute_heat_flux(MultiFab& heat_flux,
void compute_heat_flux(Vector<Real>& heat_flux_vec,
                       const MultiFab& state,
                       const Real& a,
                       const Real& b,
                       const Real& c,
                       const Real& rest,
                       const Geometry& geom,
                       const amrex::Box& domain,
                       const int n_ensembles,
                       const int n_particles) {


    // compute heat flux
    // convective: flux[0] = (ke + pe) * v
    // conductive: flux[1] = (stress) * v
    MultiFab heat_flux(state.boxArray(),state.DistributionMap(),3,0);
    for (MFIter mfi(state); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real const>&  state_fab     = state.array(mfi);
        const Array4<Real      >&  heat_flux_fab = heat_flux.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real r_r      = state_fab(i,j,k,0);
            Real r_l      = state_fab(i-1,j,k,0);
            Real v_i      = state_fab(i,j,k,1);
            Real pe_r     = (1./2.)*a*r_r*r_r + (1./3.)*b*r_r*r_r*r_r + (1./4.)*c*r_r*r_r*r_r*r_r;
            Real pe_l     = (1./2.)*a*r_l*r_l + (1./3.)*b*r_l*r_l*r_l + (1./4.)*c*r_l*r_l*r_l*r_l;
            Real pe       = 0.5*(pe_l + pe_r);
            Real ke       = 0.5*v_i*v_i;

            // convective
            heat_flux_fab(i,j,k,0) = v_i*(pe+ke);

            // conductive
            Real dx_r     = state_fab(i,j,k,0) + rest;
            // dx_r          = rest;
            Real dx_l     = state_fab(i-1,j,k,0) + rest;
            // dx_l          = rest;
            Real fbond_r  = a*r_r + b*r_r*r_r + c*r_r*r_r*r_r;
            Real fbond_l  = a*r_l + b*r_l*r_l + c*r_l*r_l*r_l; 
            
            // Real v_ip      = state_fab(i+1,j,k,1);
            // Real v_in      = state_fab(i-1,j,k,1);

//            heat_flux_fab(i,j,k,1) = 0.5 * fbond_r * dx_r * (v_i + v_j);
            // heat_flux_fab(i,j,k,1)  = (rest/2.0) * (v_i + v_in) * fbond_l; // i-1 --- i
            // heat_flux_fab(i,j,k,1) += (rest/2.0) * (v_i + v_ip) * fbond_r; // i --- i+1
            heat_flux_fab(i,j,k,1) = 0.5*(fbond_r * dx_r - fbond_l * dx_l)*v_i;

            heat_flux_fab(i,j,k,2) = heat_flux_fab(i,j,k,0) + heat_flux_fab(i,j,k,1);

        });
    }
    heat_flux.FillBoundary(geom.periodicity());
    Gpu::HostVector<Real> heat_flux_avg(n_ensembles);
    heat_flux_avg = sumToLine(heat_flux,2,1,domain,1);
    for (int n=0; n<n_ensembles; ++n) {
        heat_flux_vec[n] = heat_flux_avg[n]/Real(n_particles);
    }

}

void ComputePhiFromState(MultiFab& phi,
                         const Real& r_eq,
                         const Real& p_eq,
                         const Real& e_eq,
                         const Real& R_00,
                         const Real& R_01,
                         const Real& R_02,
                         const Real& R_10,
                         const Real& R_11,
                         const Real& R_12,
                         const Real& R_20,
                         const Real& R_21,
                         const Real& R_22) {

    // convert state to state-state_eq
    phi.plus(-r_eq,0,1,0);
    phi.plus(-p_eq,1,1,0);
    phi.plus(-e_eq,2,1,0);

    // phi = R(state-state_eq)
    for (MFIter mfi(phi); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<Real>& phi_fab = phi.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real R0 = R_00*phi_fab(i,j,k,0) + R_01*phi_fab(i,j,k,1) + R_02*phi_fab(i,j,k,2);
            Real R1 = R_10*phi_fab(i,j,k,0) + R_11*phi_fab(i,j,k,1) + R_12*phi_fab(i,j,k,2);
            Real R2 = R_20*phi_fab(i,j,k,0) + R_21*phi_fab(i,j,k,1) + R_22*phi_fab(i,j,k,2);

            phi_fab(i,j,k,0) = R0;
            phi_fab(i,j,k,1) = R1;
            phi_fab(i,j,k,2) = R2;

        });
    }

}

void ComputeCalphaalpha(MultiFab& C_alphaalpha,
                        const MultiFab& phi,
                        const MultiFab& phi0) {

    for (MFIter mfi(C_alphaalpha); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        const Array4<      Real>& C_alphaalpha_fab = C_alphaalpha.array(mfi);
        const Array4<const Real>& phi_fab          = phi.array(mfi);
        const Array4<const Real>& phi0_fab         = phi0.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            C_alphaalpha_fab(i,j,k,0) = phi_fab(i,j,k,0) * phi0_fab(0,j,k,0);
            C_alphaalpha_fab(i,j,k,1) = phi_fab(i,j,k,1) * phi0_fab(0,j,k,1);
            C_alphaalpha_fab(i,j,k,2) = phi_fab(i,j,k,2) * phi0_fab(0,j,k,2);
        });
    }

}


