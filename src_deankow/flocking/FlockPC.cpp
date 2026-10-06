#include <cmath>

#include <AMReX_Random.H>

#include "FlockPC.H"

using namespace amrex;

void
FlockPC::InitParticles (FlockingParams const& fp, int nq)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto dx  = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const Real dthq = flock_two_pi/nq;

    // integral of f over a cell: midpoint in x, y, nq-point rule in theta
    auto cell_weight = [&] (int i, int j) -> Real {
        const Real x = plo[0] + (i+0.5)*dx[0];
        const Real y = plo[1] + (j+0.5)*dx[1];
        Real w = 0.0;
        for (int q = 0; q < nq; ++q) { w += flocking_ic(x, y, (q+0.5)*dthq, fp); }
        return w*dthq*dx[0]*dx[1];
    };

    Real total = 0.0;
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) { total += cell_weight(i, j); }
        }
    }
    ParallelDescriptor::ReduceRealSum(total);

    const int my_cpu = ParallelDescriptor::MyProc();
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto& ptile = GetParticles(lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
            for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                const int n = static_cast<int>(fp.num_part*cell_weight(i, j)/total + amrex::Random());
                for (int ip = 0; ip < n; ++ip) {
                    const Real x = plo[0] + (i + amrex::Random())*dx[0];
                    const Real y = plo[1] + (j + amrex::Random())*dx[1];
                    // theta from f(x, y, .) by rejection, bound from nq samples
                    Real fmax = 0.0;
                    for (int q = 0; q < nq; ++q) { fmax = amrex::max(fmax, flocking_ic(x, y, (q+0.5)*dthq, fp)); }
                    fmax *= 1.1;
                    Real th = flock_two_pi*amrex::Random();
                    if (fmax > 0.0) {
                        while (amrex::Random()*fmax > flocking_ic(x, y, th, fp)) {
                            th = flock_two_pi*amrex::Random();
                        }
                    }
                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = my_cpu;
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.rdata(FlockRealIdx::theta) = th;
                    ptile.push_back(p);
                }
            }
        }
    }

    Redistribute();
    amrex::Print() << "Initialized " << TotalNumberOfParticles() << " particles (num_part = "
                   << fp.num_part << ")\n";
}

void
FlockPC::Advance (Real dt, FlockingParams const& fp)
{
    const int lev = 0;
    const Real v  = fp.speed;
    const Real sd = std::sqrt(2.0*fp.diff_coeff*dt);

    for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
        auto& aos = pti.GetArrayOfStructs();
        ParticleType* pstruct = aos().data();
        const int np = aos.numParticles();
        amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int ip, RandomEngine const& engine) noexcept
        {
            ParticleType& p = pstruct[ip];
            Real th = p.rdata(FlockRealIdx::theta);
            p.pos(0) += v*std::cos(th)*dt;
            p.pos(1) += v*std::sin(th)*dt;
            th += sd*amrex::RandomNormal(0.0, 1.0, engine);
            th -= flock_two_pi*std::floor(th/flock_two_pi);
            if (th >= flock_two_pi) { th = 0.0; }   // guard against round-off
            p.rdata(FlockRealIdx::theta) = th;
        });
    }

    // periodic wrap in x and y and move particles to their new boxes
    Redistribute();
}

void
FlockPC::DepositHistogram (MultiFab& hist, int ntheta) const
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto dx  = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const Real dth = flock_two_pi/ntheta;

    hist.setVal(0.0);
    for (ParConstIterType pti(*this, lev); pti.isValid(); ++pti) {
        const auto& aos = pti.GetArrayOfStructs();
        const ParticleType* pstruct = aos().data();
        const int np = aos.numParticles();
        auto const& h = hist.array(pti);
        const Box vbx = pti.validbox();
        const auto lo = lbound(vbx);
        const auto hi = ubound(vbx);
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int ip) noexcept
        {
            const ParticleType& p = pstruct[ip];
            int i = static_cast<int>(std::floor((p.pos(0) - plo[0])/dx[0]));
            int j = static_cast<int>(std::floor((p.pos(1) - plo[1])/dx[1]));
            int k = static_cast<int>(p.rdata(FlockRealIdx::theta)/dth);
            i = amrex::min(amrex::max(i, lo.x), hi.x);
            j = amrex::min(amrex::max(j, lo.y), hi.y);
            k = amrex::min(amrex::max(k, 0), ntheta-1);
            amrex::Gpu::Atomic::AddNoRet(&h(i,j,0,k), Real(1.0));
        });
    }

    // counts -> density normalized to integrate to 1, then the theta moments
    const Real ntot = static_cast<Real>(TotalNumberOfParticles());
    const Real scale = 1.0/(ntot*dx[0]*dx[1]*dth);
    for (MFIter mfi(hist); mfi.isValid(); ++mfi) {
        auto const& h = hist.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            Real rho = 0.0, px = 0.0, py = 0.0;
            for (int k = 0; k < ntheta; ++k) {
                const Real phik = h(i,j,kk,k)*scale;
                h(i,j,kk,k) = phik;
                const Real th = (k+0.5)*dth;
                rho += phik*dth;
                px  += phik*std::cos(th)*dth;
                py  += phik*std::sin(th)*dth;
            }
            h(i,j,kk,ntheta)   = rho;
            h(i,j,kk,ntheta+1) = px;
            h(i,j,kk,ntheta+2) = py;
        });
    }
}
