#include <cmath>

#include <AMReX_GpuContainers.H>
#include <AMReX_ParticleReduce.H>
#include <AMReX_Random.H>

#include "InertialPC.H"

using namespace amrex;

void
InertialPC::InitParticles (InertialParams const& ip)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto dx  = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();

    auto cell_weight = [&] (int i, int j) -> Real {
        return inertial_ic(plo[0] + (i+0.5)*dx[0], plo[1] + (j+0.5)*dx[1], ip);
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
                const int n = static_cast<int>(ip.num_part*cell_weight(i, j)/total + amrex::Random());
                for (int ip_ = 0; ip_ < n; ++ip_) {
                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = my_cpu;
                    p.pos(0) = plo[0] + (i + amrex::Random())*dx[0];
                    p.pos(1) = plo[1] + (j + amrex::Random())*dx[1];
                    ptile.push_back(p);
                }
            }
        }
    }

    Redistribute();
    amrex::Print() << "Initialized " << TotalNumberOfParticles() << " particles (num_part = "
                   << ip.num_part << ")\n";
}

InertialPC::MoveStats
InertialPC::Advance (Real dt, InertialParams const& ip)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto dx = geom.CellSizeArray();
    const Real vmax = ip.vmax;
    // keep v strictly inside the domain so the particle stays in a valid cell
    const Real vin = vmax*(1.0 - 1.e-12);
    const Real sd = std::sqrt(2.0*ip.diff_coeff*dt);
    // BAOAB: exact Ornstein-Uhlenbeck factor and noise amplitude of the O step
    const int baoab = ip.particle_integrator;
    const Real c_ou = std::exp(-ip.gamma*dt);
    const Real s_ou = std::sqrt(inertial_T(ip)*(1.0 - c_ou*c_ou));
    const InertialParams ip_loc = ip;

    Gpu::DeviceVector<Long> counts(3, 0);
    Long* cnt = counts.data();

    for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
        auto& aos = pti.GetArrayOfStructs();
        ParticleType* pstruct = aos().data();
        const int np = aos.numParticles();
        amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int n, RandomEngine const& engine) noexcept
        {
            ParticleType& p = pstruct[n];
            const Real x = p.pos(0);
            const Real v = p.pos(1);
            Real mx, mv;
            if (baoab) {
                // B A O A B: half kick, half drift, exact OU step for friction
                // and noise, half drift, half kick
                Real vv = v - 0.5*dt*inertial_dU(x, ip_loc);
                Real xx = x + 0.5*dt*vv;
                vv = c_ou*vv + s_ou*amrex::RandomNormal(0.0, 1.0, engine);
                xx += 0.5*dt*vv;
                vv -= 0.5*dt*inertial_dU(xx, ip_loc);
                mx = xx - x;
                mv = vv - v;
            } else {
                // Euler-Maruyama with the start-of-step x and v
                mx = v*dt;
                mv = (-ip_loc.gamma*v - inertial_dU(x, ip_loc))*dt
                   + sd*amrex::RandomNormal(0.0, 1.0, engine);
            }
            if (std::abs(mx) > dx[0]) {
                mx = std::copysign(dx[0], mx);
                Gpu::Atomic::AddNoRet(&cnt[0], Long(1));
            }
            if (std::abs(mv) > dx[1]) {
                mv = std::copysign(dx[1], mv);
                Gpu::Atomic::AddNoRet(&cnt[1], Long(1));
            }
            Real vn = v + mv;
            if (vn > vmax || vn < -vmax) {
                vn = (vn > vmax) ? 2.0*vmax - vn : -2.0*vmax - vn;
                Gpu::Atomic::AddNoRet(&cnt[2], Long(1));
            }
            p.pos(0) = x + mx;
            p.pos(1) = amrex::min(amrex::max(vn, -vin), vin);
        });
    }

    // periodic wrap in x and move particles to their new boxes
    Redistribute();

    Vector<Long> h(3);
    Gpu::copy(Gpu::deviceToHost, counts.begin(), counts.end(), h.begin());
    ParallelDescriptor::ReduceLongSum(h.data(), 3);
    MoveStats st;
    st.trunc_x = h[0];
    st.trunc_v = h[1];
    st.reflect = h[2];
    return st;
}

void
InertialPC::DepositHistogram (MultiFab& hist) const
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto dx  = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();

    hist.setVal(0.0);
    for (ParConstIterType pti(*this, lev); pti.isValid(); ++pti) {
        const auto& aos = pti.GetArrayOfStructs();
        const ParticleType* pstruct = aos().data();
        const int np = aos.numParticles();
        auto const& h = hist.array(pti);
        const Box vbx = pti.validbox();
        const auto lo = lbound(vbx);
        const auto hi = ubound(vbx);
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int n) noexcept
        {
            const ParticleType& p = pstruct[n];
            int i = static_cast<int>(std::floor((p.pos(0) - plo[0])/dx[0]));
            int j = static_cast<int>(std::floor((p.pos(1) - plo[1])/dx[1]));
            i = amrex::min(amrex::max(i, lo.x), hi.x);
            j = amrex::min(amrex::max(j, lo.y), hi.y);
            Gpu::Atomic::AddNoRet(&h(i,j,0,0), Real(1.0));
        });
    }
    const Real ntot = static_cast<Real>(TotalNumberOfParticles());
    hist.mult(1.0/(ntot*dx[0]*dx[1]), 0, 1);
}

void
InertialPC::Moments (InertialParams const& ip, Real m[4]) const
{
    const InertialParams ip_loc = ip;
    using PType = InertialPC::ParticleType;
    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_ops;
    auto r = amrex::ParticleReduce<ReduceData<Real, Real, Real, Real>>(*this,
        [=] AMREX_GPU_DEVICE (const PType& p) -> GpuTuple<Real, Real, Real, Real>
        {
            const Real d = inertial_dx(p.pos(0), ip_loc.pot_x0, ip_loc);
            const Real v = p.pos(1);
            return {d, v, d*d, v*v};
        }, reduce_ops);
    // ParticleReduce is local to this rank
    Real s[4] = {amrex::get<0>(r), amrex::get<1>(r), amrex::get<2>(r), amrex::get<3>(r)};
    ParallelDescriptor::ReduceRealSum(s, 4);
    const Real n = static_cast<Real>(TotalNumberOfParticles());
    for (int c = 0; c < 4; ++c) { m[c] = s[c]/n; }
}
