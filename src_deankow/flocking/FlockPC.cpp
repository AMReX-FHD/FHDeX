#include <cmath>

#include <AMReX_NeighborList.H>
#include <AMReX_ParticleReduce.H>
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
                    p.rdata(FlockRealIdx::qx) = 0.0;
                    p.rdata(FlockRealIdx::qy) = 0.0;
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
    const Real eps = fp.align_eps;
    const Real sd = std::sqrt(2.0*fp.diff_coeff*dt);

    for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
        auto& aos = pti.GetArrayOfStructs();
        ParticleType* pstruct = aos().data();
        const int np = aos.numParticles();
        amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int ip, RandomEngine const& engine) noexcept
        {
            ParticleType& p = pstruct[ip];
            Real th = p.rdata(FlockRealIdx::theta);
            const Real c = std::cos(th), s = std::sin(th);
            p.pos(0) += v*c*dt;
            p.pos(1) += v*s*dt;
            if (eps != 0.0) {
                th += eps*(c*p.rdata(FlockRealIdx::qy) - s*p.rdata(FlockRealIdx::qx))*dt;
            }
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
FlockPC::DepositPolarizationCIC (MultiFab& pol) const
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto dx  = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    AMREX_ALWAYS_ASSERT(pol.nGrow() >= 1 && pol.nComp() >= 2);

    pol.setVal(0.0);
    for (ParConstIterType pti(*this, lev); pti.isValid(); ++pti) {
        const auto& aos = pti.GetArrayOfStructs();
        const ParticleType* pstruct = aos().data();
        const int np = aos.numParticles();
        auto const& r = pol.array(pti);
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int ip) noexcept
        {
            const ParticleType& p = pstruct[ip];
            const Real lx = (p.pos(0) - plo[0])/dx[0] - 0.5;
            const Real ly = (p.pos(1) - plo[1])/dx[1] - 0.5;
            const int i = static_cast<int>(std::floor(lx));
            const int j = static_cast<int>(std::floor(ly));
            const Real fx = lx - i;
            const Real fy = ly - j;
            const Real th = p.rdata(FlockRealIdx::theta);
            const Real w[2] = {std::cos(th), std::sin(th)};
            for (int n = 0; n < 2; ++n) {
                amrex::Gpu::Atomic::AddNoRet(&r(i  ,j  ,0,n), w[n]*(1.0-fx)*(1.0-fy));
                amrex::Gpu::Atomic::AddNoRet(&r(i+1,j  ,0,n), w[n]*     fx *(1.0-fy));
                amrex::Gpu::Atomic::AddNoRet(&r(i  ,j+1,0,n), w[n]*(1.0-fx)*     fy );
                amrex::Gpu::Atomic::AddNoRet(&r(i+1,j+1,0,n), w[n]*     fx *     fy );
            }
        });
    }
    pol.SumBoundary(geom.periodicity());
    const Real ntot = static_cast<Real>(TotalNumberOfParticles());
    pol.mult(1.0/(ntot*dx[0]*dx[1]), 0, 2);
    pol.FillBoundary(geom.periodicity());
}

void
FlockPC::InterpolateCIC (MultiFab const& q)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto dx  = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    AMREX_ALWAYS_ASSERT(q.nGrow() >= 1 && q.nComp() >= 2);

    for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
        auto& aos = pti.GetArrayOfStructs();
        ParticleType* pstruct = aos().data();
        const int np = aos.numParticles();
        auto const& r = q.const_array(pti);
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int ip) noexcept
        {
            ParticleType& p = pstruct[ip];
            const Real lx = (p.pos(0) - plo[0])/dx[0] - 0.5;
            const Real ly = (p.pos(1) - plo[1])/dx[1] - 0.5;
            const int i = static_cast<int>(std::floor(lx));
            const int j = static_cast<int>(std::floor(ly));
            const Real fx = lx - i;
            const Real fy = ly - j;
            for (int n = 0; n < 2; ++n) {
                p.rdata(FlockRealIdx::qx + n) = (1.0-fx)*(1.0-fy)*r(i,j,0,n) + fx*(1.0-fy)*r(i+1,j,0,n)
                                              + (1.0-fx)*fy*r(i,j+1,0,n) + fx*fy*r(i+1,j+1,0,n);
            }
        });
    }
}

void
FlockPC::CompareWithCIC (MultiFab const& q, Real& maxdiff, Real& rmsdiff, Real& rmsval) const
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto dx  = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    maxdiff = 0.0;
    Real s2 = 0.0, v2 = 0.0;
    for (ParConstIterType pti(*this, lev); pti.isValid(); ++pti) {
        const auto& aos = pti.GetArrayOfStructs();
        const ParticleType* pstruct = aos().data();
        const int np = aos.numParticles();
        auto const& r = q.const_array(pti);
        for (int ip = 0; ip < np; ++ip) {
            const ParticleType& p = pstruct[ip];
            const Real lx = (p.pos(0) - plo[0])/dx[0] - 0.5;
            const Real ly = (p.pos(1) - plo[1])/dx[1] - 0.5;
            const int i = static_cast<int>(std::floor(lx));
            const int j = static_cast<int>(std::floor(ly));
            const Real fx = lx - i;
            const Real fy = ly - j;
            for (int n = 0; n < 2; ++n) {
                const Real cic = (1.0-fx)*(1.0-fy)*r(i,j,0,n) + fx*(1.0-fy)*r(i+1,j,0,n)
                               + (1.0-fx)*fy*r(i,j+1,0,n) + fx*fy*r(i+1,j+1,0,n);
                const Real val = p.rdata(FlockRealIdx::qx + n);
                const Real d = cic - val;
                maxdiff = amrex::max(maxdiff, std::abs(d));
                s2 += d*d;
                v2 += val*val;
            }
        }
    }
    ParallelDescriptor::ReduceRealMax(maxdiff);
    ParallelDescriptor::ReduceRealSum(s2);
    ParallelDescriptor::ReduceRealSum(v2);
    const Real n = static_cast<Real>(TotalNumberOfParticles());
    rmsdiff = std::sqrt(s2/n);
    rmsval  = std::sqrt(v2/n);
}

Real
FlockPC::MaxTurningRate (FlockingParams const& fp) const
{
    using PType = FlockPC::ParticleType;
    Real q2 = amrex::ReduceMax(*this, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
    {
        return p.rdata(FlockRealIdx::qx)*p.rdata(FlockRealIdx::qx)
             + p.rdata(FlockRealIdx::qy)*p.rdata(FlockRealIdx::qy);
    });
    return fp.align_eps*std::sqrt(q2);
}

void
FlockPC::ComputeAlignPairs (FlockingParams const& fp)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const Real lx = geom.ProbLength(0);
    const Real ly = geom.ProbLength(1);
    const Real rc = flock_align_cutoff(fp);
    const Real rc2 = rc*rc;
    const Real rinv = 1.0/fp.align_R;
    const Real hmin = amrex::min(geom.CellSize(0), geom.CellSize(1));
    if (rc > m_neighbor_cells*hmin) {
        Abort("FlockPC::ComputeAlignPairs: the kernel cutoff exceeds the neighbour range");
    }
    const Real ninv = 1.0/static_cast<Real>(TotalNumberOfParticles());

    fillNeighbors();

    for (ParIterType pti(*this, lev); pti.isValid(); ++pti) {
        auto& ptile = ParticlesAt(lev, pti);
        auto& aos = ptile.GetArrayOfStructs();
        const int np = ptile.numRealParticles();
        ParticleType* pstruct = aos().data();

        auto check_pair = [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p1, const ParticleType& p2) noexcept
        {
            Real dxij = p2.pos(0) - p1.pos(0);
            Real dyij = p2.pos(1) - p1.pos(1);
            dxij -= lx*std::floor(dxij/lx + 0.5);
            dyij -= ly*std::floor(dyij/ly + 0.5);
            return (dxij*dxij + dyij*dyij < rc2);
        };
        amrex::NeighborList<ParticleType> nlist;
        Box bx = pti.tilebox();
        bx.grow(m_neighbor_cells);
        nlist.build(ptile, bx, geom, check_pair, m_neighbor_cells);
        auto ndata = nlist.data();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            // the particle itself, G(0) = 1
            Real sx = std::cos(p.rdata(FlockRealIdx::theta));
            Real sy = std::sin(p.rdata(FlockRealIdx::theta));
            for (auto const& q : ndata.getNeighbors(i)) {
                Real dxij = q.pos(0) - p.pos(0);
                Real dyij = q.pos(1) - p.pos(1);
                dxij -= lx*std::floor(dxij/lx + 0.5);
                dyij -= ly*std::floor(dyij/ly + 0.5);
                const Real g = flock_align_G(std::sqrt(dxij*dxij + dyij*dyij)*rinv);
                sx += g*std::cos(q.rdata(FlockRealIdx::theta));
                sy += g*std::sin(q.rdata(FlockRealIdx::theta));
            }
            p.rdata(FlockRealIdx::qx) = sx*ninv;
            p.rdata(FlockRealIdx::qy) = sy*ninv;
        });
    }

    clearNeighbors();
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
