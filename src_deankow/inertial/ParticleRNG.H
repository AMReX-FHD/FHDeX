#ifndef INERTIAL_PARTICLE_RNG_H_
#define INERTIAL_PARTICLE_RNG_H_

#include <cmath>
#include <cstdint>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

// Counter-based random numbers for particles: a draw depends only on
// (seed, particle id, cpu, step, draw index), so a particle gets the same
// increments wherever (on whichever patch or rank) it is advanced, and runs do
// not depend on the domain decomposition.

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::uint64_t prng_mix (std::uint64_t x) noexcept
{
    // splitmix64 finalizer
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
std::uint64_t prng_key (std::uint64_t seed, std::int64_t id, int cpu, std::int64_t step,
                        std::uint64_t draw) noexcept
{
    std::uint64_t h = prng_mix(seed ^ 0x5851f42d4c957f2dULL);
    h = prng_mix(h ^ static_cast<std::uint64_t>(id));
    h = prng_mix(h ^ (static_cast<std::uint64_t>(static_cast<std::uint32_t>(cpu)) << 1));
    h = prng_mix(h ^ static_cast<std::uint64_t>(step));
    return prng_mix(h ^ draw);
}

// uniform in (0, 1]
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real prng_uniform (std::uint64_t k) noexcept
{
    return (static_cast<amrex::Real>(k >> 11) + 1.0) * (1.0/9007199254740992.0);
}

// standard normal number n (0, 1, ...) for this particle and step, by Box-Muller
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real prng_normal (std::uint64_t seed, std::int64_t id, int cpu, std::int64_t step,
                         int n) noexcept
{
    const amrex::Real u1 = prng_uniform(prng_key(seed, id, cpu, step, 2*static_cast<std::uint64_t>(n)));
    const amrex::Real u2 = prng_uniform(prng_key(seed, id, cpu, step, 2*static_cast<std::uint64_t>(n) + 1));
    return std::sqrt(-2.0*std::log(u1)) * std::cos(6.283185307179586476925*u2);
}

#endif
