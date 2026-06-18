#ifndef RESERVOIRSTAG_K_H_
#define RESERVOIRSTAG_K_H_

#include "compressible_functions_stag.H"
#include "compressible_functions.H"
#include "common_functions.H"
#include "rng_functions.H"
#include <math.h>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
random_cross_vel(const amrex::Real& a, const amrex::Real& vT, bool tangential, amrex::Real Vtangential,
                 amrex::RandomEngine const& engine)
{
    constexpr amrex::Real sqrtPI = 1.77245385091;
    constexpr amrex::Real small_number = 1.0e-13;

    amrex::Real v;
    if (!tangential) { // normal velocity crossing
        amrex::Real z;
        if (a > small_number) { // a > 0
            bool stop = false;
            while (!stop) {
                amrex::Real u = amrex::Random(engine);
                if ((a*sqrtPI)/(a*sqrtPI + 1.0 + a*a) > u) {
                    z = -1.0*std::abs(amrex::RandomNormal(0.0,sqrt(0.5),engine));
                    stop = true;
                }
                else if ((a*sqrtPI + 1.0)/(a*sqrtPI + 1.0 + a*a) > u) {
                    z = -1.0*sqrt(-1.0*std::log(amrex::Random(engine)));
                    stop = true;
                }
                else {
                    z = a*(1.0-sqrt(amrex::Random(engine)));
                    if (exp(1.0*z*z) > amrex::Random(engine)) stop = true;
                }
            }
            v = (a-z)*vT;
        }
        else if (a < -1.0*small_number) { // a < 0
            bool stop = false;
            while (!stop) {
                z = -1.0*sqrt(a*a - std::log(amrex::Random(engine)));
                if (a > z*(1.0 - amrex::Random(engine))) stop = true;
            }
            v = (a-z)*vT;
        }
        else { // a = 0
            v = vT*sqrt(-1.0*std::log(amrex::Random(engine)));
        }
    }
    else { // tangential velocity crossing
        v = Vtangential + vT*amrex::RandomNormal(0.0,sqrt(0.5),engine);
    }
    return v;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void
poisson_process_reservoir(const amrex::GpuArray<amrex::Real,MAX_SPECIES>& mass,
                          const amrex::GpuArray<amrex::Real,MAX_SPECIES>& rhoYk,
                          const amrex::Real& T, // temperature
                          const amrex::Real& Vnormal, // normal velocity
                          const int& nspec,
                          const amrex::Real& area,
                          const amrex::Real& kB,
                          const amrex::Real& dt,
                          amrex::Real& mass_cross,
                          amrex::GpuArray<amrex::Real,3>& mom_cross,
                          amrex::Real& en_cross,
                          amrex::GpuArray<amrex::Real,MAX_SPECIES>& spec_mass_cross,
                          int dim,         // dimension
                          amrex::Real VT1, // tangential velocity 1
                          amrex::Real VT2, // tangential velocity 2
                          amrex::RandomEngine const& engine)

{
    constexpr amrex::Real small_number = 1.0e-13;
    constexpr amrex::Real sqrtPI = 1.77245385091;

    mass_cross = 0.0;
    en_cross   = 0.0;
    for (int n=0;n<3;++n) {
        mom_cross[n] = 0.0;
    }
    for (int n=0;n<nspecies;++n) {
        spec_mass_cross[n] = 0.0;
    }

    amrex::GpuArray<amrex::Real,MAX_SPECIES> particle_rate;
    amrex::GpuArray<amrex::Real,MAX_SPECIES> timer;
    amrex::GpuArray<amrex::Real,MAX_SPECIES> thermal_speed;
    amrex::GpuArray<amrex::Real,MAX_SPECIES> speed_ratio;
    amrex::GpuArray<bool,MAX_SPECIES> stop_flag;

    for (int n=0;n<nspec;++n) {
        thermal_speed[n] = sqrt(2.0*kB*T/mass[n]);
        speed_ratio[n] = Vnormal/thermal_speed[n];
        amrex::Real num_density = rhoYk[n]/mass[n];
        if (std::abs(speed_ratio[n]) > small_number) {
            particle_rate[n] = (1.0/(2.0*sqrtPI))*num_density*area*thermal_speed[n]*
                (exp(-1.0*speed_ratio[n]*speed_ratio[n]) + speed_ratio[n]*sqrtPI*(1.0 + erf(speed_ratio[n])));
            if (particle_rate[n] < 0.0) {
                amrex::Abort("negative particle rate encountered");
            }
        }
        else {
            particle_rate[n] = (1.0/(2.0*sqrtPI))*num_density*area*thermal_speed[n];
            if (particle_rate[n] < 0.0) {
                amrex::Abort("negative particle rate encountered");
            }
        }
        timer[n] = 0.0;
        stop_flag[n] = false;
    }

    // Run Poisson Process
    bool stop = false;
    while(!stop) {

        // add a particle for each species (if stop[n] = false)
        amrex::Real vel_rand;
        for (int n=0;n<nspec;++n) {
            if (!stop_flag[n]) {

                // normal component
                mass_cross += mass[n]; // add a particle crossing for total mass
                spec_mass_cross[n] += mass[n]; // add a particle crossing for species mass
                vel_rand = random_cross_vel(speed_ratio[n],thermal_speed[n],false,0.0,engine); // get random velocity
                mom_cross[0] += mass[n]*vel_rand; // add a particle crossing for total momentum
                en_cross  += 0.5*mass[n]*(vel_rand)*(vel_rand); // add a particle crossing for total energy

                // tangential component 1
                if (dim >= 2) {
                    vel_rand = random_cross_vel(speed_ratio[n],thermal_speed[n],true,VT1,engine); // get random velocity
                    mom_cross[1] += mass[n]*vel_rand; // add a particle crossing for total momentum
                    en_cross  += 0.5*mass[n]*(vel_rand)*(vel_rand); // add a particle crossing for total energy
                }

                // tangential component 2
                if (dim == 3) {
                    vel_rand = random_cross_vel(speed_ratio[n],thermal_speed[n],true,VT2,engine); // get random velocity
                    mom_cross[2] += mass[n]*vel_rand; // add a particle crossing for total momentum
                    en_cross  += 0.5*mass[n]*(vel_rand)*(vel_rand); // add a particle crossing for total energy
                }
            }
        }

        stop = true;
        for (int n=0;n<nspec;++n) { // update dt for each species
            if (!stop_flag[n]) {
                amrex::Real rn   = amrex::Random(engine);
                timer[n]         = timer[n] - std::log(rn)/(particle_rate[n]);
                if (timer[n] >= dt) stop_flag[n] = true;
            }
            stop = stop and stop_flag[n];
        }
    }
}

#endif
