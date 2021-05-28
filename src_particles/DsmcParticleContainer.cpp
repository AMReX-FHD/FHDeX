#include "DsmcParticleContainer.H"

//#include "particle_functions_K.H"
#include "paramplane_functions_K.H"
#include <math.h>


FhdParticleContainer::FhdParticleContainer(const Geometry & geom,
                              const DistributionMapping & dmap,
                              const BoxArray            & ba,
                              int ncells)
    : NeighborParticleContainer<FHD_realData::count, FHD_intData::count> (geom, dmap, ba, ncells)
{
    BL_PROFILE_VAR("FhdParticleContainer()",FhdParticleContainer);

    realParticles = 0;
    simParticles = 0;

    for(int i=0;i<nspecies;i++) {
        
        properties[i].mass = mass[i];
        properties[i].radius = diameter[i]/2.0;
        properties[i].Neff = particle_neff; // From DSMC, this will be set to 1 for electolyte calcs
        properties[i].R = k_B/properties[i].mass; //used a lot in kinetic stats cals, bu not otherwise necessary for electrolytes
        properties[i].T = T_init[i];

        if (particle_count[i] >= 0) {

            properties[i].total = particle_count[i];
            properties[i].n0 = particle_neff*properties[i].total/domainVol;
            
            Print() << "Species " << i << " count " << properties[i].total << "\n";
            Print() << "Species " << i << " n0 " << properties[i].total << "\n";
        }
        else {
            properties[i].total = (int)amrex::Math::ceil(particle_n0[i]*domainVol/particle_neff);
            properties[i].n0 = properties[i].total/domainVol;

            Print() << "Species " << i << " count " << properties[i].total << "\n";
            Print() << "Species " << i << " n0 adjusted to " << properties[i].n0 << "\n";
        }

        Print() << "Species " << i << "\n";
        Print() << "Mass " << properties[i].mass << "\n";
        Print() << "Radius " << properties[i].radius << "\n";

        realParticles = realParticles + properties[i].total*particle_neff;
        simParticles = simParticles + properties[i].total;

    }

    totalCollisionCells = n_cells[0]*n_cells[1]*n_cells[2];
    domainVol = (prob_hi[0] - prob_lo[0])*(prob_hi[1] - prob_lo[1])*(prob_hi[2] - prob_lo[2]);

    Print() << "Total real particles: " << realParticles << "\n";
    Print() << "Total sim particles: " << simParticles << "\n";

    Print() << "Collision cells: " << totalCollisionCells << "\n";
    Print() << "Sim particles per cell: " << simParticles/totalCollisionCells << "\n";


}


