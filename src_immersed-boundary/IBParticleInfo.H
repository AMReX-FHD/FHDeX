#ifndef _IB_ParticleInfo_H_
#define _IB_ParticleInfo_H_

#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>

using namespace amrex;


/********************************************************************************
 *                                                                              *
 * PARTICLE INFO structures used to communicate particle data from              *
 * IBParticleContainer to various other utilities (such as IBCore).             *
 *                                                                              *
 * NOTE: these data structures need to match the fortran versions in:           *
 *       ib_fort_utils.f90 / module ib_fort_utils                               *
 *                                                                              *
 *******************************************************************************/




struct IBP_info {
    RealVect pos;     // Particle position
    RealVect vel;     // Particle velocity
    IntVect  index;   // Index (of pos) in the grid

    Real radius;      // Particle radius
    int  id;          // Particle ID (local initial MPI rank)
    int  cpu;         // MPI rank (where the particle was created)

    int  real;        // Is the particle real (not a neighbor)
                      // NOTE: not a bool because of c/fortran-interface

    std::pair<int, int> asPairIndex() const {
        std::pair<int, int> pindex(id, cpu);
        return pindex;
    }

    // only compare id and cpu
    bool operator== (const IBP_info & other) {
        return ( (this->id == other.id) && (this->cpu == other.cpu) );
    }
};



struct IBM_info {
    RealVect pos;     // Particle position
    RealVect vel;     // Particle velocity
    RealVect force;   // Marker forces

    int  id;          // Particle ID (local initial MPI rank)
    int  cpu;         // MPI rank (where the particle was created)

    int  real;        // Is the particle real (not a neighbor)
                      // NOTE: not a bool because of c/fortran-interface

    std::pair<int, int> asPairIndex() const {
        std::pair<int, int> pindex(id, cpu);
        return pindex;
    }

    // only compare id and cpu
    bool operator== (const IBM_info & other) {
        return ( (this->id == other.id) && (this->cpu == other.cpu) );
    }
};


#endif
