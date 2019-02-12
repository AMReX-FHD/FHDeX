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
    RealVect pos;
    RealVect vel;
    IntVect  index;

    Real radius;
    int  id;
    int  cpu;
};



#endif
