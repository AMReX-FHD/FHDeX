#ifndef _DsmcParticleContainer_H_
#define _DsmcParticleContainer_H_


#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Particles.H>
#include <AMReX_Periodicity.H>
#include <AMReX_NeighborParticles.H>

#include <common_namespace.H>

#include <species.H>
#include "rng_functions.H"
#include "paramPlane.H"
#include "particle_functions_F.H"


using namespace amrex;
using namespace std;


// IBM => Immmersed Boundary Marker
struct DSMC_realData {
    //Analogous to particle realData (p.m_data)
    enum {
        radius = 0,
        velx,
        vely,
        velz,
        forcex,
        forcey,
        forcez,
        pred_posx,
        pred_posy,
        pred_posz,
        pred_velx,
        pred_vely,
        pred_velz,
        pred_forcex,
        pred_forcey,
        pred_forcez,
        vx,
        vy,
        vz,
        fx,
        fy,
        fz,
        ux,
        uy,
        uz,
        mass,
        R,
        q,
        accelFactor,
        dragFactor,
        ox,
        oy,
        oz,
        ax,
        ay,
        az,
        travelTime,
        diffAv,
        stepCount,
        multi,
        dryDiff,
        wetDiff,
        totalDiff,
        sigma,
        eepsilon,
        potential,
	p3m_radius,
        count    // Awesome little trick! (only works if first field is 0)
    };
};



struct DSMC_intData {
    //Analogous to particle intData (p.m_data)
    enum {
        sorted,
        i,
        j,
        k,
        species,
        count
    };
};

class DsmcParticleContainer
    : public amrex::NeighborParticleContainer<DSMC_realData::count, DSMC_intData::count>
{
public:

    using FhdParIter = ParIter<DSMC_realData::count, DSMC_intData::count>;

    DsmcParticleContainer(const Geometry            & geom, 
                                  const DistributionMapping & dmap,
                                  const BoxArray            & ba,
                                  int ncells);

    void InitParticles(species* particleInfo, const Real* dxp);

    void WriteParticlesAscii(string asciiName);

    void DoRFD(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF, const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           const MultiFab& cellCenters,
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const paramPlane* paramPlaneList, const int paramPlaneCount, int sw);

    void computeForcesNL(const MultiFab& charge, const MultiFab& coords, const Real* dx);

    void MoveParticlesDSMC(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount,Real time, int* flux);

    void MoveIons(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF, const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const MultiFab& mobility,
                                           const paramPlane* paramPlaneList, const int paramPlaneCount, int sw);


    void SpreadIons(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF, const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                                           const MultiFab& charge,
                                           const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                                           const MultiFab& cellCenters,
                                           std::array<MultiFab, AMREX_SPACEDIM>& source,
                                           std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                                           const paramPlane* paramPlaneList, const int paramPlaneCount, int sw);

    void SyncMembrane(double* spec3xPos, double* spec3yPos, double* spec3zPos, double* spec3xForce, double* spec3yForce, double* spec3zForce, int length, int step, const species* particleInfo);

    void collectFields(const Real dt, const Real* dxPotential, const MultiFab& RealCenterCoords, const Geometry geomF, MultiFab& charge, MultiFab& chargeTemp, MultiFab& mass, MultiFab& massTemp);

    void InitCollisionCells(              MultiFab& collisionPairs,
                                          MultiFab& collisionFactor, 
                                          MultiFab& cellVols, const species particleInfo, const Real delt);

    void CollideParticles(              MultiFab& collisionPairs,
                                          MultiFab& collisionFactor, 
                                          MultiFab& cellVols, const species particleInfo, const Real delt);


    void InitializeFields(    MultiFab& particleInstant,
                              MultiFab& cellVols, const species particleInfo);

    void EvaluateStats(
                              MultiFab& particleInstant,
                              MultiFab& particleMeans,
                              MultiFab& particleVars,
                              MultiFab& cellVols, species particleInfo, const Real delt, int steps);

  

    int numWrongCell();


    void ReBin();

    static bool use_neighbor_list;
    static bool sort_neighbor_list;



protected:

        void UpdateFortranStructures();

        void UpdateCellVectors();

        // used to store vectors of particle indices on a cell-by-cell basis
//        std::map<int, amrex::BaseFab<std::vector<int> > > m_cell_vectors;
        std::map<int, std::vector<std::vector<int> > > m_cell_vectors;

        // primitive data for passing the above into Fortran
        std::map<int, amrex::BaseFab<int> > m_vector_size;
        std::map<int, amrex::BaseFab<int*> > m_vector_ptrs;

        // used to keep track of whether we need to rebuild the above or not
        bool m_vectors_initialized = false;
        amrex::BoxArray::RefID m_BARef;
        amrex::DistributionMapping::RefID m_DMRef;


private:

//    void correctCellVectors(int old_index, int new_index, int grid, const ParticleType& p);

    int n_list;

    struct CHECK_PAIR {
        bool operator()(const ParticleType & p1, const ParticleType & p2) const {
            return true;
        }
    };
};


#endif
