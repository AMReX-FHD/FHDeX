#ifndef _IBParticleContainer_H_
#define _IBParticleContainer_H_


#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Particles.H>
#include <AMReX_Periodicity.H>
#include <AMReX_NeighborParticles.H>
#include <MFUtil.H>
#include <IBParticleInfo.H>
#include <IBCore.H>


using namespace amrex;
// using namespace std;

class mfix_level;
class IBCore;


//TODO: Check members lables with "Audit"

// IBP => Immmersed Boundary Particle: fairly simple immersed boundaries
// (analogous to MFIX Particles)

struct IBP_realData {
    //Analogous to particle realData (p.m_data)
    // TODO: orientation?
    enum {
        radius = 0,
        volume,    // .................................... TODO: Audit
        mass,
        density,   // .................................... TODO: Audit
        oneOverI,  // inverse moment of intertia
        velx,
        vely,
        velz,
        omegax,
        omegay,
        omegaz,
        dragx,    // ..................................... TODO: Audit
        dragy,    // ..................................... TODO: Audit
        dragz,    // ..................................... TODO: Audit
        count     // Awesome little trick!
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "radius",
            "volume",
            "mass",
            "density",
            "oneOverI",
            "velx",
            "vely",
            "velz",
            "omegax",
            "omegay",
            "omegaz",
            "dragx",
            "dragy",
            "dragz"
        };
    };
};



struct IBP_intData {
    //Analogous to particle intData (p.m_data)
    enum {
        phase = 0,  // ................................... TODO: Audit
        state,      // ................................... TODO: Audit
        count
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "phase",
            "state"
        };
    };
};



// TODO: final is used here because IBParticleContainer extends a template
// class into a non-template class: any further specialization might need
// different template parameters



class IBParIter final 
    : public amrex::ParIter<IBP_realData::count, IBP_intData::count, 0, 0>
{

public:

    using ParIter<IBP_realData::count, IBP_intData::count, 0, 0>::ParIter;
};



class IBParticleContainer final
    : public amrex::NeighborParticleContainer<IBP_realData::count, IBP_intData::count>
{

public:

    using NeighborParticleContainer<IBP_realData::count, IBP_intData::count>
        ::NeighborParticleContainer;

    // Get number of particles
    int NumberOfParticles(IBParIter & pti){ return pti.GetArrayOfStructs().size();};

    IBParticleContainer(amrex::AmrCore * amr_core, int n_nbhd);
    IBParticleContainer(const Geometry & geom, const DistributionMapping & dmap,
                        const BoxArray & ba, int n_nbhd);

    // TODO: Audit
    virtual ~IBParticleContainer() {};



    /****************************************************************************
    * Aux data types                                                            *
    ****************************************************************************/

    // Used by PICMultiDeposition when iterating over particles (TODO: Audit)
    //using MyConstParIter = ParConstIter<IBP_realData::count, IBP_intData::count, 0, 0>;

    // Data types representing slightly different representations of particle
    // data (for convenience and avoiding auto &).
    //using RealType      = typename ParticleType::RealType;
    //using ParticleLevel = std::map<std::pair<int, int>, ParticleTileType>;
    // using AoS = ArrayOfStructs<IBP_realData::count, IBP_intData::count>;

    // (ID, init CPU) tuple: unique to each particle
    using ParticleIndex = std::pair<int, int>;

    // Particle Position, Velocity, and Radius
    using PPVR = struct {
        RealVect pos;
        RealVect vel;
        double   rad;
    };


    /****************************************************************************
    * Initialization                                                            *
    ****************************************************************************/

    void InitList(int lev,
                  const Vector<RealVect> & pos,
                  const Vector<Real> & r,
                  const Vector<Real> & rho);



    /****************************************************************************
     * Immersed-Boundary Representation                                         *
     ***************************************************************************/

    const Vector<RealVect> & MarkerPositions(int lev, const ParticleIndex & pindex) const {
        // Do not use std::map::operator[] as this would discard the const qualifier
        return marker_positions[lev].at(pindex);
    };

    void FillMarkerPositions(int lev, int n_marker);

    void SpreadMarkers(int lev, const ParticleIndex & pindex, const Vector<RealVect> & f_in,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_out);

    void InterpolateMarkers(int lev, const ParticleIndex & pindex, Vector<RealVect> & f_out,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_in);



    /****************************************************************************
    * Physics                                                                   *
    ****************************************************************************/

    void InterpolateParticleForces(int lev,
            const std::array<MultiFab, AMREX_SPACEDIM> & force, const IBCore & ib_core,
            std::map<ParticleIndex, std::array<Real, AMREX_SPACEDIM>> & particle_forces
        );

    void MoveIBParticles(int lev, Real dt,
            const std::map<ParticleIndex, std::array<Real, AMREX_SPACEDIM>> & particle_forces
        );


    // TODO: do we still need this?
    void AllocData();
    // TODO: do we still need this?
    void AllocateArrays(int lev, int nghost);



    /****************************************************************************
    * Data Regrid / Data Transfer                                               *
    ****************************************************************************/

    //void RegridArrays(int lev, BoxArray& new_grids, DistributionMapping& new_dmap);
    // TODO: do we still need this?
    void CopyFluidData(int lev, const mfix_level & mf_lev);

    int get_nghost() const {return nghost;};



    /****************************************************************************
     * Output particle data                                                     *
     ***************************************************************************/

    void PrintParticleData(int lev);

    void LocalIBParticleInfo(Vector<IBP_info> & part_info, int lev, PairIndex index,
                             bool unique = false);
    Vector<IBP_info> LocalIBParticleInfo(int lev, PairIndex index, bool unique = false);

    void NeighborIBParticleInfo(Vector<IBP_info> & npart_info, int lev, PairIndex index,
                                bool unique = false);
    Vector<IBP_info> NeighborIBParticleInfo(int lev, PairIndex index, bool unique = false);

    void IBParticleInfo(Vector<IBP_info> & part_info, int lev, PairIndex index,
                        bool unique = false);
    Vector<IBP_info> IBParticleInfo(int lev, PairIndex index, bool unique = false);



private:


    ///
    /// We add a particle to the neighbor list if it's center is within 3 particle radii
    ///
    struct {
        bool operator()(const ParticleType& p1, const ParticleType& p2) const
        {
            return AMREX_D_TERM(   (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0)) ,
                                 + (p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)) ,
                                 + (p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)) )
                <= 9.0*p1.rdata(0)*p1.rdata(0);

        }
    } CheckPair;


    /***************************************************************************
     *                                                                         *
     * IB Marker and interpolation data                                        *
     *                                                                         *
     **************************************************************************/

    // Marker positions per particle and level
    Vector<std::map<ParticleIndex, Vector<RealVect>>> marker_positions;
    // Posotions on faces
    Vector<std::array<MultiFab, AMREX_SPACEDIM>>      face_coords;


    static bool use_neighbor_list;
    static bool sort_neighbor_list;

    void InitInternals(int ngrow);

    void ReadStaticParameters();


    // TODO: this might bto be used anymore:
    amrex::AmrCore * m_amr_core;
    int nghost;
};


#endif
