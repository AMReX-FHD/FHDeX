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
#include <IBMarkerContainer.H>
#include <IBCore.H>


using namespace amrex;

class mfix_level;
class IBCore;



//TODO: Check members lables with "Audit"


// IBP => Immmersed Boundary Particle
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
    int NumberOfParticles(IBParIter & pti){ return pti.GetArrayOfStructs().numParticles();};

    IBParticleContainer(AmrCore * amr_core, int n_nbhd);
    IBParticleContainer(const Geometry & geom, const DistributionMapping & dmap,
                        const BoxArray & ba, int n_nbhd);

    // TODO: Audit
    virtual ~IBParticleContainer() {};



    /****************************************************************************
    * Auxiliary data types                                                      *
    ****************************************************************************/

    // Data types representing slightly different representations of particle
    // data (for convenience and avoiding auto &).
    // using RealType      = typename ParticleType::RealType;
    // using ParticleLevel = std::map<std::pair<int, int>, ParticleTileType>;
    // using AoS = ArrayOfStructs<IBP_realData::count, IBP_intData::count>;

    // (ID, initial CPU) tuple: unique to each particle
    using ParticleIndex = std::pair<int, int>;
    using IBMarkerMap   = std::map<ParticleIndex, Vector<RealVect>>;

    // Particle Position, Velocity, Orientation, and Radius (abbriated particle info struct)
    using PPVOR = struct {
        RealVect pos;
        RealVect vel;
        RealVect ori;
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

    const Vector<RealVect> & MarkerVelocities(int lev, const ParticleIndex & pindex) const {
        // Do not use std::map::operator[] as this would discard the const qualifier
        return marker_velocities[lev].at(pindex);
    };

    const Vector<RealVect> & MarkerForces(int lev, const ParticleIndex & pindex) const {
        // Do not use std::map::operator[] as this would discard the const qualifier
        return marker_forces[lev].at(pindex);
    };


    void FillMarkerPositions(int lev, int n_marker);

    void SpreadMarkers(int lev, const ParticleIndex & pindex,
                       const Vector<RealVect> & f_in,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;

    void SpreadMarkers(int lev, const ParticleIndex & pindex,
                       const Vector<RealVect> & f_in,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_out,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const;

    void SpreadKernel(const Box& bx, std::array<MultiFab, AMREX_SPACEDIM> & f_out,
                      std::array<MultiFab, AMREX_SPACEDIM> & f_weights,
                      const Vector<std::array<MultiFab, AMREX_SPACEDIM>> & face_coords,
                      const Real* pos, const Real* v_spread,
                      const Real* dx, int* nghost, MFIter* mfi) const;

    void InvInterpolateMarkers(int lev, const ParticleIndex & pindex,
                               const Vector<RealVect> & f_in,
                               std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;

    void InterpolateMarkers(int lev, const ParticleIndex & pindex,
                            Vector<RealVect> & f_out,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const;

    void InterpolateMarkers(int lev, const ParticleIndex & pindex,
                            Vector<RealVect> & f_out,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_in) const;

    void InvSpreadMarkers(int lev, const ParticleIndex & pindex,
                          Vector<RealVect> & f_out,
                          const std::array<MultiFab, AMREX_SPACEDIM> & f_in) const;


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


    /****************************************************************************
     * Output particle data                                                     *
     ***************************************************************************/

    int get_nghost() const {return nghost;};

    void PrintParticleData(int lev);

    void LocalIBParticleInfo(Vector<IBP_info> & part_info, int lev, PairIndex index,
                             bool unique = false) const;
    Vector<IBP_info> LocalIBParticleInfo(int lev, PairIndex index,
                                         bool unique = false) const;
    Vector<IBP_info> LocalIBParticleInfo(int lev, bool unique = false) const;


    void NeighborIBParticleInfo(Vector<IBP_info> & npart_info, int lev, PairIndex index,
                                bool unique = false) const ;
    Vector<IBP_info> NeighborIBParticleInfo(int lev, PairIndex index,
                                            bool unique = false) const;
    Vector<IBP_info> NeighborIBParticleInfo(int lev, bool unique = false) const;


    void IBParticleInfo(Vector<IBP_info> & part_info, int lev, PairIndex index,
                        bool unique = false) const;
    Vector<IBP_info> IBParticleInfo(int lev, PairIndex index,
                                    bool unique = false) const;
    Vector<IBP_info> IBParticleInfo(int lev, bool unique = false) const;

    const Vector<std::array<MultiFab, AMREX_SPACEDIM>> & get_face_coords() {
        return  face_coords;
    }

 //   const Vector<RealVect> get_pcenter_coords() {
 //       return  pos(0);

private:


    ///
    /// We add a particle to the neighbor list if it's center is within 3 particle radii
    ///
    struct {
        bool operator() (const ParticleType & p1, const ParticleType & p2) const {
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
    Vector<IBMarkerMap> marker_positions;
    Vector<IBMarkerMap> marker_velocities;
    Vector<IBMarkerMap> marker_forces;

    // Positions on faces
    Vector<std::array<MultiFab, AMREX_SPACEDIM>> face_coords;


    static bool use_neighbor_list;
    static bool sort_neighbor_list;

    void InitInternals(int ngrow);

    void ReadStaticParameters();

    // TODO: this might bto be used anymore:
    AmrCore * m_amr_core;
    int nghost;
};


#endif
