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
};

struct IBP_intData {
    //Analogous to particle intData (p.m_data)
    enum {
        phase = 0,  // ................................... TODO: Audit
        state,      // ................................... TODO: Audit
        count
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

    // Used by PICMultiDeposition when iterating over particles (TODO: Audit)
    //using MyConstParIter = ParConstIter<IBP_realData::count, IBP_intData::count, 0, 0>;

    // Data types representing slightly different representations of particle
    // data (for convenience and avoiding auto &).
    //using RealType      = typename ParticleType::RealType;
    //using ParticleLevel = std::map<std::pair<int, int>, ParticleTileType>;
    // using AoS = ArrayOfStructs<IBP_realData::count, IBP_intData::count>;

    // Get number of particles
    int NumberOfParticles(IBParIter & pti){ return pti.GetArrayOfStructs().size(); };

    IBParticleContainer(amrex::AmrCore * amr_core, int n_nbhd);
    IBParticleContainer(const Geometry & geom, const DistributionMapping & dmap,
                        const BoxArray & ba, int n_nbhd);

    // TODO: Audit
    virtual ~IBParticleContainer() {};


    /****************************************************************************
    * Initialization                                                            *
    ****************************************************************************/

    void InitList(int lev,
                  const Vector<RealVect> & pos,
                  const Vector<Real> & r,
                  const Vector<Real> & rho);


    /****************************************************************************
    * Physics                                                                   *
    ****************************************************************************/

    // (ID, init CPU) tuple: unique to each particle
    using ParticleIndex = std::pair<int, int>;

    void InterpolateParticleForces(
            const std::array<MultiFab, AMREX_SPACEDIM> & force, const IBCore & ib_core, int lev,
            std::map<ParticleIndex, std::array<Real, AMREX_SPACEDIM>> particle_forces
        );

    void MoveIBParticles(int lev, Real dt,
            std::map<ParticleIndex, std::array<Real, AMREX_SPACEDIM>> particle_forces
        );


    // TODO: do we still need this?
    void AllocData();
    // TODO: do we still need this?
    void AllocateArrays(int lev, int nghost);


    /****************************************************************************
    * Data Regrid / Data Transfer                                               *
    ****************************************************************************/

    //void RegridArrays(int lev, BoxArray& new_grids, DistributionMapping& new_dmap);
    void CopyFluidData(int lev, const mfix_level & mf_lev);

    int get_nghost() const {return nghost;};


    /****************************************************************************
     * Output particle data                                                     *
     ***************************************************************************/

    void PrintParticleData(int lev);

    void LocalIBParticleInfo(Vector<IBP_info> & part_info, int lev, PairIndex index);
    Vector<IBP_info> LocalIBParticleInfo(int lev, PairIndex index);

    void NeighborIBParticleInfo(Vector<IBP_info> & npart_info, int lev, PairIndex index);
    Vector<IBP_info> NeighborIBParticleInfo(int lev, PairIndex index);

    void IBParticleInfo(Vector<IBP_info> & part_info, int lev, PairIndex index);
    Vector<IBP_info> IBParticleInfo(int lev, PairIndex index);



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


    static bool use_neighbor_list;
    static bool sort_neighbor_list;

    void InitInternals();

    void ReadStaticParameters();


    // TODO: this might bto be used anymore:
    amrex::AmrCore * m_amr_core;
    int nghost;
};


#endif
