#ifndef _DSMCParticleContainer_H_
#define _DSMCParticleContainer_H_

#include "paramPlane.H"
#include <AMReX_NeighborParticles.H>
#include <AMReX_MultiFab.H>
#include <common_namespace.H>

using namespace std;
using namespace common;

struct FHD_realData {
    enum {
        velx = 0,
        vely,
        velz,
        boostx,
        boosty,
        boostz,
        R,
        timeFrac,
        omega,
        lambda,
        radius,
        mass,
        travelTime,
        count
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "velx",
            "vely",
            "velz",
            "boostx",
            "boosty",
            "boostz",
            "R",
            "timeFrac",
            "omega",
            "lambda",
            "radius",
            "mass",
            "travelTime"
        };
    };
};



struct FHD_intData {
    enum {
        sorted,
        i,
        j,
        k,
        species,
        newSpecies,
        fluxRec,
        count
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "sorted",
            "i",
            "j",
            "k",
            "species",
            "newSpecies",
            "fluxRec"
        };
    };
};

typedef struct {
    int type;
    long total;
    double mass;
    double radius;
    double n0;
    double R;
    double Neff;
    double partVol;
    double part2cellVol;
} dsmcSpecies;

typedef struct {
    double alpha;
    double csx;
} dsmcInterSpecies;

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
IntVect getPartCell (P const& p,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                         const Box& domain) noexcept
{
    IntVect iv = IntVect(
        AMREX_D_DECL(int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                     int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                     int(amrex::Math::floor((p.pos(2)-plo[2])*dxi[2]))));
    //iv -= domain.smallEnd();
    return iv;
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
IntVect getPartCell2 (P const& p,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                         const Box& domain) noexcept
{
    IntVect iv = IntVect(
        AMREX_D_DECL(int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                     int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                     int(amrex::Math::floor((p.pos(2)-plo[2])*dxi[2]))));
    cout << "pos: " << p.pos(0) << ", plo: " << plo[0] << ", dxi: " << dxi[0] << endl;
    return iv;
}

/*template <typename P>*/
/*AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE*/
/*void getCellList (int* cellList, int* permArray, int* offArray, IntVect* cell, int spec) noexcept*/
/*{*/
/*    */
/*}*/

struct getBin
{
    GpuArray<Real,AMREX_SPACEDIM> plo;
    GpuArray<Real,AMREX_SPACEDIM> dxi;
    Box domain;
    int ncells;

    template <typename ParticleType>
    AMREX_GPU_HOST_DEVICE
    unsigned int operator() (const ParticleType& p) const noexcept
    {
        //Box tbx;
        auto iv = getPartCell(p, plo, dxi, domain);
        //auto tid = getTileIndex(iv, box, true, bin_size, tbx);
        long imap = domain.index(iv);
        int spec = p.idata(FHD_intData::species);
        int map = imap;
        if(spec > 0)
        {
            for(int i=0;i<spec;i++)
            {
                map += ncells;
            }
        }
        //Print() << "map: " << map << ", spec: " << spec << endl;
        return static_cast<unsigned int>(map);
    }
};


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int mapBin(IntVect iv, int spec, Box domain)
{
    auto cellLo = domain.smallEnd();
    auto cellHi = domain.bigEnd();
    int ncells = (cellHi[0]-cellLo[0]+1)*(cellHi[1]-cellLo[1]+1)*(cellHi[2]-cellLo[2]+1);

    long imap = domain.index(iv);
    int map = imap;
    if(spec > 0)
    {
        for(int i=0;i<spec;i++)
        {
            map += ncells;
        }
    }
    return map;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int getBinSize(unsigned int* offArray, IntVect iv, int spec, Box domain)
{
    int bin = mapBin(iv, spec, domain);
    return offArray[bin+1]-offArray[bin];
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
unsigned int* getCellList(unsigned int* permArray, unsigned int* offArray, IntVect iv, int spec, Box domain)
{
    int bin = mapBin(iv, spec, domain);
    return &permArray[offArray[bin]];
}



class FhdParticleContainer
    : public amrex::NeighborParticleContainer<FHD_realData::count, FHD_intData::count>
{
public:

    using FhdParIter = ParIter<FHD_realData::count, FHD_intData::count>;

    FhdParticleContainer(const Geometry & geom,
        const DistributionMapping & dmap,
        const BoxArray & ba,
        int ncells);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// DSMC Routines
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void InitParticles(Real & dt);
    void ReInitParticles();
    void InitCollisionCells();

    void Source(const Real dt, paramPlane* paramPlaneList, const int paramPlaneCount, MultiFab& mfcuInst);
    void SourcePhonons(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Main DSMC routines
    void SortParticles();
    void SortParticlesDB();

    void CalcSelections(Real dt);
    void CollideParticles(Real dt);
    void CollideParticles2(Real dt);

    void MoveParticlesCPP(const Real dt, paramPlane* paramPlaneList, const int paramPlaneCount);
    void MovePhononsCPP(const Real dt, paramPlane* paramPlaneList, const int paramPlaneCount, const int step, const int istep, iMultiFab& bCell);
    void externalForce(const Real dt);
    void updateTimeStep(const Geometry& geom, Real& dt);
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MultiFab mfselect;
    MultiFab mfvrmax;
    MultiFab mfphi;

    MultiFab mfCollisions;
    int expectedCollisions[MAX_SPECIES], countedCollisions[MAX_SPECIES];
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Outputs
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void EvaluateStats(MultiFab& cuInst,
        MultiFab& cuMeans,
        MultiFab& cuVars,
        MultiFab& primInst,
        MultiFab& primMeans,
        MultiFab& primVars,
        MultiFab& mfcvlInst,
        MultiFab& mfcvlMeans,
        MultiFab& mfQMeans,
        MultiFab& coVars,
        MultiFab& spatialCross1D,
        const int steps,
        Real time);

   void EvaluateStatsPhonon(MultiFab& cuInst,
        MultiFab& cuMeans,
        MultiFab& cuVars,
        const int steps,
        Real time);

    void OutputParticles();

    void zeroCells();


  //Real Tg0[MAX_SPECIES];
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Scalars
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    GpuArray<dsmcSpecies, MAX_SPECIES> properties;
    GpuArray<dsmcInterSpecies, MAX_SPECIES*MAX_SPECIES> interproperties;

    long realParticles;
    long simParticles;

    long totalCollisionCells; // might only be helpful for statistics
    Real collisionCellVol; // assumed constant for all collision cells
    Real ocollisionCellVol;
    Real domainVol;
    Real domSize[3];
    //Real dx[3];

    int getSpeciesIndex(int species1, int species2);

    // Evaluating radial distribution function
    Real g0_Ma_Ahmadi(int species1, int species2, Real phi1, Real phi2);
    const Real chi_max = 500.0; // value is arbitrary (change as needed or maybe as input?)
    const Real phi_max = 0.643;

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constants
    Real pi_usr = 4.0*atan(1.0);

protected:
    // used to store vectors of particle indices on a cell-by-cell basis
    std::map<int, std::vector<vector<int> > > m_cell_vectors[MAX_SPECIES];
    //std::vector<std::map<int, std::vector<Gpu::ManagedVector<int> > >> m_cell_vectors;
    //std::vector<std::vector<Gpu::ManagedVector<int>>>  m_cell_vectors[MAX_SPECIES];

    DenseBins<ParticleType> m_bins;
};


#endif