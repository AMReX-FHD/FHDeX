#ifndef _FhdParticleContainer_H_
#define _FhdParticleContainer_H_

#include "IBParticleInfo.H"
#include "IBMarkerContainerBase.H"
#include "species.H"

#include "paramPlane.H"
//#include "paramplane_functions_F.H"

//#include "particle_functions_F.H"
#include "particle_functions.H"

using namespace std;

// IBM => Immmersed Boundary Marker
struct FHD_realData {
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
        spring,
        omega,
        lambda,
        count    // Awesome little trick! (only works if first field is 0)
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "radius",
            "velx",
            "vely",
            "velz",
            "forcex",
            "forcey",
            "forcez",
            "pred_posx",
            "pred_posy",
            "pred_posz",
            "pred_velx",
            "pred_vely",
            "pred_velz",
            "pred_forcex",
            "pred_forcey",
            "pred_forcez",
            "vx",
            "vy",
            "vz",
            "fx",
            "fy",
            "fz",
            "ux",
            "uy",
            "uz",
            "mass",
            "R",
            "q",
            "accelFactor",
            "dragFactor",
            "ox",
            "oy",
            "oz",
            "ax",
            "ay",
            "az",
            "travelTime",
            "diffAv",
            "stepCount",
            "multi",
            "dryDiff",
            "wetDiff",
            "totalDiff",
            "sigma",
            "eepsilon",
            "potential",
            "p3m_radius",
            "spring",
            "omega",
            "lambda"
        };
    };
};

struct FHD_intData {
    //Analogous to particle intData (p.m_data)
    enum {
        sorted,
        i,
        j,
        k,
        species,
        visible,
        pinned,
        newSpecies,
        count
    };

    static Vector<std::string> names() {
        return Vector<std::string> {
            "sorted",
            "i",
            "j",
            "k",
            "species",
            "visible",
            "pinned",
            "new_species"
        };
    };
};

struct FHD_realDataLW {
    //Analogous to particle realData (p.m_data)
    enum {
        velx,
        vely,
        velz,
        forcex,
        forcey,
        forcez,
        posx,
        posy,
        posz,
        count    // Awesome little trick! (only works if first field is 0)
    };
};

typedef struct {
    Real x;
    Real y;
    Real z;
} Triplet;

class FhdParIter
    : public IBMarIterBase<FHD_realData::count, FHD_intData::count>
{
public:
    using IBMarIterBase<FHD_realData::count, FHD_intData::count>::IBMarIterBase;
};

class FhdParticleContainer
    : public IBMarkerContainerBase<FHD_realData, FHD_intData>
{
public:
    using IBMarkerContainerBase<FHD_realData, FHD_intData>
        ::IBMarkerContainerBase;

    using MyConstIBMarIter = IBMarkerContainerBase<FHD_realData, FHD_intData>
        ::MyConstIBMarIter;

    FhdParticleContainer(AmrCore * amr_core, int n_nbhd);
    FhdParticleContainer(const Geometry & geom, const Geometry & geomF, const DistributionMapping & dmap,
                         const BoxArray & ba, const BoxArray & bc, int n_nbhd, int ngF);

    virtual ~FhdParticleContainer() {};

    // (ID, initial CPU) tuple: unique to each particle
    using MarkerIndex = typename IBMarkerContainerBase<FHD_realData, FHD_intData>
        ::PairIndex;

    void InitParticles(species* particleInfo, const Real* dxp);
    void InitParticlesFromFile(species* particleInfo, const Real* dxp);

    void ReInitParticles();//species* particleInfo, const Real* dxp, Real * posX, Real * posY, Real * posZ, Real * charge, Real * sigma, Real * epsilon, int * speciesV, Real * diffdry, Real * diffwet, Real * difftotal);

    void PostRestart();
    void WriteParticlesAscii(string asciiName);
    void forceFunction(Real dt);
    void velNorm();
    void pinForce();
    void potentialDistribution(long totalParticles, const int step, const species* particleInfo);
    void fillMobilityMatrix(int id, int comp);
    void clearMobilityMatrix();
    void writeMat();

    void DoRFD(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF,
               const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
               const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
               const MultiFab& cellCenters,
               std::array<MultiFab, AMREX_SPACEDIM>& source,
               std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
               const paramPlane* paramPlaneList, const int paramPlaneCount, int sw);

    void PrintParticles();
    void SetPosition(int id, Real x, Real y, Real z);
    void SetVel(int id, Real x, Real y, Real z);
    void SetForce(int id, Real x, Real y, Real z);
    void writeVel(int id);
    void MeanSqrCalc(int lev, int step);
    void BuildCorrectionTable(const Real* dx, int setMeasureFinal);

    void DoRFDbase(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF,
                   const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                   const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                   const MultiFab& cellCenters,
                   std::array<MultiFab, AMREX_SPACEDIM>& source,
                   std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                   const paramPlane* paramPlaneList, const int paramPlaneCount, int sw);

    void computeForcesNL(const MultiFab& charge, const MultiFab& coords, const Real* dx);
    void computeForcesNLGPU(const MultiFab& charge, const MultiFab& coords, const Real* dx);
    void computeForcesCoulombGPU(long totalParticles);
    void MoveParticlesDSMC(const Real dt, const paramPlane* paramPlaneList, const int paramPlaneCount, Real time, int* flux);
    void MoveIonsCPP(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF,
                     const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                     const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                     std::array<MultiFab, AMREX_SPACEDIM>& source,
                     std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                     paramPlane* paramPlaneList, const int paramPlaneCount, int sw);

    void SpreadIons(const Real dt, const Real* dxFluid, const Real* dxE, const Geometry geomF,
                    const std::array<MultiFab, AMREX_SPACEDIM>& umac, const std::array<MultiFab, AMREX_SPACEDIM>& efield,
                    const MultiFab& charge,
                    const std::array<MultiFab, AMREX_SPACEDIM>& RealFaceCoords,
                    const MultiFab& cellCenters,
                    std::array<MultiFab, AMREX_SPACEDIM>& source,
                    std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp,
                    const paramPlane* paramPlaneList, const int paramPlaneCount, int sw);

    void SpreadIonsGPU(const Real* dxFluid, const Real* dxE, const Geometry geomF,
                       const std::array<MultiFab, AMREX_SPACEDIM>& umac,
                       const std::array<MultiFab, AMREX_SPACEDIM>& coords,
                       std::array<MultiFab, AMREX_SPACEDIM>& efield,
                       std::array<MultiFab, AMREX_SPACEDIM>& source,
                       std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp);

    void SpreadIonsGPU(const Real* dxFluid, const Geometry geomF,
                       const std::array<MultiFab, AMREX_SPACEDIM>& umac,
                       const std::array<MultiFab, AMREX_SPACEDIM>& coords,
                       std::array<MultiFab, AMREX_SPACEDIM>& source,
                       std::array<MultiFab, AMREX_SPACEDIM>& sourceTemp);

    //void SyncMembrane(double* spec3xPos, double* spec3yPos, double* spec3zPos, double* spec3xForce, double* spec3yForce, double* spec3zForce, int length, int step, const species* particleInfo);

    void RadialDistribution(long totalParticles, const int step, const species* particleInfo);
    void CartesianDistribution(long totalParticles, const int step, const species* particleInfo);
    void collectFields(const Real dt, const Real* dxPotential, const MultiFab& RealCenterCoords, const Geometry geomF, MultiFab& charge, MultiFab& chargeTemp, MultiFab& mass, MultiFab& massTemp);
    void collectFieldsGPU(const Real dt, const Real* dxPotential, const MultiFab& RealCenterCoords, const Geometry geomP, MultiFab& charge, MultiFab& chargeTemp, MultiFab& mass, MultiFab& massTemp);
    void InitCollisionCells(MultiFab& collisionPairs,
                            MultiFab& collisionFactor,
                            MultiFab& cellVols, const species particleInfo, const Real delt);
    void CollideParticles(MultiFab& collisionPairs,
                          MultiFab& collisionFactor,
                          MultiFab& cellVols, const species particleInfo, const Real delt);
    void InitializeFields(MultiFab& particleInstant,
                          MultiFab& cellVols, const species particleInfo);
    void EvaluateStats(MultiFab& particleInstant,
                       MultiFab& particleMeans, species particleInfo, const Real delt, int steps);
    void invertMatrix();
    void TwoParticleCorrelation();
    void GetAllParticlePositions(Real* posx, Real* posy, Real* posz, int totalParticles);

    /****************************************************************************
     *                                                                          *
     * Access marker (amrex particle) data                                      *
     *                                                                          *
     ***************************************************************************/

    const IntVector & get_neighbor_list(int lev, const PairIndex & index) const {
        return neighbor_list[lev].at(index);
    };

    ParticleType * get_neighbors(int lev, const PairIndex & index) const {
        return (ParticleType *) neighbors[lev].at(index).GetArrayOfStructs().dataPtr();;
    };

    int numWrongCell();
    void ReBin();

    static bool use_neighbor_list;
    static bool sort_neighbor_list;

    int totalBins;
    int doRedist;
    Real *nearestN;

    Real *meanRadialDistribution   ;
    Real *meanRadialDistribution_pp;
    Real *meanRadialDistribution_pm;
    Real *meanRadialDistribution_mm;

    Real *binVolRadial;
    int radialStatsCount;

    Real *meanXDistribution   ;
    Real *meanXDistribution_pp;
    Real *meanXDistribution_pm;
    Real *meanXDistribution_mm;

    Real *meanYDistribution   ;
    Real *meanYDistribution_pp;
    Real *meanYDistribution_pm;
    Real *meanYDistribution_mm;

    Real *meanZDistribution   ;
    Real *meanZDistribution_pp;
    Real *meanZDistribution_pm;
    Real *meanZDistribution_mm;

    Real binVolCartesian;
    int cartesianStatsCount;

    int  threepmBins = 50;
    int  threepmCurrentBin = 1;
    int  threepmCurrentSample = 1;
    int  threepmSamples = 1000;
    Real threepmVals[50];
    Real threepmMax[50];
    Real threepmMin[50];
    Real threepmPoints[50];
    Real threepmRange = 5;

    int bottomListLength;
    int topListLength;

    Triplet* bottomList;
    Triplet* topList;

    //protected:

    // used to store vectors of particle indices on a cell-by-cell basis
    //    std::map<int, amrex::BaseFab<std::vector<int> > > m_cell_vectors;
    std::map<int, std::vector<std::vector<int> > > m_cell_vectors;

    // primitive data for passing the above into Fortran
    std::map<int, amrex::BaseFab<int> > m_vector_size;
    std::map<int, amrex::BaseFab<int*> > m_vector_ptrs;

    // used to keep track of whether we need to rebuild the above or not
    bool m_vectors_initialized = false;
    amrex::BoxArray::RefID m_BARef;
    amrex::DistributionMapping::RefID m_DMRef;

    //private:

    //void correctCellVectors(int old_index, int new_index, int grid, const ParticleType& p);

    int n_list;

    struct CHECK_PAIR {
        AMREX_GPU_HOST_DEVICE AMREX_INLINE
        bool operator()(const ParticleType & p1, const ParticleType & p2) const {
            return true;
        }
    };
};

#endif
