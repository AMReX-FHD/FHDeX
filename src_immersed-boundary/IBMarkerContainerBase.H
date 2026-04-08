#ifndef _IBMarkerContainerBase_H_
#define _IBMarkerContainerBase_H_


#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_Particles.H>
#include <AMReX_Periodicity.H>
#include <AMReX_NeighborParticles.H>
#include <IBParticleInfo.H>
#include <common_namespace.H>


using namespace amrex;


template <int NStructReal, int NStructInt>
class IBMarIterBase
    : public ParIter<NStructReal, NStructInt, 0, 0>
{

public:

    using ContainerType = ParticleContainer<NStructReal, NStructInt, 0, 0>;

    IBMarIterBase (ContainerType & pc, int level)
        : ParIter<NStructReal, NStructInt, 0, 0>(pc,level)
        {}

    IBMarIterBase (ContainerType & pc, int level, MFItInfo& info)
        : ParIter<NStructReal, NStructInt, 0, 0>(pc,level,info)
        {}
};



template <typename StructReal, typename StructInt>
class IBMarkerContainerBase
    : public NeighborParticleContainer<StructReal::count, StructInt::count>
{

public:

    using NeighborParticleContainer<StructReal::count, StructInt::count>
          ::NeighborParticleContainer;

    using MyConstIBMarIter = ParConstIter<StructReal::count, StructInt::count, 0, 0>;
    using MyIBMarIter      = IBMarIterBase<StructReal::count, StructInt::count>;

    using ParticleType = typename NeighborParticleContainer<StructReal::count, StructInt::count>::ParticleType;
    using PairIndex = typename NeighborParticleContainer<StructReal::count, StructInt::count>::PairIndex;
    using AoS = typename NeighborParticleContainer<StructReal::count, StructInt::count>::AoS;

    // indexing tiles (box index, local tile index)
    using TileIndex = std::pair<int, int>;
    // indexing particles and corresponding spot in the neighbor list
    using MarkerListIndex = std::pair<int, int>;


    IBMarkerContainerBase(AmrCore * amr_core, int n_nbhd);
    IBMarkerContainerBase(const Geometry & geom, const DistributionMapping & dmap,
                          const BoxArray & ba, int n_nbhd);

    IBMarkerContainerBase(const Geometry & geom, const Geometry & geomF, const DistributionMapping & dmap,
                          const BoxArray & ba, const BoxArray & baF, int n_nbhd, int ngF);

    virtual ~IBMarkerContainerBase() {};


    // (ID, initial CPU) tuple: unique to each particle
    using MarkerIndex = std::pair<int, int>;


    // Get number of particles
    int NumberOfMarkers(IBMarIterBase<StructReal::count, StructInt::count> & pti){
        return pti.GetArrayOfStructs().numParticles();
    };

    // Get the next particle ID -- can be used to find the number of particles
    // initialized on this GPU already
    int NextID() {
        return ParticleType::the_next_id;
    }

    const Gpu::ManagedDeviceVector<int> & getNumIDs() {
        return num_ids;
    }

    int getNumIDs(int rank) {
        return num_ids[rank];
    }

    const Gpu::ManagedDeviceVector<int> & getCPUOffset() {
        return cpu_offset;
    }

    int getCPUOffset(int cpu) {
        return cpu_offset[cpu];
    }

    int getTotalNumIDs() {
        return total_num_ids;
    }

    // Get total number of markers accross all ranks.
    // TODO: clean up this version
    int getTotalMarkers(){
        return totalMarkers;
    };

    // Load bonds from data file.
    void loadBonds(int totalP, char* filename);


    /****************************************************************************
     *                                                                          *
     * Marker Spreading and Interpolation                                       *
     *                                                                          *
     ***************************************************************************/

    void initRankLists(int totalM);
    void loadPinMatrix(int totalP, char* filename);
    void pinnedParticleInversion();

    int getTotalPinnedMarkers();

    //___________________________________________________________________________
    // Static spreading operators applied to raw data arrays
    static void SpreadMarkers(
                const Vector<RealVect> & f_in,
                const Vector<RealVect> & f_pos,
                      std::array<MultiFab, AMREX_SPACEDIM> & f_out,
                      std::array<MultiFab, AMREX_SPACEDIM> & f_weights,
                const std::array<MultiFab, AMREX_SPACEDIM> & coords,
                const Real* dx, int ghost
            );

    static void SpreadMarkers(
                const Vector<RealVect> & f_in,
                const Vector<RealVect> & f_pos,
                const Box & tile_box,
                      std::array<      FArrayBox *, AMREX_SPACEDIM> & f_out,
                      std::array<      FArrayBox *, AMREX_SPACEDIM> & f_weights,
                const std::array<const FArrayBox *, AMREX_SPACEDIM> & coords,
                const Real* dx, int ghost
            );

    static void SpreadKernel(const Box& bx,
                      std::array<      FArrayBox *, AMREX_SPACEDIM> & f_out,
                      std::array<      FArrayBox *, AMREX_SPACEDIM> & f_weights,
                      const std::array<const FArrayBox *, AMREX_SPACEDIM> & coords,
                      GpuArray<Real, AMREX_SPACEDIM>const& pos,
                      GpuArray<Real, AMREX_SPACEDIM>const& v_spread,
                      const Real* dx,
                      int* nghost);

    static void SpreadKernelGpu(const AoS& aos,
                      const Box& bx,
                      std::array<      FArrayBox *, AMREX_SPACEDIM> & f_out,
                      std::array<      FArrayBox *, AMREX_SPACEDIM> & f_weights,
                      const std::array<const FArrayBox *, AMREX_SPACEDIM> & coords,
                      const Real* dx,
                      int* nghost);
    //---------------------------------------------------------------------------



    //___________________________________________________________________________
    // Spreading operators applied to internal force fields (in amrex::Particles)
    void SpreadMarkers(int lev,
                       std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;

    void SpreadMarkers(int lev,
                             std::array<MultiFab, AMREX_SPACEDIM> & f_out,
                       const std::array<MultiFab, AMREX_SPACEDIM> & coords,
                       const Real* dx, int ghost) const;

    void SpreadMarkersGpu(int lev,
            std::array<MultiFab, AMREX_SPACEDIM> & f_out,
            const std::array<MultiFab, AMREX_SPACEDIM> & coords,
            const Real * dx, int ghost);

    void SpreadPredictor(int lev,
                         std::array<MultiFab, AMREX_SPACEDIM> & f_out) const;
    //---------------------------------------------------------------------------



    //___________________________________________________________________________
    // Static interpolation operators applied to raw data arrays
    static void InterpolateMarkers(
                      Vector<RealVect> & f_out,
                const Vector<RealVect> & f_pos,
                const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
                const std::array<MultiFab, AMREX_SPACEDIM> & f_weights,
                const std::array<MultiFab, AMREX_SPACEDIM> & coords,
                const Real * dx, const int ghost
            );

    static void InterpolateMarkers(
                      Vector<RealVect> & f_out,
                const Vector<RealVect> & f_pos,
                const Box & tile_box,
                const std::array<const FArrayBox *, AMREX_SPACEDIM> & f_in,
                const std::array<const FArrayBox *, AMREX_SPACEDIM> & f_weights,
                const std::array<const FArrayBox *, AMREX_SPACEDIM> & coords,
                const Real * dx, const int ghost
            );

    static void InterpolateKernelGpu(AoS& aos,
                      const Box& bx,
                      const std::array<const FArrayBox *, AMREX_SPACEDIM> & f_in,
                      const std::array<const FArrayBox *, AMREX_SPACEDIM> & f_weights,
                      const std::array<const FArrayBox *, AMREX_SPACEDIM> & coords,
                      const Real* dx,
                      const int* nghost,
                      int& check);

    //---------------------------------------------------------------------------



    void InterpolateMarkers(
            int lev,
            const Real * dx,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
            const std::array<MultiFab, AMREX_SPACEDIM> & coords
        );

    void InterpolateMarkersGpu(
            int lev,
            const Real * dx,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
            const std::array<MultiFab, AMREX_SPACEDIM> & coords);

   void InterpolateMarkersGpu(
            int lev,
            const Real * dx,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
            const std::array<MultiFab, AMREX_SPACEDIM> & coords, Real& checkR);
    //___________________________________________________________________________
    // Interpolation operators applied to internal force fields (in amrex::Particles)
    void InterpolateMarkers(int lev,
                            const std::array<MultiFab, AMREX_SPACEDIM> & f_in);

    void InterpolatePredictor(int lev,
                              const std::array<MultiFab, AMREX_SPACEDIM> & f_in);
    //---------------------------------------------------------------------------



    /****************************************************************************
     *                                                                          *
     * Physics                                                                  *
     *                                                                          *
     ***************************************************************************/


    void MoveMarkers(int lev, Real dt);
    void MovePredictor(int lev, Real dt);
    void ResetMarkers(int lev);
    void ResetPredictor(int lev);



    /****************************************************************************
     *                                                                          *
     * Stochastic drift                                                         *
     *                                                                          *
     ***************************************************************************/

    void RFD(int lev, const Real* dx,
                   std::array<MultiFab, AMREX_SPACEDIM> & f_out,
             const std::array<MultiFab, AMREX_SPACEDIM> & coords);



    /****************************************************************************
     *                                                                          *
     * Access marker (amrex particle) data                                      *
     *                                                                          *
     ***************************************************************************/

    void UpdatePIDMap();

    void PullDown(int lev, Real * list, int element, int totalParticles);

    void PullDown(int lev, Vector<Real> & list, int element);

    void PullDownInt(int lev, int * list, int element, int totalParticles);

    void PullDownInt(int lev, Vector<int> & list, int element);

    void PushUpAdd(int lev, Real * list, int element, int totalParticles);

    int get_nghost() const {return nghost;};

    void PrintMarkerData(int lev) const;

    const Vector<std::array<MultiFab, AMREX_SPACEDIM>> & get_face_coords() const {
        return  face_coords;
    }
    const std::array<MultiFab, AMREX_SPACEDIM> & get_face_coords(int lev) const {
        return  face_coords[lev];
    }

    std::array<MultiFab, AMREX_SPACEDIM> f_weights;


protected:

    int nghost;

    void InitInternals(int ngrow);
    void ReadStaticParameters();

    Real integrate_es(Real beta_in, int w_in);

    AmrCore * m_amr_core;

    Gpu::ManagedDeviceVector<int> ranksIdSorted;
    Gpu::ManagedDeviceVector<int> idsRankSorted;
    Gpu::ManagedDeviceVector<int> rankTotals;

    Gpu::ManagedDeviceVector<int> head_index;
    Gpu::ManagedDeviceVector<int> num_bond;
    Gpu::ManagedDeviceVector<int> bond_atom;
    Gpu::ManagedDeviceVector<amrex::Real> bond_coeff1;
    Gpu::ManagedDeviceVector<amrex::Real> bond_coeff2;

    Gpu::ManagedDeviceVector<amrex::Real> pinMatrix;

    int totalMarkers;
    int totalPinnedMarkers;

    // TODO: why is this static inline?
    static inline Gpu::ManagedDeviceVector<Real> norm_es;

private:

    // Positions on faces
    Vector<std::array<MultiFab, AMREX_SPACEDIM>> face_coords;

    // Number of paricle IDs for each rank:
    Gpu::ManagedDeviceVector<int> num_ids;
    // Offset of each rank in the particle list
    Gpu::ManagedDeviceVector<int> cpu_offset;
    // Number of markers -- TODO: deal with duplicaton with totalMarkers
    int total_num_ids;
};


#include "IBMarkerContainerBaseI.H"


#endif
