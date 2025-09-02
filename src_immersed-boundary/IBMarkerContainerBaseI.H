#ifndef _IBMarkerContainerBaseI_H_
#define _IBMarkerContainerBaseI_H_



#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)

#include <common_functions.H>

#include <IBMarkerContainer.H>
#include <ib_functions_F.H>
#include "rng_functions.H"

#include "kernel_functions_K.H"

#include "particle_functions.H"

#include <iostream>
#include <fstream>

using namespace amrex;

#include <type_traits>

//Returns true or false depending on whether the particle has an integer member called 'visible'
template <typename T1, typename = void>
struct has_visible : std::false_type{};

template <typename T1>
struct has_visible<T1, decltype((void)T1::visible, void())> : std::true_type {};

#define HAS_VISIBLE has_visible<StructInt>::value

//Returns true or false depending on whether the particle has an integer member called 'visible'
template <typename T2, typename = void>
struct has_pinned : std::false_type{};

template <typename T2>
struct has_pinned<T2, decltype((void)T2::pinned, void())> : std::true_type {};

#define HAS_PINNED has_pinned<StructInt>::value


// Needs to be specified once a concrete class is created
// bool IBMarkerContainer::use_neighbor_list  {true};
// bool IBMarkerContainer::sort_neighbor_list {false};



template <typename StructReal, typename StructInt>
IBMarkerContainerBase<StructReal, StructInt>::IBMarkerContainerBase(
            const Geometry & geom,
            const DistributionMapping & dmap,
            const BoxArray & ba,
            int n_nbhd
        ) : NeighborParticleContainer<StructReal::count, StructInt::count>(
            geom, dmap, ba, n_nbhd
        ),
    nghost(n_nbhd)
{
    InitInternals(n_nbhd);


    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(ba, dmap, 1, n_nbhd);
        f_weights[d].setVal(-1.0);
    }


}

template <typename StructReal, typename StructInt>
IBMarkerContainerBase<StructReal, StructInt>::IBMarkerContainerBase(
            const Geometry & geom,
            const Geometry & geomF,
            const DistributionMapping & dmap,
            const BoxArray & ba,
            const BoxArray & baF,
            int n_nbhd,
            int ngF
        ) : NeighborParticleContainer<StructReal::count, StructInt::count>(
            geom, dmap, ba, n_nbhd
        ),
    nghost(n_nbhd)
{
    InitInternals(n_nbhd);

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(baF, dmap, 1, ngF);
        f_weights[d].setVal(-1.0);
    }


}



template <typename StructReal, typename StructInt>
IBMarkerContainerBase<StructReal, StructInt>::IBMarkerContainerBase(
            AmrCore * amr_core, int n_nbhd
        ) : NeighborParticleContainer<StructReal::count, StructReal::count>(
            amr_core->GetParGDB(), n_nbhd
        ),
    m_amr_core(amr_core),
    nghost(n_nbhd)
{
    InitInternals(n_nbhd);
}



template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::MoveMarkers(int lev, Real dt) {

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        for (int i = 0; i < np; ++ i) {
            ParticleType & part = particles[i];

            for (int d=0; d<AMREX_SPACEDIM; ++d)
                part.pos(d) += dt * part.rdata(StructReal::velx + d);
        }
    }
}

template <typename StructReal, typename StructInt>
int IBMarkerContainerBase<StructReal, StructInt>::getTotalPinnedMarkers() {

    return totalPinnedMarkers;
}

template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::RFD(
            int lev, const Real * dx,
                  std::array<MultiFab, AMREX_SPACEDIM> & f_out,
            const std::array<MultiFab, AMREX_SPACEDIM> & coords)
{
    // timer for profiling
    BL_PROFILE_VAR("RFD()",RFD);

    Real delta = 1e-4*dx[0];
    Real deltaover2 = delta/2.0;
    Real force = variance_coef_mom*k_B*T_init[0]/delta;
    const Geometry & geom = this->Geom(0);

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        const int np = this->GetParticles(lev)[index].numRealParticles();
        auto& plev = this->GetParticles(lev);
        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();

        amrex::ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
        {
            ParticleType & part = particles[i];
            int pinned = 0;
            if(HAS_PINNED)
            {
                pinned = part.idata(StructInt::pinned);
            }
            if(pinned==0)
            {
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    // using pred_pos to temporarily store the random direction -
                    // hopefully this doesn't cause problems elsewhere
//                part.rdata(StructReal::pred_posx + d) =  get_particle_normal_func();
                    part.rdata(StructReal::pred_posx + d) = amrex::RandomNormal(0.,1.,engine);

                    // move particle to q^n + (dt/2) * Wtilde
                    part.pos(d) += deltaover2 * part.rdata(StructReal::pred_posx + d);
                    // scale by k_B*T/delta
                    // NOTE: this assumed that the pred_pos* fields are
                    // displacement vectors wrt pos(*)

                    part.rdata(StructReal::forcex + d) =
                        force * part.rdata(StructReal::pred_posx + d);

                }
            }
        });
    }

    // spread this term into the velocity forcing
    SpreadMarkersGpu(lev, f_out, coords, dx, 1);

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        const int np = this->GetParticles(lev)[index].numRealParticles();
        auto& plev = this->GetParticles(lev);
        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();

        AMREX_FOR_1D( np, i,
        {
            ParticleType & part = particles[i];
            int pinned = 0;
            if(HAS_PINNED)
            {
                pinned = part.idata(StructInt::pinned);
            }
            if(pinned==0)
            {
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                // move particle to q^n - (dt/2) * Wtilde (by subtracting dt *
                // Wtilde)
                part.pos(d) -= delta * part.rdata(StructReal::pred_posx + d);
                // scale by -k_B*T/delta
                // NOTE: this assumed that the pred_pos* fields are
                // displacement vectors wrt pos(*)

                part.rdata(StructReal::forcex + d) =
                    -force * part.rdata(StructReal::pred_posx + d);

            }
            }
        });
    }

    // spread this term into the velocity forcing
    SpreadMarkersGpu(lev, f_out, coords, dx, 1);

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        const int np = this->GetParticles(lev)[index].numRealParticles();
        auto& plev = this->GetParticles(lev);
        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();

        AMREX_FOR_1D( np, i,
        {
            ParticleType & part = particles[i];
            int pinned = 0;
            if(HAS_PINNED)
            {
                pinned = part.idata(StructInt::pinned);
            }
            if(pinned==0)
            {
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                // move particle back to original position q^n (by adding dt/2
                // * Wtilde)
                part.pos(d) += deltaover2 * part.rdata(StructReal::pred_posx + d);
                // scale by 0
                part.rdata(StructReal::forcex + d) = 0;
            }
            }
        });

    }

    f_out[0].SumBoundary(geom.periodicity());
    f_out[1].SumBoundary(geom.periodicity());
    f_out[2].SumBoundary(geom.periodicity());
}


template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::initRankLists(int totalM)
{
    // timer for profiling
    BL_PROFILE_VAR("initRankLists()",initRankLists);

    totalMarkers = totalM;

    ranksIdSorted.resize(totalMarkers, 0);
    idsRankSorted.resize(totalMarkers, 0);
    rankTotals.resize(ParallelDescriptor::NProcs(), 0);

    const int lev = 0;

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        const int np = this->GetParticles(lev)[index].numRealParticles();
        auto& plev = this->GetParticles(lev);
        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();

        // THIS DOES NOT COMPILE WITH CUDA?
        for(int i=0; i<np;i++)
        {
            ParticleType & part = particles[i];

            ranksIdSorted[part.id()-1] = part.cpu();
            rankTotals[ParallelDescriptor::MyProc()]++;
        }
    }

    // ParallelDescriptor::ReduceIntSum(idsRankSorted);
    // ParallelDescriptor::ReduceIntSum(rankTotals);

    for (int i=0; i<totalMarkers; ++i) {
        int temp = idsRankSorted[i];
        ParallelDescriptor::ReduceIntSum(temp);
        idsRankSorted[i] = temp;
    }

    for (int i=0; i<ParallelDescriptor::NProcs(); ++i) {
        int temp = rankTotals[i];
        ParallelDescriptor::ReduceIntSum(temp);
        rankTotals[i] = temp;
    }

}


template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::loadPinMatrix(int totalP, char* filename)
{
    // timer for profiling
    BL_PROFILE_VAR("initRankLists()",initRankLists);

    totalPinnedMarkers = totalP;

    pinMatrix.resize(pow(3*totalP,2),0);

    std::ifstream ifs("invOut", std::ios::binary);

    //std::ifstream matFile("invOut");

    for(int i=0; i<pow(3*totalP,2); i++)
    {
        double read;
        ifs.read(reinterpret_cast<char*>( &read ), sizeof read);
        pinMatrix[i] = read;

        //matFile >> pinMatrix[i];
        if(i%1000000==0)
        {
            Print() << i << " of " << pow(3*totalP,2) << ": " << pinMatrix[i] << std::endl;
        }

    }

    ifs.close();
    //matFile.close();

}

// load bond info into two giant arrays, one for bonded particle IDs and one for bond type
template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::loadBonds(int totalP, char* filename)
{
    // timer for profiling
    BL_PROFILE_VAR("loadBonds()",loadBonds);

    int atom1, atom2, atom_current;
    Real k_bond, x0_bond;

    head_index.resize(totalP,0);
    num_bond.resize(totalP,0);
    bond_atom.resize(common::nbonds,0);
    bond_coeff1.resize(common::nbonds,0);
    bond_coeff2.resize(common::nbonds,0);

    int bcount = 0;
    atom_current = -1;

    std::ifstream bondFile(filename); // bond file is generated so that atom2 >= atom1
    std::string line;
    while (std::getline(bondFile, line)) {
    //for(int i_bond=0; i_bond < common::nbonds; i_bond++) {

        std::istringstream ss(line);

        ss >> atom1;
        ss >> atom2;
        ss >> k_bond;
        ss >> x0_bond;

        //Print() << atom1 << " " << atom2 << " " << k_bond << std::endl;

        if(atom1>atom_current)
        {
            atom_current = atom1;
            head_index[atom_current] = bcount;
        }

        bond_coeff1[bcount] = k_bond;
        bond_coeff2[bcount] = x0_bond;
        bond_atom[bcount] = atom2;
        num_bond[atom1] = num_bond[atom1]+1;

        bcount++;
    }

    if (bcount != common::nbonds) {
        Abort("Missing bond info!");
    } else {
        Print() << "Loaded " << bcount << " bonds." << std::endl;
    }

    bondFile.close();

    for (int i=0; i<totalP; i++) {
        Print() << i << " of " << totalP << ": " << head_index[i] << " " << num_bond[i] << std::endl;
    }

}

template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::MovePredictor(int lev, Real dt) {

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());

        AoS & particles = this->GetParticles(lev).at(index).GetArrayOfStructs();
        long np = this->GetParticles(lev).at(index).numParticles();

        for (int i = 0; i < np; ++ i) {
            ParticleType & part = particles[i];

            for (int d=0; d<AMREX_SPACEDIM; ++d)
                part.rdata(StructReal::pred_posx + d) +=
                    dt * part.rdata(StructReal::pred_velx + d);
        }
    }
}

//Find a more efficient way to do this...
template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::pinnedParticleInversion()
{
   // Print() << "inv1, " << totalMarkers;
    Real velx[totalMarkers];
    Real vely[totalMarkers];
    Real velz[totalMarkers];
    int  pinned[totalMarkers];
    Real  forcex[totalMarkers];
    Real  forcey[totalMarkers];
    Real  forcez[totalMarkers];

    PullDown(0, velx, StructReal::velx, totalMarkers);
    PullDown(0, vely, StructReal::vely, totalMarkers);
    PullDown(0, velz, StructReal::velz, totalMarkers);
    PullDownInt(0, pinned, StructInt::pinned, totalMarkers);
  //  Print() << "inv2\n";

    Real rhs[totalPinnedMarkers*3];
    Real lhs[totalPinnedMarkers*3];
    int k = 0;
    Real rhsNorm = 0;
    for(int i=0;i<totalMarkers;i++)
    {
        if(pinned[i] == 1)
        {
            rhs[k] = -velx[i];
            rhs[k+1] = -vely[i];
            rhs[k+2] = -velz[i];

           // Print() << velx[i] << ", " << vely[i] << ", " << velz[i] << std::endl;
          //  Print() << rhs[k] << ", " << rhs[k+1] << ", " << rhs[k+2] << std::endl;
            rhsNorm = rhsNorm + rhs[k]*rhs[k] + rhs[k+1]*rhs[k+1] + rhs[k+2]*rhs[k+2];

            k = k+3;
        }
    }
    //Print() << "rhsNorm: " << sqrt(rhsNorm) << std::endl;
    //Print() << 3*totalPinnedMarkers << std::endl;
    for(int i=0;i<(3*totalPinnedMarkers);i++)
    {
        lhs[i]=0;
        for(int j=0;j<(3*totalPinnedMarkers);j++)
        {
            lhs[i] = lhs[i] + pinMatrix[i*3*totalPinnedMarkers + j]*rhs[j];
            //Print() << i << ", " << j << ", " << i*3*totalPinnedMarkers + j << ", " << pinMatrix[i*totalPinnedMarkers + j] << std::endl;

        }
    }

    k=0;
    for(int i=0;i<totalMarkers;i++)
    {
        if(pinned[i] == 1)
        {
            forcex[i] = lhs[k];
            forcey[i] = lhs[k+1];
            forcez[i] = lhs[k+2];

            k = k+3;

            //Print() << forcex[i] << ", " << forcey[i] << ", " << forcez[i] << std::endl;
        }

    }

    int lev=0;
    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        const int np = this->GetParticles(lev)[index].numRealParticles();
        auto& plev = this->GetParticles(lev);
        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();

//            AMREX_FOR_1D( np, i,
            for(int i=0; i<np; i++)
            {
                ParticleType & part = particles[i];
                if(part.idata(StructInt::pinned) == 1)
                {
                    part.rdata(StructReal::forcex) = forcex[part.id()-1];
                    part.rdata(StructReal::forcey) = forcey[part.id()-1];
                    part.rdata(StructReal::forcez) = forcez[part.id()-1];

                    //Print() << part.id() << ", " << part.rdata(StructReal::forcex) << ", " << part.rdata(StructReal::forcey) << ", " << part.rdata(StructReal::forcez) << std::endl;
                }

            }
//            );

    }
}


template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::PullDown(
            int lev, Real * list, int element, int totalParticles)
{
    // timer for profiling
    BL_PROFILE_VAR("PullDown()",PullDown);

    Real temp;

    for (int i = 0; i < totalParticles; ++i) {
        //std::cout << i << " of " << totalParticles << std::endl;
        list[i] = 0;
    }


    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        const int np = this->GetParticles(lev)[index].numRealParticles();
        auto& plev = this->GetParticles(lev);
        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();

        //std::cout << "Proc " << ParallelDescriptor::MyProc() << " doing " << np << " particles\n";

        if (element >= 0) {

            //std::cout << "Proc " << ParallelDescriptor::MyProc() << " 1a start\n";
            AMREX_FOR_1D( np, i,
            {
                ParticleType & part = particles[i];
                int id = part.id();
                //TODO: Why "radius"?!
                list[id-1] = particles[i].rdata(StructReal::radius + element);
            });
            //std::cout << "Proc " << ParallelDescriptor::MyProc() << " 1a finish\n";

        } else {

            //std::cout << "Proc " << ParallelDescriptor::MyProc() << " 1b start\n";
            AMREX_FOR_1D( np, i,
            {
                ParticleType & part = particles[i];
                int id = part.id();
                list[id-1] = part.pos((-element)-1);

            //  std::cout << "Proc " << ParallelDescriptor::MyProc() << " storing " << part.pos((-element)-1)
            //                               << " from " << i << " in position " << id-1 << std::endl;
            });
            //std::cout << "Proc " << ParallelDescriptor::MyProc() << " 1b finish\n";
        }

    }
    //std::cout << "Proc " << ParallelDescriptor::MyProc() << " 2\n";
    //I'm sure there is an array version of this but this will do for now.
    // There is: void amrex::ParallelDescriptor::ReduceRealSum (Vector< std::reference_wrapper< Real > > &&)
    // the commented code below should do the trick (but it has not been tested)
    for (int i=0; i<totalParticles; ++i) {
        temp = list[i];
        ParallelDescriptor::ReduceRealSum(temp);
        list[i] = temp;
    }
    // Vectorized alternative (warning: untested!):
    // // Create a vector wrapping the `list` array:
    // Vector<std::reference_wrapper<Real>> tmp(totalParticles);
    // for (int i=0; i<totalParticles; ++i) tmp[i] = & list[i];

    // // Vectorized parallel reduce:
    // ParallelDescriptor::ReduceRealSum(std::move(tmp)); // Ask Weiqun why tmp needs to be &&?
    //std::cout << "Proc " << ParallelDescriptor::MyProc() << " 3\n";

}


template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::UpdatePIDMap() {

    BL_PROFILE_VAR("UpdatePIDMap()", UpdatePIDMap);

    int rank = ParallelDescriptor::MyProc();
    int size = ParallelDescriptor::NProcs();

    // Number of IDs associated with this rank
    num_ids[rank] = this->NextID() - 1;

    ParallelDescriptor::ReduceIntSum(
        Vector<std::reference_wrapper<int>>(std::begin(num_ids), std::end(num_ids))
    );

    int offset = 0;
    for (int i=0; i < size; i++) {
        cpu_offset[i] = offset;
        offset += num_ids[i];
    }

    total_num_ids = offset;
}


template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::PullDown(
            int lev, Vector<Real> & list, int element
) {
    // timer for profiling
    BL_PROFILE_VAR("PullDown()", PullDown);

//    Gpu::ManagedDeviceVector<int> cpu_offset;
//    Gpu::ManagedDeviceVector<int> cpu_offset;

    int * cpu_offset_ptr = cpu_offset.data();
    Real * list_ptr      = list.data();

    for (int i = 0; i < total_num_ids; ++i) { list[i] = 0;}

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {
        PairIndex index(pti.index(), pti.LocalTileIndex());
        const int np = this->GetParticles(lev)[index].numRealParticles();
        auto & plev  = this->GetParticles(lev);
        auto & ptile = plev[index];
        auto & aos   = ptile.GetArrayOfStructs();
        ParticleType * particles = aos().dataPtr();

        if (element >= 0) {
            AMREX_FOR_1D(np, i, {
                ParticleType & part = particles[i];
                int id  = part.id();
                int cpu = part.cpu();

                // Particle IDs start at 1, CPUs at 0 -- urgh!
                list_ptr[id-1 + cpu_offset_ptr[cpu]] = particles[i].rdata(element);
            });
        } else {
            AMREX_FOR_1D(np, i, {
                ParticleType & part = particles[i];
                int id  = part.id();
                int cpu = part.cpu();

                // Particle IDs start at 1, CPUs at 0 -- urgh!
                list_ptr[id-1 + cpu_offset_ptr[cpu]] = part.pos((-element)-1);
            });
        }

    }

    ParallelDescriptor::ReduceRealSum(
        Vector<std::reference_wrapper<Real>>(begin(list), end(list))
    );
}


template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::PullDownInt(
            int lev, int * list, int element, int totalParticles)
{
    // timer for profiling
    BL_PROFILE_VAR("PullDownInt()",PullDownInt);

    int temp;

    for (int i = 0; i < totalParticles; ++i) {
        //std::cout << i << " of " << totalParticles << std::endl;
        list[i] = 0;
    }


    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        const int np = this->GetParticles(lev)[index].numRealParticles();
        auto& plev = this->GetParticles(lev);
        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();

        if (element >= 0) {

            AMREX_FOR_1D( np, i,
            {
                ParticleType & part = particles[i];
                int id = part.id();
                //std::cout << id << " of " << totalParticles << std::endl;
                // TODO: Why "sorted"?
                list[id-1] = particles[i].idata(StructInt::sorted + element);
            });

        }else
        {
            AMREX_FOR_1D( np, i,
            {
                ParticleType & part = particles[i];
                int id = part.id();
                //std::cout << id << " of " << totalParticles << std::endl;
                list[id-1] = particles[i].cpu();
            });
        }
    }

    //I'm sure there is an array version of this but this will do for now.
    // There is: void amrex::ParallelDescriptor::ReduceRealSum (Vector< std::reference_wrapper< Real > > &&)
    // the commented code below should do the trick (but it has not been tested)
    for (int i=0; i<totalParticles; ++i) {
        temp = list[i];
        ParallelDescriptor::ReduceIntSum(temp);
        list[i] = temp;
    }
    // Vectorized alternative (warning: untested!):
    // // Create a vector wrapping the `list` array:
    // Vector<std::reference_wrapper<Real>> tmp(totalParticles);
    // for (int i=0; i<totalParticles; ++i) tmp[i] = & list[i];

    // // Vectorized parallel reduce:
    // ParallelDescriptor::ReduceRealSum(std::move(tmp)); // Ask Weiqun why tmp needs to be &&?
}


template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::PullDownInt(
            int lev, Vector<int> & list, int element
) {
    // timer for profiling
    BL_PROFILE_VAR("PullDownInt()", PullDownInt);

    int * cpu_offset_ptr = cpu_offset.data();
    int * list_ptr       = list.data();

    for (int i = 0; i < total_num_ids; ++i) { list[i] = 0;}

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {
        PairIndex index(pti.index(), pti.LocalTileIndex());
        const int np = this->GetParticles(lev)[index].numRealParticles();
        auto & plev  = this->GetParticles(lev);
        auto & ptile = plev[index];
        auto & aos   = ptile.GetArrayOfStructs();
        ParticleType * particles = aos().dataPtr();

        if (element >= 0) {
            AMREX_FOR_1D(np, i, {
                ParticleType & part = particles[i];
                int id  = part.id();
                int cpu = part.cpu();

                // Particle IDs start at 1, CPUs at 0 -- urgh!
                list_ptr[id-1 + cpu_offset_ptr[cpu]] = particles[i].idata(element);
            });
        } else {
            AMREX_FOR_1D(np, i, {
                ParticleType & part = particles[i];
                int id  = part.id();
                int cpu = part.cpu();

                // Particle IDs start at 1, CPUs at 0 -- urgh!
                if (element == -1) {
                    list_ptr[id-1 + cpu_offset_ptr[cpu]] = part.id();
                } else if (element == -2) {
                    list_ptr[id-1 + cpu_offset_ptr[cpu]] = part.cpu();
                } else{
                    list_ptr[id-1 + cpu_offset_ptr[cpu]] = -1;
                }
            });
        }

    }

    ParallelDescriptor::ReduceIntSum(
        Vector<std::reference_wrapper<int>>(begin(list), end(list))
    );
}

template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::PushUpAdd(
            int lev, Real * list, int element, int totalParticles)
{
    // timer for profiling
    BL_PROFILE_VAR("PushUpAdd()",PushUpAdd);

    Real temp;

    //I'm sure there is an array version of this but this will do for now.
    // There is: void amrex::ParallelDescriptor::ReduceRealSum (Vector< std::reference_wrapper< Real > > &&)
    // the commented code below should do the trick (but it has not been tested)
    for(int i=0; i<totalParticles; ++i) {
        temp = list[i];
        ParallelDescriptor::ReduceRealSum(temp);
        list[i] = temp;
    }
    // Vectorized alternative (warning: untested!):
    // // Create a vector wrapping the `list` array:
    // Vector<std::reference_wrapper<Real>> tmp(totalParticles);
    // for (int i=0; i<totalParticles; ++i) tmp[i] = & list[i];

    // // Vectorized parallel reduce:
    // ParallelDescriptor::ReduceRealSum(std::move(tmp)); // Ask Weiqun why tmp needs to be &&?


    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        const int np = this->GetParticles(lev)[index].numRealParticles();
        auto& plev = this->GetParticles(lev);
        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();

        if(element >= 0) {

            AMREX_FOR_1D( np, i,
            {
                ParticleType & part = particles[i];
                int id = part.id();

                particles[i].rdata(element) = particles[i].rdata(element) + list[id-1];
            });

        } else {

            AMREX_FOR_1D( np, i,
            {
                ParticleType & part = particles[i];
                int id = part.id();

                particles[i].pos((-element)-1) = particles[i].pos((-element)-1) + list[id-1];
            });
        }
    }
}



template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::ResetMarkers(int lev)
{
    // timer for profiling
    BL_PROFILE_VAR("ResetMarkers()",ResetMarkers);

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        const int np = this->GetParticles(lev)[index].numRealParticles();
        auto& plev = this->GetParticles(lev);
        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        ParticleType* particles = aos().dataPtr();

        AMREX_FOR_1D( np, i,
        {
            ParticleType & part = particles[i];

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                part.rdata(StructReal::velx + d)   = 0.;
                part.rdata(StructReal::forcex + d) = 0.;
            }
        });
    }
}



template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::ResetPredictor(int lev) {

    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        TileIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = this->GetParticles(lev).at(index);
        long np = this->GetParticles(lev).at(index).numParticles();

        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                part.rdata(StructReal::pred_posx + d)   = 0.;
                part.rdata(StructReal::pred_velx + d)   = 0.;
                part.rdata(StructReal::pred_forcex + d) = 0.;
            }
        }
    }
}



template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::SpreadMarkers(
            const Vector<RealVect> & f_in,
            const Vector<RealVect> & f_pos,
                  std::array<MultiFab, AMREX_SPACEDIM> & f_out,
                  std::array<MultiFab, AMREX_SPACEDIM> & f_weights,
            const std::array<MultiFab, AMREX_SPACEDIM> & coords,
            const Real * dx, const int ghost
        ) {

    for (MFIter mfi(f_out[0]); mfi.isValid(); ++mfi) {

        // Convert stagered boxes into cell centered
        Box tile_box = enclosedCells(mfi.tilebox());

        std::array<      FArrayBox *, AMREX_SPACEDIM> f_out_fab;
        std::array<      FArrayBox *, AMREX_SPACEDIM> f_weights_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> coords_fab;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            f_out_fab[d]     = & f_out[d][mfi];
            f_weights_fab[d] = & f_weights[d][mfi];
            coords_fab[d]    = & coords[d][mfi];
        }

        SpreadMarkers(f_in, f_pos, tile_box,
                      f_out_fab, f_weights_fab, coords_fab, dx, ghost);
    }
}

//template <typename StructReal, typename StructInt>
//void IBMarkerContainerBase<StructReal, StructInt>::SpreadMarkersCC(
//            const Vector<RealVect> & f_in,
//            const Vector<RealVect> & f_pos,
//                  std::array<MultiFab, AMREX_SPACEDIM> & f_out,
//                  std::array<MultiFab, AMREX_SPACEDIM> & f_weights,
//            const std::array<MultiFab, AMREX_SPACEDIM> & coords,
//            const Real * dx, const int ghost
//        ) {

//    for (MFIter mfi(f_out[0]); mfi.isValid(); ++mfi) {

//        // Convert stagered boxes into cell centered
//        Box tile_box = enclosedCells(mfi.tilebox());

//        std::array<      FArrayBox *, AMREX_SPACEDIM> f_out_fab;
//        std::array<      FArrayBox *, AMREX_SPACEDIM> f_weights_fab;
//        std::array<const FArrayBox *, AMREX_SPACEDIM> coords_fab;
//        for (int d=0; d<AMREX_SPACEDIM; ++d) {
//            f_out_fab[d]     = & f_out[d][mfi];
//            f_weights_fab[d] = & f_weights[d][mfi];
//            coords_fab[d]    = & coords[d][mfi];
//        }

//        SpreadMarkers(f_in, f_pos, tile_box,
//                      f_out_fab, f_weights_fab, coords_fab, dx, ghost);
//    }
//}



template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::SpreadMarkers(
            const Vector<RealVect> & f_in,
            const Vector<RealVect> & f_pos,
            const Box & tile_box,
                  std::array<      FArrayBox *, AMREX_SPACEDIM> & f_out,
                  std::array<      FArrayBox *, AMREX_SPACEDIM> & f_weights,
            const std::array<const FArrayBox *, AMREX_SPACEDIM> & coords,
            const Real* dx, int ghost)
{
    // timer for profiling
    BL_PROFILE_VAR("SpreadMarkers()",SpreadMarkers);

    const int n_marker = f_pos.size();

    Vector<int> pkernel_fluid_in(1);
    for (int i=0; i<1; ++i) {
        pkernel_fluid_in[i] = pkernel_fluid[i];
    }

    // Convert stagered boxes into cell centered
    Box bx = enclosedCells(f_out[0]->box());

    spread_markers(BL_TO_FORTRAN_BOX(bx),
                   BL_TO_FORTRAN_BOX(tile_box),
                   BL_TO_FORTRAN_ANYD(* f_out[0]),
#if (AMREX_SPACEDIM > 1)
                   BL_TO_FORTRAN_ANYD(* f_out[1]),
#endif
#if (AMREX_SPACEDIM > 2)
                   BL_TO_FORTRAN_ANYD(* f_out[2]),
#endif
                   BL_TO_FORTRAN_ANYD(* f_weights[0]),
#if (AMREX_SPACEDIM > 1)
                   BL_TO_FORTRAN_ANYD(* f_weights[1]),
#endif
#if (AMREX_SPACEDIM > 2)
                   BL_TO_FORTRAN_ANYD(* f_weights[2]),
#endif
                   BL_TO_FORTRAN_ANYD(* coords[0]),
#if (AMREX_SPACEDIM > 1)
                   BL_TO_FORTRAN_ANYD(* coords[1]),
#endif
#if (AMREX_SPACEDIM > 2)
                   BL_TO_FORTRAN_ANYD(* coords[2]),
#endif
                   f_pos.dataPtr(),
                   f_in.dataPtr(),
                   & n_marker,
                   dx, & ghost,
                   pkernel_fluid_in.dataPtr());

}

template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::SpreadKernelGpu(const AoS& aos,
         const Box& bx,
         std::array<     FArrayBox *, AMREX_SPACEDIM> & f_out,
         std::array<     FArrayBox *, AMREX_SPACEDIM> & f_weights,
         const std::array<const FArrayBox *, AMREX_SPACEDIM> & coords,
         const Real* dx,
         int* nghost)
{
    // timer for profiling
    BL_PROFILE_VAR("SpreadKernelGpu()",SpreadKernelGpu);

    GpuArray<Real, AMREX_SPACEDIM> invdx;
    Real invvol = 1.0;

    for (int i=0; i<AMREX_SPACEDIM; ++i)
        invdx[i]=1.0/dx[i];

    for (int i=0; i<AMREX_SPACEDIM; ++i)
        invvol *= invdx[i];

    const int spaceDim = AMREX_SPACEDIM;

    GpuArray<int, 3> bx_lo = {bx.loVect()[0], bx.loVect()[1], bx.loVect()[2]};
    GpuArray<int, 3> bx_hi = {bx.hiVect()[0], bx.hiVect()[1], bx.hiVect()[2]};

    Array4<Real> const& fout_x = f_out[0]->array();
    Array4<Real> const& fout_y = f_out[1]->array();
    Array4<Real> const& fout_z = f_out[2]->array();
    Array4<Real> const& fweights_x = f_weights[0]->array();
    Array4<Real> const& fweights_y = f_weights[1]->array();
    Array4<Real> const& fweights_z = f_weights[2]->array();
    Array4<const Real> const& coords_x = coords[0]->array();
    Array4<const Real> const& coords_y = coords[1]->array();
    Array4<const Real> const& coords_z = coords[2]->array();

    const auto Np = aos.numParticles();
    const auto pstruct = aos().dataPtr();
//    const int pk = pkernel_fluid;
    int ng = *nghost;
    Real* norm_ptr = norm_es.data();

    AMREX_FOR_1D( Np, ip,
    {
        const ParticleType& p = pstruct[ip];
        //Real norm_es = integral(eskernel_beta[p.idata(StructInt::species)-1], eskernel_fluid[p.idata(StructInt::species)-1]);

        int lo_dim[3];
        int hi_dim[3];
        int gs;

        if(pkernel_fluid[p.idata(StructInt::species)-1] == 3)
        {
            gs = 2;
        }
        else if (pkernel_fluid[p.idata(StructInt::species)-1] == 4)
        {
            gs = 3;
        }
        else if (pkernel_fluid[p.idata(StructInt::species)-1] == 1)
        {
            gs = 1;
        }
        else if (pkernel_fluid[p.idata(StructInt::species)-1] == 6)
        {
            gs = 4;
        }
        else if (eskernel_fluid[p.idata(StructInt::species)-1] > 0)
        {
            gs = static_cast<int>(floor(eskernel_fluid[p.idata(StructInt::species)-1]*0.5)+1);
        }
        //Print() << "gs for spreading is " << gs << "\n";

        if (ng == 0)
        {
            int loc = static_cast<int>(p.pos(0) * invdx[0] - gs);
            lo_dim[0] = (loc >= bx_lo[0]) ? loc : bx_lo[0];
            loc = static_cast<int>(p.pos(0) * invdx[0] + gs);
            hi_dim[0] = (loc <= bx_hi[0]) ? loc : bx_hi[0];
            loc = static_cast<int>(p.pos(1) * invdx[1] - gs);
            lo_dim[1] = (loc >= bx_lo[1]) ? loc : bx_lo[1];
            loc = static_cast<int>(p.pos(1) * invdx[1] + gs);
            hi_dim[1] = (loc <= bx_hi[1]) ? loc : bx_hi[1];
            loc = static_cast<int>(p.pos(2) * invdx[2] - gs);
            lo_dim[2] = (loc >= bx_lo[2]) ? loc : bx_lo[2];
            loc = static_cast<int>(p.pos(2) * invdx[2] + gs);
            hi_dim[2] = (loc <= bx_hi[2]) ? loc : bx_hi[2];
        }
        else
        {
            lo_dim[0] = static_cast<int>(p.pos(0) * invdx[0] - gs);
            hi_dim[0] = static_cast<int>(p.pos(0) * invdx[0] + gs);
            lo_dim[1] = static_cast<int>(p.pos(1) * invdx[1] - gs);
            hi_dim[1] = static_cast<int>(p.pos(1) * invdx[1] + gs);
            lo_dim[2] = static_cast<int>(p.pos(2) * invdx[2] - gs);
            hi_dim[2] = static_cast<int>(p.pos(2) * invdx[2] + gs);
        }

        int vis = 1;
        if(HAS_VISIBLE)
        {
            //Print() << "Has visible, "<< p.idata(StructInt::visible) << std::endl;
            if(p.idata(StructInt::visible) != 1)
            {
                vis = 0;
            }
        }

        //std::cout << "normalization factor " << 1 << " are " << norm_es[1] << std::endl;
        // x-components
        if(vis == 1)
        {
            for (int k = lo_dim[2] ; k < hi_dim[2]; ++k) {
                for (int j = lo_dim[1] ; j < hi_dim[1]; ++j) {
                    for (int i = lo_dim[0] ; i < hi_dim[0]+1; ++i) {
                      Real pos_grid[3];
                      pos_grid[0]=(p.pos(0)-coords_x(i,j,k,0))*invdx[0];
                      pos_grid[1]=(p.pos(1)-coords_x(i,j,k,1))*invdx[1];
                      pos_grid[2]=(p.pos(2)-coords_x(i,j,k,2))*invdx[2];
                      Real pos_grid_norm = std::sqrt(pos_grid[0]*pos_grid[0]+pos_grid[1]*pos_grid[1]+pos_grid[2]*pos_grid[2]);

                      Real weight = 1.0;
                      for (int ii=0; ii<spaceDim; ++ii)
                      {
                         if(pkernel_fluid[p.idata(StructInt::species)-1] == 3)
                         {
                            weight*=Kernel3P()(pos_grid[ii]);
                         }
                         else if (pkernel_fluid[p.idata(StructInt::species)-1] == 4)
                         {
                            weight*=Kernel4P()(pos_grid[ii]);
                         }
                         else if (pkernel_fluid[p.idata(StructInt::species)-1] == 1)
                         {
                            weight*=Kernel1P()(pos_grid[ii]);
                         }
                         else if (pkernel_fluid[p.idata(StructInt::species)-1] == 6)
                         {
                            weight*=Kernel6P()(pos_grid[ii]);
                         }
                         else if (eskernel_fluid[p.idata(StructInt::species)-1] > 0)
                         {
                            weight*=KernelES()(pos_grid[ii],eskernel_beta[p.idata(StructInt::species)-1], eskernel_fluid[p.idata(StructInt::species)-1])/norm_ptr[p.idata(StructInt::species)-1];
                         }

                      }

                      //std::cout << "Adding " << p.rdata(StructReal::forcex + 0) << ", " << weight*invvol << std::endl;
                      amrex::Gpu::Atomic::Add(&fout_x(i,j,k), p.rdata(StructReal::forcex + 0)*weight*invvol);
                      amrex::Gpu::Atomic::Add(&fweights_x(i,j,k), weight);
                    }
                }
            }

            // y-components
            for (int k = lo_dim[2] ; k < hi_dim[2]; ++k) {
                for (int j = lo_dim[1] ; j < hi_dim[1]+1; ++j) {
                    for (int i = lo_dim[0] ; i < hi_dim[0]; ++i) {
                      Real pos_grid[3];
                      pos_grid[0]=(p.pos(0)-coords_y(i,j,k,0))*invdx[0];
                      pos_grid[1]=(p.pos(1)-coords_y(i,j,k,1))*invdx[1];
                      pos_grid[2]=(p.pos(2)-coords_y(i,j,k,2))*invdx[2];
                      Real pos_grid_norm = std::sqrt(pos_grid[0]*pos_grid[0]+pos_grid[1]*pos_grid[1]+pos_grid[2]*pos_grid[2]);

                      Real weight = 1.0;
                      for (int ii=0; ii<spaceDim; ++ii)
                      {
                         if(pkernel_fluid[p.idata(StructInt::species)-1] == 3)
                         {
                            weight*=Kernel3P()(pos_grid[ii]);
                         }
                         else if (pkernel_fluid[p.idata(StructInt::species)-1] == 4)
                         {
                            weight*=Kernel4P()(pos_grid[ii]);
                         }
                         else if (pkernel_fluid[p.idata(StructInt::species)-1] == 1)
                         {
                            weight*=Kernel1P()(pos_grid[ii]);
                         }
                         else if (pkernel_fluid[p.idata(StructInt::species)-1] == 6)
                         {
                            weight*=Kernel6P()(pos_grid[ii]);
                         }
                         else if (eskernel_fluid[p.idata(StructInt::species)-1] > 0)
                         {
                            weight*=KernelES()(pos_grid[ii],eskernel_beta[p.idata(StructInt::species)-1], eskernel_fluid[p.idata(StructInt::species)-1])/norm_ptr[p.idata(StructInt::species)-1];
                         }
                      }



                      amrex::Gpu::Atomic::Add(&fout_y(i,j,k), p.rdata(StructReal::forcex + 1)*weight*invvol);
                      amrex::Gpu::Atomic::Add(&fweights_y(i,j,k), weight);
                    }
                }
            }

            // z-components
            for (int k = lo_dim[2] ; k < hi_dim[2]+1; ++k) {
                for (int j = lo_dim[1] ; j < hi_dim[1]; ++j) {
                    for (int i = lo_dim[0] ; i < hi_dim[0]; ++i) {
                      Real pos_grid[3];
                      pos_grid[0]=(p.pos(0)-coords_z(i,j,k,0))*invdx[0];
                      pos_grid[1]=(p.pos(1)-coords_z(i,j,k,1))*invdx[1];
                      pos_grid[2]=(p.pos(2)-coords_z(i,j,k,2))*invdx[2];
                      Real pos_grid_norm = std::sqrt(pos_grid[0]*pos_grid[0]+pos_grid[1]*pos_grid[1]+pos_grid[2]*pos_grid[2]);

                      Real weight = 1.0;
                      for (int ii=0; ii<spaceDim; ++ii)
                      {
                         if(pkernel_fluid[p.idata(StructInt::species)-1] == 3)
                         {
                            weight*=Kernel3P()(pos_grid[ii]);
                         }
                         else if (pkernel_fluid[p.idata(StructInt::species)-1] == 4)
                         {
                            weight*=Kernel4P()(pos_grid[ii]);
                         }
                         else if (pkernel_fluid[p.idata(StructInt::species)-1] == 1)
                         {
                            weight*=Kernel1P()(pos_grid[ii]);
                         }
                         else if (pkernel_fluid[p.idata(StructInt::species)-1] == 6)
                         {
                            weight*=Kernel6P()(pos_grid[ii]);
                         }
                         else if (eskernel_fluid[p.idata(StructInt::species)-1] > 0)
                         {
                            weight*=KernelES()(pos_grid[ii],eskernel_beta[p.idata(StructInt::species)-1], eskernel_fluid[p.idata(StructInt::species)-1])/norm_ptr[p.idata(StructInt::species)-1];
                         }
                      }



                      amrex::Gpu::Atomic::Add(&fout_z(i,j,k), p.rdata(StructReal::forcex + 2)*weight*invvol);
                      amrex::Gpu::Atomic::Add(&fweights_z(i,j,k), weight);
                    }
                }
            }
        }
    });
}

template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::SpreadMarkersGpu(
            int lev,
            std::array<MultiFab, AMREX_SPACEDIM> & f_out,
            const std::array<MultiFab, AMREX_SPACEDIM> & coords,
            const Real * dx, int ghost)
{
    // timer for profiling
    BL_PROFILE_VAR("SpreadMarkersGpu()",SpreadMarkersGpu);

    for (MyConstIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();

        auto& particle_tile = this->GetParticles(lev).at(std::make_pair(grid_id,tile_id));
        const auto& particles = particle_tile.GetArrayOfStructs();

        Box tile_box = enclosedCells(f_out[0][pti].box());

        std::array<      FArrayBox *, AMREX_SPACEDIM> f_out_fab;
        std::array<      FArrayBox *, AMREX_SPACEDIM> f_weights_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> coords_fab;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            f_out_fab[d]     = & f_out[d][pti];
            f_weights_fab[d] = & f_weights[d][pti];
            coords_fab[d]    = & coords[d][pti];
        }

        SpreadKernelGpu(particles, tile_box, f_out_fab, f_weights_fab, coords_fab, dx, & ghost);
    }
}

template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::InterpolateMarkersGpu(
            int lev,
            const Real * dx,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
            const std::array<MultiFab, AMREX_SPACEDIM> & coords,
            Real& checkR)
{
    // timer for profiling
    BL_PROFILE_VAR("InterpolateMarkersGpu()",InterpolateMarkersGpu);

    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_in[d].boxArray(), f_in[d].DistributionMap(),
                            1, f_in[d].nGrow());
        f_weights[d].setVal(-1.);
    }
    int rejected_proc = 0;
    int check = 0;

    for (MyConstIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        const int grid_id = pti.index();
        const int tile_id = pti.LocalTileIndex();

        auto& particle_tile = this->GetParticles(lev).at(std::make_pair(grid_id,tile_id));
        auto& particles = particle_tile.GetArrayOfStructs();

        Box tile_box = enclosedCells(f_in[0][pti].box());

        std::array<const FArrayBox *, AMREX_SPACEDIM> f_in_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> f_weights_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> coords_fab;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            f_in_fab[d]     = & f_in[d][pti];
            f_weights_fab[d] = & f_weights[d][pti];
            coords_fab[d]    = & coords[d][pti];
        }

        //Print() << "Here1!\n";
        int gs = 0;
        InterpolateKernelGpu(particles, tile_box, f_in_fab, f_weights_fab, coords_fab, dx, &gs, check);

        //Print() << "Here2!" << check << "\n";
        rejected_proc += check;
    }

    checkR = rejected_proc;
}

template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::InterpolateMarkersGpu(
            int lev,
            const Real * dx,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
            const std::array<MultiFab, AMREX_SPACEDIM> & coords)
{

    Real check;
    InterpolateMarkersGpu(lev, dx, f_in, coords,check);

    if (check)
    {
       amrex::Print() << "Particles out of bound in InterpolateKernelGpu" << std::endl;
    }
}

template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::SpreadMarkers(
            int lev,
                  std::array<MultiFab, AMREX_SPACEDIM> & f_out,
            const std::array<MultiFab, AMREX_SPACEDIM> & coords,
            const Real * dx, int ghost
        ) const {


    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_out[d].boxArray(), f_out[d].DistributionMap(),
                            1, f_out[d].nGrow());
        f_weights[d].setVal(0.);
    }


    // for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {
    for (MyConstIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        // Marker (non-neighbor particle) data for current tile
        TileIndex index(pti.index(), pti.LocalTileIndex());
        const auto & particle_data = this->GetParticles(lev).at(index);


        //_______________________________________________________________________
        // Fill vector of marker positions and forces (for current level)
        long np = this->GetParticles(lev).at(index).numParticles();

        Vector<RealVect> marker_positions(np), marker_forces(np);

        const AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            const ParticleType & part = particles[i];

            RealVect ppos, pfor;
            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                ppos[d] = part.pos(d);
                pfor[d] = part.rdata(StructReal::forcex + d);
            }

            marker_positions[i] = ppos;
            marker_forces[i]    = pfor;
        }

        //_______________________________________________________________________
        // Spread the non-neighbor particles (markers)

        // Convert stagered boxes into cell centered
        Box tile_box = enclosedCells(f_out[0][pti].box());

        std::array<      FArrayBox *, AMREX_SPACEDIM> f_out_fab;
        std::array<      FArrayBox *, AMREX_SPACEDIM> f_weights_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> coords_fab;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            f_out_fab[d]     = & f_out[d][pti];
            f_weights_fab[d] = & f_weights[d][pti];
            coords_fab[d]    = & coords[d][pti];
        }

        SpreadMarkers(marker_forces, marker_positions, tile_box,
                      f_out_fab, f_weights_fab, coords_fab, dx, ghost);
    }
}



template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::SpreadMarkers(
            int lev,
            std::array<MultiFab, AMREX_SPACEDIM> & f_out
        ) const {

    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = this->Geom(0);
    const Real     *   dx = geom.CellSize();

    SpreadMarkers(lev, f_out, face_coords[lev], dx, 0);
}



template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::SpreadPredictor(
            int lev,
            std::array<MultiFab, AMREX_SPACEDIM> & f_out
        ) const {


    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = this->Geom(0);
    const Real     *   dx = geom.CellSize();

    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_out[d].boxArray(), f_out[d].DistributionMap(),
                            1, f_out[d].nGrow());
        f_weights[d].setVal(0.);
    }


    // for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {
    for (MyConstIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        // Marker (non-neighbor particle) data for current tile
        TileIndex index(pti.index(), pti.LocalTileIndex());
        const auto & particle_data = this->GetParticles(lev).at(index);

        //_______________________________________________________________________
        // Fill vector of marker positions and forces (for current level)
        long np = this->GetParticles(lev).at(index).numParticles();

        Vector<RealVect> marker_positions(np), marker_forces(np);

        const AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            const ParticleType & part = particles[i];

            RealVect ppos, pfor;

            for (int d=0; d<AMREX_SPACEDIM; ++d) {
                ppos[d] = part.pos(d) + part.rdata(StructReal::pred_posx + d);
                pfor[d] = part.rdata(StructReal::pred_forcex + d);
            }

            marker_positions[i] = ppos;
            marker_forces[i]    = pfor;
        }

        //_______________________________________________________________________
        // Spread the non-neighbor particles (markers)

        // Convert stagered boxes into cell centered
        Box tile_box = enclosedCells(f_out[0][pti].box());

        std::array<      FArrayBox *, AMREX_SPACEDIM> f_out_fab;
        std::array<      FArrayBox *, AMREX_SPACEDIM> f_weights_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> coords_fab;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            f_out_fab[d]     = & f_out[d][pti];
            f_weights_fab[d] = & f_weights[d][pti];
            coords_fab[d]    = & face_coords[lev][d][pti];
        }

        SpreadMarkers(marker_forces, marker_positions, tile_box,
                      f_out_fab, f_weights_fab, coords_fab, dx, 0);
    }
}



template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::InterpolateMarkers(
                  Vector<RealVect> & f_out,
            const Vector<RealVect> & f_pos,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_weights,
            const std::array<MultiFab, AMREX_SPACEDIM> & coords,
            const Real * dx, const int ghost
        ) {



    for (MFIter mfi(f_in[0]); mfi.isValid(); ++mfi) {

        // Convert stagered boxes into cell centered
        Box tile_box = enclosedCells(mfi.tilebox());

        std::array<const FArrayBox *, AMREX_SPACEDIM> f_in_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> f_weights_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> coords_fab;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            f_in_fab[d]      = & f_in[d][mfi];
            f_weights_fab[d] = & f_weights[d][mfi];
            coords_fab[d]    = & coords[d][mfi];
        }

        InterpolateMarkers(f_out, f_pos, tile_box,
                           f_in_fab, f_weights_fab, coords_fab, dx, ghost);
    }
}



template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::InterpolateMarkers(
                  Vector<RealVect> & f_out,
            const Vector<RealVect> & f_pos,
            const Box & tile_box,
            const std::array<const FArrayBox *, AMREX_SPACEDIM> & f_in,
            const std::array<const FArrayBox *, AMREX_SPACEDIM> & f_weights,
            const std::array<const FArrayBox *, AMREX_SPACEDIM> & coords,
            const Real * dx, const int ghost)
{
    // timer for profiling
    BL_PROFILE_VAR("InterpolateMarkers()",InterpolateMarkers);

    const int n_marker = f_pos.size();

    // Convert stagered boxes into cell centered
    Box bx = enclosedCells(f_in[0]->box());

    Vector<int> pkernel_fluid_in(1);
    for (int i=0; i<1; ++i) {
        pkernel_fluid_in[i] = pkernel_fluid[i];
    }

    interpolate_markers(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_BOX(tile_box),
                       BL_TO_FORTRAN_ANYD(* f_in[0]),
#if (AMREX_SPACEDIM > 1)
                       BL_TO_FORTRAN_ANYD(* f_in[1]),
#endif
#if (AMREX_SPACEDIM > 2)
                       BL_TO_FORTRAN_ANYD(* f_in[2]),
#endif
                       BL_TO_FORTRAN_ANYD(* f_weights[0]),
#if (AMREX_SPACEDIM > 1)
                       BL_TO_FORTRAN_ANYD(* f_weights[1]),
#endif
#if (AMREX_SPACEDIM > 2)
                       BL_TO_FORTRAN_ANYD(* f_weights[2]),
#endif
                       BL_TO_FORTRAN_ANYD(* coords[0]),
#if (AMREX_SPACEDIM > 1)
                       BL_TO_FORTRAN_ANYD(* coords[1]),
#endif
#if (AMREX_SPACEDIM > 2)
                       BL_TO_FORTRAN_ANYD(* coords[2]),
#endif
                       f_pos.dataPtr(),
                       f_out.dataPtr(),
                       & n_marker,
                       dx,
                       pkernel_fluid_in.dataPtr() );

}

template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::InterpolateKernelGpu(AoS& aos,
         const Box& bx,
         const std::array<const FArrayBox *, AMREX_SPACEDIM> & f_in,
         const std::array<const FArrayBox *, AMREX_SPACEDIM> & f_weights,
         const std::array<const FArrayBox *, AMREX_SPACEDIM> & coords,
         const Real* dx,
         const int* nghost,
         int& check)
{
    // timer for profiling
    BL_PROFILE_VAR("InterpolateKernelGpu()",InterpolateKernelGpu);

    GpuArray<Real, AMREX_SPACEDIM> invdx;
    //Gpu::ManagedDeviceVector<int> checkg(1, check);

    Real invvol = 1.0;

    for (int i=0; i<AMREX_SPACEDIM; ++i)
        invdx[i]=1.0/dx[i];

    for (int i=0; i<AMREX_SPACEDIM; ++i)
        invvol *= invdx[i];

    const int spaceDim = AMREX_SPACEDIM;

    GpuArray<int, 3> bx_lo = {bx.loVect()[0], bx.loVect()[1], bx.loVect()[2]};
    GpuArray<int, 3> bx_hi = {bx.hiVect()[0], bx.hiVect()[1], bx.hiVect()[2]};
    Array4<const Real> const& fin_x = f_in[0]->array();
    Array4<const Real> const& fin_y = f_in[1]->array();
    Array4<const Real> const& fin_z = f_in[2]->array();
    Array4<const Real> const& fweights_x = f_weights[0]->array();
    Array4<const Real> const& fweights_y = f_weights[1]->array();
    Array4<const Real> const& fweights_z = f_weights[2]->array();
    Array4<const Real> const& coords_x = coords[0]->array();
    Array4<const Real> const& coords_y = coords[1]->array();
    Array4<const Real> const& coords_z = coords[2]->array();

    const auto Np = aos.numParticles();
    const auto pstruct = aos().dataPtr();
    int ng = nghost[0];

    check = 0;
    int vcheck = 0;

    Gpu::DeviceScalar<Real> check_gpu(check);
    Real* pcheck = check_gpu.dataPtr();

    Gpu::DeviceScalar<Real> vcheck_gpu(vcheck);
    Real* pvcheck = vcheck_gpu.dataPtr();

    Real* norm_ptr = norm_es.data();

    AMREX_FOR_1D( Np, ip,
    {
        ParticleType& p = pstruct[ip];
        //Real norm_es = integral(eskernel_beta[p.idata(StructInt::species)-1], eskernel_fluid[p.idata(StructInt::species)-1]);

        int vis = 1;

        if(HAS_VISIBLE)
        {
            //Print() << "Has visible, "<< p.idata(StructInt::visible) << std::endl;
            if(p.idata(StructInt::visible) != 1)
            {
                vis = 0;
            }
        }

        int lo_dim[3];
        int hi_dim[3];
        int gs;
        int checkg = 0;

        if(pkernel_fluid[p.idata(StructInt::species)-1] == 3)
        {
            gs = 2;
        }
        else if (pkernel_fluid[p.idata(StructInt::species)-1] == 4)
        {
            gs = 3;
        }
        else if (pkernel_fluid[p.idata(StructInt::species)-1] == 1)
        {
            gs = 1;
        }
        else if (pkernel_fluid[p.idata(StructInt::species)-1] == 6)
        {
            gs = 4;
        }
        else if (eskernel_fluid[p.idata(StructInt::species)-1] > 0)
        {
            gs = static_cast<int>(floor(eskernel_fluid[p.idata(StructInt::species)-1]*0.5)+1);
        }
        //Print() << "gs for interpolating is " << gs << "\n";

//        if (ng == 0) //this condition isn't in the Fortan version? -DRL
//        {
//            int loc = static_cast<int>(p.pos(0) * invdx[0] - gs);
//            lo_dim[0] = (loc >= bx_lo[0]) ? loc : bx_lo[0];
//            if (loc<bx_lo[0]) checkg[0] = 1;

//            loc = static_cast<int>(p.pos(0) * invdx[0] + gs);
//            hi_dim[0] = (loc <= bx_hi[0]) ? loc : bx_hi[0];
//            if (loc>bx_hi[0]) checkg[0] = 1;

//            loc = static_cast<int>(p.pos(1) * invdx[1] - gs);
//            lo_dim[1] = (loc >= bx_lo[1]) ? loc : bx_lo[1];
//            if (loc<bx_lo[1]) checkg[0] = 1;

//            loc = static_cast<int>(p.pos(1) * invdx[1] + gs);
//            hi_dim[1] = (loc <= bx_hi[1]) ? loc : bx_hi[1];
//            if (loc>bx_hi[1]) checkg[0] = 1;

//            loc = static_cast<int>(p.pos(2) * invdx[2] - gs);
//            lo_dim[2] = (loc >= bx_lo[2]) ? loc : bx_lo[2];
//            if (loc<bx_lo[2]) checkg[0] = 1;

//            loc = static_cast<int>(p.pos(2) * invdx[2] + gs);
//            hi_dim[2] = (loc <= bx_hi[2]) ? loc : bx_hi[2];
//            if (loc>bx_hi[2]) checkg[0] = 1;
//        }
//        else

        checkg = 0;

        {
            lo_dim[0] = static_cast<int>(p.pos(0) * invdx[0] - gs);
            if (lo_dim[0]<bx_lo[0]) checkg = 1;
            hi_dim[0] = static_cast<int>(p.pos(0) * invdx[0] + gs);
            if (hi_dim[0]>bx_hi[0]) checkg = 1;
            lo_dim[1] = static_cast<int>(p.pos(1) * invdx[1] - gs);
            if (lo_dim[1]<bx_lo[1]) checkg = 1;
            hi_dim[1] = static_cast<int>(p.pos(1) * invdx[1] + gs);
            if (hi_dim[1]>bx_hi[1]) checkg = 1;
            lo_dim[2] = static_cast<int>(p.pos(2) * invdx[2] - gs);
            if (lo_dim[2]<bx_lo[2]) checkg = 1;
            hi_dim[2] = static_cast<int>(p.pos(2) * invdx[2] + gs);
            if (hi_dim[2]>bx_hi[2]) checkg = 1;
        }

        if(checkg == 0 && vis == 1)
        //if(true)
        {

        p.rdata(StructReal::velx + 0) = 0;
        p.rdata(StructReal::velx + 1) = 0;
        p.rdata(StructReal::velx + 2) = 0;
        //std::cout << "SUCCESS\n";
        double wSum = 0;

        //std::cout << "normalization factor " << 1 << " are " << norm_es[1] << std::endl;
        // x-components
        for (int k = lo_dim[2] ; k < hi_dim[2]; ++k) {
            for (int j = lo_dim[1] ; j < hi_dim[1]; ++j) {
                for (int i = lo_dim[0] ; i < hi_dim[0]+1; ++i) {

                  Real pos_grid[3];
                  pos_grid[0]=(p.pos(0)-coords_x(i,j,k,0))*invdx[0];
                  pos_grid[1]=(p.pos(1)-coords_x(i,j,k,1))*invdx[1];
                  pos_grid[2]=(p.pos(2)-coords_x(i,j,k,2))*invdx[2];
                  Real pos_grid_norm = std::sqrt(pos_grid[0]*pos_grid[0]+pos_grid[1]*pos_grid[1]+pos_grid[2]*pos_grid[2]);

                  Real weight = 1.0;
                  for (int ii=0; ii<spaceDim; ++ii)
                  {
                     if(pkernel_fluid[p.idata(StructInt::species)-1] == 3)
                     {
                        weight*=Kernel3P()(pos_grid[ii]);
                     }
                     else if (pkernel_fluid[p.idata(StructInt::species)-1] == 4)
                     {
                        weight*=Kernel4P()(pos_grid[ii]);
                     }
                     else if (pkernel_fluid[p.idata(StructInt::species)-1] == 1)
                     {
                        weight*=Kernel1P()(pos_grid[ii]);
                     }
                     else if (pkernel_fluid[p.idata(StructInt::species)-1] == 6)
                     {
                        weight*=Kernel6P()(pos_grid[ii]);
                     }
                     else if (eskernel_fluid[p.idata(StructInt::species)-1] > 0)
                     {
                        weight*=KernelES()(pos_grid[ii],eskernel_beta[p.idata(StructInt::species)-1], eskernel_fluid[p.idata(StructInt::species)-1])/norm_ptr[p.idata(StructInt::species)-1];
                        //Print() << "weight is " << weight << "\n";
                     }
                  }

                  Real wfrac = 1.0;
                  if (fweights_x(i,j,k) > 0)
                  {
                     wfrac = weight/fweights_x(i,j,k);
                  }
                  else
                  {
                     wfrac = 1.0;
                  }

                  wSum += weight;

                  p.rdata(StructReal::velx + 0) += fin_x(i,j,k)*wfrac*weight;

                 // Print() << i << ", " << j << ", " << k << ", " << coords_x(i,j,k,0) << ", "<< coords_x(i,j,k,1) << ", "<< coords_x(i,j,k,2) << ", " << fin_x(i,j,k) << ", " << weight << "\n";
                }
            }
        }
//        Print() << "Total weight: " << wSum << "\n";

        // y-components
        for (int k = lo_dim[2] ; k < hi_dim[2]; ++k) {
            for (int j = lo_dim[1] ; j < hi_dim[1]+1; ++j) {
                for (int i = lo_dim[0] ; i < hi_dim[0]; ++i) {
                  Real pos_grid[3];
                  pos_grid[0]=(p.pos(0)-coords_y(i,j,k,0))*invdx[0];
                  pos_grid[1]=(p.pos(1)-coords_y(i,j,k,1))*invdx[1];
                  pos_grid[2]=(p.pos(2)-coords_y(i,j,k,2))*invdx[2];
                  Real pos_grid_norm = std::sqrt(pos_grid[0]*pos_grid[0]+pos_grid[1]*pos_grid[1]+pos_grid[2]*pos_grid[2]);

                  Real weight = 1.0;
                  for (int ii=0; ii<spaceDim; ++ii)
                  {
                     if(pkernel_fluid[p.idata(StructInt::species)-1] == 3)
                     {
                        weight*=Kernel3P()(pos_grid[ii]);
                     }
                     else if (pkernel_fluid[p.idata(StructInt::species)-1] == 4)
                     {
                        weight*=Kernel4P()(pos_grid[ii]);
                     }
                     else if (pkernel_fluid[p.idata(StructInt::species)-1] == 1)
                     {
                        weight*=Kernel1P()(pos_grid[ii]);
                     }
                     else if (pkernel_fluid[p.idata(StructInt::species)-1] == 6)
                     {
                        weight*=Kernel6P()(pos_grid[ii]);
                     }
                     else if (eskernel_fluid[p.idata(StructInt::species)-1] > 0)
                     {
                        weight*=KernelES()(pos_grid[ii],eskernel_beta[p.idata(StructInt::species)-1], eskernel_fluid[p.idata(StructInt::species)-1])/norm_ptr[p.idata(StructInt::species)-1];
                     }
                  }

                  Real wfrac = 1.0;
                  if (fweights_y(i,j,k) > 0)
                  {
                     wfrac = weight/fweights_y(i,j,k);
                  }
                  else
                  {
                     wfrac = 1.0;
                  }


                  p.rdata(StructReal::velx + 1) += fin_y(i,j,k)*wfrac*weight;
                }
            }
        }

        // z-components
        for (int k = lo_dim[2] ; k < hi_dim[2]+1; ++k) {
            for (int j = lo_dim[1] ; j < hi_dim[1]; ++j) {
                for (int i = lo_dim[0] ; i < hi_dim[0]; ++i) {
                  Real pos_grid[3];
                  pos_grid[0]=(p.pos(0)-coords_z(i,j,k,0))*invdx[0];
                  pos_grid[1]=(p.pos(1)-coords_z(i,j,k,1))*invdx[1];
                  pos_grid[2]=(p.pos(2)-coords_z(i,j,k,2))*invdx[2];
                  Real pos_grid_norm = std::sqrt(pos_grid[0]*pos_grid[0]+pos_grid[1]*pos_grid[1]+pos_grid[2]*pos_grid[2]);

                  Real weight = 1.0;
                  for (int ii=0; ii<spaceDim; ++ii)
                  {
                     if(pkernel_fluid[p.idata(StructInt::species)-1] == 3)
                     {
                        weight*=Kernel3P()(pos_grid[ii]);
                     }
                     else if (pkernel_fluid[p.idata(StructInt::species)-1] == 4)
                     {
                        weight*=Kernel4P()(pos_grid[ii]);
                     }
                     else if (pkernel_fluid[p.idata(StructInt::species)-1] == 1)
                     {
                        weight*=Kernel1P()(pos_grid[ii]);
                     }
                     else if (pkernel_fluid[p.idata(StructInt::species)-1] == 6)
                     {
                        weight*=Kernel6P()(pos_grid[ii]);
                     }
                     else if (eskernel_fluid[p.idata(StructInt::species)-1] > 0)
                     {
                        weight*=KernelES()(pos_grid[ii],eskernel_beta[p.idata(StructInt::species)-1], eskernel_fluid[p.idata(StructInt::species)-1])/norm_ptr[p.idata(StructInt::species)-1];
                     }
                  }

                  Real wfrac = 1.0;
                  if (fweights_z(i,j,k) > 0)
                  {
                     wfrac = weight/fweights_z(i,j,k);
                  }
                  else
                  {
                     wfrac = 1.0;
                  }

                  p.rdata(StructReal::velx + 2) += fin_z(i,j,k)*wfrac*weight;
                }
            }
        }
        }else if(checkg == 1)
        {
           pcheck[0]++;
           //std::cout << "Fail\n";
        }else
        {
           pvcheck[0]++;
           //Print() << "Ignoring invisible!\n";
        }
    });

    check = check_gpu.dataValue();
    vcheck = vcheck_gpu.dataValue();
    //std::cout << "HERE! " << check << std::endl;

}

template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::InterpolateMarkers(
            int lev,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_in
        ) {


    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = this->Geom(0);
    const Real     *   dx = geom.CellSize();

    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_in[d].boxArray(), f_in[d].DistributionMap(),
                            1, f_in[d].nGrow());
        f_weights[d].setVal(-1.);
    }


    // for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {
    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        // Marker (non-neighbor particle) data for current tile
        TileIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = this->GetParticles(lev).at(index);

        //_______________________________________________________________________
        // Fill vector of marker positions and forces (for current level)
        long np = this->GetParticles(lev).at(index).numParticles();

        Vector<RealVect> marker_positions(np), marker_forces(np);

        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            const ParticleType & part = particles[i];

            RealVect ppos, pfor;
            for (int d=0; d<AMREX_SPACEDIM; ++d)
                ppos[d] = part.pos(d);

            // pfor should be (0, .., 0)

            marker_positions[i] = ppos;
            marker_forces[i]    = pfor;
        }


        //_______________________________________________________________________
        // Interpolate the non-neighbor particles (markers)

        // Convert stagered boxes into cell centered
        Box tile_box = enclosedCells(f_in[0][pti].box());

        std::array<const FArrayBox *, AMREX_SPACEDIM> f_in_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> f_weights_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> coords_fab;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            f_in_fab[d]      = & f_in[d][pti];
            f_weights_fab[d] = & f_weights[d][pti];
            coords_fab[d]    = & face_coords[lev][d][pti];
        }

        InterpolateMarkers(marker_forces, marker_positions, tile_box,
                           f_in_fab, f_weights_fab, coords_fab, dx, 0);

        // Add interpolated markers back to the particles (markers)
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];

            for (int d=0; d<AMREX_SPACEDIM; ++d)
                part.rdata(StructReal::velx + d) += marker_forces[i][d];
        }
    }
}

template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::InterpolateMarkers(
            int lev,
            const Real * dx,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
            const std::array<MultiFab, AMREX_SPACEDIM> & coords
        ) {


    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = this->Geom(0);
//    const Real     *   dx = geom.CellSize();

    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_in[d].boxArray(), f_in[d].DistributionMap(),
                            1, f_in[d].nGrow());
        f_weights[d].setVal(-1.);
    }


    // for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {
    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        // Marker (non-neighbor particle) data for current tile
        TileIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = this->GetParticles(lev).at(index);

        //_______________________________________________________________________
        // Fill vector of marker positions and vels (for current level)
        long np = this->GetParticles(lev).at(index).numParticles();

        Vector<RealVect> marker_positions(np), marker_vels(np);

        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            const ParticleType & part = particles[i];

            RealVect ppos, pfor;
            for (int d=0; d<AMREX_SPACEDIM; ++d)
                ppos[d] = part.pos(d);

            // pfor should be (0, .., 0)

            marker_positions[i] = ppos;
            marker_vels[i]    = pfor;
        }


        //_______________________________________________________________________
        // Interpolate the non-neighbor particles (markers)

        // Convert stagered boxes into cell centered
        Box tile_box = enclosedCells(f_in[0][pti].box());

        std::array<const FArrayBox *, AMREX_SPACEDIM> f_in_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> f_weights_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> coords_fab;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            f_in_fab[d]      = & f_in[d][pti];
            f_weights_fab[d] = & f_weights[d][pti];
            coords_fab[d]    = & coords[d][pti];
        }

        InterpolateMarkers(marker_vels, marker_positions, tile_box,
                           f_in_fab, f_weights_fab, coords_fab, dx, 0);

        // Add interpolated markers back to the particles (markers)
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];

            for (int d=0; d<AMREX_SPACEDIM; ++d)
                part.rdata(StructReal::velx + d) = marker_vels[i][d];
        }
    }
}



template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::InterpolatePredictor(
            int lev,
            const std::array<MultiFab, AMREX_SPACEDIM> & f_in
        ) {


    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = this->Geom(0);
    const Real     *   dx = geom.CellSize();

    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_in[d].boxArray(), f_in[d].DistributionMap(),
                            1, f_in[d].nGrow());
        f_weights[d].setVal(-1.);
    }


    // for (MFIter pti = MakeMFIter(lev, false); pti.isValid(); ++pti) {
    for (MyIBMarIter pti(* this, lev); pti.isValid(); ++pti) {

        // Marker (non-neighbor particle) data for current tile
        TileIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = this->GetParticles(lev).at(index);

        //_______________________________________________________________________
        // Fill vector of marker positions and forces (for current level)
        long np = this->GetParticles(lev).at(index).numParticles();

        Vector<RealVect> marker_positions(np), marker_forces(np);

        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            const ParticleType & part = particles[i];

            RealVect ppos, pfor;

            for (int d=0; d<AMREX_SPACEDIM; ++d)
                ppos[d] = part.pos(d) + part.rdata(StructReal::pred_posx + d);

            // pfor should be (0, .., 0)

            marker_positions[i] = ppos;
            marker_forces[i]    = pfor;
        }


        //_______________________________________________________________________
        // Interpolate the non-neighbor particles (markers)

        // Convert stagered boxes into cell centered
        Box tile_box = enclosedCells(f_in[0][pti].box());

        std::array<const FArrayBox *, AMREX_SPACEDIM> f_in_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> f_weights_fab;
        std::array<const FArrayBox *, AMREX_SPACEDIM> coords_fab;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            f_in_fab[d]      = & f_in[d][pti];
            f_weights_fab[d] = & f_weights[d][pti];
            coords_fab[d]    = & face_coords[lev][d][pti];
        }

        InterpolateMarkers(marker_forces, marker_positions, tile_box,
                           f_in_fab, f_weights_fab, coords_fab, dx, 0);

        // Add interpolated markers back to the particles (markers)
        for (int i=0; i<np; ++i) {
            ParticleType & part = particles[i];

            for (int d=0; d<AMREX_SPACEDIM; ++d)
                part.rdata(StructReal::pred_velx + d) += marker_forces[i][d];
        }
    }
}



template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::PrintMarkerData(int lev) const {

    // Inverse cell-size vector => max is used for determining IBParticle
    // radius in units of cell size
    Vector<Real> inv_dx = {AMREX_D_DECL(this->Geom(lev).InvCellSize(0),
                                        this->Geom(lev).InvCellSize(1),
                                        this->Geom(lev).InvCellSize(2)   )};

    // Find max inv_dx (in case we have an anisotropic grid)
    Real mx_inv_dx = * std::max_element(inv_dx.begin(), inv_dx.end());


    amrex::AllPrintToFile("ib_marker_data") << "Particles on each box:" << std::endl;


    long local_count = 0;

    // ParIter skips tiles without particles => Iterate over MultiFab instead
    // of ParticleIter. Note also that AmrexParticleContainer uses strange
    // tiling => don't turn off tiling (particles are stored in tile)
    for(MFIter pti = this->MakeMFIter(lev, this->do_tiling); pti.isValid(); ++pti) {
        // MuliFabs are indexed using a pair: (BoxArray index, tile index):
        TileIndex index(pti.index(), pti.LocalTileIndex());

        // Neighbours are stored as raw data (see below)
        int ng = this->neighbors[lev].at(index).size();

        auto & particle_data = this->GetParticles(lev).at(index);
        long np = this->GetParticles(lev).at(index).numParticles();

        local_count += np;

        // Print current box info
        AllPrintToFile("ib_marker_data") << "Box:"         << pti.index()
                                         << " "            << pti.tilebox()
                                         << ", count: "    << np
                                         << ", nb count: " << ng
                                         << std::endl;

        // Print IBMarker
        AllPrintToFile("ib_marker_data") << " * IBMarkers:" << std::endl;

        const AoS & particles = particle_data.GetArrayOfStructs();
        for(int i = 0; i < np; i++){
            const ParticleType & part = particles[i];

            AllPrintToFile("ib_marker_data") << "   +--> " << part << std::endl;
        }

        // Print neighbour IBParticles
        AllPrintToFile("ib_marker_data") << " * Grown IBParticles:" << std::endl;

        // TODO: HAXOR!!! This should be fixed ASAP: if I understand this
        // correctly, the neighbor data contains the particle data as a binary
        // array (char). By casting to ParticleType, what we're doing is
        // interpreting the data in neighbours[index] as valid particle data.
        // Also we stride the neighbors[index] array in units of
        // sizeof(ParticleData). All of this is a little too dangerous for my
        // taste: never hide what you're doing from your compiler!!!
        const ParticleType * nbhd_data = (ParticleType *) this->neighbors[lev].at(index).GetArrayOfStructs().dataPtr();
        for(int i = 0; i < ng; i++){
            const ParticleType & part = nbhd_data[i];

            AllPrintToFile("ib_marker_data") << "   +--> " << part << std::endl;
        }
    }

    AllPrintToFile("ib_marker_data") << "Total for this process: "
                                       << local_count << std::endl << std::endl;
}



template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::InitInternals(int ngrow) {

    ReadStaticParameters();

    this->SetVerbose(0);

    // Turn off certain components for ghost particle communication
    // Field numbers: {0, 1, 2} => {x, y, z} particle coordinates
    //      => 3 corresponds to the start of IBM_realData
    for (int i=3; i < StructReal::count + 3; ++i)
        this->setRealCommComp(i,  true);

    // Field numbers: {0, 1} => {ID, CPU}
    //      => 2 corresponds to the start of IBM_intData
    // We _do_ want the the neighbour particles to have ID and cpu init data.
    for (int i = 2; i < StructInt::count + 2; ++i)
        this->setIntCommComp(i, true);


    //for (int d=0; d<AMREX_SPACEDIM; ++d) {
    //    this->setRealCommComp(StructReal::pred_forcex + d, false);
    //    this->setRealCommComp(StructReal::forcex + d, false);
    //}

    this->setEnableInverse(true);


    /****************************************************************************
     *                                                                          *
     * Fill auxiallry data used by interpolation                                *
     *   -> face_coords: the face-centered coordinates used by the fluid grids  *
     *                                                                          *
     ***************************************************************************/

    // Number of particle IDs for each rank:
    num_ids.resize(ParallelDescriptor::NProcs(), 0);
    // Offset of each rank in the particle list
    cpu_offset.resize(ParallelDescriptor::NProcs(), 0);

    // TODO: this is only assuming 1 fluid level (level 0)
    int lev = 0;

    face_coords.resize(lev + 1);
    const BoxArray & ba            = this->ParticleBoxArray(lev);
    const DistributionMapping & dm = this->ParticleDistributionMap(lev);

    // Determine the minimal ngrow (which depends on the kernel order): This
    // needs to be AT LEAST 1 larger than the `gs` value used in
    // [spread,interpolate]_kernel
    int ngc_min_fluid = * (
        std::max_element(
            common::pkernel_fluid.begin(),
            common::pkernel_fluid.begin() + common::nspecies
        )
    ) == 3 ? 3 : 5; // GALEN - FLUID KERNEL
    int ngc_min_es = * (
        std::max_element(
            common::pkernel_es.begin(),
            common::pkernel_es.begin() + common::nspecies
        )
    ) == 3 ? 3 : 5;

    // Let's be on the safe side and use the larger out of the
    // fluid/electro-static kernels
    int ngrow_min    = std::max(ngc_min_es, ngc_min_fluid);
    int ngrow_coords = std::max(ngrow, ngrow_min);

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        const BoxArray ba_fc = convert(ba, nodal_flag_dir[d]);
        face_coords[lev][d].define(ba_fc, dm, AMREX_SPACEDIM, ngrow_coords);
    }

    const Geometry & geom = this->Geom(lev);
    FindFaceCoords(face_coords[lev], geom);

    norm_es.resize(common::nspecies);
    for (int i=0; i<common::nspecies; ++i) {
        norm_es[i] = integrate_es(
            common::eskernel_beta[i], common::eskernel_fluid[i]
        );
    }
}


template <typename StructReal, typename StructInt>
void IBMarkerContainerBase<StructReal, StructInt>::ReadStaticParameters() {
    static bool initialized = false;

    if (!initialized) {
        ParmParse pp("particles");

//#ifndef AMREX_USE_GPU
//        // AMReX default is false => enable by default
//        this->do_tiling = true;
//        // Allow user to overwrite
//        pp.query("do_tiling",  this->do_tiling);
//#endif

        // If tiling is enabled, make sure that the tile size is at least the
        // number of ghost cells (otherwise strange things happen)
        if (this->do_tiling)
            this->tile_size = IntVect{AMREX_D_DECL(max_grid_size[0],
                                                   max_grid_size[1],
                                                   max_grid_size[2])};
            //this->tile_size = IntVect{AMREX_D_DECL(nghost, nghost, nghost)};

        // User can overwrite
        Vector<int> ts(BL_SPACEDIM);
        if (pp.queryarr("tile_size", ts))
            this->tile_size = IntVect(ts);

        // pp.query("use_neighbor_list", use_neighbor_list);
        // pp.query("sort_neighbor_list", sort_neighbor_list);

        initialized = true;
    }
}


template <typename StructReal, typename StructInt>
Real IBMarkerContainerBase<StructReal, StructInt>::integrate_es(Real beta_in, int w_in) {
    Real beta = beta_in;
    Real w = static_cast<Real>(w_in);
    int nbin = 100;
    Real dxi = w/nbin;
    Real integral = 0.0;

    // mid-point rule for integration
    for (int i = 0; i < nbin; i++) {
        integral += dxi*std::exp(beta*(std::sqrt(1 - 4*(-w*0.5+i*dxi+0.5*dxi)*(-w*0.5+i*dxi+0.5*dxi)/w/w) - 1));
    }
    //std::cout << "denominator is " << integral << std::endl;
    return integral;
}


#endif
