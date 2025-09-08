#include <AMReX.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <AMReX_VisMF.H>  // amrex::VisMF::Write(MultiFab)

#include <common_functions.H>

#include <IBParticleContainer.H>
#include <ib_functions_F.H>
#include <MFUtil.H>



using namespace amrex;

bool IBParticleContainer::use_neighbor_list  {true};
bool IBParticleContainer::sort_neighbor_list {false};



IBParticleContainer::IBParticleContainer(const Geometry & geom,
                                         const DistributionMapping & dmap,
                                         const BoxArray & ba,
                                         int n_nbhd)
    : NeighborParticleContainer<IBP_realData::count, IBP_intData::count>(
            geom, dmap, ba, n_nbhd
        ),
      nghost(n_nbhd)
{
    InitInternals(n_nbhd);
}



IBParticleContainer::IBParticleContainer(AmrCore * amr_core, int n_nbhd)
    : NeighborParticleContainer<IBP_realData::count, IBP_intData::count>(
            amr_core->GetParGDB(), n_nbhd
        ),
      m_amr_core(amr_core),
      nghost(n_nbhd)
{
    InitInternals(n_nbhd);
}



void IBParticleContainer::InitList(int lev,
                                   const Vector<RealVect> & pos,
                                   const Vector<Real> & r,
                                   const Vector<Real> & rho) {

    // Inverse cell-size vector => used for determining index corresponding to
    // IBParticle position (pos)
    RealVect inv_dx = RealVect(
        AMREX_D_DECL(
            Geom(lev).InvCellSize(0),
            Geom(lev).InvCellSize(1),
            Geom(lev).InvCellSize(2)
        )
    );


    int total_np = 0;

    // This uses the particle tile size. Note that the default is to tile so if
    // we remove the true and don't explicitly add false it will still tile.
    for (MFIter mfi = MakeMFIter(lev, true); mfi.isValid(); ++mfi) {

        // This is particles per grid so we reset to 0
        int pcount = 0;

        // Current tile box
        const Box & tile_box = mfi.tilebox();

        // Create a particle container for this grid and add the
        // immersed-boundary particles to it if the particle's position (pos[i])
        // is within current tile box.
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto & particles = GetParticles(lev)[std::make_pair(grid_id,tile_id)];

        for(int i = 0; i < pos.size(); i++) {
            // IntVect representing particle's position in the tile_box grid.
            RealVect pos_grid = pos[i];
            pos_grid *= inv_dx;
            IntVect pos_ind = IntVect(
                AMREX_D_DECL(
                    (int) pos_grid[0],
                    (int) pos_grid[1],
                    (int) pos_grid[2]
                )
            );

            // Add particle at position pos iff it's vector index is contained
            // within tile_box.
            if(tile_box.contains(pos_ind)) {
                pcount ++;

                ParticleType p_new;

                // Set id and cpu for this particle
                p_new.id()  = ParticleType::NextID();
                p_new.cpu() = ParallelDescriptor::MyProc();

                // Set particle position
                p_new.pos(0) = pos[i][0];
                p_new.pos(1) = pos[i][1];
                p_new.pos(2) = pos[i][2];

                // Initialize particle velocity (and angular velocity) as well
                // as drag to 0
                p_new.rdata(IBP_realData::velx)   = 0.;
                p_new.rdata(IBP_realData::vely)   = 0.;
                p_new.rdata(IBP_realData::velz)   = 0.;

                p_new.rdata(IBP_realData::omegax) = 0.;
                p_new.rdata(IBP_realData::omegay) = 0.;
                p_new.rdata(IBP_realData::omegaz) = 0.;

                p_new.rdata(IBP_realData::dragx)  = 0.;
                p_new.rdata(IBP_realData::dragy)  = 0.;
                p_new.rdata(IBP_realData::dragz)  = 0.;

                // Calculate particle volume, mass, moment of inertia from
                // radius and density
                int pstate = 1;
                Real pradius = r[i], pdensity = rho[i];
                Real pvol, pmass, pomoi, pomega;
                // TODO: this only works for spherical particles
                set_particle_properties(& pstate, & pradius, & pdensity,
                                        & pvol,   & pmass,   & pomoi,     & pomega );

                // TODO: Audit
                p_new.idata(IBP_intData::phase)     = 1;
                p_new.idata(IBP_intData::state)     = pstate;

                p_new.rdata(IBP_realData::volume)   = pvol;
                p_new.rdata(IBP_realData::density)  = pdensity;
                p_new.rdata(IBP_realData::mass)     = pmass;
                p_new.rdata(IBP_realData::oneOverI) = pomoi;
                p_new.rdata(IBP_realData::radius)   = pradius;

                // Add to the data structure
                particles.push_back(p_new);
            }
        }

        const int np = pcount;
        total_np += np;
    }

    ParallelDescriptor::ReduceIntSum(total_np,ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "Total number of generated particles: " << total_np << std::endl;

    // We shouldn't need this if the particles are tiled with one tile per
    // grid, but otherwise we do need this to move particles from tile 0 to the
    // correct tile.
    Redistribute();
}



void IBParticleContainer::InitInternals(int ngrow) {
    ReadStaticParameters();

    this->SetVerbose(0);

    // Turn off certain components for ghost particle communication
    // Field numbers: {0, 1, 2} => {x, y, z} particle coordinates
    //      => 3 corresponds to the start of IBP_realData
    setRealCommComp(4, false);   // IBP_realData.volume
    setRealCommComp(5, false);   // IBP_realData.mass
    setRealCommComp(6, false);   // IBP_realData.density
    setRealCommComp(7, false);   // IBP_realData.oneOverI
    setRealCommComp(14, false);  // IBP_realData.dragx
    setRealCommComp(15, false);  // IBP_realData.dragy
    setRealCommComp(16, false);  // IBP_realData.dragz

    // Field numbers: {0, 1} => {ID, CPU}
    //      => 2 corresponds to the start of IBP_intData
    // We _do_ want the the neighbour particles to have ID and cpu init data.
    //setIntCommComp(0, false);  // IBP_intData.phase
    //setIntCommComp(1, false);  // IBP_intData.state
    setIntCommComp(3, false);    // IBP_intData.phase


    /****************************************************************************
     *                                                                          *
     * Fill auxiallry data used by interpolsation                               *
     *   -> face_coords: the face-centered coordinates used by the fluid grids  *
     *                                                                          *
     ***************************************************************************/

    // TODO: this is only assuming 1 fluid level (level 0)
    int lev = 0;

    face_coords.resize(lev + 1);
    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        const BoxArray ba_fc = convert(ba, nodal_flag_dir[d]);
        face_coords[lev][d].define(ba_fc, dm, AMREX_SPACEDIM, ngrow);
    }

    const Geometry & geom = Geom(lev);
    FindFaceCoords(face_coords[lev], geom);
}



void IBParticleContainer::ReadStaticParameters() {
    static bool initialized = false;

    if (!initialized) {
        ParmParse pp("particles");

 //       // AMReX default is false => enable by default
 //       do_tiling = true;
 //       // Allow user to overwrite
 //       pp.query("do_tiling",  do_tiling);

        // If tiling is enabled, make sure that the tile size is at least the
        // number of ghost cells (otherwise strange things happen)
        if (do_tiling)
            tile_size = IntVect{AMREX_D_DECL(nghost, nghost, nghost)};
        // User can overwrite
        Vector<int> ts(BL_SPACEDIM);
        if (pp.queryarr("tile_size", ts))
            tile_size = IntVect(ts);

        pp.query("use_neighbor_list", use_neighbor_list);
        pp.query("sort_neighbor_list", sort_neighbor_list);

        initialized = true;
    }
}



void IBParticleContainer::FillMarkerPositions(int lev, int n_marker) {

    /****************************************************************************
     *                                                                          *
     * Collect the list of particles and neighbor particles                     *
     *                                                                          *
     ***************************************************************************/

    fillNeighbors();

    std::map<ParticleIndex, PPVOR> ib_ppvor;

    //___________________________________________________________________________
    // Add particles to list
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {
        // MuliFabs are indexed using a pair: (BoxArray index, tile index):
        PairIndex index(pti.index(), pti.LocalTileIndex());

        auto & particle_data = GetParticles(lev)[index];
        long np = GetParticles(lev)[index].numParticles();

        // Iterate over local particle data
        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];
            ParticleIndex pindex(part.id(), part.cpu());

            // Add to list if it's not already there (don't forget the
            // critical, just in case std::vector::push_back is not atomic).
            // TODO: this implementation needs work to increase performance.
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                auto part_it = ib_ppvor.find(pindex);
                if (part_it == ib_ppvor.end()) {
                    // Position of IBParticle
                    RealVect pos = RealVect(
                        AMREX_D_DECL(part.pos(0), part.pos(1), part.pos(2))
                    );

                    // Velocity of IBParticle
                    RealVect vel = RealVect(
                        AMREX_D_DECL(part.rdata(IBP_realData::velx),
                                     part.rdata(IBP_realData::vely),
                                     part.rdata(IBP_realData::velz)   )
                    );
                   // HACK orientation hard coded in
                    RealVect ori = RealVect(
                        AMREX_D_DECL(0.,0.,1. ));

                    // Radius
                    Real r = part.rdata(IBP_realData::radius);

                    PPVOR ppvor{pos, vel, ori, r};
                    ib_ppvor[pindex] = ppvor;
                }
            }
        }

        long ng = neighbors[lev][index].size();

        // Iterate over neighbor particle data
        ParticleType * nbhd_data = (ParticleType *) neighbors[lev][index].GetArrayOfStructs().dataPtr();
        for (int i = 0; i < ng; ++i) {
            ParticleType & part = nbhd_data[i];
            ParticleIndex pindex(part.id(), part.cpu());

            // Add to list if it's not already there (don't forget the
            // critical, just in case std::vector::push_back is not atomic).
            // TODO: this implementation needs work to increase performance.
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                auto part_it = ib_ppvor.find(pindex);
                if (part_it == ib_ppvor.end()) {
                    // Position of IBParticle
                    RealVect pos = RealVect(
                        AMREX_D_DECL(part.pos(0), part.pos(1), part.pos(2))
                    );

                    // Velocity of IBParticle
                    RealVect vel = RealVect(
                        AMREX_D_DECL(part.rdata(IBP_realData::velx),
                                     part.rdata(IBP_realData::vely),
                                     part.rdata(IBP_realData::velz)   )
                    );
                   // HACK orientation hard coded in

                    RealVect ori = RealVect(
                        AMREX_D_DECL(0.,0.,1. ));

                    // Radius
                    Real r = part.rdata(IBP_realData::radius);

                    PPVOR ppvor{pos, vel, ori, r};
                    ib_ppvor[pindex] = ppvor;
                }
            }
        }
    }


    /****************************************************************************
     *                                                                          *
     * Fill Markert list                                                        *
     *                                                                          *
     ***************************************************************************/

    //___________________________________________________________________________
    // Based on the paper:
    // >*Distributing many points on a sphere*, E. B. Saff, A. B. J. Kuijlaars,
    // *The Mathematical Intelligencer*, **19** (1997)

    //___________________________________________________________________________
    // Ensure that the marker lists have enough levels, and clear previous ones
    if (marker_positions.size() <= lev) {
        marker_positions.resize(lev+1);
        marker_velocities.resize(lev+1);
        marker_forces.resize(lev+1);
    }
    marker_positions[lev].clear();
    marker_velocities[lev].clear();
    marker_forces[lev].clear();


    //___________________________________________________________________________
    // Compute marker coordinates for each particle
    double inv_sqrt_n = 1./std::sqrt(n_marker);

    for (const auto & elt : ib_ppvor) {
        // elt = (particle ID, particle data)
        //        ^^ first ^^, ^^ second  ^^

        //_______________________________________________________________________
        // Create blank marker list, and access particle data
        // ... initialized to (0..0) by default constructor
        marker_positions[lev][elt.first].resize(n_marker);
        marker_velocities[lev][elt.first].resize(n_marker);
        marker_forces[lev][elt.first].resize(n_marker);


        double   r     = elt.second.rad*0.8; // HACK: put markers slightly inside
        RealVect pos_0 = elt.second.pos;

        //_______________________________________________________________________
        // Fill marker using Saff spiral
        double phi = 0.;
        for (int i=0; i<n_marker; ++i) {

            // Compute polar coordinates of marker positions
            double ck    = -1. + (2.*i)/(n_marker-1);
            double theta = std::acos(ck);

            if ( (i==0) || (i==n_marker-1) ) phi = 0;
            else phi = std::fmod(phi + 3.6*inv_sqrt_n/std::sqrt(1-ck*ck), 2*M_PI);

            // Convert to cartesian coordinates
            RealVect pos;
#if   (AMREX_SPACEDIM == 2)
            pos[0] = pos_0[0] + r*std::sin(theta);
            pos[1] = pos_0[1] + r*std::cos(theta);
#elif (AMREX_SPACEDIM == 3)
            pos[0] = pos_0[0] + r*std::sin(theta)*std::cos(phi);
            pos[1] = pos_0[1] + r*std::sin(theta)*std::sin(phi);
            pos[2] = pos_0[2] + r*std::cos(theta);
#endif

            // Add to list
            marker_positions[lev][elt.first][i] = pos;
        }
    }

}



void IBParticleContainer::SpreadMarkers(int lev, const ParticleIndex & pindex,
        const Vector<RealVect> & f_in, std::array<MultiFab, AMREX_SPACEDIM> & f_out,
        std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const {


    //___________________________________________________________________________
    // Don't do anything if pindex isn't on this rank
    auto part_it = marker_positions[lev].find(pindex);
    if (part_it == marker_positions[lev].end())
        return;


    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = Geom(0);
    const Real     *   dx = geom.CellSize();

    const int n_marker = marker_positions[lev].at(pindex).size();


    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBParticleContainer is on a differnt grid
    // than the grid which we're spreading to

    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);

    MultiFab dummy(ba, dm, 1, f_out[0].nGrow());

    int ghost = 0;

    for (MFIter mfi(dummy); mfi.isValid(); ++mfi) {

        const Box & bx       = mfi.growntilebox();
        //const Box & tile_box = mfi.tilebox();
        const Box & tile_box = mfi.growntilebox(); // HACK, use tilebox going forward
/*
        spread_markers(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_BOX(tile_box),
                       BL_TO_FORTRAN_ANYD(f_out[0][mfi]),
#if (AMREX_SPACEDIM > 1)
                       BL_TO_FORTRAN_ANYD(f_out[1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                       BL_TO_FORTRAN_ANYD(f_out[2][mfi]),
#endif
                       BL_TO_FORTRAN_ANYD(f_weights[0][mfi]),
#if (AMREX_SPACEDIM > 1)
                       BL_TO_FORTRAN_ANYD(f_weights[1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                       BL_TO_FORTRAN_ANYD(f_weights[2][mfi]),
#endif
                       BL_TO_FORTRAN_ANYD(face_coords[lev][0][mfi]),
#if (AMREX_SPACEDIM > 1)
                       BL_TO_FORTRAN_ANYD(face_coords[lev][1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                       BL_TO_FORTRAN_ANYD(face_coords[lev][2][mfi]),
#endif
                       marker_positions[lev].at(pindex).dataPtr(),
                       f_in.dataPtr(),
                       & n_marker,
                       dx, & ghost);
*/
         Real pos[AMREX_SPACEDIM];
         Real v_spread[AMREX_SPACEDIM];

         for (int i=0; i<n_marker; ++i) {
             for (int j=0; j<AMREX_SPACEDIM; ++j)
                 pos[j] = marker_positions[lev].at(pindex)[j][i];
                 if (ghost == 0){
                         if (pos[0] < tile_box.loVect()[0]*dx[0]) continue;
                         if (pos[0] >= (tile_box.hiVect()[0]+1)*dx[0]) continue;

                         if (pos[1] < tile_box.loVect()[1]*dx[1]) continue;
                         if (pos[1] >= (tile_box.hiVect()[1]+1)*dx[1]) continue;

                         if (pos[2] < tile_box.loVect()[2]*dx[2]) continue;
                         if (pos[2] >= (tile_box.hiVect()[2]+1)*dx[2]) continue;
                 }
                 for (int j=0; j<AMREX_SPACEDIM; ++j)
                         v_spread[j] = f_in[j][i];

                 SpreadKernel(bx, f_out, f_weights, face_coords, pos, v_spread, dx, & ghost, & mfi);
         }

    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real kernel_3p(Real r_in)
{
    Real kernel_3p;

    //internal variables
    Real r;
    Real r1;
    Real r2;

    //initialize r
    r = r_in;
    r1 = amrex::Math::abs(r_in);
    r2 = r1*r1;

    if (r1 <= 0.5){
       kernel_3p = (1+std::sqrt(1-(3*r2)))/3.0;
    }else if (r1 <= 1.5){
       kernel_3p = (5.0-(3.0*r1)-std::sqrt(-3*(1-r1)*(1-r1)+1))/6.0;
    }else{
       kernel_3p = 0.0;
    }
    return kernel_3p;
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real kernel_6p(Real r_in)
{
    //The 6-point kernel function, based on the paper:
    //>*A Gaussian-like immersed-boundary kernel with three continuous derivatives and improved translational invariance*
    //>Yuanxun Bao, Jason Kaye, Charles S. Peskin, *Journal of Computational Physics* **316**, 139 (2016)
    //>https://dx.doi.org/10.1016/j.jcp.2016.04.024
    //Note also: https://github.com/stochasticHydroTools/IBMethod/blob/master/IBKernels/Kernels.c because the paper above has mistakes (but the repo's OK)

    Real kernel_6p;

    //internal parameters
    Real K = 59.0/60 - std::sqrt(29.0)/20.0;
    Real sgn = (3.0/2 - K)/amrex::Math::abs(3.0/2 - K);

    //pre-computed ratios
    Real inv16 = 1.0/16.0;
    Real inv8 = 1.0/8.0;
    Real inv12 = 1.0/12.0;
    Real inv4 = 1.0/4.0;
    Real inv6 = 1.0/6.0;
    Real rat58 = 5.0/8.0;

    //internal variables
    Real r = r_in;

    //compute kernel function
    if (r <= -3){
       kernel_6p = 0.0;
    }else if (r <= -2){
       r += 3;
       kernel_6p = phi1(r);
    }else if (r <= -1){
       r += 2;
       kernel_6p = -3*phi1(r) - inv16 + inv8*(K+(r*r)) + inv12*((3*K-1)*r+(r*r*r));
    }else if (r <= 0){
       r += 1;
       kernel_6p = 2*phi1(r) + inv4 + inv6*((4-3*K)*r-(r*r*r));
    }else if (r <= 1){
       kernel_6p = 2*phi1(r) + rat58 - inv4*(K+(r*r));
    }else if (r <= 2){
       r -= 1;
       kernel_6p = -3*phi1(r) + inv4 - inv6*((4-3*K)*r-(r*r*r));
    }else if (r <= 3){
       r -= 2;
       kernel_6p = phi1(r) - inv16 + inv8*(K+(r*r)) - inv12*((3*K-1)*r+(r*r*r));
    }else{
       kernel_6p = 0.0;
    }

    return kernel_6p;
}


void IBParticleContainer::SpreadKernel(const Box& bx, std::array<MultiFab, AMREX_SPACEDIM> & f_out, std::array<MultiFab, AMREX_SPACEDIM> & f_weights, const Vector<std::array<MultiFab, AMREX_SPACEDIM>> & face_coords, const Real* pos, const Real* v_spread, const Real* dx, int* nghost, MFIter* mfi) const {

    int i, j, k, ilo, ihi, jlo, jhi, klo, khi, gs;
    Real (*kernel_ptr) (Real );
    if(pkernel_fluid[0] == 3) {
        kernel_ptr = &kernel_3p;
        gs = 2;
    }
    else {
        kernel_ptr = &kernel_6p;
        gs = 4;
    }

    Real invdx[AMREX_SPACEDIM];
    Real invvol = 1.0;

    for (int i=0; i<AMREX_SPACEDIM; ++i)
        invdx[i]=1.0/dx[i];

    for (int i=0; i<AMREX_SPACEDIM; ++i)
        invvol *= invdx[i];

    if (nghost == 0){
        int loc = amrex::Math::floor(pos[0] * invdx[0] - gs);
        ilo = amrex::max(bx.loVect()[0], loc);
        loc = amrex::Math::floor(pos[0] * invdx[0] + gs);
        ihi = amrex::max(bx.hiVect()[0], loc);
        loc = amrex::Math::floor(pos[1] * invdx[1] - gs);
        jlo = amrex::max(bx.loVect()[1], loc);
        loc = amrex::Math::floor(pos[1] * invdx[1] + gs);
        jhi = amrex::max(bx.hiVect()[1], loc);
        loc = amrex::Math::floor(pos[2] * invdx[2] - gs);
        klo = amrex::max(bx.loVect()[2], loc);
        loc = amrex::Math::floor(pos[2] * invdx[2] + gs);
        khi = amrex::max(bx.hiVect()[2], loc);
    }else{
        ilo = amrex::Math::floor(pos[0] * invdx[0] - gs);
        ihi = amrex::Math::floor(pos[0] * invdx[0] + gs);
        jlo = amrex::Math::floor(pos[1] * invdx[1] - gs);
        jhi = amrex::Math::floor(pos[1] * invdx[1] + gs);
        klo = amrex::Math::floor(pos[2] * invdx[2] - gs);
        khi = amrex::Math::floor(pos[2] * invdx[2] + gs);
    }

    IntVect scalx_lo(ilo,jlo,klo);
    IntVect scalx_hi(ihi+1,jhi,khi);
    const Box bx_x(scalx_lo, scalx_hi);
    IntVect scaly_lo(ilo,jlo,klo);
    IntVect scaly_hi(ihi,jhi+1,khi);
    const Box bx_y(scaly_lo, scaly_hi);
    IntVect scalz_lo(ilo,jlo,klo);
    IntVect scalz_hi(ihi,jhi,khi+1);
    const Box bx_z(scalz_lo, scalz_hi);


    Array4<Real> const& fout_x = f_out[0].array(*mfi);
    Array4<Real> const& fout_y = f_out[1].array(*mfi);
    Array4<Real> const& fout_z = f_out[2].array(*mfi);
    Array4<Real> const& fweights_x = f_weights[0].array(*mfi);
    Array4<Real> const& fweights_y = f_weights[1].array(*mfi);
    Array4<Real> const& fweights_z = f_weights[2].array(*mfi);
    Array4<const Real> const& faceCoordsx_x = face_coords[0][0].array(*mfi);
    Array4<const Real> const& faceCoordsx_y = face_coords[0][1].array(*mfi);
    Array4<const Real> const& faceCoordsx_z = face_coords[0][2].array(*mfi);
    Array4<const Real> const& faceCoordsy_x = face_coords[1][0].array(*mfi);
    Array4<const Real> const& faceCoordsy_y = face_coords[1][1].array(*mfi);
    Array4<const Real> const& faceCoordsy_z = face_coords[1][2].array(*mfi);
    Array4<const Real> const& faceCoordsz_x = face_coords[2][0].array(*mfi);
    Array4<const Real> const& faceCoordsz_y = face_coords[2][1].array(*mfi);
    Array4<const Real> const& faceCoordsz_z = face_coords[2][2].array(*mfi);

    const int spaceDim = AMREX_SPACEDIM;

    // x-components
    amrex::ParallelFor(bx_x,
      [invdx,invvol,v_spread,fout_x,fweights_x,faceCoordsx_x,faceCoordsx_y,faceCoordsx_z,spaceDim,kernel_ptr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
         Real pos_grid[3];
         pos_grid[0]=(pos_grid[0]-faceCoordsx_x(i,j,k))*invdx[0];
         pos_grid[1]=(pos_grid[1]-faceCoordsx_y(i,j,k))*invdx[1];
         pos_grid[2]=(pos_grid[2]-faceCoordsx_z(i,j,k))*invdx[2];

         Real weight = 1.0;
         for (int i=0; i<spaceDim; ++i)
             weight*=(*kernel_ptr)(pos_grid[i]);

         fout_x(i,j,k) += v_spread[0]*weight*invvol;
         fweights_x(i,j,k) += weight;
    });

    // y-components
    amrex::ParallelFor(bx_y,
      [invdx,invvol,v_spread,fout_y,fweights_y,faceCoordsy_x,faceCoordsy_y,faceCoordsy_z,spaceDim,kernel_ptr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
         Real pos_grid[3];
         pos_grid[0]=(pos_grid[0]-faceCoordsy_x(i,j,k))*invdx[0];
         pos_grid[1]=(pos_grid[1]-faceCoordsy_y(i,j,k))*invdx[1];
         pos_grid[2]=(pos_grid[2]-faceCoordsy_z(i,j,k))*invdx[2];

         Real weight = 1.0;
         for (int i=0; i<spaceDim; ++i)
             weight*=(*kernel_ptr)(pos_grid[i]);

         fout_y(i,j,k) += v_spread[1]*weight*invvol;
         fweights_y(i,j,k) += weight;
    });

    // z-components
    amrex::ParallelFor(bx_z,
      [invdx,invvol,v_spread,fout_z,fweights_z,faceCoordsz_x,faceCoordsz_y,faceCoordsz_z,spaceDim,kernel_ptr]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
         Real pos_grid[3];
         pos_grid[0]=(pos_grid[0]-faceCoordsz_x(i,j,k))*invdx[0];
         pos_grid[1]=(pos_grid[1]-faceCoordsz_y(i,j,k))*invdx[1];
         pos_grid[2]=(pos_grid[2]-faceCoordsz_z(i,j,k))*invdx[2];

         Real weight = 1.0;
         for (int i=0; i<spaceDim; ++i)
             weight*=(*kernel_ptr)(pos_grid[i]);

         fout_z(i,j,k) += v_spread[0]*weight*invvol;
         fweights_z(i,j,k) += weight;
    });
}


void IBParticleContainer::SpreadMarkers(int lev, const ParticleIndex & pindex,
        const Vector<RealVect> & f_in, std::array<MultiFab, AMREX_SPACEDIM> & f_out) const {

    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_out[d].boxArray(), f_out[d].DistributionMap(),
                            1, f_out[d].nGrow());
        f_weights[d].setVal(0.);
    }

    SpreadMarkers(lev, pindex, f_in, f_out, f_weights);

}



void IBParticleContainer::InvInterpolateMarkers(int lev, const ParticleIndex & pindex,
        const Vector<RealVect> & f_in, std::array<MultiFab, AMREX_SPACEDIM> & f_out) const {


    //___________________________________________________________________________
    // Don't do anything if pindex isn't on this rank
    auto part_it = marker_positions[lev].find(pindex);
    if (part_it == marker_positions[lev].end())
        return;


    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = Geom(0);
    const Real     *   dx = geom.CellSize();

    const int n_marker = marker_positions[lev].at(pindex).size();


    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBParticleContainer is on a differnt grid
    // than the grid which we're inverse interpolating to

    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);

    MultiFab dummy(ba, dm, 1, f_out[0].nGrow());

    Vector<int> pkernel_fluid_in(1);
    for (int i=0; i<1; ++i) {
        pkernel_fluid_in[i] = pkernel_fluid[i];
    }

    for (MFIter mfi(dummy); mfi.isValid(); ++mfi) {

        const Box & bx       = mfi.growntilebox();
        // const Box & tile_box = mfi.tilebox();
        const Box & tile_box = mfi.growntilebox(); // HACK, use tilebox going forward

        inv_interpolate_markers(BL_TO_FORTRAN_BOX(bx),
                                BL_TO_FORTRAN_BOX(tile_box),
                                BL_TO_FORTRAN_ANYD(f_out[0][mfi]),
#if (AMREX_SPACEDIM > 1)
                                BL_TO_FORTRAN_ANYD(f_out[1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                                BL_TO_FORTRAN_ANYD(f_out[2][mfi]),
#endif
                                BL_TO_FORTRAN_ANYD(face_coords[lev][0][mfi]),
#if (AMREX_SPACEDIM > 1)
                                BL_TO_FORTRAN_ANYD(face_coords[lev][1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                                BL_TO_FORTRAN_ANYD(face_coords[lev][2][mfi]),
#endif
                                marker_positions[lev].at(pindex).dataPtr(),
                                f_in.dataPtr(),
                                & n_marker,
                                dx,
                                pkernel_fluid_in.dataPtr() );
    }
}



void IBParticleContainer::InterpolateMarkers(int lev, const ParticleIndex & pindex,
        Vector<RealVect> & f_out, const std::array<MultiFab, AMREX_SPACEDIM> & f_in,
        const std::array<MultiFab, AMREX_SPACEDIM> & f_weights) const {


    //___________________________________________________________________________
    // Don't do anything if pindex isn't on this rank
    auto part_it = marker_positions[lev].find(pindex);
    if (part_it == marker_positions[lev].end())
        return;


    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = Geom(0);
    const Real     *   dx = geom.CellSize();

    const int n_marker = marker_positions[lev].at(pindex).size();


    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBParticleContainer is on a differnt grid
    // than the grid which we're interpolating from

    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);

    MultiFab dummy(ba, dm, 1, f_in[0].nGrow());

    Vector<int> pkernel_fluid_in(1);
    for (int i=0; i<1; ++i) {
        pkernel_fluid_in[i] = pkernel_fluid[i];
    }

    for (MFIter mfi(dummy); mfi.isValid(); ++mfi) {

        Box bx       = mfi.growntilebox();
        Box tile_box = mfi.tilebox();

        interpolate_markers(BL_TO_FORTRAN_BOX(bx),
                            BL_TO_FORTRAN_BOX(tile_box),
                            BL_TO_FORTRAN_ANYD(f_in[0][mfi]),
#if (AMREX_SPACEDIM > 1)
                            BL_TO_FORTRAN_ANYD(f_in[1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                            BL_TO_FORTRAN_ANYD(f_in[2][mfi]),
#endif
                            BL_TO_FORTRAN_ANYD(f_weights[0][mfi]),
#if (AMREX_SPACEDIM > 1)
                            BL_TO_FORTRAN_ANYD(f_weights[1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                            BL_TO_FORTRAN_ANYD(f_weights[2][mfi]),
#endif
                            BL_TO_FORTRAN_ANYD(face_coords[lev][0][mfi]),
#if (AMREX_SPACEDIM > 1)
                            BL_TO_FORTRAN_ANYD(face_coords[lev][1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                            BL_TO_FORTRAN_ANYD(face_coords[lev][2][mfi]),
#endif
                            marker_positions[lev].at(pindex).dataPtr(),
                            f_out.dataPtr(),
                            & n_marker,
                            dx,
                            pkernel_fluid_in.dataPtr() );
    }
}



void IBParticleContainer::InterpolateMarkers(int lev, const ParticleIndex & pindex,
        Vector<RealVect> & f_out, const std::array<MultiFab, AMREX_SPACEDIM> & f_in) const {

    //___________________________________________________________________________
    // We don't need these spreading weights => create a dummy MF
    std::array<MultiFab, AMREX_SPACEDIM> f_weights;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        f_weights[d].define(f_in[d].boxArray(), f_in[d].DistributionMap(),
                            1, f_in[d].nGrow());
        f_weights[d].setVal(-1.); // Set to <0 to guarantee that weights are ignored
    }

    InterpolateMarkers(lev, pindex, f_out, f_in, f_weights);

}



void IBParticleContainer::InvSpreadMarkers(int lev, const ParticleIndex & pindex,
        Vector<RealVect> & f_out, const std::array<MultiFab, AMREX_SPACEDIM> & f_in) const {


    //___________________________________________________________________________
    // Don't do anything if pindex isn't on this rank
    auto part_it = marker_positions[lev].find(pindex);
    if (part_it == marker_positions[lev].end())
        return;


    //___________________________________________________________________________
    // Geometry data
    const Geometry & geom = Geom(0);
    const Real     *   dx = geom.CellSize();

    const int n_marker = marker_positions[lev].at(pindex).size();


    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBParticleContainer is on a differnt grid
    // than the grid which we're inverse spreading from

    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);

    MultiFab dummy(ba, dm, 1, f_in[0].nGrow());

    Vector<int> pkernel_fluid_in(1);
    for (int i=0; i<1; ++i) {
        pkernel_fluid_in[i] = pkernel_fluid[i];
    }

    for (MFIter mfi(dummy); mfi.isValid(); ++mfi) {

        const Box & bx       = mfi.growntilebox();
        // const Box & tile_box = mfi.tilebox();
        const Box & tile_box = mfi.growntilebox(); // HACK, use tilebox going forward

        inv_spread_markers(BL_TO_FORTRAN_BOX(bx),
                           BL_TO_FORTRAN_BOX(tile_box),
                           BL_TO_FORTRAN_ANYD(f_in[0][mfi]),
#if (AMREX_SPACEDIM > 1)
                           BL_TO_FORTRAN_ANYD(f_in[1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                           BL_TO_FORTRAN_ANYD(f_in[2][mfi]),
#endif
                           BL_TO_FORTRAN_ANYD(face_coords[lev][0][mfi]),
#if (AMREX_SPACEDIM > 1)
                           BL_TO_FORTRAN_ANYD(face_coords[lev][1][mfi]),
#endif
#if (AMREX_SPACEDIM > 2)
                           BL_TO_FORTRAN_ANYD(face_coords[lev][2][mfi]),
#endif
                           marker_positions[lev].at(pindex).dataPtr(),
                           f_out.dataPtr(),
                           & n_marker,
                           dx,
                        pkernel_fluid_in.dataPtr() );
    }
}



void IBParticleContainer::InterpolateParticleForces(int lev,
        const std::array<MultiFab, AMREX_SPACEDIM> & force, const IBCore & ib_core,
        std::map<ParticleIndex, std::array<Real, AMREX_SPACEDIM>> & particle_forces
    ) {


    fillNeighbors();


    /****************************************************************************
     *                                                                          *
     * Collect the list of particles which are also contained in the IBCore     *
     *                                                                          *
     ***************************************************************************/

    Vector<ParticleIndex> index_list;
    std::map<ParticleIndex, int> ibm_index_map;

    //___________________________________________________________________________
    // Add particles to list
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {
        // MuliFabs are indexed using a pair: (BoxArray index, tile index):
        PairIndex index(pti.index(), pti.LocalTileIndex());

        auto & particle_data = GetParticles(lev)[index];
        long np = GetParticles(lev)[index].numParticles();

        // Iterate over local particle data
        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];
            ParticleIndex pindex(part.id(), part.cpu());

            int ibm_index = ib_core.get_IBMIndex(pindex);
            if (ibm_index == -1) continue;

            // Add to list if it's not already there (don't forget the
            // critical, just in case std::vector::push_back is not atomic).
            // TODO: this implementation needs work to increase performance.
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                auto search = std::find(std::begin(index_list), std::end(index_list), pindex);
                if (search == std::end(index_list)) {
                    index_list.push_back(pindex);
                    ibm_index_map[pindex] = ibm_index;
                }
            }
        }


        long ng = neighbors[lev][index].size();

        // Iterate over neighbor particle data
        // ParticleType * nbhd_data = (ParticleType *) neighbors[lev][index].dataPtr();
        ParticleVector & nbhd_data = GetNeighbors(lev, pti.index(), pti.LocalTileIndex()).GetArrayOfStructs()();
        for (int i = 0; i < ng; ++i) {
            ParticleType & part = nbhd_data[i];
            ParticleIndex pindex(part.id(), part.cpu());

            int ibm_index = ib_core.get_IBMIndex(pindex);
            if (ibm_index == -1) continue;

            // Add to list if it's not already there (don't forget the
            // critical, just in case std::vector::push_back is not atomic).
            // TODO: this implementation needs work to increase performance.
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                auto search = std::find(std::begin(index_list), std::end(index_list), pindex);
                if (search == std::end(index_list)) {
                    index_list.push_back(pindex);
                    ibm_index_map[pindex] = ibm_index;
                }
            }
        }
    }


    /****************************************************************************
     *                                                                          *
     * Ensure that the force data `force_ibm` has enough ghost cells            *
     *                                                                          *
     ***************************************************************************/

    std::array<MultiFab, AMREX_SPACEDIM> force_buffer;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        force_buffer[d].define(
                force[d].boxArray(), force[d].DistributionMap(),
                1, get_nghost()
            );

        force_buffer[d].setVal(0.);
        MultiFab::Copy(force_buffer[d], force[d], 0, 0, 1, 0);
        force_buffer[d].FillBoundary(Geom(lev).periodicity());
    }


    // Iterate over cell-centered MultiFab `dummy` as reference for face-centered data
    BoxArray force_grids_cc = convert(force[0].boxArray(), IntVect::TheCellVector());
    MultiFab dummy(force_grids_cc, force[0].DistributionMap(), 1, get_nghost());


    std::map<ParticleIndex, std::array<FArrayBox, AMREX_SPACEDIM>> pforce_buffer;
    for (const ParticleIndex & pindex : index_list) {

        int index_ibm = ibm_index_map[pindex];
        Box pbox_cc   = ib_core.get_IBMBox(index_ibm);

        std::array<FArrayBox, AMREX_SPACEDIM> pforce;
        for (int d=0; d<AMREX_SPACEDIM; ++d) {
            Box pbox_face = convert(pbox_cc, nodal_flag_dir[d]);

            // FArrayBox has deleted its copy and copy-assignement
            // constructors. We'll therefore use a particular feature of
            // std::map::operator[] => if pindex isn't already in the map, the
            // default FArrayBox constructor is invoked.
            pforce_buffer[pindex][d].resize(pbox_face);
            pforce_buffer[pindex][d].setVal(0.);
        }
    }


#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(dummy, true); mfi.isValid(); ++ mfi) {
        const Box & tile_box = mfi.growntilebox();

        for (const ParticleIndex & pindex : index_list) {
            int index_ibm = ibm_index_map[pindex];

            Box pbox_cc     = ib_core.get_IBMBox(index_ibm);
            Box work_region = tile_box & pbox_cc;

            if (work_region.ok()) {
                for (int d=0; d<AMREX_SPACEDIM; ++d) {
                    Box wregion_fc = convert(work_region, nodal_flag_dir[d]);

                    pforce_buffer[pindex][d].copy(
                            force_buffer[d][mfi],
                            wregion_fc, 0,
                            wregion_fc, 0, 1);
                }
            }
        }
    }




    /****************************************************************************
     *                                                                          *
     * Compute Hydrodynamic Forces                                              *
     *                                                                          *
     ***************************************************************************/

    for (const ParticleIndex & pindex : index_list) {
        std::array<Real, AMREX_SPACEDIM> f_trans;
        int index_ibm = ibm_index_map[pindex];

        ib_core.InterpolateForce(pforce_buffer[pindex], lev, index_ibm, pindex, f_trans);
        particle_forces[pindex] = f_trans;
    }
}



void IBParticleContainer::MoveIBParticles(int lev, Real dt,
        const std::map<ParticleIndex, std::array<Real, AMREX_SPACEDIM>> & particle_forces) {


    for (IBParIter pti(*this, lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = GetParticles(lev)[index];
        long np = GetParticles(lev)[index].numParticles();

        AoS & particles = particle_data.GetArrayOfStructs();
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];
            ParticleIndex pindex(part.id(), part.cpu());

            // std::map::operator[] requires non-const particle_forces => use at()
            std::array<Real, AMREX_SPACEDIM> f = particle_forces.at(pindex);
            Real mass = part.rdata(IBP_realData::mass);

            // Standard Euler update. TODO: RK or Verlet?
            part.rdata(IBP_realData::velx) += dt * f[0] / mass;
            part.rdata(IBP_realData::vely) += dt * f[1] / mass;
            part.rdata(IBP_realData::velz) += dt * f[2] / mass;

            part.pos(0) += dt *  part.rdata(IBP_realData::velx);
            part.pos(1) += dt *  part.rdata(IBP_realData::vely);
            part.pos(2) += dt *  part.rdata(IBP_realData::velz);
        }
    }
}



void IBParticleContainer::PrintParticleData(int lev) {

    // Inverse cell-size vector => max is used for determining IBParticle
    // radius in units of cell size
    Vector<Real> inv_dx = {
            AMREX_D_DECL(
                Geom(lev).InvCellSize(0),
                Geom(lev).InvCellSize(1),
                Geom(lev).InvCellSize(2)
            )
        };
    // Find max inv_dx (in case we have an anisotropic grid)
    Real mx_inv_dx = * std::max_element(inv_dx.begin(), inv_dx.end());

    amrex::AllPrintToFile("ib_particle_data") << "Particles on each box:" << std::endl;

    fillNeighbors();

    long local_count = 0;

    // ParIter skips tiles without particles => Iterate over MultiFab instead
    // of ParticleIter. Note also that AmrexParticleContainer uses wired tiling =>
    // turn tiling off
    for(MFIter pti = MakeMFIter(lev, true); pti.isValid(); ++pti) {
        // MuliFabs are indexed using a pair: (BoxArray index, tile index):
        PairIndex index(pti.index(), pti.LocalTileIndex());

        // Neighbours are stored as raw data (see below)
        int ng = neighbors[lev][index].size();

        auto & particle_data = GetParticles(lev)[index];
        long np = GetParticles(lev)[index].numParticles();

        local_count += np;

        // Print current box info
        AllPrintToFile("ib_particle_data") << "Box:"         << pti.index()
                                           << " "            << pti.tilebox()
                                           << ", count: "    << np
                                           << ", nb count: " << ng
                                           << std::endl;

        // Print IBParticle
        AllPrintToFile("ib_particle_data") << " * IBPartcies:" << std::endl;

        //AoS & particles = pti.GetArrayOfStructs();
        AoS & particles = particle_data.GetArrayOfStructs();
        for(int i = 0; i < np; i++){
            ParticleType & part = particles[i];
            Real r              = part.rdata(IBP_realData::radius);

            int r_ncx = (int) (r * mx_inv_dx);

            AllPrintToFile("ib_particle_data") << "   +- " << part << std::endl;
            AllPrintToFile("ib_particle_data") << "   +---> Radius [NCells]: " << r_ncx << std::endl;
        }

        // Print neighbour IBParticles
        AllPrintToFile("ib_particle_data") << " * Grown IBParticles:" << std::endl;

        // TODO: HAXOR!!! This should be fixed ASAP: if I understand this
        // correctly, the neighbor data contains the particle data as a binary
        // array (char). By casting to ParticleType, what we're doing is
        // interpreting the data in neighbours[index] as valid particle data.
        // Also we stride the neighbors[index] array in units of
        // sizeof(ParticleData). All of this is a little too dangerous for my
        // taste: never hide what you're doing from your compiler!!!
        ParticleType * nbhd_data = (ParticleType *) neighbors[lev][index].GetArrayOfStructs().dataPtr();
        for (int i = 0; i < ng; i++){
            ParticleType & part = nbhd_data[i];
            Real r              = part.rdata(IBP_realData::radius);

            int r_ncx = (int) (r * mx_inv_dx);

            AllPrintToFile("ib_particle_data") << "   +- " << part << std::endl;
            AllPrintToFile("ib_particle_data") << "   +---> Radius [NCells]: " << r_ncx << std::endl;
        }
    }

    AllPrintToFile("ib_particle_data") << "Total for this process: "
                                       << local_count << std::endl << std::endl;
}



void IBParticleContainer::LocalIBParticleInfo(Vector<IBP_info> & info,
                                              int lev, PairIndex index,
                                              bool unique) const {

    // Inverse cell-size vector => used for determining index corresponding to
    // IBParticle position (pos)
    RealVect inv_dx = RealVect(
            AMREX_D_DECL(
                Geom(lev).InvCellSize(0),
                Geom(lev).InvCellSize(1),
                Geom(lev).InvCellSize(2)
            )
        );


    auto & particle_data = GetParticles(lev).at(index);
    long np = GetParticles(lev).at(index).numParticles();

    // Iterate over local particle data
    const AoS & particles = particle_data.GetArrayOfStructs();
    for(int i = 0; i < np; i++){
        const ParticleType & part = particles[i];

        // Position of IBParticle
        RealVect pos = RealVect(
                AMREX_D_DECL(part.pos(0), part.pos(1), part.pos(2))
            );

        // Velocity of IBParticle
        RealVect vel = RealVect(
                AMREX_D_DECL(part.rdata(IBP_realData::velx),
                             part.rdata(IBP_realData::vely),
                             part.rdata(IBP_realData::velz)   )
            );
        // HACK HARD CODE IN POSITION
        RealVect ori = RealVect(
                AMREX_D_DECL(0.,0.,1. )
            );

        // Position of IBParticle on grid
        RealVect pos_grid = pos * inv_dx;
        IntVect  pos_ind  = IntVect(
                AMREX_D_DECL( (int) pos_grid[0], (int) pos_grid[1], (int) pos_grid[2] )
            );

        // Radius
        Real r = part.rdata(IBP_realData::radius);

        // Construct info struct
        IBP_info part_info;
        part_info.pos    = pos;
        part_info.vel    = vel;
        part_info.ori    = ori;
        part_info.index  = pos_ind;
        part_info.radius = r;
        part_info.id     = part.id();
        part_info.cpu    = part.cpu();
        part_info.real   = 1; // 1 => real (non-neighbor particle)

        // Add to list

        if (unique) {
            // If in unique-mode: Don't add unless `part_info` is not already in `info`
            const auto & search = std::find(std::begin(info), std::end(info), part_info);
            if ( search == std::end(info)) {
#ifdef _OPENMP
#pragma omp critical
#endif
                { info.push_back(part_info); }
            }
        } else {
#ifdef _OPENMP
#pragma omp critical
#endif
            { info.push_back(part_info); }
        }
    }
}



Vector<IBP_info> IBParticleContainer::LocalIBParticleInfo(int lev, PairIndex index,
                                                          bool unique) const {

    // Allocate Particle Info vector
    Vector<IBP_info> info;

    //___________________________________________________________________________
    // Fill Particle Info vector with local (non-neighbour) data
    LocalIBParticleInfo(info, lev, index, unique);


    return info;
}



Vector<IBP_info> IBParticleContainer::LocalIBParticleInfo(int lev, bool unique) const {

    // Allocate Particle Info vector
    Vector<IBP_info> info;

    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBParticleContainer is on a differnt grid
    // than the grid which we're searching for particles (this should usually
    // be fine though)

    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);

    MultiFab dummy(ba, dm, 1, 1);

    //___________________________________________________________________________
    // Iterate over `dummy` looking for particles. NOTE: use the
    // IBParticleContainer tile size
    for (MFIter mfi(dummy, tile_size); mfi.isValid(); ++mfi){
        PairIndex index(mfi.index(), mfi.LocalTileIndex());
        LocalIBParticleInfo(info, lev, index, true);
    }


    return info;
}



void IBParticleContainer::NeighborIBParticleInfo(Vector<IBP_info> & info,
                                                 int lev, PairIndex index,
                                                 bool unique) const {

    RealVect inv_dx = RealVect(
            AMREX_D_DECL(
                Geom(lev).InvCellSize(0),
                Geom(lev).InvCellSize(1),
                Geom(lev).InvCellSize(2)
            )
        );

    int ng = neighbors[lev].at(index).size();

    // Iterate over neighbour particles:
    // TODO: HAXOR!!! This should be fixed ASAP: if I understand this correctly,
    // the neighbor data contains the particle data as a binary array (char). By
    // casting to ParticleType, what we're doing is interpreting the data in
    // neighbours[index] as valid particle data. Also we stride the
    // neighbors[index] array in units of sizeof(ParticleData). All of this is a
    // little too dangerous for my taste: never hide what you're doing from your
    // compiler!!!
    const ParticleType * nbhd_data = (ParticleType *) neighbors[lev].at(index).GetArrayOfStructs().dataPtr();
    for(int i = 0; i < ng; i++){
        const ParticleType & part = nbhd_data[i];

        // Position of neighbour IBParticle
        RealVect pos = RealVect(
                AMREX_D_DECL(part.pos(0), part.pos(1), part.pos(2))
            );

        // Velocity of IBParticle
        RealVect vel = RealVect(
                AMREX_D_DECL(part.rdata(IBP_realData::velx),
                             part.rdata(IBP_realData::vely),
                             part.rdata(IBP_realData::velz)   )
            );

        // Position of neighbour IBParticle on grid
        RealVect pos_grid = pos * inv_dx;
        IntVect  pos_ind  = IntVect(
                AMREX_D_DECL( (int) pos_grid[0], (int) pos_grid[1], (int) pos_grid[2] )
            );
        // HACK HARD CODE IN POSITION
        RealVect ori = RealVect(
                AMREX_D_DECL(0. ,0. ,1. )
            );

        // Radius
        Real r = part.rdata(IBP_realData::radius);

        // Construct info struct
        IBP_info part_info;
        part_info.pos    = pos;
        part_info.vel    = vel;
        part_info.ori    = ori;
        part_info.index  = pos_ind;
        part_info.radius = r;
        part_info.id     = part.id();
        part_info.cpu    = part.cpu();
        part_info.real   = 0; // 0 => neighbor particle

        // Add to list


        if (unique) {
            // If in unique-mode: Don't add unless `part_info` is not already in `info`
            const auto & search = std::find(std::begin(info), std::end(info), part_info);
            if ( search == std::end(info)) {
#ifdef _OPENMP
#pragma omp critical
#endif
                { info.push_back(part_info); }
            }
        } else {
#ifdef _OPENMP
#pragma omp critical
#endif
            { info.push_back(part_info); }
        }
    }
}



Vector<IBP_info> IBParticleContainer::NeighborIBParticleInfo(int lev, PairIndex index,
                                                             bool unique) const {

    // Allocate Particle Info vector
    Vector<IBP_info> info;

    //___________________________________________________________________________
    // Fill Particle Info vector with neighbour data
    NeighborIBParticleInfo(info, lev, index, unique);


    return info;
}



Vector<IBP_info> IBParticleContainer::NeighborIBParticleInfo(int lev, bool unique) const {

    // Allocate Particle Info vector
    Vector<IBP_info> info;

    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBParticleContainer is on a differnt grid
    // than the grid which we're searching for particles (this should usually
    // be fine though)

    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);

    MultiFab dummy(ba, dm, 1, 1);

    //___________________________________________________________________________
    // Iterate over `dummy` looking for particles. NOTE: use the
    // IBParticleContainer tile size
    for (MFIter mfi(dummy, tile_size); mfi.isValid(); ++mfi){
        PairIndex index(mfi.index(), mfi.LocalTileIndex());
        NeighborIBParticleInfo(info, lev, index, true);
    }


    return info;
}



void IBParticleContainer::IBParticleInfo(Vector<IBP_info> & info, int lev, PairIndex index,
                                         bool unique) const {

    //___________________________________________________________________________
    // Fill Particle Info vector with local (non-neighbour) and neighbour data
       LocalIBParticleInfo(info, lev, index, unique);
    NeighborIBParticleInfo(info, lev, index, unique);
}



Vector<IBP_info> IBParticleContainer::IBParticleInfo(int lev, PairIndex index,
                                                     bool unique) const {

    // Allocate Particle Info vector
    Vector<IBP_info> info;

    //___________________________________________________________________________
    // Fill Particle Info vector with local (non-neighbour) and neighbour data
       LocalIBParticleInfo(info, lev, index, unique);
    NeighborIBParticleInfo(info, lev, index, unique);


    return info;
}



Vector<IBP_info> IBParticleContainer::IBParticleInfo(int lev, bool unique) const {

    // Allocate Particle Info vector
    Vector<IBP_info> info;

    //___________________________________________________________________________
    // Cell-centered MultiFab used as a reference for iterating over data
    // WARNING: this will break if IBParticleContainer is on a differnt grid
    // than the grid which we're searching for particles (this should usually
    // be fine though)

    const BoxArray & ba            = ParticleBoxArray(lev);
    const DistributionMapping & dm = ParticleDistributionMap(lev);

    MultiFab dummy(ba, dm, 1, 1);

    //___________________________________________________________________________
    // Iterate over `dummy` looking for particles. NOTE: use the
    // IBParticleContainer tile size
    for (MFIter mfi(dummy, tile_size); mfi.isValid(); ++mfi){
        PairIndex index(mfi.index(), mfi.LocalTileIndex());
        IBParticleInfo(info, lev, index, true);
    }


    return info;
}
