#include "compressible_functions_stag.H"
#include "compressible_functions.H"
#include "reservoirStag_K.H"
#include "common_functions.H"
#include "rng_functions.H"
#include <math.h>

///////////////////////////////////////////////////////////////////////////
// compute momentum and fluxes at the reservoir-FHD interface /////////////
///////////////////////////////////////////////////////////////////////////
void
ComputeFluxMomReservoir(const MultiFab& cons0_in, const MultiFab& prim0_in,
                        const std::array<MultiFab, AMREX_SPACEDIM>& vel0,
                        std::array<MultiFab, AMREX_SPACEDIM>& cumom_res,
                        std::array<MultiFab, AMREX_SPACEDIM>& faceflux_res,
                        const amrex::Geometry& geom,
                        const amrex::Real dt)
{
    BL_PROFILE_VAR("ComputeFluxMomReservoir()",ComputeFluxMomReservoir);

    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    Box dom(geom.Domain());

    GpuArray<Real,MAX_SPECIES> mass;
    for (int l=0;l<nspecies;++l) {
        mass[l] = molmass[l]/avogadro;
    }

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        cumom_res[d].setVal(0.0);
        faceflux_res[d].setVal(0.0);
    }

    Real area, vol;

    // Reservoir in LO X
    if (bc_mass_lo[0] == 4) {

        area = dx[1]*dx[2];
        vol  = dx[0]*dx[1]*dx[2];

        // face-based flux (mass and energy) and normal momentum /////
        //////////////////////////////////////////////////////////////

        // domain grown nodally based on faceflux_res[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_res[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xlo = amrex::bdryNode(dom_x, Orientation(0, Orientation::low));

        for (MFIter mfi(faceflux_res[0]); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xlo;

            const Array4<const Real> cons0   = cons0_in.array(mfi);
            const Array4<const Real> prim0   = prim0_in.array(mfi);
            const Array4<const Real>& xvel0  = (vel0[0]).array(mfi);
            const Array4<const Real>& yvel0  = (vel0[1]).array(mfi);
            const Array4<const Real>& zvel0  = (vel0[2]).array(mfi);

            const Array4<Real>& xmom  = (cumom_res[0]).array(mfi);
            const Array4<Real>& xflux = (faceflux_res[0]).array(mfi);

            if (b.ok()) {
                amrex::ParallelForRNG(b, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
                {

                    ////////////////// from reservoir //////////////////

                    ////////////////////////////////////////////////////
                    // number of particles in FHD cell (for correction)
                    Real N = 0.0;
                    GpuArray<Real,MAX_SPECIES> N_i;
                    for (int n=0;n<nspecies;++n) {
                        N_i[n] = vol*(cons0(i,j,k,5+n)/mass[n]);
                        N += N_i[n];
                    }
                    ////////////////////////////////////////////////////

                    Real T = prim0(i-1,j,k,4);
                    Real Vx, Vy, Vz;
                    Vx = 0.0; // normal
                    Vy = 0.0;
                    Vz = 0.0;
                    Real mass_cross; // total mass crossing at the reservoir interface
                    GpuArray<Real,3> mom_cross; // total momentum crossing at the reservoir interface
                    Real en_cross; // total energy crossing at the reservoir interface
                    GpuArray<Real,MAX_SPECIES> spec_mass_cross;
                    GpuArray<Real,MAX_SPECIES> rhoYk;
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i-1,j,k,5+n);
                    }
                    if (do_1D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vy,Vz,engine);
                    }
                    else if (do_2D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vy,Vz,engine);
                    }
                    else {
                        poisson_process_reservoir(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vy,Vz,engine);
                    }
                    xflux(i,j,k,0) += (1.0 - (1.0/(12.0*N)))*mass_cross/(dt*area); // update mass flux
                    xflux(i,j,k,4) += (1.0 + (1.0/( 4.0*N)))*en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) += (1.0 - (1.0/(12.0*N)))*spec_mass_cross[n]/(dt*area); // update species flux
                    }

                    if (do_1D) {
                        xmom(i,j,k) += mom_cross[0]/vol;
                    }
                    else if (do_2D) {
                        xmom(i,j,k) += mom_cross[0]/vol;
                        xflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
                    }
                    else {
                        xmom(i,j,k) += mom_cross[0]/vol;
                        xflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
                        xflux(i,j,k,2) += mom_cross[2]/dt/area; // tangential flux
                    }

                    ////////////////////////////////////////////////////

                    ////////////////// to reservoir ////////////////////

                    T = prim0(i,j,k,4);
                    Vx = -0.5*(xvel0(i,j,k)+xvel0(i+1,j,k));
                    Vy =  0.5*(yvel0(i,j,k)+yvel0(i,j+1,k));
                    Vz =  0.5*(zvel0(i,j,k)+zvel0(i,j,k+1));
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i,j,k,5+n);
                    }
                    if (do_1D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vy,Vz,engine);
                    }
                    else if (do_2D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vy,Vz,engine);
                    }
                    else {
                        poisson_process_reservoir(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vy,Vz,engine);
                    }
                    xflux(i,j,k,0) -= mass_cross/(dt*area); // update mass flux
                    xflux(i,j,k,4) -= en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) -= spec_mass_cross[n]/(dt*area); // update species flux
                    }

                    if (do_1D) {
                        xmom(i,j,k) -= mom_cross[0]/vol;
                    }
                    else if (do_2D) {
                        xmom(i,j,k) -= mom_cross[0]/vol;
                        xflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
                    }
                    else {
                        xmom(i,j,k) -= mom_cross[0]/vol;
                        xflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
                        xflux(i,j,k,2) -= mom_cross[2]/dt/area; // tangential flux
                    }

                    ////////////////////////////////////////////////////
                });
            }
        }
    }

    // Reservoir in HI X
    if (bc_mass_hi[0] == 4) {

        area = dx[1]*dx[2];
        vol  = dx[0]*dx[1]*dx[2];

        // face-based flux (mass and energy) and normal momentum /////
        //////////////////////////////////////////////////////////////

        // domain grown nodally based on faceflux_res[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_res[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xhi = amrex::bdryNode(dom_x, Orientation(0, Orientation::high));

        for (MFIter mfi(faceflux_res[0]); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xhi;

            const Array4<const Real> cons0   = cons0_in.array(mfi);
            const Array4<const Real> prim0   = prim0_in.array(mfi);
            const Array4<const Real>& xvel0  = (vel0[0]).array(mfi);
            const Array4<const Real>& yvel0  = (vel0[1]).array(mfi);
            const Array4<const Real>& zvel0  = (vel0[2]).array(mfi);

            const Array4<Real>& xmom  = (cumom_res[0]).array(mfi);
            const Array4<Real>& xflux = (faceflux_res[0]).array(mfi);

            if (b.ok()) {
                amrex::ParallelForRNG(b, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
                {

                    ////////////////// from reservoir //////////////////

                    ////////////////////////////////////////////////////
                    // number of particles in FHD cell (for correction)
                    Real N = 0.0;
                    GpuArray<Real,MAX_SPECIES> N_i;
                    for (int n=0;n<nspecies;++n) {
                        N_i[n] = vol*(cons0(i-1,j,k,5+n)/mass[n]);
                        N += N_i[n];
                    }
                    ////////////////////////////////////////////////////

                    Real T = prim0(i,j,k,4);
                    Real Vx, Vy, Vz;
                    Vx = 0.0; // normal
                    Vy = 0.0;
                    Vz = 0.0;
                    Real mass_cross; // total mass crossing at the reservoir interface
                    GpuArray<Real,3> mom_cross; // total momentum crossing at the reservoir interface
                    Real en_cross; // total energy crossing at the reservoir interface
                    GpuArray<Real,MAX_SPECIES> spec_mass_cross;
                    GpuArray<Real,MAX_SPECIES> rhoYk;
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i,j,k,5+n);
                    }
                    if (do_1D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vy,Vz,engine);
                    }
                    else if (do_2D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vy,Vz,engine);
                    }
                    else {
                        poisson_process_reservoir(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vy,Vz,engine);
                    }
                    xflux(i,j,k,0) -= (1.0 - (1.0/(12.0*N)))*mass_cross/(dt*area); // update mass flux
                    xflux(i,j,k,4) -= (1.0 + (1.0/( 4.0*N)))*en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) -= (1.0 - (1.0/(12.0*N)))*spec_mass_cross[n]/(dt*area); // update species flux
                    }

                    if (do_1D) {
                        xmom(i,j,k) -= mom_cross[0]/vol;
                    }
                    else if (do_2D) {
                        xmom(i,j,k) -= mom_cross[0]/vol;
                        xflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
                    }
                    else {
                        xmom(i,j,k) -= mom_cross[0]/vol;
                        xflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
                        xflux(i,j,k,2) -= mom_cross[2]/dt/area; // tangential flux
                    }

                    ////////////////////////////////////////////////////

                    ////////////////// to reservoir ////////////////////

                    T = prim0(i-1,j,k,4);
                    Vx = 0.5*(xvel0(i-1,j,k)+xvel0(i,j,k));
                    Vy = 0.5*(yvel0(i-1,j,k)+yvel0(i-1,j+1,k));
                    Vz = 0.5*(zvel0(i-1,j,k)+zvel0(i-1,j,k+1));
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i-1,j,k,5+n);
                    }
                    if (do_1D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vy,Vz,engine);
                    }
                    else if (do_2D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vy,Vz,engine);
                    }
                    else {
                        poisson_process_reservoir(mass,rhoYk,T,Vx,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vy,Vz,engine);
                    }
                    xflux(i,j,k,0) += mass_cross/(dt*area); // update mass flux
                    xflux(i,j,k,4) += en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) += spec_mass_cross[n]/(dt*area);
                    }

                    if (do_1D) {
                        xmom(i,j,k) += mom_cross[0]/vol;  // add momentum for particles going other way
                    }
                    else if (do_2D) {
                        xmom(i,j,k) += mom_cross[0]/vol;  // add momentum for particles going other way
                        xflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
                    }
                    else {
                        xmom(i,j,k) += mom_cross[0]/vol;  // add momentum for particles going other way
                        xflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
                        xflux(i,j,k,2) += mom_cross[2]/dt/area; // tangential flux
                    }

                    ////////////////////////////////////////////////////

                });
            }
        }
    }

    // Reservoir in LO Y
    if (bc_mass_lo[1] == 4) {

        area = dx[0]*dx[2];
        vol  = dx[0]*dx[1]*dx[2];

        // face-based flux (mass and energy) and normal momentum /////
        //////////////////////////////////////////////////////////////

        // domain grown nodally based on faceflux_res[0] nodality (x)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_res[1].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_ylo = amrex::bdryNode(dom_y, Orientation(1, Orientation::low));

        for (MFIter mfi(faceflux_res[1]); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_ylo;

            const Array4<const Real> cons0   = cons0_in.array(mfi);
            const Array4<const Real> prim0   = prim0_in.array(mfi);
            const Array4<const Real>& xvel0  = (vel0[0]).array(mfi);
            const Array4<const Real>& yvel0  = (vel0[1]).array(mfi);
            const Array4<const Real>& zvel0  = (vel0[2]).array(mfi);

            const Array4<Real>& ymom  = (cumom_res[1]).array(mfi);
            const Array4<Real>& yflux = (faceflux_res[1]).array(mfi);

            if (b.ok()) {
                amrex::ParallelForRNG(b, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
                {

                    ////////////////// from reservoir //////////////////

                    ////////////////////////////////////////////////////
                    // number of particles in FHD cell (for correction)
                    Real N = 0.0;
                    GpuArray<Real,MAX_SPECIES> N_i;
                    for (int n=0;n<nspecies;++n) {
                        N_i[n] = vol*(cons0(i,j,k,5+n)/mass[n]);
                        N += N_i[n];
                    }
                    ////////////////////////////////////////////////////

                    Real T = prim0(i,j-1,k,4);
                    Real Vx, Vy, Vz;
                    Vx = 0.0;
                    Vy = 0.0; // normal
                    Vz = 0.0;
                    Real mass_cross; // total mass crossing at the reservoir interface
                    GpuArray<Real,3> mom_cross; // total momentum crossing at the reservoir interface
                    Real en_cross; // total energy crossing at the reservoir interface
                    GpuArray<Real,MAX_SPECIES> spec_mass_cross;
                    GpuArray<Real,MAX_SPECIES> rhoYk;
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i,j-1,k,5+n);
                    }
                    if (do_1D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vy,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vx,Vz,engine);
                    }
                    else if (do_2D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vy,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vx,Vz,engine);
                    }
                    else {
                        poisson_process_reservoir(mass,rhoYk,T,Vy,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vx,Vz,engine);
                    }
                    yflux(i,j,k,0) += (1.0 - (1.0/(12.0*N)))*mass_cross/(dt*area); // update mass flux
                    yflux(i,j,k,4) += (1.0 + (1.0/( 4.0*N)))*en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        yflux(i,j,k,5+n) += (1.0 - (1.0/(12.0*N)))*spec_mass_cross[n]/(dt*area); // update species flux
                    }

                    if (do_1D) {
                        ymom(i,j,k) += mom_cross[0]/vol;
                    }
                    else if (do_2D) {
                        ymom(i,j,k) += mom_cross[0]/vol;
                        yflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
                    }
                    else {
                        ymom(i,j,k) += mom_cross[0]/vol;
                        yflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
                        yflux(i,j,k,2) += mom_cross[2]/dt/area; // tangential flux
                    }

                    ////////////////////////////////////////////////////

                    ////////////////// to reservoir ////////////////////

                    T = prim0(i,j,k,4);
                    Vx =  0.5*(xvel0(i,j,k)+xvel0(i+1,j,k));
                    Vy =  -0.5*(yvel0(i,j,k)+yvel0(i,j+1,k));
                    Vz =  0.5*(zvel0(i,j,k)+zvel0(i,j,k+1));
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i,j,k,5+n);
                    }
                    if (do_1D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vy,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vx,Vz,engine);
                    }
                    else if (do_2D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vy,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vx,Vz,engine);
                    }
                    else {
                        poisson_process_reservoir(mass,rhoYk,T,Vy,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vx,Vz,engine);
                    }
                    yflux(i,j,k,0) -= mass_cross/(dt*area); // update mass flux
                    yflux(i,j,k,4) -= en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        yflux(i,j,k,5+n) -= spec_mass_cross[n]/(dt*area); // update species flux
                    }

                    if (do_1D) {
                        ymom(i,j,k) -= mom_cross[0]/vol;
                    }
                    else if (do_2D) {
                        ymom(i,j,k) -= mom_cross[0]/vol;
                        yflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
                    }
                    else {
                        ymom(i,j,k) -= mom_cross[0]/vol;
                        yflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
                        yflux(i,j,k,2) -= mom_cross[2]/dt/area; // tangential flux
                    }

                    ////////////////////////////////////////////////////
                });
            }
        }
    }

    // Reservoir in HI Y
    if (bc_mass_hi[1] == 4) {

        area = dx[0]*dx[2];
        vol  = dx[0]*dx[1]*dx[2];

        // face-based flux (mass and energy) and normal momentum /////
        //////////////////////////////////////////////////////////////

        // domain grown nodally based on faceflux_res[0] nodality (x)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_res[1].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yhi = amrex::bdryNode(dom_y, Orientation(1, Orientation::high));

        for (MFIter mfi(faceflux_res[1]); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yhi;

            const Array4<const Real> cons0   = cons0_in.array(mfi);
            const Array4<const Real> prim0   = prim0_in.array(mfi);
            const Array4<const Real>& xvel0  = (vel0[0]).array(mfi);
            const Array4<const Real>& yvel0  = (vel0[1]).array(mfi);
            const Array4<const Real>& zvel0  = (vel0[2]).array(mfi);

            const Array4<Real>& ymom  = (cumom_res[1]).array(mfi);
            const Array4<Real>& yflux = (faceflux_res[1]).array(mfi);

            if (b.ok()) {
                amrex::ParallelForRNG(b, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
                {

                    ////////////////// from reservoir //////////////////

                    ////////////////////////////////////////////////////
                    // number of particles in FHD cell (for correction)
                    Real N = 0.0;
                    GpuArray<Real,MAX_SPECIES> N_i;
                    for (int n=0;n<nspecies;++n) {
                        N_i[n] = vol*(cons0(i,j-1,k,5+n)/mass[n]);
                        N += N_i[n];
                    }
                    ////////////////////////////////////////////////////

                    Real T = prim0(i,j,k,4);
                    Real Vx, Vy, Vz;
                    Vx = 0.0;
                    Vy = 0.0; // normal
                    Vz = 0.0;
                    Real mass_cross; // total mass crossing at the reservoir interface
                    GpuArray<Real,3> mom_cross; // total momentum crossing at the reservoir interface
                    Real en_cross; // total energy crossing at the reservoir interface
                    GpuArray<Real,MAX_SPECIES> spec_mass_cross;
                    GpuArray<Real,MAX_SPECIES> rhoYk;
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i,j,k,5+n);
                    }
                    if (do_1D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vy,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vx,Vz,engine);
                    }
                    else if (do_2D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vy,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vx,Vz,engine);
                    }
                    else {
                        poisson_process_reservoir(mass,rhoYk,T,Vy,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vx,Vz,engine);
                    }
                    yflux(i,j,k,0) -= (1.0 - (1.0/(12.0*N)))*mass_cross/(dt*area); // update mass flux
                    yflux(i,j,k,4) -= (1.0 + (1.0/( 4.0*N)))*en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        yflux(i,j,k,5+n) -= (1.0 - (1.0/(12.0*N)))*spec_mass_cross[n]/(dt*area); // update species flux
                    }

                    if (do_1D) {
                        ymom(i,j,k) -= mom_cross[0]/vol;
                    }
                    else if (do_2D) {
                        ymom(i,j,k) -= mom_cross[0]/vol;
                        yflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
                    }
                    else {
                        ymom(i,j,k) -= mom_cross[0]/vol;
                        yflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
                        yflux(i,j,k,2) -= mom_cross[2]/dt/area; // tangential flux
                    }

                    ////////////////////////////////////////////////////

                    ////////////////// to reservoir ////////////////////

                    T = prim0(i,j-1,k,4);
                    Vx = 0.5*(xvel0(i,j-1,k)+xvel0(i+1,j-1,k));
                    Vy = 0.5*(yvel0(i,j-1,k)+yvel0(i,j,k));
                    Vz = 0.5*(zvel0(i,j-1,k)+zvel0(i,j-1,k+1));
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i,j-1,k,5+n);
                    }
                    if (do_1D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vy,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vx,Vz,engine);
                    }
                    else if (do_2D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vy,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vx,Vz,engine);
                    }
                    else {
                        poisson_process_reservoir(mass,rhoYk,T,Vy,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vx,Vz,engine);
                    }
                    yflux(i,j,k,0) += mass_cross/(dt*area); // update mass flux
                    yflux(i,j,k,4) += en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        yflux(i,j,k,5+n) += spec_mass_cross[n]/(dt*area);
                    }

                    if (do_1D) {
                        ymom(i,j,k) += mom_cross[0]/vol;  // add momentum for particles going other way
                    }
                    else if (do_2D) {
                        ymom(i,j,k) += mom_cross[0]/vol;  // add momentum for particles going other way
                        yflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
                    }
                    else {
                        ymom(i,j,k) += mom_cross[0]/vol;  // add momentum for particles going other way
                        yflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
                        yflux(i,j,k,2) += mom_cross[2]/dt/area; // tangential flux
                    }

                    ////////////////////////////////////////////////////

                });
            }
        }
    }

    // Reservoir in LO Z
    if (bc_mass_lo[2] == 4) {

        area = dx[0]*dx[1];
        vol  = dx[0]*dx[1]*dx[2];

        // face-based flux (mass and energy) and normal momentum /////
        //////////////////////////////////////////////////////////////

        // domain grown nodally based on faceflux_res[1] nodality (y)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_res[2].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zlo = amrex::bdryNode(dom_z, Orientation(2, Orientation::low));

        for (MFIter mfi(faceflux_res[2]); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zlo;

            const Array4<const Real> cons0   = cons0_in.array(mfi);
            const Array4<const Real> prim0   = prim0_in.array(mfi);
            const Array4<const Real>& xvel0  = (vel0[0]).array(mfi);
            const Array4<const Real>& yvel0  = (vel0[1]).array(mfi);
            const Array4<const Real>& zvel0  = (vel0[2]).array(mfi);

            const Array4<Real>& zmom  = (cumom_res[2]).array(mfi);
            const Array4<Real>& zflux = (faceflux_res[2]).array(mfi);

            if (b.ok()) {
                amrex::ParallelForRNG(b, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
                {

                    ////////////////// from reservoir //////////////////

                    ////////////////////////////////////////////////////
                    // number of particles in FHD cell (for correction)
                    Real N = 0.0;
                    GpuArray<Real,MAX_SPECIES> N_i;
                    for (int n=0;n<nspecies;++n) {
                        N_i[n] = vol*(cons0(i,j,k,5+n)/mass[n]);
                        N += N_i[n];
                    }
                    ////////////////////////////////////////////////////

                    Real T = prim0(i,j,k-1,4);
                    Real Vx, Vy, Vz;
                    Vx = 0.0;
                    Vy = 0.0;
                    Vz = 0.0; // normal
                    Real mass_cross; // total mass crossing at the reservoir interface
                    GpuArray<Real,3> mom_cross; // total momentum crossing at the reservoir interface
                    Real en_cross; // total energy crossing at the reservoir interface
                    GpuArray<Real,MAX_SPECIES> spec_mass_cross;
                    GpuArray<Real,MAX_SPECIES> rhoYk;
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i,j,k-1,5+n);
                    }
                    if (do_1D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vz,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vx,Vy,engine);
                    }
                    else if (do_2D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vz,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vx,Vy,engine);
                    }
                    else {
                        poisson_process_reservoir(mass,rhoYk,T,Vz,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vx,Vy,engine);
                    }
                    zflux(i,j,k,0) += (1.0 - (1.0/(12.0*N)))*mass_cross/(dt*area); // update mass flux
                    zflux(i,j,k,4) += (1.0 + (1.0/( 4.0*N)))*en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        zflux(i,j,k,5+n) += (1.0 - (1.0/(12.0*N)))*spec_mass_cross[n]/(dt*area); // update species flux
                    }

                    if (do_1D) {
                        zmom(i,j,k) += mom_cross[0]/vol;
                    }
                    else if (do_2D) {
                        zmom(i,j,k) += mom_cross[0]/vol;
                        zflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
                    }
                    else {
                        zmom(i,j,k) += mom_cross[0]/vol;
                        zflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
                        zflux(i,j,k,2) += mom_cross[2]/dt/area; // tangential flux
                    }

                    ////////////////////////////////////////////////////

                    ////////////////// to reservoir ////////////////////

                    T = prim0(i,j,k,4);
                    Vx =  0.5*(xvel0(i,j,k)+xvel0(i+1,j,k));
                    Vy =  0.5*(yvel0(i,j,k)+yvel0(i,j+1,k));
                    Vz =  -0.5*(zvel0(i,j,k)+zvel0(i,j,k+1));
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i,j,k,5+n);
                    }
                    if (do_1D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vz,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vx,Vy,engine);
                    }
                    else if (do_2D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vz,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vx,Vy,engine);
                    }
                    else {
                        poisson_process_reservoir(mass,rhoYk,T,Vz,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vx,Vy,engine);
                    }
                    zflux(i,j,k,0) -= mass_cross/(dt*area); // update mass flux
                    zflux(i,j,k,4) -= en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        zflux(i,j,k,5+n) -= spec_mass_cross[n]/(dt*area); // update species flux
                    }

                    if (do_1D) {
                        zmom(i,j,k) -= mom_cross[0]/vol;
                    }
                    else if (do_2D) {
                        zmom(i,j,k) -= mom_cross[0]/vol;
                        zflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
                    }
                    else {
                        zmom(i,j,k) -= mom_cross[0]/vol;
                        zflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
                        zflux(i,j,k,2) -= mom_cross[2]/dt/area; // tangential flux
                    }

                    ////////////////////////////////////////////////////
                });
            }
        }
    }

    // Reservoir in HI Z
    if (bc_mass_hi[2] == 4) {

        area = dx[0]*dx[1];
        vol  = dx[0]*dx[1]*dx[2];

        // face-based flux (mass and energy) and normal momentum /////
        //////////////////////////////////////////////////////////////

        // domain grown nodally based on faceflux_res[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_res[2].ixType());

        // this is the z-lo domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zhi = amrex::bdryNode(dom_z, Orientation(2, Orientation::high));

        for (MFIter mfi(faceflux_res[2]); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zhi;

            const Array4<const Real> cons0   = cons0_in.array(mfi);
            const Array4<const Real> prim0   = prim0_in.array(mfi);
            const Array4<const Real>& xvel0  = (vel0[0]).array(mfi);
            const Array4<const Real>& yvel0  = (vel0[1]).array(mfi);
            const Array4<const Real>& zvel0  = (vel0[2]).array(mfi);

            const Array4<Real>& zmom  = (cumom_res[2]).array(mfi);
            const Array4<Real>& zflux = (faceflux_res[2]).array(mfi);

            if (b.ok()) {
                amrex::ParallelForRNG(b, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
                {

                    ////////////////// from reservoir //////////////////

                    ////////////////////////////////////////////////////
                    // number of particles in FHD cell (for correction)
                    Real N = 0.0;
                    GpuArray<Real,MAX_SPECIES> N_i;
                    for (int n=0;n<nspecies;++n) {
                        N_i[n] = vol*(cons0(i,j,k-1,5+n)/mass[n]);
                        N += N_i[n];
                    }
                    ////////////////////////////////////////////////////

                    Real T = prim0(i,j,k,4);
                    Real Vx, Vy, Vz;
                    Vx = 0.0;
                    Vy = 0.0;
                    Vz = 0.0; // normal
                    Real mass_cross; // total mass crossing at the reservoir interface
                    GpuArray<Real,3> mom_cross; // total momentum crossing at the reservoir interface
                    Real en_cross; // total energy crossing at the reservoir interface
                    GpuArray<Real,MAX_SPECIES> spec_mass_cross;
                    GpuArray<Real,MAX_SPECIES> rhoYk;
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i,j,k,5+n);
                    }
                    if (do_1D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vz,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vx,Vy,engine);
                    }
                    else if (do_2D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vz,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vx,Vy,engine);
                    }
                    else {
                        poisson_process_reservoir(mass,rhoYk,T,Vz,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vx,Vy,engine);
                    }
                    zflux(i,j,k,0) -= (1.0 - (1.0/(12.0*N)))*mass_cross/(dt*area); // update mass flux
                    zflux(i,j,k,4) -= (1.0 + (1.0/( 4.0*N)))*en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        zflux(i,j,k,5+n) -= (1.0 - (1.0/(12.0*N)))*spec_mass_cross[n]/(dt*area); // update species flux
                    }

                    if (do_1D) {
                        zmom(i,j,k) -= mom_cross[0]/vol;
                    }
                    else if (do_2D) {
                        zmom(i,j,k) -= mom_cross[0]/vol;
                        zflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
                    }
                    else {
                        zmom(i,j,k) -= mom_cross[0]/vol;
                        zflux(i,j,k,1) -= mom_cross[1]/dt/area; // tangential flux
                        zflux(i,j,k,2) -= mom_cross[2]/dt/area; // tangential flux
                    }

                    ////////////////////////////////////////////////////

                    ////////////////// to reservoir ////////////////////

                    T = prim0(i,j,k-1,4);
                    Vx = 0.5*(xvel0(i,j,k-1)+xvel0(i+1,j,k-1));
                    Vy = 0.5*(yvel0(i,j,k-1)+yvel0(i,j+1,k-1));
                    Vz = 0.5*(zvel0(i,j,k-1)+zvel0(i,j,k));
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i,j,k-1,5+n);
                    }
                    if (do_1D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vz,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,1,Vx,Vy,engine);
                    }
                    else if (do_2D) {
                        poisson_process_reservoir(mass,rhoYk,T,Vz,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,2,Vx,Vy,engine);
                    }
                    else {
                        poisson_process_reservoir(mass,rhoYk,T,Vz,nspecies,area,k_B,dt,
                                              mass_cross,mom_cross,en_cross,spec_mass_cross,3,Vx,Vy,engine);
                    }
                    zflux(i,j,k,0) += mass_cross/(dt*area); // update mass flux
                    zflux(i,j,k,4) += en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        zflux(i,j,k,5+n) += spec_mass_cross[n]/(dt*area);
                    }

                    if (do_1D) {
                        zmom(i,j,k) += mom_cross[0]/vol;  // add momentum for particles going other way
                    }
                    else if (do_2D) {
                        zmom(i,j,k) += mom_cross[0]/vol;  // add momentum for particles going other way
                        zflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
                    }
                    else {
                        zmom(i,j,k) += mom_cross[0]/vol;  // add momentum for particles going other way
                        zflux(i,j,k,1) += mom_cross[1]/dt/area; // tangential flux
                        zflux(i,j,k,2) += mom_cross[2]/dt/area; // tangential flux
                    }

                    ////////////////////////////////////////////////////

                });
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////
// reset fluxes at the reservoir-FHD interface ////////////////////////////
///////////////////////////////////////////////////////////////////////////
void
ResetReservoirFluxes(const std::array<MultiFab, AMREX_SPACEDIM>& faceflux_res,
                     std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                     std::array< MultiFab, 2 >& edgeflux_x,
                     std::array< MultiFab, 2 >& edgeflux_y,
                     std::array< MultiFab, 2 >& edgeflux_z,
                     const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("ResetReservoirFluxes()",ResetReservoirFluxes);

    // Reservoir in LO X
    if (bc_mass_lo[0] == 4) {

        /////////////// reset energy and mass fluxes ////////////////////////////

        // domain grown nodally based on faceflux_res[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xlo = amrex::bdryNode(dom_x, Orientation(0, Orientation::low));

        for (MFIter mfi(faceflux[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xlo;

            const Array4<Real>& xflux           = (faceflux[0]).array(mfi);
            const Array4<const Real>& xflux_res = (faceflux_res[0]).array(mfi);

            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    xflux(i,j,k,0) = xflux_res(i,j,k,0); // mass flux
                    xflux(i,j,k,4) = xflux_res(i,j,k,4); // energy flux
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) = xflux_res(i,j,k,5+n); // species flux
                    }
                });
            }
        }

        ///////////////////////////////////////////////////////////////////////////

        /////////////// reset tangential momentum fluxes //////////////////////////

        if (do_1D) {
        }
        else { // works both for 2D and 3D

            // domain grown nodally based on edgeflux_x[0] nodality (xy)
            const Box& dom_xy = amrex::convert(geom.Domain(), edgeflux_x[0].ixType());

            // this is the x-lo domain boundary box (xy nodality)
            // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
            const Box& dom_xy_xlo = amrex::bdryNode(dom_xy, Orientation(0, Orientation::low));

            for (MFIter mfi(edgeflux_x[0]); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.fabbox();
                const Box& b = bx & dom_xy_xlo;
                Array4<Real> const& edgex_v = (edgeflux_x[0]).array(mfi);
                const Array4<const Real>& xflux_res = (faceflux_res[0]).array(mfi);
                if (b.ok()) {
                    amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        edgex_v(i,j,k) = 0.5*(xflux_res(i,j,k,1) + xflux_res(i,j-1,k,1));
                    });
                }
            }

            if ((!do_1D) and (!do_2D)) { // z-momentum

                // domain grown nodally based on edgeflux_x[1] nodality (xz)
                const Box& dom_xz = amrex::convert(geom.Domain(), edgeflux_x[1].ixType());

                // this is the x-lo domain boundary box (xz nodality)
                // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
                const Box& dom_xz_xlo = amrex::bdryNode(dom_xz, Orientation(0, Orientation::low));

                for (MFIter mfi(edgeflux_x[1]); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.fabbox();
                    const Box& b = bx & dom_xz_xlo;
                    Array4<Real> const& edgex_w = (edgeflux_x[1]).array(mfi);
                    const Array4<const Real>& xflux_res = (faceflux_res[0]).array(mfi);
                    if (b.ok()) {
                        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            edgex_w(i,j,k) = 0.5*(xflux_res(i,j,k,2) + xflux_res(i,j,k-1,2));
                        });
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////

    }

    // Reservoir in HI X
    if (bc_mass_hi[0] == 4) {

        /////////////// reset energy and mass fluxes ////////////////////////////

        // domain grown nodally based on faceflux_res[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xhi = amrex::bdryNode(dom_x, Orientation(0, Orientation::high));

        for (MFIter mfi(faceflux[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xhi;

            const Array4<Real>& xflux           = (faceflux[0]).array(mfi);
            const Array4<const Real>& xflux_res = (faceflux_res[0]).array(mfi);

            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    xflux(i,j,k,0) = xflux_res(i,j,k,0); // mass flux
                    xflux(i,j,k,4) = xflux_res(i,j,k,4); // energy flux
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) = xflux_res(i,j,k,5+n); // species flux
                    }
                });
            }
        }

        ///////////////////////////////////////////////////////////////////////////

        /////////////// reset tangential momentum fluxes //////////////////////////

        if (do_1D) {
        }
        else { // works both for 2D and 3D

            // domain grown nodally based on edgeflux_x[0] nodality (xy)
            const Box& dom_xy = amrex::convert(geom.Domain(), edgeflux_x[0].ixType());

            // this is the x-lo domain boundary box (xy nodality)
            // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
            const Box& dom_xy_xhi = amrex::bdryNode(dom_xy, Orientation(0, Orientation::high));

            for (MFIter mfi(edgeflux_x[0]); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.fabbox();
                const Box& b = bx & dom_xy_xhi;
                Array4<Real> const& edgex_v = (edgeflux_x[0]).array(mfi);
                const Array4<const Real>& xflux_res = (faceflux_res[0]).array(mfi);
                if (b.ok()) {
                    amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        edgex_v(i,j,k) = 0.5*(xflux_res(i,j,k,1) + xflux_res(i,j-1,k,1));
                    });
                }
            }

            if ((!do_1D) and (!do_2D)) { // z-momentum

                // domain grown nodally based on edgeflux_x[1] nodality (xz)
                const Box& dom_xz = amrex::convert(geom.Domain(), edgeflux_x[1].ixType());

                // this is the x-lo domain boundary box (xz nodality)
                // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
                const Box& dom_xz_xhi = amrex::bdryNode(dom_xz, Orientation(0, Orientation::high));

                for (MFIter mfi(edgeflux_x[1]); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.fabbox();
                    const Box& b = bx & dom_xz_xhi;
                    Array4<Real> const& edgex_w = (edgeflux_x[1]).array(mfi);
                    const Array4<const Real>& xflux_res = (faceflux_res[0]).array(mfi);
                    if (b.ok()) {
                        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            edgex_w(i,j,k) = 0.5*(xflux_res(i,j,k,2) + xflux_res(i,j,k-1,2));
                        });
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////
    }

    // Reservoir in LO Y
    if (bc_mass_lo[1] == 4) {

        // domain grown nodally based on faceflux[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux[1].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_ylo = amrex::bdryNode(dom_y, Orientation(1, Orientation::low));

        for (MFIter mfi(faceflux[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_ylo;

            const Array4<Real>& yflux           = (faceflux[1]).array(mfi);
            const Array4<const Real>& yflux_res = (faceflux_res[1]).array(mfi);

            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    yflux(i,j,k,0) = yflux_res(i,j,k,0); // mass flux
                    yflux(i,j,k,4) = yflux_res(i,j,k,4); // energy flux
                    for (int n=0;n<nspecies;++n) {
                        yflux(i,j,k,5+n) = yflux_res(i,j,k,5+n); // species flux
                    }
                });
            }
        }

        ///////////////////////////////////////////////////////////////////////////

        /////////////// reset tangential momentum fluxes //////////////////////////

        if (do_1D) {
        }
        else { // works both for 2D and 3D

            // domain grown nodally based on edgeflux_x[0] nodality (xy)
            const Box& dom_xy = amrex::convert(geom.Domain(), edgeflux_y[0].ixType());

            // this is the x-lo domain boundary box (xy nodality)
            // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
            const Box& dom_xy_ylo = amrex::bdryNode(dom_xy, Orientation(1, Orientation::low));

            for (MFIter mfi(edgeflux_y[0]); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.fabbox();
                const Box& b = bx & dom_xy_ylo;
                Array4<Real> const& edgey_u = (edgeflux_y[0]).array(mfi);
                const Array4<const Real>& yflux_res = (faceflux_res[1]).array(mfi);
                if (b.ok()) {
                    amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        edgey_u(i,j,k) = 0.5*(yflux_res(i,j,k,1) + yflux_res(i-1,j,k,1));
                    });
                }
            }

            if ((!do_1D) and (!do_2D)) { // z-momentum

                // domain grown nodally based on edgeflux_x[1] nodality (xz)
                const Box& dom_yz = amrex::convert(geom.Domain(), edgeflux_y[1].ixType());

                // this is the x-lo domain boundary box (xz nodality)
                // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
                const Box& dom_yz_ylo = amrex::bdryNode(dom_yz, Orientation(1, Orientation::low));

                for (MFIter mfi(edgeflux_y[1]); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.fabbox();
                    const Box& b = bx & dom_yz_ylo;
                    Array4<Real> const& edgey_w = (edgeflux_y[1]).array(mfi);
                    const Array4<const Real>& yflux_res = (faceflux_res[1]).array(mfi);
                    if (b.ok()) {
                        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            edgey_w(i,j,k) = 0.5*(yflux_res(i,j,k,2) + yflux_res(i,j,k-1,2));
                        });
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////
    }


    // Reservoir in HI Y
    if (bc_mass_hi[1] == 4) {

        // domain grown nodally based on faceflux[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux[1].ixType());

        // this is the y-hi domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yhi = amrex::bdryNode(dom_y, Orientation(1, Orientation::high));

        for (MFIter mfi(faceflux[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yhi;

            const Array4<Real>& yflux           = (faceflux[1]).array(mfi);
            const Array4<const Real>& yflux_res = (faceflux_res[1]).array(mfi);

            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    yflux(i,j,k,0) = yflux_res(i,j,k,0); // mass flux
                    yflux(i,j,k,4) = yflux_res(i,j,k,4); // energy flux
                    for (int n=0;n<nspecies;++n) {
                        yflux(i,j,k,5+n) = yflux_res(i,j,k,5+n); // species flux
                    }
                });
            }
        }

        ///////////////////////////////////////////////////////////////////////////

        /////////////// reset tangential momentum fluxes //////////////////////////

        if (do_1D) {
        }
        else { // works both for 2D and 3D

            // domain grown nodally based on edgeflux_x[0] nodality (xy)
            const Box& dom_xy = amrex::convert(geom.Domain(), edgeflux_y[0].ixType());

            // this is the x-lo domain boundary box (xy nodality)
            // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
            const Box& dom_xy_yhi = amrex::bdryNode(dom_xy, Orientation(1, Orientation::high));

            for (MFIter mfi(edgeflux_y[0]); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.fabbox();
                const Box& b = bx & dom_xy_yhi;
                Array4<Real> const& edgey_u = (edgeflux_y[0]).array(mfi);
                const Array4<const Real>& yflux_res = (faceflux_res[1]).array(mfi);
                if (b.ok()) {
                    amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        edgey_u(i,j,k) = 0.5*(yflux_res(i,j,k,1) + yflux_res(i-1,j,k,1));
                    });
                }
            }

            if ((!do_1D) and (!do_2D)) { // z-momentum

                // domain grown nodally based on edgeflux_x[1] nodality (xz)
                const Box& dom_yz = amrex::convert(geom.Domain(), edgeflux_y[1].ixType());

                // this is the x-lo domain boundary box (xz nodality)
                // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
                const Box& dom_yz_yhi = amrex::bdryNode(dom_yz, Orientation(1, Orientation::high));

                for (MFIter mfi(edgeflux_y[1]); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.fabbox();
                    const Box& b = bx & dom_yz_yhi;
                    Array4<Real> const& edgey_w = (edgeflux_y[1]).array(mfi);
                    const Array4<const Real>& yflux_res = (faceflux_res[1]).array(mfi);
                    if (b.ok()) {
                        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            edgey_w(i,j,k) = 0.5*(yflux_res(i,j,k,2) + yflux_res(i,j,k-1,2));
                        });
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////
    }

    // Reservoir in LO Z
    if (bc_mass_lo[2] == 4) {

        // domain grown nodally based on faceflux[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux[2].ixType());

        // this is the z-lo domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zlo = amrex::bdryNode(dom_z, Orientation(2, Orientation::low));

        for (MFIter mfi(faceflux[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zlo;

            const Array4<Real>& zflux           = (faceflux[2]).array(mfi);
            const Array4<const Real>& zflux_res = (faceflux_res[2]).array(mfi);

            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    zflux(i,j,k,0) = zflux_res(i,j,k,0); // mass flux
                    zflux(i,j,k,4) = zflux_res(i,j,k,4); // energy flux
                    for (int n=0;n<nspecies;++n) {
                        zflux(i,j,k,5+n) = zflux_res(i,j,k,5+n); // species flux
                    }
                });
            }
        }

        ///////////////////////////////////////////////////////////////////////////

        /////////////// reset tangential momentum fluxes //////////////////////////

        if (do_1D) {
        }
        else { // works both for 2D and 3D

            // domain grown nodally based on edgeflux_x[0] nodality (xy)
            const Box& dom_xz = amrex::convert(geom.Domain(), edgeflux_z[0].ixType());

            // this is the x-lo domain boundary box (xy nodality)
            // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
            const Box& dom_xz_zlo = amrex::bdryNode(dom_xz, Orientation(2, Orientation::low));

            for (MFIter mfi(edgeflux_z[0]); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.fabbox();
                const Box& b = bx & dom_xz_zlo;
                Array4<Real> const& edgez_u = (edgeflux_z[0]).array(mfi);
                const Array4<const Real>& zflux_res = (faceflux_res[2]).array(mfi);
                if (b.ok()) {
                    amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        edgez_u(i,j,k) = 0.5*(zflux_res(i,j,k,1) + zflux_res(i-1,j,k,1));
                    });
                }
            }

            if ((!do_1D) and (!do_2D)) { // z-momentum

                // domain grown nodally based on edgeflux_x[1] nodality (xz)
                const Box& dom_yz = amrex::convert(geom.Domain(), edgeflux_z[1].ixType());

                // this is the x-lo domain boundary box (xz nodality)
                // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
                const Box& dom_yz_zlo = amrex::bdryNode(dom_yz, Orientation(2, Orientation::low));

                for (MFIter mfi(edgeflux_z[1]); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.fabbox();
                    const Box& b = bx & dom_yz_zlo;
                    Array4<Real> const& edgez_v = (edgeflux_z[1]).array(mfi);
                    const Array4<const Real>& zflux_res = (faceflux_res[2]).array(mfi);
                    if (b.ok()) {
                        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            edgez_v(i,j,k) = 0.5*(zflux_res(i,j,k,2) + zflux_res(i,j-1,k,2));
                        });
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////
    }

    // Reservoir in HI Z
    if (bc_mass_hi[2] == 4) {

        // domain grown nodally based on faceflux[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux[2].ixType());

        // this is the z-hi domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zhi = amrex::bdryNode(dom_z, Orientation(2, Orientation::high));

        for (MFIter mfi(faceflux[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zhi;

            const Array4<Real>& zflux           = (faceflux[2]).array(mfi);
            const Array4<const Real>& zflux_res = (faceflux_res[2]).array(mfi);

            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    zflux(i,j,k,0) = zflux_res(i,j,k,0); // mass flux
                    zflux(i,j,k,4) = zflux_res(i,j,k,4); // energy flux
                    for (int n=0;n<nspecies;++n) {
                        zflux(i,j,k,5+n) = zflux_res(i,j,k,5+n); // species flux
                    }
                });
            }
        }

        ///////////////////////////////////////////////////////////////////////////

        /////////////// reset tangential momentum fluxes //////////////////////////

        if (do_1D) {
        }
        else { // works both for 2D and 3D

            // domain grown nodally based on edgeflux_x[0] nodality (xy)
            const Box& dom_xz = amrex::convert(geom.Domain(), edgeflux_z[0].ixType());

            // this is the x-lo domain boundary box (xy nodality)
            // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
            const Box& dom_xz_zhi = amrex::bdryNode(dom_xz, Orientation(2, Orientation::high));

            for (MFIter mfi(edgeflux_z[0]); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.fabbox();
                const Box& b = bx & dom_xz_zhi;
                Array4<Real> const& edgez_u = (edgeflux_z[0]).array(mfi);
                const Array4<const Real>& zflux_res = (faceflux_res[2]).array(mfi);
                if (b.ok()) {
                    amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        edgez_u(i,j,k) = 0.5*(zflux_res(i,j,k,1) + zflux_res(i-1,j,k,1));
                    });
                }
            }

            if ((!do_1D) and (!do_2D)) { // z-momentum

                // domain grown nodally based on edgeflux_x[1] nodality (xz)
                const Box& dom_yz = amrex::convert(geom.Domain(), edgeflux_z[1].ixType());

                // this is the x-lo domain boundary box (xz nodality)
                // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
                const Box& dom_yz_zhi = amrex::bdryNode(dom_yz, Orientation(2, Orientation::high));

                for (MFIter mfi(edgeflux_z[1]); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.fabbox();
                    const Box& b = bx & dom_yz_zhi;
                    Array4<Real> const& edgez_v = (edgeflux_z[1]).array(mfi);
                    const Array4<const Real>& zflux_res = (faceflux_res[2]).array(mfi);
                    if (b.ok()) {
                        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                            edgez_v(i,j,k) = 0.5*(zflux_res(i,j,k,2) + zflux_res(i,j-1,k,2));
                        });
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////
    }

}

///////////////////////////////////////////////////////////////////////////
// reset momentum at the reservoir-FHD interface //////////////////////////
///////////////////////////////////////////////////////////////////////////
void
ResetReservoirMom(std::array<MultiFab, AMREX_SPACEDIM>& cumom,
                  const std::array<MultiFab, AMREX_SPACEDIM>& cumom_res,
                  const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("ResetReservoirMom()",ResetReservoirMom);

    // Reservoir in LO X
    if (bc_mass_lo[0] == 4) {

        // domain grown nodally based on cumom_res[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), cumom[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xlo = amrex::bdryNode(dom_x, Orientation(0, Orientation::low));

        for (MFIter mfi(cumom[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xlo;

            const Array4<Real>& xmom           = (cumom[0]).array(mfi);
            const Array4<const Real>& xmom_res = (cumom_res[0]).array(mfi);

            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    xmom(i,j,k) = xmom_res(i,j,k);
                });
            }
        }
    }

    // Reservoir in HI X
    if (bc_mass_hi[0] == 4) {

        // domain grown nodally based on cumom_res[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), cumom[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xhi = amrex::bdryNode(dom_x, Orientation(0, Orientation::high));

        for (MFIter mfi(cumom[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xhi;

            const Array4<Real>& xmom           = (cumom[0]).array(mfi);
            const Array4<const Real>& xmom_res = (cumom_res[0]).array(mfi);

            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    xmom(i,j,k) = xmom_res(i,j,k);
                });
            }
        }
    }

    // Reservoir in LO Y
    if (bc_mass_lo[1] == 4) {

        // domain grown nodally based on cumom[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), cumom[1].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_ylo = amrex::bdryNode(dom_y, Orientation(1, Orientation::low));

        for (MFIter mfi(cumom[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_ylo;

            const Array4<Real>& ymom           = (cumom[1]).array(mfi);
            const Array4<const Real>& ymom_res = (cumom_res[1]).array(mfi);

            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ymom(i,j,k) = ymom_res(i,j,k);
                });
            }
        }
    }

    // Reservoir in HI Y
    if (bc_mass_hi[1] == 4) {

        // domain grown nodally based on cumom[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), cumom[1].ixType());

        // this is the y-hi domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yhi = amrex::bdryNode(dom_y, Orientation(1, Orientation::high));

        for (MFIter mfi(cumom[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yhi;

            const Array4<Real>& ymom           = (cumom[1]).array(mfi);
            const Array4<const Real>& ymom_res = (cumom_res[1]).array(mfi);

            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    ymom(i,j,k) = ymom_res(i,j,k);
                });
            }
        }
    }

    // Reservoir in LO Z
    if (bc_mass_lo[2] == 4) {

        // domain grown nodally based on cumom[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), cumom[2].ixType());

        // this is the z-lo domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zlo = amrex::bdryNode(dom_z, Orientation(2, Orientation::low));

        for (MFIter mfi(cumom[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zlo;

            const Array4<Real>& zmom           = (cumom[2]).array(mfi);
            const Array4<const Real>& zmom_res = (cumom_res[2]).array(mfi);

            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    zmom(i,j,k) = zmom_res(i,j,k);
                });
            }
        }
    }

    // Reservoir in HI Z
    if (bc_mass_hi[2] == 4) {

        // domain grown nodally based on cumom[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), cumom[2].ixType());

        // this is the z-hi domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zhi = amrex::bdryNode(dom_z, Orientation(2, Orientation::high));

        for (MFIter mfi(cumom[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zhi;

            const Array4<Real>& zmom           = (cumom[2]).array(mfi);
            const Array4<const Real>& zmom_res = (cumom_res[2]).array(mfi);

            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    zmom(i,j,k) = zmom_res(i,j,k);
                });
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////
// Reflux conserved qtys at the cell next to reservoir ////////////////////
// not used in the current implementation /////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//void
//ReFluxCons(MultiFab& cu, const MultiFab& cu0,
//           const std::array<MultiFab, AMREX_SPACEDIM>& faceflux_res,
//           const std::array<MultiFab, AMREX_SPACEDIM>& faceflux_cont,
//           const amrex::Geometry& geom,
//           const amrex::Real dt)
//{
//    BL_PROFILE_VAR("ReFluxCons()",ReFluxCons);
//
//    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
//    Box dom(geom.Domain());
//
//    for ( MFIter mfi(cu0); mfi.isValid(); ++mfi) {
//
//        const Box& bx = mfi.validbox();
//
//        AMREX_D_TERM(Array4<Real const> const& xflux_res = faceflux_res[0].array(mfi);,
//                     Array4<Real const> const& yflux_res = faceflux_res[1].array(mfi);,
//                     Array4<Real const> const& zflux_res = faceflux_res[2].array(mfi););
//        AMREX_D_TERM(Array4<Real const> const& xflux_cont = faceflux_cont[0].array(mfi);,
//                     Array4<Real const> const& yflux_cont = faceflux_cont[1].array(mfi);,
//                     Array4<Real const> const& zflux_cont = faceflux_cont[2].array(mfi););
//
//        const Array4<Real>& cons             = cu.array(mfi);
//        const Array4<const Real>& cons0      = cu0.array(mfi);
//
//        // Reservoir in LO X
//        if ((bc_mass_lo[0] == 4) and (bx.smallEnd(0) <= dom.smallEnd(0))) {
//            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//            {
//                if (i == dom.smallEnd(0)) {
//                    cons(i,j,k,0) = cons0(i,j,k,0)
//                      - (dt/dx[0])*(xflux_cont(i,j,k,0) - xflux_res(i,j,k,0)); // correct density
//                    cons(i,j,k,4) = cons0(i,j,k,4)
//                      - (dt/dx[0])*(xflux_cont(i,j,k,4) - xflux_res(i,j,k,4)); // correct en. density
//                    for (int n=0;n<nspecies;++n) {
//                        cons(i,j,k,n+5) = cons0(i,j,k,n+5)
//                      - (dt/dx[0])*(xflux_cont(i,j,k,n+5) - xflux_res(i,j,k,n+5)); // correct species
//                    }
//                }
//            });
//        }
//
//        // Reservoir in HI X
//        if ((bc_mass_hi[0] == 4) and (bx.bigEnd(0) >= dom.bigEnd(0))) {
//            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//            {
//                if (i == dom.bigEnd(0)) {
//                    cons(i,j,k,0) = cons0(i,j,k,0)
//                      + (dt/dx[0])*(xflux_cont(i+1,j,k,0) - xflux_res(i+1,j,k,0)); // correct density
//                    cons(i,j,k,4) = cons0(i,j,k,4)
//                      + (dt/dx[0])*(xflux_cont(i+1,j,k,4) - xflux_res(i+1,j,k,4)); // correct en. density
//                    for (int n=0;n<nspecies;++n) {
//                        cons(i,j,k,n+5) = cons0(i,j,k,n+5)
//                      + (dt/dx[0])*(xflux_cont(i+1,j,k,n+5) - xflux_res(i+1,j,k,n+5)); // correct species
//                    }
//                }
//            });
//        }
//    }
//}

