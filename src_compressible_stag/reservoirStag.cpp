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

    Real N_A = Runiv/k_B; // Avagadro's number
    Real vol = dx[0]*dx[1]*dx[2];
    Real PI = 4.0*atan(1.0);
    Real sqrtPI = sqrt(PI);
    
    GpuArray<Real,MAX_SPECIES> mass;
    for (int l=0;l<nspecies;++l) {
        mass[l] = molmass[l]/(N_A);
    }

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        cumom_res[d].setVal(0.0);
        faceflux_res[d].setVal(0.0);
    }


    // Reservoir in LO X
    if (bc_mass_lo[0] == 4) {

        Real area = dx[1]*dx[2];
        Real vol  = dx[0]*dx[1]*dx[2];

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
            
            const Array4<Real>& xmom  = (cumom_res[0]).array(mfi);
            const Array4<Real>& xflux = (faceflux_res[0]).array(mfi);

            if (b.ok()) {
                amrex::ParallelForRNG(b, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
                {
                    // from reservoir
                    Real T = prim0(i-1,j,k,4);
                    Real V = 0.5*(xvel0(i-1,j,k)+xvel0(i,j,k));
                    V = 0.0;
                    Real mass_cross; // total mass crossing at the reservoir interface
                    Real mom_cross; // total momentum crossing at the reservoir interface
                    Real en_cross; // total energy crossing at the reservoir interface
                    GpuArray<Real,MAX_SPECIES> spec_mass_cross;
                    GpuArray<Real,MAX_SPECIES> rhoYk;
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i-1,j,k,5+n);
                    }
                    poisson_process_reservoir(mass,rhoYk,T,V,nspecies,area,k_B,dt,mass_cross,mom_cross,en_cross,spec_mass_cross,engine);
                    xflux(i,j,k,0) += mass_cross/(dt*area); // update mass flux
                    xflux(i,j,k,4) += en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) += spec_mass_cross[n]/(dt*area); // update species flux
                    }
                    xmom(i,j,k)    += mom_cross/vol; // set face momentum
                    //xmom(i,j,k)    += mass_cross/(dt*area); // set face momentum

                    // to reservoir
                    T = prim0(i,j,k,4);
                    V = -0.5*(xvel0(i,j,k)+xvel0(i+1,j,k));
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i,j,k,5+n);
                    }
                    poisson_process_reservoir(mass,rhoYk,T,V,nspecies,area,k_B,dt,mass_cross,mom_cross,en_cross,spec_mass_cross,engine);
                    xflux(i,j,k,0) -= mass_cross/(dt*area); // update mass flux
                    xflux(i,j,k,4) -= en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) -= spec_mass_cross[n]/(dt*area); // update species flux
                    }
                    xmom(i,j,k)    -= mom_cross/vol; // subtract momentum for particles going other way
                    //xmom(i,j,k)    -= mass_cross/(dt*area); // subtract momentum for particles going other way
                });
            }
        }
    }
    
    // Reservoir in HI X
    if (bc_mass_hi[0] == 4) {

        Real area = dx[1]*dx[2];
        Real vol  = dx[0]*dx[1]*dx[2];

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
            
            const Array4<Real>& xmom  = (cumom_res[0]).array(mfi);
            const Array4<Real>& xflux = (faceflux_res[0]).array(mfi);

            if (b.ok()) {
                amrex::ParallelForRNG(b, [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
                {
                    // from reservoir
                    Real T = prim0(i,j,k,4);
                    Real V = -0.5*(xvel0(i,j,k)+xvel0(i+1,j,k));
                    V = 0.0;
                    Real mass_cross; // total mass crossing at the reservoir interface
                    Real mom_cross; // total momentum crossing at the reservoir interface
                    Real en_cross; // total energy crossing at the reservoir interface
                    GpuArray<Real,MAX_SPECIES> spec_mass_cross;
                    GpuArray<Real,MAX_SPECIES> rhoYk;
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i,j,k,5+n);
                    }
                    poisson_process_reservoir(mass,rhoYk,T,V,nspecies,area,k_B,dt,mass_cross,mom_cross,en_cross,spec_mass_cross,engine);
                    xflux(i,j,k,0) -= mass_cross/(dt*area); // update mass flux
                    xflux(i,j,k,4) -= en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) -= spec_mass_cross[n]/(dt*area); // update species flux
                    }
                    xmom(i,j,k)    -= mom_cross/vol;
                    //xmom(i,j,k)    -= mass_cross/(dt*area);

                    // to reservoir
                    T = prim0(i-1,j,k,4);
                    V = 0.5*(xvel0(i-1,j,k)+xvel0(i,j,k));
                    for (int n=0;n<nspecies;++n) {
                        rhoYk[n] = cons0(i-1,j,k,5+n);
                    }
                    poisson_process_reservoir(mass,rhoYk,T,V,nspecies,area,k_B,dt,mass_cross,mom_cross,en_cross,spec_mass_cross,engine);
                    xflux(i,j,k,0) += mass_cross/(dt*area); // update mass flux
                    xflux(i,j,k,4) += en_cross/(dt*area); // update energy flux
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) += spec_mass_cross[n]/(dt*area);
                    }
                    xmom(i,j,k)    += mom_cross/vol;  // add momentum for particles going other way
                    //xmom(i,j,k)    += mass_cross/(dt*area);  // add momentum for particles going other way
                });
            }
        }
    }

//    faceflux_res[0].OverrideSync(geom.periodicity());
//    faceflux_res[1].OverrideSync(geom.periodicity());
//    faceflux_res[2].OverrideSync(geom.periodicity());
//    cumom_res[0].OverrideSync(geom.periodicity());
//    cumom_res[1].OverrideSync(geom.periodicity());
//    cumom_res[2].OverrideSync(geom.periodicity());
}

///////////////////////////////////////////////////////////////////////////
// reset fluxes at the reservoir-FHD interface ////////////////////////////
///////////////////////////////////////////////////////////////////////////
void 
ResetReservoirFluxes(std::array<MultiFab, AMREX_SPACEDIM>& faceflux,
                     const std::array<MultiFab, AMREX_SPACEDIM>& faceflux_res,
                     const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("ResetReservoirFluxes()",ResetReservoirFluxes);

    // Reservoir in LO X
    if (bc_mass_lo[0] == 4) {

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
    }

    // Reservoir in HI X
    if (bc_mass_hi[0] == 4) {

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
///////////////////////////////////////////////////////////////////////////
void 
ReFluxCons(MultiFab& cu, const MultiFab& cu0,
           const std::array<MultiFab, AMREX_SPACEDIM>& faceflux_res,
           const std::array<MultiFab, AMREX_SPACEDIM>& faceflux_cont,
           const amrex::Geometry& geom,
           const amrex::Real dt)
{
    BL_PROFILE_VAR("ReFluxCons()",ReFluxCons);
    
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    Box dom(geom.Domain());
    
    for ( MFIter mfi(cu0); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();

        AMREX_D_TERM(Array4<Real const> const& xflux_res = faceflux_res[0].array(mfi);,
                     Array4<Real const> const& yflux_res = faceflux_res[1].array(mfi);,
                     Array4<Real const> const& zflux_res = faceflux_res[2].array(mfi););
        AMREX_D_TERM(Array4<Real const> const& xflux_cont = faceflux_cont[0].array(mfi);,
                     Array4<Real const> const& yflux_cont = faceflux_cont[1].array(mfi);,
                     Array4<Real const> const& zflux_cont = faceflux_cont[2].array(mfi););

        const Array4<Real>& cons             = cu.array(mfi);
        const Array4<const Real>& cons0      = cu0.array(mfi);

        // Reservoir in LO X
        if ((bc_mass_lo[0] == 4) and (bx.smallEnd(0) <= dom.smallEnd(0))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == dom.smallEnd(0)) {
                    cons(i,j,k,0) = cons0(i,j,k,0) 
                      - (dt/dx[0])*(xflux_cont(i,j,k,0) - xflux_res(i,j,k,0)); // correct density
                    cons(i,j,k,4) = cons0(i,j,k,4) 
                      - (dt/dx[0])*(xflux_cont(i,j,k,4) - xflux_res(i,j,k,4)); // correct en. density
                    for (int n=0;n<nspecies;++n) {
                        cons(i,j,k,n+5) = cons0(i,j,k,n+5) 
                      - (dt/dx[0])*(xflux_cont(i,j,k,n+5) - xflux_res(i,j,k,n+5)); // correct species
                    }
                }
            });
        }

        // Reservoir in HI X
        if ((bc_mass_hi[0] == 4) and (bx.bigEnd(0) >= dom.bigEnd(0))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == dom.bigEnd(0)) {
                    cons(i,j,k,0) = cons0(i,j,k,0) 
                      + (dt/dx[0])*(xflux_cont(i+1,j,k,0) - xflux_res(i+1,j,k,0)); // correct density
                    cons(i,j,k,4) = cons0(i,j,k,4) 
                      + (dt/dx[0])*(xflux_cont(i+1,j,k,4) - xflux_res(i+1,j,k,4)); // correct en. density
                    for (int n=0;n<nspecies;++n) {
                        cons(i,j,k,n+5) = cons0(i,j,k,n+5) 
                      + (dt/dx[0])*(xflux_cont(i+1,j,k,n+5) - xflux_res(i+1,j,k,n+5)); // correct species
                    }
                }
            });
        }
    }
}

