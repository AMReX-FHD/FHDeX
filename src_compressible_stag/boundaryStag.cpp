#include "compressible_functions.H"
#include "compressible_functions_stag.H"
#include "common_functions.H"

void SetupBCStag() {
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_vel_lo[i] == -1) {
            bc_mass_lo[i] = -1;
            bc_mass_hi[i] = -1;
            bc_therm_lo[i] = -1;
            bc_therm_hi[i] = -1;
        }
    }
}

void SetupCWallStag() {

    Real sumx, sumy;

    // Compute Xk or Yk at the wall, depending on which is defined
    // For reservoirs, also compute pressure in the reservoir depending on reservoir t, rho, Yk
    // X walls
    if ((bc_mass_lo[0] == 2) or (bc_mass_lo[0] >= 3)) {
       sumx = 0.;
       sumy = 0.;
       for (int ns=0; ns<nspecies; ++ns) {
          sumx = sumx + bc_Xk_x_lo[ns];
          sumy = sumy + bc_Yk_x_lo[ns];
       }
       if (amrex::Math::abs(sumx-1) < 1.e-10) {
           GetMassfrac(bc_Xk_x_lo,bc_Yk_x_lo);;
       } else if (amrex::Math::abs(sumy-1) < 1.e-10) {
           GetMolfrac(bc_Yk_x_lo,bc_Xk_x_lo);
       }
       else {
           Abort("SetupCWallStag: lo-x; mass or mole fractions do not sum to 1");
       }
    }

    if (bc_mass_lo[0] >= 3) {
        // set reservoir pressure equal to inital ambient pressure (for no flow)
        // if t_lo/hi is positive, compute rho_lo/hi (default)
        // if rho_lo/hi is positive, rewrite t_lo/hi (from input script)

        if (p_lo[0] <= 0.0) { // set reservoir pressure to ambient if not specified
            GpuArray<Real,MAX_SPECIES> massvec;
            for (int ns=0;ns<nspecies;++ns) massvec[ns] = rhobar[ns];
            GetPressureGas(p_lo[0],massvec,rho0,T_init[0]);
        }

        if (rho_lo[0] <= 0.0) { // specify reservoir density if not specified
            GetDensity(p_lo[0],rho_lo[0],t_lo[0],bc_Yk_x_lo);
        }
        else if (t_lo[0] <= 0.0) { // specify reservoir temperature if not specified
            Real molmix = 0.;
            for (int n=0; n<nspecies; ++n) {
                molmix += bc_Yk_x_lo[n]/molmass[n];
            }
            molmix = 1./molmix;
            t_lo[0] = p_lo[0]*(molmix/Runiv)/rho_lo[0];
        }
    }

    if ((bc_mass_hi[0] == 2) or (bc_mass_hi[0] >= 3)) {
       sumx = 0.;
       sumy = 0.;
       for (int ns=0; ns<nspecies; ++ns) {
          sumx = sumx + bc_Xk_x_hi[ns];
          sumy = sumy + bc_Yk_x_hi[ns];
       }
       if (amrex::Math::abs(sumx-1) < 1.e-10) {
          GetMassfrac(bc_Xk_x_hi,bc_Yk_x_hi);
       } else if (amrex::Math::abs(sumy-1) < 1.e-10) {
          GetMolfrac(bc_Yk_x_hi,bc_Xk_x_hi);
       } else {
           Abort("SetupCWallStag: hi-x; mass or mole fractions do not sum to 1");
       }
    }

    if (bc_mass_hi[0] >= 3) {
        // set reservoir pressure equal to inital ambient pressure (for no flow)
        // if t_lo/hi is positive, compute rho_lo/hi (default)
        // if rho_lo/hi is positive, rewrite t_lo/hi (from input script)

        if (p_hi[0] <= 0.0) {  // set reservoir pressure to ambient if not specified
            GpuArray<Real,MAX_SPECIES> massvec;
            for (int ns=0;ns<nspecies;++ns) massvec[ns] = rhobar[ns];
            GetPressureGas(p_hi[0],massvec,rho0,T_init[0]);
        }

        if (rho_hi[0] <= 0.0) { // specify reservoir density  if not specified
            GetDensity(p_hi[0],rho_hi[0],t_hi[0],bc_Yk_x_hi);
        }
        else if (t_hi[0] <= 0.0) { // specify reservoir temperature if not specified
            Real molmix = 0.;
            for (int n=0; n<nspecies; ++n) {
                molmix += bc_Yk_x_hi[n]/molmass[n];
            }
            molmix = 1./molmix;
            t_hi[0] = p_hi[0]*(molmix/Runiv)/rho_hi[0];
        }
    }

    // Y walls
    if ((bc_mass_lo[1] == 2) or (bc_mass_lo[1] >= 3)) {
       sumx = 0.;
       sumy = 0.;
       for (int ns=0; ns<nspecies; ++ns) {
          sumx = sumx + bc_Xk_y_lo[ns];
          sumy = sumy + bc_Yk_y_lo[ns];
       }
       if (amrex::Math::abs(sumx-1) < 1.e-10) {
          GetMassfrac(bc_Xk_y_lo,bc_Yk_y_lo);
       } else if (amrex::Math::abs(sumy-1) < 1.e-10) {
          GetMolfrac(bc_Yk_y_lo,bc_Xk_y_lo);
       } else {
           Abort("SetupCWallStag: lo-y; mass or mole fractions do not sum to 1");
       }
    }

    if (bc_mass_lo[1] >= 3) {
        // set reservoir pressure equal to inital ambient pressure (for no flow)
        // if t_lo/hi is positive, compute rho_lo/hi (default)
        // if rho_lo/hi is positive, rewrite t_lo/hi (from input script)

        if (p_lo[1] <= 0.0) {  // set reservoir pressure to ambient if not specified
            GpuArray<Real,MAX_SPECIES> massvec;
            for (int ns=0;ns<nspecies;++ns) massvec[ns] = rhobar[ns];
            GetPressureGas(p_lo[1],massvec,rho0,T_init[0]);
        }

        if (rho_lo[1] <= 0.0) { // specify reservoir density  if not specified
            GetDensity(p_lo[1],rho_lo[1],t_lo[1],bc_Yk_y_lo);
        }
        else if (t_lo[1] <= 0.0) { // specify reservoir temperature if not specified
            Real molmix = 0.;
            for (int n=0; n<nspecies; ++n) {
                molmix += bc_Yk_y_lo[n]/molmass[n];
            }
            molmix = 1./molmix;
            t_lo[1] = p_lo[1]*(molmix/Runiv)/rho_lo[1];
        }
    }

    if ((bc_mass_hi[1] == 2) or (bc_mass_hi[1] >= 3)) {
       sumx = 0.;
       sumy = 0.;
       for (int ns=0; ns<nspecies; ++ns) {
          sumx = sumx + bc_Xk_y_hi[ns];
          sumy = sumy + bc_Yk_y_hi[ns];
       }
       if (amrex::Math::abs(sumx-1) < 1.e-10) {
          GetMassfrac(bc_Xk_y_hi,bc_Yk_y_hi);
       } else if (amrex::Math::abs(sumy-1) < 1.e-10) {
          GetMolfrac(bc_Yk_y_hi,bc_Xk_y_hi);
       } else {
           Abort("SetupCWallStag: hi-y; mass or mole fractions do not sum to 1");
       }
    }

    if (bc_mass_hi[1] >= 3) {
        // set reservoir pressure equal to inital ambient pressure (for no flow)
        // if t_lo/hi is positive, compute rho_lo/hi (default)
        // if rho_lo/hi is positive, rewrite t_lo/hi (from input script)

        if (p_hi[1] <= 0.0) {  // set reservoir pressure to ambient if not specified
            GpuArray<Real,MAX_SPECIES> massvec;
            for (int ns=0;ns<nspecies;++ns) massvec[ns] = rhobar[ns];
            GetPressureGas(p_hi[1],massvec,rho0,T_init[0]);
        }

        if (rho_hi[1] <= 0.0) { // specify reservoir density  if not specified
            GetDensity(p_hi[1],rho_hi[1],t_hi[1],bc_Yk_y_hi);
        }
        else if (t_hi[1] <= 0.0) { // specify reservoir temperature if not specified
            Real molmix = 0.;
            for (int n=0; n<nspecies; ++n) {
                molmix += bc_Yk_y_hi[n]/molmass[n];
            }
            molmix = 1./molmix;
            t_hi[1] = p_hi[1]*(molmix/Runiv)/rho_hi[1];
        }
    }

    // Z walls
    if ((bc_mass_lo[2] == 2) or (bc_mass_lo[2] >= 3)) {
       sumx = 0.;
       sumy = 0.;
       for (int ns=0; ns<nspecies; ++ns) {
          sumx = sumx + bc_Xk_z_lo[ns];
          sumy = sumy + bc_Yk_z_lo[ns];
       }
       if (amrex::Math::abs(sumx-1) < 1.e-10) {
          GetMassfrac(bc_Xk_z_lo,bc_Yk_z_lo);
       } else if (amrex::Math::abs(sumy-1) < 1.e-10) {
          GetMolfrac(bc_Yk_z_lo,bc_Xk_z_lo);
       } else {
           Abort("SetupCWallStag: lo-z; mass or mole fractions do not sum to 1");
       }
    }

    if (bc_mass_lo[2] >= 3) {
        // set reservoir pressure equal to inital ambient pressure (for no flow)
        // if t_lo/hi is positive, compute rho_lo/hi (default)
        // if rho_lo/hi is positive, rewrite t_lo/hi (from input script)

        if (p_lo[2] <= 0.0) {  // set reservoir pressure to ambient if not specified
            GpuArray<Real,MAX_SPECIES> massvec;
            for (int ns=0;ns<nspecies;++ns) massvec[ns] = rhobar[ns];
            GetPressureGas(p_lo[2],massvec,rho0,T_init[0]);
        }

        if (rho_lo[2] <= 0.0) { // specify reservoir density  if not specified
            GetDensity(p_lo[2],rho_lo[2],t_lo[2],bc_Yk_z_lo);
        }
        else if (t_lo[2] <= 0.0) { // specify reservoir temperature if not specified
            Real molmix = 0.;
            for (int n=0; n<nspecies; ++n) {
                molmix += bc_Yk_z_lo[n]/molmass[n];
            }
            molmix = 1./molmix;
            t_lo[2] = p_lo[2]*(molmix/Runiv)/rho_lo[2];
        }
    }

    if ((bc_mass_hi[2] == 2) or (bc_mass_hi[2] >= 3)) {
       sumx = 0.;
       sumy = 0.;
       for (int ns=0; ns<nspecies; ++ns) {
          sumx = sumx + bc_Xk_z_hi[ns];
          sumy = sumy + bc_Yk_z_hi[ns];
       }
       if (amrex::Math::abs(sumx-1) < 1.e-10) {
          GetMassfrac(bc_Xk_z_hi,bc_Yk_z_hi);
       } else if (amrex::Math::abs(sumy-1) < 1.e-10) {
          GetMolfrac(bc_Yk_z_hi,bc_Xk_z_hi);
       } else {
           Abort("SetupCWallStag: hi-z; mass or mole fractions do not sum to 1");
       }
    }

    if (bc_mass_hi[2] >= 3) {
        // set reservoir pressure equal to inital ambient pressure (for no flow)
        // if t_lo/hi is positive, compute rho_lo/hi (default)
        // if rho_lo/hi is positive, rewrite t_lo/hi (from input script)

        if (p_hi[2] <= 0.0) {  // set reservoir pressure to ambient if not specified
            GpuArray<Real,MAX_SPECIES> massvec;
            for (int ns=0;ns<nspecies;++ns) massvec[ns] = rhobar[ns];
            GetPressureGas(p_hi[2],massvec,rho0,T_init[0]);
        }

        if (rho_hi[2] <= 0.0) { // specify reservoir density  if not specified
            GetDensity(p_hi[2],rho_hi[2],t_hi[2],bc_Yk_z_hi);
        }
        else if (t_hi[2] <= 0.0) { // specify reservoir temperature if not specified
            Real molmix = 0.;
            for (int n=0; n<nspecies; ++n) {
                molmix += bc_Yk_z_hi[n]/molmass[n];
            }
            molmix = 1./molmix;
            t_hi[2] = p_hi[2]*(molmix/Runiv)/rho_hi[2];
        }
    }
}

// Set boundary and ghost cells for staggered compressible code based on BCs
void setBCStag(MultiFab& prim_in, MultiFab& cons_in,
                 std::array< MultiFab, AMREX_SPACEDIM >& cumom_in,
                 std::array< MultiFab, AMREX_SPACEDIM >& vel_in,
                 const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("setBCStag()",setBCStag);

    int ng_c = cons_in.nGrow();
    int ng_p = prim_in.nGrow();
    if (ng_c != ng_p) {
        Abort("setBC: prim and cons need the same number of ghost cells");
    }
    int ng_m = cumom_in[0].nGrow();
    int ng_v = vel_in[0].nGrow();
    if (ng_m != ng_v) {
        Abort("setBC: momentum and velocity need the same number of ghost cells");
    }

    if (membrane_cell >= 0) { // set adiabatic slip BC at the membrane
        BCMem(prim_in, cons_in, cumom_in, vel_in, geom);
    }

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        BCMassTempPress(prim_in, cons_in, geom, i);
        BCMomNormal(cumom_in[i], vel_in[i], cons_in, geom, i);
        BCMomTrans(cumom_in[i], vel_in[i], geom, i);
    }
    BCRhoRhoE(cons_in, prim_in, cumom_in, geom);
}

// set species and total density flux to zero for wall boundary conditions
// set the diffusive momentum flux to zero in the reservoir cells
void BCWallReservoirFluxStag(std::array< MultiFab, AMREX_SPACEDIM >& faceflux,
                              std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in,
                             const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("BCWallReservoirFluxStag()",BCWallReservoirFluxStag);

    // Wall BC: set species and total density flux to zero for wall boundary conditions
    // LO X
    if (bc_mass_lo[0] == 1) {

        // domain grown nodally based on faceflux[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xlo = amrex::bdryNode(dom_x, Orientation(0, Orientation::low));

        for (MFIter mfi(faceflux[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xlo;
            Array4<Real> const& flux = (faceflux[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) = 0.;
                    }
                    // density
                    flux(i,j,k,0) = 0.;
                    // Dufour
                    flux(i,j,k,nvars+3) = 0.;
                });
            }
        }
    }
    // HI X
    if (bc_mass_hi[0] == 1) {

        // domain grown nodally based on faceflux[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux[0].ixType());

        // this is the x-hi domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xhi = amrex::bdryNode(dom_x, Orientation(0, Orientation::high));

        for (MFIter mfi(faceflux[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xhi;
            Array4<Real> const& flux = (faceflux[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) = 0.;
                    }
                    // density
                    flux(i,j,k,0) = 0.;
                    // Dufour
                    flux(i,j,k,nvars+3) = 0.;
                });
            }
        }
    }
    // LO Y
    if (bc_mass_lo[1] == 1) {

        // domain grown nodally based on faceflux[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux[1].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_ylo = amrex::bdryNode(dom_y, Orientation(1, Orientation::low));

        for (MFIter mfi(faceflux[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_ylo;
            Array4<Real> const& flux = (faceflux[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) = 0.;
                    }
                    // density
                    flux(i,j,k,0) = 0.;
                    // Dufour
                    flux(i,j,k,nvars+3) = 0.;
                });
            }
        }
    }
    // HI Y
    if (bc_mass_hi[1] == 1) {

        // domain grown nodally based on faceflux[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux[1].ixType());

        // this is the y-hi domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yhi = amrex::bdryNode(dom_y, Orientation(1, Orientation::high));

        for (MFIter mfi(faceflux[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yhi;
            Array4<Real> const& flux = (faceflux[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) = 0.;
                    }
                    // density
                    flux(i,j,k,0) = 0.;
                    // Dufour
                    flux(i,j,k,nvars+3) = 0.;
                });
            }
        }
    }
    // LO Z
    if (bc_mass_lo[2] == 1) {

        // domain grown nodally based on faceflux[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux[2].ixType());

        // this is the z-lo domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zlo = amrex::bdryNode(dom_z, Orientation(2, Orientation::low));

        for (MFIter mfi(faceflux[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zlo;
            Array4<Real> const& flux = (faceflux[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) = 0.;
                    }
                    // density
                    flux(i,j,k,0) = 0.;
                    // Dufour
                    flux(i,j,k,nvars+3) = 0.;
                });
            }
        }
    }
    // HI Z
    if (bc_mass_hi[2] == 1) {

        // domain grown nodally based on faceflux[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux[2].ixType());

        // this is the z-hi domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zhi = amrex::bdryNode(dom_z, Orientation(2, Orientation::high));

        for (MFIter mfi(faceflux[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zhi;
            Array4<Real> const& flux = (faceflux[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) = 0.;
                    }
                    // density
                    flux(i,j,k,0) = 0.;
                    // Dufour
                    flux(i,j,k,nvars+3) = 0.;
                });
            }
        }
    }

    // RESERVOIR BC: set the diffusive momentum flux to zero in the reservoir cells
    Box dom(geom.Domain());
    for ( MFIter mfi(cenflux_in[0]); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(1);

        const Array4<Real>& cenx = cenflux_in[0].array(mfi);
        const Array4<Real>& ceny = cenflux_in[1].array(mfi);
        const Array4<Real>& cenz = cenflux_in[2].array(mfi);

        // LO X
        if ((bc_mass_lo[0] == 3) && (bx.smallEnd(0) < dom.smallEnd(0))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) cenx(i,j,k) = 0.0;
            });
        }
        if ((bc_mass_lo[0] == 4) && (bx.smallEnd(0) < dom.smallEnd(0))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) cenx(i,j,k) = cenx(i+1,j,k);
            });
        }

        // HI X
        if ((bc_mass_hi[0] == 3) && (bx.bigEnd(0) > dom.bigEnd(0))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > dom.bigEnd(0)) cenx(i,j,k) = 0.0;
            });
        }
        if ((bc_mass_hi[0] == 4) && (bx.bigEnd(0) > dom.bigEnd(0))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > dom.bigEnd(0)) cenx(i,j,k) = cenx(i-1,j,k);
            });
        }

        // LO Y
        if ((bc_mass_lo[1] == 3) && (bx.smallEnd(1) < dom.smallEnd(1))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) ceny(i,j,k) = 0.0;
            });
        }
        if ((bc_mass_lo[1] == 4) && (bx.smallEnd(1) < dom.smallEnd(1))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) ceny(i,j,k) = ceny(i,j+1,k);
            });
        }

        // HI Y
        if ((bc_mass_hi[1] == 3) && (bx.bigEnd(1) > dom.bigEnd(1))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j > dom.bigEnd(1)) ceny(i,j,k) = 0.0;
            });
        }
        if ((bc_mass_hi[1] == 4) && (bx.bigEnd(1) > dom.bigEnd(1))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j > dom.bigEnd(1)) ceny(i,j,k) = ceny(i,j-1,k);
            });
        }

        // LO Z
        if ((bc_mass_lo[2] == 3) && (bx.smallEnd(2) < dom.smallEnd(2))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) cenz(i,j,k) = 0.0;
            });
        }
        if ((bc_mass_lo[2] == 4) && (bx.smallEnd(2) < dom.smallEnd(2))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) cenz(i,j,k) = cenz(i,j,k+1);;
            });
        }

        // HI Z
        if ((bc_mass_hi[2] == 3) && (bx.bigEnd(2) > dom.bigEnd(2))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k > dom.bigEnd(2)) cenz(i,j,k) = 0.0;
            });
        }
        if ((bc_mass_hi[2] == 4) && (bx.bigEnd(2) > dom.bigEnd(2))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k > dom.bigEnd(2)) cenz(i,j,k) = cenz(i,j,k-1);;
            });
        }
    }

}

// Set adiabatic slip boundary condition at the membrane
void BCMem(MultiFab& prim_in, MultiFab& cons_in,
           std::array< MultiFab, AMREX_SPACEDIM >& cumom_in,
           std::array< MultiFab, AMREX_SPACEDIM >& vel_in,
           const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("BCMem()",BCMem);

    Box dom(geom.Domain());
    int ng_p = prim_in.nGrow();

    // first set adiabatic temperature and pressure, and a wall
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng_p);
        int lo = bx.smallEnd(0) + ng_p;
        int hi = bx.bigEnd(0) - ng_p ;

        const Array4<Real>& prim = prim_in.array(mfi);

        // membrane at the left end (cell to the right of the membrane)
        if (lo == membrane_cell) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < lo) {
                    prim(i,j,k,4) = prim(2*lo-i-1,j,k,4);
                    prim(i,j,k,5) = prim(2*lo-i-1,j,k,5);
                    for (int n=6; n<nprimvars; ++n) {
                        prim(i,j,k,n) = prim(2*lo-i-1,j,k,n);
                    }
                }
            });
        }

        // membrane at the right end (cell to the left of the membrane)
        else if (hi == membrane_cell - 1) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > hi) {
                    prim(i,j,k,4) = prim(2*hi-i+1,j,k,4);
                    prim(i,j,k,5) = prim(2*hi-i+1,j,k,5);
                    for (int n=6; n<nprimvars; ++n) {
                        prim(i,j,k,n) = prim(2*hi-i+1,j,k,n);
                    }
                }
            });
        }
    }

    // next set normal velocity and momentum (for zero normal flux)
    for ( MFIter mfi(vel_in[0]); mfi.isValid(); ++mfi) {

        int ngv = vel_in[0].nGrow();
        const Box& bx = mfi.growntilebox(ngv);
        int lo = bx.smallEnd(0) + ngv;
        int hi = bx.bigEnd(0) - ngv;

        const Array4<Real>& vel = vel_in[0].array(mfi);
        const Array4<Real>& mom = cumom_in[0].array(mfi);

        // membrane at the left end (cell to the right of the membrane)
        if (lo == membrane_cell) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < lo) {
                    vel(i,j,k) = -vel(2*lo-i,j,k);
                    mom(i,j,k) = -mom(2*lo-i,j,k);
                }
                else if (i == lo) {
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }

        // membrane at the right end (cell to the left of the membrane)
        if (hi == membrane_cell) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > hi) {
                    vel(i,j,k) = -vel(2*hi-i,j,k);
                    mom(i,j,k) = -mom(2*hi-i,j,k);
                }
                else if (i == hi) {
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }
    }

    // next set Y tangential velocity and momentum (for slip BC)
    for ( MFIter mfi(vel_in[1]); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(vel_in[1].nGrow());

        const Array4<Real>& vel = vel_in[1].array(mfi);
        const Array4<Real>& mom = cumom_in[1].array(mfi);

        // membrane at the left end (cell to the right of the membrane)
        if (bx.smallEnd(0) == membrane_cell) {

            int lo = bx.smallEnd(0);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < bx.smallEnd(0)) {
                    vel(i,j,k) = vel(2*lo-i-1,j,k);
                    mom(i,j,k) = mom(2*lo-i-1,j,k);
                }
            });
        }

        // membrane at the right end (cell to the left of the membrane)
        else if (bx.bigEnd(0) == membrane_cell - 1) {

            int hi = bx.bigEnd(0);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > bx.bigEnd(0)) {
                    vel(i,j,k) = vel(2*hi-i+1,j,k);
                    mom(i,j,k) = mom(2*hi-i+1,j,k);
                }
            });
        }
    }

    // next set Z tangential velocity and momentum (for slip BC)
    for ( MFIter mfi(vel_in[2]); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(vel_in[2].nGrow());

        const Array4<Real>& vel = vel_in[2].array(mfi);
        const Array4<Real>& mom = cumom_in[2].array(mfi);

        // membrane at the left end (cell to the right of the membrane)
        if (bx.smallEnd(0) == membrane_cell) {

            int lo = bx.smallEnd(0);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < bx.smallEnd(0)) {
                    vel(i,j,k) = vel(2*lo-i-1,j,k);
                    mom(i,j,k) = mom(2*lo-i-1,j,k);
                }
            });
        }

        // membrane at the right end (cell to the left of the membrane)
        else if (bx.bigEnd(0) == membrane_cell - 1) {

            int hi = bx.bigEnd(0);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > bx.bigEnd(0)) {
                    vel(i,j,k) = vel(2*hi-i+1,j,k);
                    mom(i,j,k) = mom(2*hi-i+1,j,k);
                }
            });
        }
    }

    // finally set densities and energies for slip BC at the membrane
    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng_p);
        int lo = bx.smallEnd(0) + ng_p;
        int hi = bx.bigEnd(0) - ng_p ;

        const Array4<Real>& prim = prim_in.array(mfi);
        const Array4<Real>& cons = cons_in.array(mfi);
        AMREX_D_TERM(Array4<Real const> const& momx = cumom_in[0].array(mfi);,
                     Array4<Real const> const& momy = cumom_in[1].array(mfi);,
                     Array4<Real const> const& momz = cumom_in[2].array(mfi););


        // membrane at the left end (cell to the right of the membrane)
        if (lo == membrane_cell) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < lo) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec);
                    GetEnergy(intenergy,fracvec,temp);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);

                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }

        // membrane at the right end (cell to the left of the membrane)
        else if (hi == membrane_cell - 1) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > hi) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec);
                    GetEnergy(intenergy,fracvec,temp);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);

                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }
    }

}

// Set mass, pressure and temperature on ghost cells based on BCs
void BCMassTempPress(MultiFab& prim_in,MultiFab& cons_in,const amrex::Geometry& geom,int dim)
{
    BL_PROFILE_VAR("BCMassTempPress()",BCMassTempPress);

    Box dom(geom.Domain());
    int ng_p = prim_in.nGrow();

    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng_p);

        const Array4<Real>& prim = prim_in.array(mfi);
        const Array4<Real>& cons = cons_in.array(mfi);

        // LO X
        if ((dim == 0) && (bx.smallEnd(0) < dom.smallEnd(0))) {

            int lo = dom.smallEnd(0);

            // mass fractions, wall
            if ( bc_mass_lo[0] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(2*lo-i-1,j,k,n);
                        }
                    }
                });
            }

            // mass fractions, concentration (set ghost equal to wall value)
            if (bc_mass_lo[0] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_x_lo[n];
                            prim(i,j,k,6+nspecies+n) = bc_Xk_x_lo[n];
                        }
                    }
                });
            }

            // mass fractions, density, temperature and pressure in the reservoir
            if (bc_mass_lo[0] >= 3 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_x_lo[n]; // set ghost cell equal to reservoir mass fraction
                            prim(i,j,k,6+nspecies+n) = bc_Xk_x_lo[n]; // set ghost cell equal to reservoir mole fraction

                            prim(i,j,k,0) = rho_lo[0]; // set ghost cell equal to reservoir density
                            cons(i,j,k,0) = rho_lo[0]; // set ghost cell equal to reservoir density

                            prim(i,j,k,4) = t_lo[0]; // set ghost cell equal to reservoir temperature
                            prim(i,j,k,5) = p_lo[0]; // set ghost cell equal to reservoir pressure
                        }
                    }

                });
            }

            // temperature and pressure, adiabatic
            else if (bc_therm_lo[0] == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(2*lo-i-1,j,k,n);
                        }
                    }
                });
            } else if (bc_therm_lo[0] == 2) { // isothermal  (set ghost equal to wall value)
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i < dom.smallEnd(0)) {
                        prim(i,j,k,4) = t_lo[0];
                        prim(i,j,k,5) = prim(2*lo-i-1,j,k,5);
                    }
                });
            }

        } // end LO X

        // HI X
        if ((dim == 0) && (bx.bigEnd(0) > dom.bigEnd(0))) {

            int hi = dom.bigEnd(0);

            // mass fractions, wall
            if ( bc_mass_hi[0] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(2*hi-i+1,j,k,n);
                        }
                    }
                });
            }

            // mass fractions, concentration  (set ghost equal to wall value)
            if (bc_mass_hi[0] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_x_hi[n];
                            prim(i,j,k,6+nspecies+n) = bc_Xk_x_hi[n];
                        }
                    }
                });
            }

            // mass fractions, density, temperature and pressure in the reservoir
            if (bc_mass_hi[0] >= 3 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_x_hi[n]; // set ghost cell equal to reservoir mass fraction
                            prim(i,j,k,6+nspecies+n) = bc_Xk_x_hi[n]; // set ghost cell equal to reservoir mole fraction

                            prim(i,j,k,0) = rho_hi[0]; // set ghost cell equal to reservoir density
                            cons(i,j,k,0) = rho_hi[0]; // set ghost cell equal to reservoir density

                            prim(i,j,k,4) = t_hi[0]; // set ghost cell equal to reservoir temperature
                            prim(i,j,k,5) = p_hi[0]; // set ghost cell equal to reservoir pressure
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            else if (bc_therm_hi[0] == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(2*hi-i+1,j,k,n);
                        }
                    }
                });
            } else if (bc_therm_hi[0] == 2) { // isothermal  (set ghost equal to wall value)
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (i > dom.bigEnd(0)) {
                        prim(i,j,k,4) = t_hi[0];
                        prim(i,j,k,5) = prim(2*hi-i+1,j,k,5);
                    }
                });
            }

        } // end HI X

        // LO Y
        if ((dim == 1) && (bx.smallEnd(1) < dom.smallEnd(1))) {

            int lo = dom.smallEnd(1);

            // mass fractions, wall
            if ( bc_mass_lo[1] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(i,2*lo-j-1,k,n);
                        }
                    }
                });
            }

            // mass fractions, concentration (set ghost equal to wall value)
            if (bc_mass_lo[1] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_y_lo[n];
                            prim(i,j,k,6+nspecies+n) = bc_Xk_y_lo[n];
                        }
                    }
                });
            }

            // mass fractions, density, temperature and pressure in the reservoir
            if (bc_mass_lo[1] >= 3 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_y_lo[n]; // set ghost cell equal to reservoir mass fraction
                            prim(i,j,k,6+nspecies+n) = bc_Xk_y_lo[n]; // set ghost cell equal to reservoir mole fraction

                            prim(i,j,k,0) = rho_lo[1]; // set ghost cell equal to reservoir density
                            cons(i,j,k,0) = rho_lo[1]; // set ghost cell equal to reservoir density

                            prim(i,j,k,4) = t_lo[1]; // set ghost cell equal to reservoir temperature
                            prim(i,j,k,5) = p_lo[1]; // set ghost cell equal to reservoir pressure
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            else if (bc_therm_lo[1] == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(i,2*lo-j-1,k,n);
                        }
                    }
                });
            } else if (bc_therm_lo[1] == 2) { // isothermal (set ghost equal to wall value)
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j < dom.smallEnd(1)) {
                        prim(i,j,k,4) = t_lo[1];
                        prim(i,j,k,5) = prim(i,2*lo-j-1,k,5);
                    }
                });
            }

        } // end LO Y


        // HI Y
        if ((dim == 1) && (bx.bigEnd(1) > dom.bigEnd(1))) {

            int hi = dom.bigEnd(1);

            // mass fractions, wall
            if ( bc_mass_hi[1] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(i,2*hi-j+1,k,n);

                        }
                    }
                });
            }

            // mass fractions, concentration (set ghost equal to wall value)
            if (bc_mass_hi[1] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_y_hi[n];
                            prim(i,j,k,6+nspecies+n) = bc_Xk_y_hi[n];
                        }
                    }
                });
            }

            // mass fractions, density, temperature and pressure in the reservoir
            if (bc_mass_hi[1] >= 3 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_y_hi[n]; // set ghost cell equal to reservoir mass fraction
                            prim(i,j,k,6+nspecies+n) = bc_Xk_y_hi[n]; // set ghost cell equal to reservoir mole fraction

                            prim(i,j,k,0) = rho_hi[1]; // set ghost cell equal to reservoir density
                            cons(i,j,k,0) = rho_hi[1]; // set ghost cell equal to reservoir density

                            prim(i,j,k,4) = t_hi[1]; // set ghost cell equal to reservoir temperature
                            prim(i,j,k,5) = p_hi[1]; // set ghost cell equal to reservoir pressure
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            else if (bc_therm_hi[1] == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(i,2*hi-j+1,k,n);
                        }
                    }
                });
            } else if (bc_therm_hi[1] == 2) { // isothermal (set ghost equal to wall value)
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (j > dom.bigEnd(1)) {
                        prim(i,j,k,4) = t_hi[1];
                        prim(i,j,k,5) = prim(i,2*hi-j+1,k,5);
                    }
                });
            }

        } // end HI Y

        // LO Z
        if ((dim == 2) && (bx.smallEnd(2) < dom.smallEnd(2))) {

            int lo = dom.smallEnd(2);

            // mass fractions, wall
            if ( bc_mass_lo[2] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(i,j,2*lo-k-1,n);
                        }
                    }
                });
            }

            // mass fractions, concentration (set ghost equal to wall value)
            if (bc_mass_lo[2] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_z_lo[n];
                            prim(i,j,k,6+nspecies+n) = bc_Xk_z_lo[n];
                        }
                    }
                });
            }

            // mass fractions, density, temperature and pressure in the reservoir
            if (bc_mass_lo[2] >= 3 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_z_lo[n]; // set ghost cell equal to reservoir mass fraction
                            prim(i,j,k,6+nspecies+n) = bc_Xk_z_lo[n]; // set ghost cell equal to reservoir mole fraction

                            prim(i,j,k,0) = rho_lo[2]; // set ghost cell equal to reservoir density
                            cons(i,j,k,0) = rho_lo[2]; // set ghost cell equal to reservoir density

                            prim(i,j,k,4) = t_lo[2]; // set ghost cell equal to reservoir temperature
                            prim(i,j,k,5) = p_lo[2]; // set ghost cell equal to reservoir pressure
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            else if (bc_therm_lo[2] == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(i,j,2*lo-k-1,n);
                        }
                    }
                });
            } else if (bc_therm_lo[2] == 2) { // isothermal (set ghost equal to wall value)
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k < dom.smallEnd(2)) {
                        prim(i,j,k,4) = t_lo[2];
                        prim(i,j,k,5) = prim(i,j,2*lo-k-1,5);
                    }
                });
            }

        } // end LO Z


        // HI Z
        if ((dim == 2) && (bx.bigEnd(2) > dom.bigEnd(2))) {

            int hi = dom.bigEnd(2);

            // mass fractions, wall
            if ( bc_mass_hi[2] == 1 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        for (int n=6; n<nprimvars; ++n) {
                            prim(i,j,k,n) = prim(i,j,2*hi-k+1,n);

                        }
                    }
                });
            }

            // mass fractions, concentration (set ghost equal to wall value)
            if (bc_mass_hi[2] == 2 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_z_hi[n];
                            prim(i,j,k,6+nspecies+n) = bc_Xk_z_hi[n];
                        }
                    }
                });
            }

            // mass fractions, density, temperature and pressure in the reservoir
            if (bc_mass_hi[2] >= 3 && algorithm_type == 2) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        for (int n=0; n<nspecies; ++n) {
                            prim(i,j,k,6+n)          = bc_Yk_z_hi[n]; // set ghost cell equal to reservoir mass fraction
                            prim(i,j,k,6+nspecies+n) = bc_Xk_z_hi[n]; // set ghost cell equal to reservoir mole fraction

                            prim(i,j,k,0) = rho_hi[2]; // set ghost cell equal to reservoir density
                            cons(i,j,k,0) = rho_hi[2]; // set ghost cell equal to reservoir density

                            prim(i,j,k,4) = t_hi[2]; // set ghost cell equal to reservoir temperature
                            prim(i,j,k,5) = p_hi[2]; // set ghost cell equal to reservoir pressure
                        }
                    }
                });
            }

            // temperature and pressure, adiabatic
            else if (bc_therm_hi[2] == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        for (int n=4; n<6; ++n) {
                            prim(i,j,k,n) = prim(i,j,2*hi-k+1,n);
                        }
                    }
                });
            } else if (bc_therm_hi[2] == 2) { // isothermal (set ghost equal to wall value)
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (k > dom.bigEnd(2)) {
                        prim(i,j,k,4) = t_hi[2];
                        prim(i,j,k,5) = prim(i,j,2*hi-k+1,5);
                    }
                });
            }

        } // end HI Z
    }
}

// Set normal momemntum and velocity on the boundary and ghost cells of the
// staggered grid based on slip/no-slip BCs
void BCMomNormal(MultiFab& mom_in, MultiFab& vel_in, MultiFab& cons_in,
                 const amrex::Geometry& geom, int dim)
{
    BL_PROFILE_VAR("BCMomNormal()",BCMomNormal);

    Box dom(geom.Domain());
    int ng_v = vel_in.nGrow();

    for ( MFIter mfi(vel_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng_v);

        const Array4<Real>& vel  = vel_in.array(mfi);
        const Array4<Real>& mom  = mom_in.array(mfi);
        const Array4<Real>& cons = cons_in.array(mfi);

        // LO X
        if ((dim == 0) && (bc_mass_lo[0] >= 3) && (bx.smallEnd(0) <= dom.smallEnd(0))) { // reservoir

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) {
                    // set ghost velocity & momentum
                    if (bc_mass_lo[0] == 3) {
                        vel(i,j,k) = -2*mom(dom.smallEnd(0),j,k)/(cons(dom.smallEnd(0),j,k,0) + cons(dom.smallEnd(0)-1,j,k,0));
                        mom(i,j,k) = -1*mom(dom.smallEnd(0),j,k);
                    }
                    else if (bc_mass_lo[0] == 4) {
                        vel(i,j,k) = -1*vel(dom.smallEnd(0),j,k);
                        mom(i,j,k) = -1*mom(dom.smallEnd(0),j,k);
                    }
                }
                else if (i == dom.smallEnd(0)) {
                    if (bc_mass_lo[0] == 3) {
                        vel(i,j,k) = 2*mom(i,j,k)/(cons(i,j,k,0) + cons(i-1,j,k,0));
                    }
                }
            });
        }
        else if ((dim == 0) && (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) && (bx.smallEnd(0) <= dom.smallEnd(0))) { // slip/no-slip

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) {
                    // set ghost velocity & momentum (set to wall value; no extrapolate)
                    vel(i,j,k) = 0.0;
                    mom(i,j,k) = 0.0;
                }
                else if (i == dom.smallEnd(0)) {
                    // set normal velocity & momentum
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }

        // HI X
        if ((dim == 0) && (bc_mass_hi[0] >= 3) && (bx.bigEnd(0) >= dom.bigEnd(0)+1)) { //reservoir

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > dom.bigEnd(0)+1) {
                    // set ghost velocity & momentum
                    if (bc_mass_hi[0] == 3) {
                        vel(i,j,k) = -2*mom(dom.bigEnd(0)+1,j,k)/(cons(dom.bigEnd(0)+1,j,k,0) + cons(dom.bigEnd(0)+1-1,j,k,0));
                        mom(i,j,k) = -1*mom(dom.bigEnd(0)+1,j,k);
                    }
                    else if (bc_mass_hi[0] == 4) {
                        vel(i,j,k) = -1*vel(dom.bigEnd(0)+1,j,k);
                        mom(i,j,k) = -1*mom(dom.bigEnd(0)+1,j,k);
                    }
                }
                else if (i == dom.bigEnd(0)+1) {
                    if (bc_mass_hi[0] == 3) {
                        vel(i,j,k) = 2*mom(i,j,k)/(cons(i,j,k,0) + cons(i-1,j,k,0));
                    }
                }
            });
        }
        else if ((dim == 0) && (bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) && (bx.bigEnd(0) >= dom.bigEnd(0)+1)) { // slip/no-slip

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > dom.bigEnd(0)+1) {
                    // set ghost velocity & momentum (set to wall value; no extrapolate)
                    vel(i,j,k) = 0.0;
                    mom(i,j,k) = 0.0;
                }
                else if (i == dom.bigEnd(0)+1) {
                    // set normal velocity & momentum
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }

        // LO Y
        if ((dim == 1) && (bc_mass_lo[1] >= 3) && (bx.smallEnd(1) <= dom.smallEnd(1))) { // reservoir

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) {
                    // set ghost velocity & momentum
                    if (bc_mass_lo[1] == 3) {
                        vel(i,j,k) = -2*mom(i,dom.smallEnd(1),k)/(cons(i,dom.smallEnd(1),k,0) + cons(i,dom.smallEnd(1)-1,k,0));
                        mom(i,j,k) = -1*mom(i,dom.smallEnd(1),k);
                    }
                    else if (bc_mass_lo[1] == 4) {
                        vel(i,j,k) = -1*vel(i,dom.smallEnd(1),k);
                        mom(i,j,k) = -1*mom(i,dom.smallEnd(1),k);
                    }
                }
                else if (j == dom.smallEnd(1)) {
                    if (bc_mass_lo[1] == 3) {
                        vel(i,j,k) = 2*mom(i,j,k)/(cons(i,j,k,0) + cons(i,j-1,k,0));
                    }
                }
            });
        }
        else if ((dim == 1) && (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) && (bx.smallEnd(1) <= dom.smallEnd(1))) { // slip/no-slip

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) {
                    // set ghost velocity & momentum (set to wall value; no extrapolate)
                    vel(i,j,k) = 0.0;
                    mom(i,j,k) = 0.0;
                }
                else if (j == dom.smallEnd(1)) {
                    // set normal velocity & momentum
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }

        // HI Y
        if ((dim == 1) && (bc_mass_hi[1] >= 3) && (bx.bigEnd(1) >= dom.bigEnd(1)+1)) { // reservoir

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j > dom.bigEnd(1)+1) {
                    // set ghost velocity & momentum
                    if (bc_mass_hi[1] == 3) {
                        vel(i,j,k) = -2*mom(i,dom.bigEnd(1)+1,k)/(cons(i,dom.bigEnd(1)+1,k,0) + cons(i,dom.bigEnd(1)+1-1,k,0));
                        mom(i,j,k) = -1*mom(i,dom.bigEnd(1)+1,k);
                    }
                    else if (bc_mass_hi[1] == 4) {
                        vel(i,j,k) = -1*vel(i,dom.bigEnd(1)+1,k);
                        mom(i,j,k) = -1*mom(i,dom.bigEnd(1)+1,k);
                    }
                }
                else if (j == dom.bigEnd(1)+1) {
                    if (bc_mass_hi[1] == 3) {
                        vel(i,j,k) = 2*mom(i,j,k)/(cons(i,j,k,0) + cons(i,j-1,k,0));
                    }
                }
            });
        }
        else if ((dim == 1) && (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) && (bx.bigEnd(1) >= dom.bigEnd(1)+1)) { // slip/no-slip

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j > dom.bigEnd(1)+1) {
                    // set ghost velocity & momentum (set to wall value; no extrapolate)
                    vel(i,j,k) = 0.0;
                    mom(i,j,k) = 0.0;
                }
                else if (j == dom.bigEnd(1)+1) {
                    // set normal velocity & momentum
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }

        // LO Z
        if ((dim == 2) && (bc_mass_lo[2] >= 3) && (bx.smallEnd(2) <= dom.smallEnd(2))) { // reservoir

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) {
                    // set ghost velocity & momentum
                    if (bc_mass_lo[2] == 3) {
                        vel(i,j,k) = -2*mom(i,j,dom.smallEnd(2))/(cons(i,j,dom.smallEnd(2),0) + cons(i,j,dom.smallEnd(2)-1,0));
                        mom(i,j,k) = -1*mom(i,j,dom.smallEnd(2));
                    }
                    else if (bc_mass_lo[2] == 4) {
                        vel(i,j,k) = -1*vel(i,j,dom.smallEnd(2));
                        mom(i,j,k) = -1*mom(i,j,dom.smallEnd(2));
                    }
                }
                else if (k == dom.smallEnd(2)) {
                    if (bc_mass_lo[2] == 3) {
                        vel(i,j,k) = 2*mom(i,j,k)/(cons(i,j,k,0) + cons(i,j,k-1,0));
                    }
                }
            });
        }
        else if ((dim == 2) && (bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) && (bx.smallEnd(2) <= dom.smallEnd(2))) { // slip/no-slip

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) {
                    // set ghost velocity & momentum (set to wall value; no extrapolate)
                    vel(i,j,k) = 0.0;
                    mom(i,j,k) = 0.0;
                }
                else if (k == dom.smallEnd(2)) {
                    // set normal velocity & momentum
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }

        // HI Z
        if ((dim == 2) && (bc_mass_hi[2] >= 3) && (bx.bigEnd(2) >= dom.bigEnd(2)+1)) { // reservoir

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k > dom.bigEnd(2)+1) {
                    // set ghost velocity & momentum
                    if (bc_mass_hi[2] == 3) {
                        vel(i,j,k) = -2*mom(i,j,dom.bigEnd(2)+1)/(cons(i,j,dom.bigEnd(2)+1,0) + cons(i,j,dom.bigEnd(2)+1-1,0));
                        mom(i,j,k) = -1*mom(i,j,dom.bigEnd(2)+1);
                    }
                    else if (bc_mass_hi[2] == 4) {
                        vel(i,j,k) = -1*vel(i,j,dom.bigEnd(2)+1);
                        mom(i,j,k) = -1*mom(i,j,dom.bigEnd(2)+1);
                    }
                }
                else if (k == dom.bigEnd(2)+1) {
                    if (bc_mass_hi[2] == 3) {
                        vel(i,j,k) = 2*mom(i,j,k)/(cons(i,j,k,0) + cons(i,j,k-1,0));
                    }
                }
            });
        }
        else if ((dim == 2) && (bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) && (bx.bigEnd(2) >= dom.bigEnd(2)+1)) { // slip/no-slip

            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k > dom.bigEnd(2)+1) {
                    // set ghost velocity & momentum (set to wall value; no extrapolate)
                    vel(i,j,k) = 0.0;
                    mom(i,j,k) = 0.0;
                }
                else if (k == dom.bigEnd(2)+1) {
                    // set normal velocity & momentum
                    vel(i,j,k) = 0.;
                    mom(i,j,k) = 0.;
                }
            });
        }
    }
}


// Set transverse momemntum and velocity on the boundary and ghost cells of the
// staggered grid based on slip/no-slip BCs
void BCMomTrans(MultiFab& mom_in, MultiFab& vel_in,
                 const amrex::Geometry& geom, int dim)
{
    BL_PROFILE_VAR("BCMomTrans()",BCMomTrans);

    if (geom.isAllPeriodic()) {
        return;
    }

    Box dom(geom.Domain());
    int ng_v = vel_in.nGrow();

    for ( MFIter mfi(vel_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng_v);

        const Array4<Real>& vel = vel_in.array(mfi);
        const Array4<Real>& mom = mom_in.array(mfi);

        // LO X
        if ((dim != 0) && (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) && (bx.smallEnd(0) < dom.smallEnd(0))) {
            Real fac;
            if (bc_mass_lo[0] >= 3) {fac = 0.0;} // reservoir
            else if (bc_vel_lo[0] == 1) {fac = 1.0;} // slip
            else if (bc_vel_lo[0] == 2) {fac = 0.0;} // no-slip (wall value in ghost cell)
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) {
                    vel(i,j,k) = fac*vel(-i-1,j,k);
                    mom(i,j,k) = fac*mom(-i-1,j,k);
                }
            });
        }

        // HI X
        if ((dim != 0) && (bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) && (bx.bigEnd(0) > dom.bigEnd(0))) {
            Real fac;
            if (bc_mass_hi[0] >= 3) {fac = 0.0;} // reservoir
            else if (bc_vel_hi[0] == 1) {fac = 1.0;} // slip
            else if (bc_vel_hi[0] == 2) {fac = 0.0;} // no-slip (wall value in ghost cell)
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > dom.bigEnd(0)) {
                    vel(i,j,k) = fac*vel(2*dom.bigEnd(0)-i+1,j,k);
                    mom(i,j,k) = fac*mom(2*dom.bigEnd(0)-i+1,j,k);
                }
            });
        }

        // LO Y
        if ((dim != 1) && (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) && (bx.smallEnd(1) < dom.smallEnd(1))) {
            Real fac;
            if (bc_mass_lo[1] >= 3) {fac = 0.0;} // reservoir
            else if (bc_vel_lo[1] == 1) {fac = 1.0;} // slip
            else if (bc_vel_lo[1] == 2) {fac = 0.0;} // no-slip (wall value in ghost cell)
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) {
                    vel(i,j,k) = fac*vel(i,-j-1,k);
                    mom(i,j,k) = fac*mom(i,-j-1,k);
                }
            });
        }

        // HI Y
        if ((dim != 1) && (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) && (bx.bigEnd(1) > dom.bigEnd(1))) {
            Real fac;
            if (bc_mass_hi[1] >= 3) {fac = 0.0;} // reservoir
            else if (bc_vel_hi[1] == 1) {fac = 1.0;} // slip
            else if (bc_vel_hi[1] == 2) {fac = 0.0;} // no-slip (wall value in ghost cell)
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j > dom.bigEnd(1)) {
                    vel(i,j,k) = fac*vel(i,2*dom.bigEnd(1)-j+1,k);
                    mom(i,j,k) = fac*mom(i,2*dom.bigEnd(1)-j+1,k);
                }
            });
        }

        // LO Z
        if ((dim != 2) && (bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) && (bx.smallEnd(2) < dom.smallEnd(2))) {
            Real fac;
            if (bc_mass_lo[2] >= 3) {fac = 0.0;} // reservoir
            else if (bc_vel_lo[2] == 1) {fac = 1.0;} // slip
            else if (bc_vel_lo[2] == 2) {fac = 0.0;} // no-slip (wall value in ghost cell)
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) {
                    vel(i,j,k) = fac*vel(i,j,-k-1);
                    mom(i,j,k) = fac*mom(i,j,-k-1);
                }
            });
        }

        // HI Z
        if ((dim != 2) && (bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) && (bx.bigEnd(2) > dom.bigEnd(2))) {
            Real fac;
            if (bc_mass_hi[2] >= 3) {fac = 0.0;} // reservoir
            else if (bc_vel_hi[2] == 1) {fac = 1.0;} // slip
            else if (bc_vel_hi[2] == 2) {fac = 0.0;} // no-slip (wall value in ghost cell)
            amrex::ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k > dom.bigEnd(2)) {
                    vel(i,j,k) = fac*vel(i,j,2*dom.bigEnd(2)-k+1);
                    mom(i,j,k) = fac*mom(i,j,2*dom.bigEnd(2)-k+1);
                }
            });
        }
    }
}


// Set density and energy density on BCs
void BCRhoRhoE(MultiFab& cons_in, MultiFab& prim_in,
               std::array< MultiFab, AMREX_SPACEDIM >& cumom_in,
               const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("BCRhoRhoE()",BCRhoRhoE);

    Box dom(geom.Domain());
    int ng_p = prim_in.nGrow();

    for ( MFIter mfi(prim_in); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng_p);

        const Array4<Real>& prim = prim_in.array(mfi);
        const Array4<Real>& cons = cons_in.array(mfi);
        AMREX_D_TERM(Array4<Real const> const& momx = cumom_in[0].array(mfi);,
                     Array4<Real const> const& momy = cumom_in[1].array(mfi);,
                     Array4<Real const> const& momz = cumom_in[2].array(mfi););

        // LO X
        if ((bc_mass_lo[0] >= 3) && (bx.smallEnd(0) < dom.smallEnd(0))) { // reservoir

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real rho = prim(i,j,k,0);
                    Real intenergy;

                    GetEnergy(intenergy,fracvec,temp);

                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }
        else if ((bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) && (bx.smallEnd(0) < dom.smallEnd(0))) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec);
                    GetEnergy(intenergy,fracvec,temp);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);

                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }

        // HI X
        if ((bc_mass_hi[0] >= 3) && (bx.bigEnd(0) > dom.bigEnd(0))) { // reservoir

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > dom.bigEnd(0)) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real rho = prim(i,j,k,0);
                    Real intenergy;

                    GetEnergy(intenergy,fracvec,temp);

                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }
        else if ((bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) && (bx.bigEnd(0) > dom.bigEnd(0))) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > dom.bigEnd(0)) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec);
                    GetEnergy(intenergy,fracvec,temp);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);

                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }

        // LO Y
        if ((bc_mass_lo[1] >= 3) && (bx.smallEnd(1) < dom.smallEnd(1))) { // reservoir

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real rho = prim(i,j,k,0);
                    Real intenergy;

                    GetEnergy(intenergy,fracvec,temp);

                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }
        else if ((bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) && (bx.smallEnd(1) < dom.smallEnd(1))) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec);
                    GetEnergy(intenergy,fracvec,temp);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);

                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }

        // HI Y
        if ((bc_mass_hi[1] >= 3) && (bx.bigEnd(1) > dom.bigEnd(1))) { // reservoir

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j > dom.bigEnd(1)) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real rho = prim(i,j,k,0);
                    Real intenergy;

                    GetEnergy(intenergy,fracvec,temp);

                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }
        else if ((bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) && (bx.bigEnd(1) > dom.bigEnd(1))) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j > dom.bigEnd(1)) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec);
                    GetEnergy(intenergy,fracvec,temp);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);

                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }

        // LO Z
        if ((bc_mass_lo[2] >= 3) && (bx.smallEnd(2) < dom.smallEnd(2))) { // reservoir

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real rho = prim(i,j,k,0);
                    Real intenergy;

                    GetEnergy(intenergy,fracvec,temp);

                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }
        else if ((bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) && (bx.smallEnd(2) < dom.smallEnd(2))) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec);
                    GetEnergy(intenergy,fracvec,temp);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);

                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }

        // HI Z
        if ((bc_mass_hi[2] >= 3) && (bx.bigEnd(2) > dom.bigEnd(2))) { // reservoir

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k > dom.bigEnd(2)) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real rho = prim(i,j,k,0);
                    Real intenergy;

                    GetEnergy(intenergy,fracvec,temp);

                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }
        else if ((bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) && (bx.bigEnd(2) > dom.bigEnd(2))) {

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k > dom.bigEnd(2)) {

                    GpuArray<Real,MAX_SPECIES> fracvec;
                    for (int n=0; n<nspecies; ++n) {
                        fracvec[n] = prim(i,j,k,6+n);
                    }

                    Real temp = prim(i,j,k,4);
                    Real pt = prim(i,j,k,5);
                    Real rho;
                    Real intenergy;

                    GetDensity(pt,rho,temp,fracvec);
                    GetEnergy(intenergy,fracvec,temp);

                    prim(i,j,k,0) = rho;
                    cons(i,j,k,0) = rho;
                    if (algorithm_type == 2) {
                        for (int n=0; n<nspecies; ++n) {
                            cons(i,j,k,5+n) = rho*prim(i,j,k,6+n);
                        }
                    }

                    Real kinenergy = 0.;
                    kinenergy += (momx(i+1,j,k) + momx(i,j,k))*(momx(i+1,j,k) + momx(i,j,k));
                    kinenergy += (momy(i,j+1,k) + momy(i,j,k))*(momy(i,j+1,k) + momy(i,j,k));
                    kinenergy += (momz(i,j,k+1) + momz(i,j,k))*(momz(i,j,k+1) + momz(i,j,k));
                    kinenergy *= (0.125/rho);

                    cons(i,j,k,4) = rho*intenergy + kinenergy;
                }
            });
        }
    }
}

void StochFluxStag(std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in, std::array< MultiFab, AMREX_SPACEDIM>& cenflux_in,
                   std::array< MultiFab, 2 >& edgeflux_x_in, std::array< MultiFab, 2 >& edgeflux_y_in,
                   std::array< MultiFab, 2 >& edgeflux_z_in, const amrex::Geometry& geom)
{
    BL_PROFILE_VAR("StochFluxStag()",StochFluxStag);

    // First we do mass boundary conditions (species fluxes reside on faces)
    // LO X
    if (bc_mass_lo[0] == 1 || bc_mass_lo[0] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = concentration   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_lo[0] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_in[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xlo = amrex::bdryNode(dom_x, Orientation(0, Orientation::low));

        for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xlo;
            Array4<Real> const& flux = (faceflux_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) *= factor;
                    }
                    // Set Dufour as well
                    flux(i,j,k,nvars+3) *= factor;
                });
            }
        }
    }
    // HI X
    if (bc_mass_hi[0] == 1 || bc_mass_hi[0] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = concentration   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_hi[0] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_in[0].ixType());

        // this is the x-hi domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xhi = amrex::bdryNode(dom_x, Orientation(0, Orientation::high));

        for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xhi;
            Array4<Real> const& flux = (faceflux_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) *= factor;
                    }
                    // Set Dufour as well
                    flux(i,j,k,nvars+3) *= factor;
                });
            }
        }

    }
    // LO Y
    if (bc_mass_lo[1] == 1 || bc_mass_lo[1] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = concentration   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_lo[1] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_in[1].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_ylo = amrex::bdryNode(dom_y, Orientation(1, Orientation::low));

        for (MFIter mfi(faceflux_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_ylo;
            Array4<Real> const& flux = (faceflux_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) *= factor;
                    }
                    // Set Dufour as well
                    flux(i,j,k,nvars+3) *= factor;
                });
            }
        }
    }
    // HI Y
    if (bc_mass_hi[1] == 1 || bc_mass_hi[1] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = concentration   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_hi[1] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_in[1].ixType());

        // this is the y-hi domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yhi = amrex::bdryNode(dom_y, Orientation(1, Orientation::high));

        for (MFIter mfi(faceflux_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yhi;
            Array4<Real> const& flux = (faceflux_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) *= factor;
                    }
                    // Set Dufour as well
                    flux(i,j,k,nvars+3) *= factor;
                });
            }
        }

    }
    // LO Z
    if (bc_mass_lo[2] == 1 || bc_mass_lo[2] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = concentration   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_lo[2] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_in[2].ixType());

        // this is the z-lo domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zlo = amrex::bdryNode(dom_z, Orientation(2, Orientation::low));

        for (MFIter mfi(faceflux_in[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zlo;
            Array4<Real> const& flux = (faceflux_in[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) *= factor;
                    }
                    // Set Dufour as well
                    flux(i,j,k,nvars+3) *= factor;
                });
            }
        }
    }
    // HI Z
    if (bc_mass_hi[2] == 1 || bc_mass_hi[2] == 2) {

        // 1 = wall        : multiply fluxes on wall by 0
        // 2 = concentration   : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_mass_hi[2] == 1) ? 0. : sqrt(2.);

        // domain grown nodally based on faceflux_in[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_in[2].ixType());

        // this is the z-hi domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zhi = amrex::bdryNode(dom_z, Orientation(2, Orientation::high));

        for (MFIter mfi(faceflux_in[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zhi;
            Array4<Real> const& flux = (faceflux_in[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        flux(i,j,k,n+5) *= factor;
                    }
                    // Set Dufour as well
                    flux(i,j,k,nvars+3) *= factor;
                });
            }
        }

    }

    // Next we do thermal boundary conditions (energy fluxes reside on faces)
    // LO X
    if (bc_therm_lo[0] == 1 || bc_therm_lo[0] == 2) {

        // 1 = adiabatic        : multiply fluxes on wall by 0
        // 2 = isothermal       : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_therm_lo[0] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_lo[0] >= 3) factor = 1.0;

        // domain grown nodally based on faceflux_in[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_in[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xlo = amrex::bdryNode(dom_x, Orientation(0, Orientation::low));

        for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xlo;
            Array4<Real> const& flux = (faceflux_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    flux(i,j,k,nvars) *= factor;
                });
            }
        }
    }
    // HI X
    if (bc_therm_hi[0] == 1 || bc_therm_hi[0] == 2) {

        // 1 = adiabatic        : multiply fluxes on wall by 0
        // 2 = isothermal       : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_therm_hi[0] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_hi[0] >= 3) factor = 1.0;

        // domain grown nodally based on faceflux_in[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_in[0].ixType());

        // this is the x-hi domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xhi = amrex::bdryNode(dom_x, Orientation(0, Orientation::high));

        for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xhi;
            Array4<Real> const& flux = (faceflux_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    flux(i,j,k,nvars) *= factor;
                });
            }
        }

    }
    // LO Y
    if (bc_therm_lo[1] == 1 || bc_therm_lo[1] == 2) {

        // 1 = adiabatic        : multiply fluxes on wall by 0
        // 2 = isothermal       : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_therm_lo[1] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_lo[1] >= 3) factor = 1.0;

        // domain grown nodally based on faceflux_in[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_in[1].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_ylo = amrex::bdryNode(dom_y, Orientation(1, Orientation::low));

        for (MFIter mfi(faceflux_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_ylo;
            Array4<Real> const& flux = (faceflux_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    flux(i,j,k,nvars) *= factor;
                });
            }
        }
    }
    // HI Y
    if (bc_therm_hi[1] == 1 || bc_therm_hi[1] == 2) {

        // 1 = adiabatic        : multiply fluxes on wall by 0
        // 2 = isothermal       : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_therm_hi[1] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_hi[1] >= 3) factor = 1.0;

        // domain grown nodally based on faceflux_in[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_in[1].ixType());

        // this is the y-hi domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yhi = amrex::bdryNode(dom_y, Orientation(1, Orientation::high));

        for (MFIter mfi(faceflux_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yhi;
            Array4<Real> const& flux = (faceflux_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    flux(i,j,k,nvars) *= factor;
                });
            }
        }

    }
    // LO Z
    if (bc_therm_lo[2] == 1 || bc_therm_lo[2] == 2) {

        // 1 = adiabatic        : multiply fluxes on wall by 0
        // 2 = isothermal       : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_therm_lo[2] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_lo[2] >= 3) factor = 1.0;

        // domain grown nodally based on faceflux_in[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_in[2].ixType());

        // this is the z-lo domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zlo = amrex::bdryNode(dom_z, Orientation(2, Orientation::low));

        for (MFIter mfi(faceflux_in[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zlo;
            Array4<Real> const& flux = (faceflux_in[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    flux(i,j,k,nvars) *= factor;
                });
            }
        }
    }
    // HI Z
    if (bc_therm_hi[2] == 1 || bc_therm_hi[2] == 2) {

        // 1 = adiabatic        : multiply fluxes on wall by 0
        // 2 = isothermal       : multiply fluxes on wall by sqrt(2)
        Real factor = (bc_therm_hi[2] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_hi[2] >= 3) factor = 1.0;

        // domain grown nodally based on faceflux_in[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_in[2].ixType());

        // this is the z-hi domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zhi = amrex::bdryNode(dom_z, Orientation(2, Orientation::high));

        for (MFIter mfi(faceflux_in[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zhi;
            Array4<Real> const& flux = (faceflux_in[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    flux(i,j,k,nvars) *= factor;
                });
            }
        }

    }

    // Last we do velocity boundary conditions (momentum flux resides on cell centers and edges)
    // But we do only edges, becuase the walls are not at cell centers
    // LO X edge, Y- and Z- momentum fluxes
    if (bc_vel_lo[0] == 1 || bc_vel_lo[0] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_lo[0] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_lo[0] >= 3) factor = 1.0;

        ////////////////////////////////////////////////
        // domain grown nodally based on edgeflux_x_in[0] nodality (xy)
        const Box& dom_xy = amrex::convert(geom.Domain(), edgeflux_x_in[0].ixType());

        // this is the x-lo domain boundary box (xy nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xy_xlo = amrex::bdryNode(dom_xy, Orientation(0, Orientation::low));

        for (MFIter mfi(edgeflux_x_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xy_xlo;
            Array4<Real> const& edgex_v = (edgeflux_x_in[0]).array(mfi);
            Array4<Real> const& edgey_u = (edgeflux_y_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    edgex_v(i,j,k) *= factor;
                    edgey_u(i,j,k) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
        // domain grown nodally based on edgeflux_x_in[1] nodality (xz)
        const Box& dom_xz = amrex::convert(geom.Domain(), edgeflux_x_in[1].ixType());

        // this is the x-lo domain boundary box (xz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xz_xlo = amrex::bdryNode(dom_xz, Orientation(0, Orientation::low));

        for (MFIter mfi(edgeflux_x_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xz_xlo;
            Array4<Real> const& edgex_w = (edgeflux_x_in[1]).array(mfi);
            Array4<Real> const& edgez_u = (edgeflux_z_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    edgex_w(i,j,k) *= factor;
                    edgez_u(i,j,k) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
        // set viscous heating
        // domain grown nodally based on faceflux_in[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_in[0].ixType());

        // this is the x-lo domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xlo = amrex::bdryNode(dom_x, Orientation(0, Orientation::low));

        for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xlo;
            Array4<Real> const& flux = (faceflux_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    flux(i,j,k,nvars+2) *= factor;
                });
            }
        }
    }

    // HI X edge, Y- and Z- momentum fluxes
    if (bc_vel_hi[0] == 1 || bc_vel_hi[0] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_hi[0] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_hi[0] >= 3) factor = 1.0;

        ////////////////////////////////////////////////
        // domain grown nodally based on edgeflux_x_in[0] nodality (xy)
        const Box& dom_xy = amrex::convert(geom.Domain(), edgeflux_x_in[0].ixType());

        // this is the x-hi domain boundary box (xy nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xy_xhi = amrex::bdryNode(dom_xy, Orientation(0, Orientation::high));

        for (MFIter mfi(edgeflux_x_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xy_xhi;
            Array4<Real> const& edgex_v = (edgeflux_x_in[0]).array(mfi);
            Array4<Real> const& edgey_u = (edgeflux_y_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    edgex_v(i,j,k) *= factor;
                    edgey_u(i,j,k) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
        // domain grown nodally based on edgeflux_x_in[1] nodality (xz)
        const Box& dom_xz = amrex::convert(geom.Domain(), edgeflux_x_in[1].ixType());

        // this is the x-hi domain boundary box (xz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xz_xhi = amrex::bdryNode(dom_xz, Orientation(0, Orientation::high));

        for (MFIter mfi(edgeflux_x_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xz_xhi;
            Array4<Real> const& edgex_w = (edgeflux_x_in[1]).array(mfi);
            Array4<Real> const& edgez_u = (edgeflux_z_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    edgex_w(i,j,k) *= factor;
                    edgez_u(i,j,k) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
        // set viscous heating
        // domain grown nodally based on faceflux_in[0] nodality (x)
        const Box& dom_x = amrex::convert(geom.Domain(), faceflux_in[0].ixType());

        // this is the x-hi domain boundary box (x nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xhi = amrex::bdryNode(dom_x, Orientation(0, Orientation::high));

        for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xhi;
            Array4<Real> const& flux = (faceflux_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    flux(i,j,k,nvars+2) *= factor;
                });
            }
        }
    }

    // LO Y edge, X- and Z- momentum fluxes
    if (bc_vel_lo[1] == 1 || bc_vel_lo[1] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_lo[1] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_lo[1] >= 3) factor = 1.0;

        ////////////////////////////////////////////////
        // domain grown nodally based on edgeflux_y_in[0] nodality (xy)
        const Box& dom_xy = amrex::convert(geom.Domain(), edgeflux_y_in[0].ixType());

        // this is the y-lo domain boundary box (xy nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xy_ylo = amrex::bdryNode(dom_xy, Orientation(1, Orientation::low));

        for (MFIter mfi(edgeflux_y_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xy_ylo;
            Array4<Real> const& edgey_u = (edgeflux_y_in[0]).array(mfi);
            Array4<Real> const& edgex_v = (edgeflux_x_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    edgey_u(i,j,k) *= factor;
                    edgex_v(i,j,k) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
        // domain grown nodally based on edgeflux_y_in[2] nodality (yz)
        const Box& dom_yz = amrex::convert(geom.Domain(), edgeflux_y_in[1].ixType());

        // this is the y-lo domain boundary box (yz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yz_ylo = amrex::bdryNode(dom_yz, Orientation(1, Orientation::low));

        for (MFIter mfi(edgeflux_y_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yz_ylo;
            Array4<Real> const& edgey_w = (edgeflux_y_in[1]).array(mfi);
            Array4<Real> const& edgez_v = (edgeflux_z_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    edgey_w(i,j,k) *= factor;
                    edgez_v(i,j,k) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
        // set viscous heating
        // domain grown nodally based on faceflux_in[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_in[1].ixType());

        // this is the y-lo domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_ylo = amrex::bdryNode(dom_y, Orientation(1, Orientation::low));

        for (MFIter mfi(faceflux_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_ylo;
            Array4<Real> const& flux = (faceflux_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    flux(i,j,k,nvars+2) *= factor;
                });
            }
        }
    }

    // HI Y edge, X- and Z- momentum fluxes
    if (bc_vel_hi[1] == 1 || bc_vel_hi[1] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_hi[1] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_hi[1] >= 3) factor = 1.0;

        ////////////////////////////////////////////////
        // domain grown nodally based on edgeflux_y_in[0] nodality (xy)
        const Box& dom_xy = amrex::convert(geom.Domain(), edgeflux_y_in[0].ixType());

        // this is the y-hi domain boundary box (xy nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xy_yhi = amrex::bdryNode(dom_xy, Orientation(1, Orientation::high));

        for (MFIter mfi(edgeflux_y_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xy_yhi;
            Array4<Real> const& edgey_u = (edgeflux_y_in[0]).array(mfi);
            Array4<Real> const& edgex_v = (edgeflux_x_in[0]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    edgey_u(i,j,k) *= factor;
                    edgex_v(i,j,k) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
        // domain grown nodally based on edgeflux_y_in[2] nodality (yz)
        const Box& dom_yz = amrex::convert(geom.Domain(), edgeflux_y_in[1].ixType());

        // this is the y-hi domain boundary box (yz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yz_yhi = amrex::bdryNode(dom_yz, Orientation(1, Orientation::high));

        for (MFIter mfi(edgeflux_y_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yz_yhi;
            Array4<Real> const& edgey_w = (edgeflux_y_in[1]).array(mfi);
            Array4<Real> const& edgez_v = (edgeflux_z_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    edgey_w(i,j,k) *= factor;
                    edgez_v(i,j,k) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
        // set viscous heating
        // domain grown nodally based on faceflux_in[1] nodality (y)
        const Box& dom_y = amrex::convert(geom.Domain(), faceflux_in[1].ixType());

        // this is the y-hi domain boundary box (y nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yhi = amrex::bdryNode(dom_y, Orientation(1, Orientation::high));

        for (MFIter mfi(faceflux_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yhi;
            Array4<Real> const& flux = (faceflux_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    flux(i,j,k,nvars+2) *= factor;
                });
            }
        }
    }

    // LO Z edge, X- and Y- momentum fluxes
    if (bc_vel_lo[2] == 1 || bc_vel_lo[2] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_lo[2] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_lo[2] >= 3) factor = 1.0;

        ////////////////////////////////////////////////
        // domain grown nodally based on edgeflux_z_in[1] nodality (xz)
        const Box& dom_xz = amrex::convert(geom.Domain(), edgeflux_z_in[0].ixType());

        // this is the z-lo domain boundary box (xz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xz_zlo = amrex::bdryNode(dom_xz, Orientation(2, Orientation::low));

        for (MFIter mfi(edgeflux_z_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xz_zlo;
            Array4<Real> const& edgez_u = (edgeflux_z_in[0]).array(mfi);
            Array4<Real> const& edgex_w = (edgeflux_x_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    edgez_u(i,j,k) *= factor;
                    edgex_w(i,j,k) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
        // domain grown nodally based on edgeflux_z_in[2] nodality (yz)
        const Box& dom_yz = amrex::convert(geom.Domain(), edgeflux_z_in[1].ixType());

        // this is the y-lo domain boundary box (yz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yz_zlo = amrex::bdryNode(dom_yz, Orientation(2, Orientation::low));

        for (MFIter mfi(edgeflux_z_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yz_zlo;
            Array4<Real> const& edgez_v = (edgeflux_z_in[1]).array(mfi);
            Array4<Real> const& edgey_w = (edgeflux_y_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    edgez_v(i,j,k) *= factor;
                    edgey_w(i,j,k) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
        // set viscous heating
        // domain grown nodally based on faceflux_in[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_in[2].ixType());

        // this is the z-lo domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zlo = amrex::bdryNode(dom_z, Orientation(2, Orientation::low));

        for (MFIter mfi(faceflux_in[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zlo;
            Array4<Real> const& flux = (faceflux_in[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    flux(i,j,k,nvars+2) *= factor;
                });
            }
        }
    }

    // HI Z edge, X- and Y- momentum fluxes
    if (bc_vel_hi[2] == 1 || bc_vel_hi[2] == 2) {

        // 1 = slip wall   : multiply fluxes on wall by 0
        // 2 = no-slip wall: multiply fluxes on wall by sqrt(2)
        Real factor = (bc_vel_hi[2] == 1) ? 0. : sqrt(2.);
        // reservoir            : unchanged
        if (bc_mass_hi[2] >= 3) factor = 1.0;

        ////////////////////////////////////////////////
        // domain grown nodally based on edgeflux_z_in[1] nodality (xz)
        const Box& dom_xz = amrex::convert(geom.Domain(), edgeflux_z_in[0].ixType());

        // this is the z-hi domain boundary box (xz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_xz_zhi = amrex::bdryNode(dom_xz, Orientation(2, Orientation::high));

        for (MFIter mfi(edgeflux_z_in[0]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_xz_zhi;
            Array4<Real> const& edgez_u = (edgeflux_z_in[0]).array(mfi);
            Array4<Real> const& edgex_w = (edgeflux_x_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    edgez_u(i,j,k) *= factor;
                    edgex_w(i,j,k) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
        // domain grown nodally based on edgeflux_z_in[2] nodality (yz)
        const Box& dom_yz = amrex::convert(geom.Domain(), edgeflux_z_in[1].ixType());

        // this is the y-hi domain boundary box (yz nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_yz_zhi = amrex::bdryNode(dom_yz, Orientation(2, Orientation::high));

        for (MFIter mfi(edgeflux_z_in[1]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_yz_zhi;
            Array4<Real> const& edgez_v = (edgeflux_z_in[1]).array(mfi);
            Array4<Real> const& edgey_w = (edgeflux_y_in[1]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    edgez_v(i,j,k) *= factor;
                    edgey_w(i,j,k) *= factor;
                });
            }
        }

        ////////////////////////////////////////////////
        // set viscous heating
        // domain grown nodally based on faceflux_in[2] nodality (z)
        const Box& dom_z = amrex::convert(geom.Domain(), faceflux_in[2].ixType());

        // this is the z-hi domain boundary box (z nodality)
        // Orientation(dir,Orientation)  -- Orientation can be ::low or ::high
        const Box& dom_zhi = amrex::bdryNode(dom_z, Orientation(2, Orientation::high));

        for (MFIter mfi(faceflux_in[2]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();
            const Box& b = bx & dom_zhi;
            Array4<Real> const& flux = (faceflux_in[2]).array(mfi);
            if (b.ok()) {
                amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    flux(i,j,k,nvars+2) *= factor;
                });
            }
        }
    }


    //////////////////////////////////////////
    // set stochastic momentum flux to zero in
    // the reservoir cells
    //////////////////////////////////////////
    Box dom(geom.Domain());

    for ( MFIter mfi(cenflux_in[0]); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(1);

        const Array4<Real>& cenx = cenflux_in[0].array(mfi);
        const Array4<Real>& ceny = cenflux_in[1].array(mfi);
        const Array4<Real>& cenz = cenflux_in[2].array(mfi);

        // LO X
        if ((bc_mass_lo[0] == 3) && (bx.smallEnd(0) < dom.smallEnd(0))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) cenx(i,j,k) = 0.0;
            });
        }
        if ((bc_mass_lo[0] == 4) && (bx.smallEnd(0) < dom.smallEnd(0))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i < dom.smallEnd(0)) cenx(i,j,k) = cenx(i+1,j,k);
            });
        }

        // HI X
        if ((bc_mass_hi[0] == 3) && (bx.bigEnd(0) > dom.bigEnd(0))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > dom.bigEnd(0)) cenx(i,j,k) = 0.0;
            });
        }
        if ((bc_mass_hi[0] == 4) && (bx.bigEnd(0) > dom.bigEnd(0))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i > dom.bigEnd(0)) cenx(i,j,k) = cenx(i-1,j,k);
            });
        }

        // LO Y
        if ((bc_mass_lo[1] == 3) && (bx.smallEnd(1) < dom.smallEnd(1))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) ceny(i,j,k) = 0.0;
            });
        }
        if ((bc_mass_lo[1] == 4) && (bx.smallEnd(1) < dom.smallEnd(1))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < dom.smallEnd(1)) ceny(i,j,k) = ceny(i,j+1,k);
            });
        }

        // HI Y
        if ((bc_mass_hi[1] == 3) && (bx.bigEnd(1) > dom.bigEnd(1))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j > dom.bigEnd(1)) ceny(i,j,k) = 0.0;
            });
        }
        if ((bc_mass_hi[1] == 4) && (bx.bigEnd(1) > dom.bigEnd(1))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j > dom.bigEnd(1)) ceny(i,j,k) = ceny(i,j-1,k);
            });
        }

        // LO Z
        if ((bc_mass_lo[2] == 3) && (bx.smallEnd(2) < dom.smallEnd(2))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) cenz(i,j,k) = 0.0;
            });
        }
        if ((bc_mass_lo[2] == 4) && (bx.smallEnd(2) < dom.smallEnd(2))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k < dom.smallEnd(2)) cenz(i,j,k) = cenz(i,j,k+1);
            });
        }

        // HI Z
        if ((bc_mass_hi[2] == 3) && (bx.bigEnd(2) > dom.bigEnd(2))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k > dom.bigEnd(2)) cenz(i,j,k) = 0.0;
            });
        }
        if ((bc_mass_hi[2] == 4) && (bx.bigEnd(2) > dom.bigEnd(2))) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k > dom.bigEnd(2)) cenz(i,j,k) = cenz(i,j,k-1);
            });
        }
    }

}

void StochFluxMem(std::array<MultiFab, AMREX_SPACEDIM>& faceflux_in, std::array< MultiFab, 2 >& edgeflux_x_in,
                   std::array< MultiFab, 2 >& edgeflux_y_in, std::array< MultiFab, 2 >& edgeflux_z_in)

{

    BL_PROFILE_VAR("StochFluxMem()",StochFluxMem);

    // The membrane is an adiabatic wall -- setup the stochastic heat and species fluxes to zero
    for (MFIter mfi(faceflux_in[0]); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();
        const Array4<Real>& xflux = faceflux_in[0].array(mfi);

        if (bx.smallEnd(0) == membrane_cell) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == membrane_cell) {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) = 0.;
                    }
                    // heat
                    xflux(i,j,k,nvars+0) = 0.; // stochastic heating (adiabatic wall)
                    xflux(i,j,k,nvars+1) = 0.; // stochastic viscous heating (normal velocity zero at membrane)
                    xflux(i,j,k,nvars+2) = 0.; // stochastic viscous heating (slip BC)
                    xflux(i,j,k,nvars+3) = 0.; // stochastic dufour
                }
            });
        }
        else if (bx.bigEnd(0) == membrane_cell) {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i == membrane_cell) {
                    // species
                    for (int n=0;n<nspecies;++n) {
                        xflux(i,j,k,5+n) = 0.;
                    }
                    // heat
                    xflux(i,j,k,nvars+0) = 0.; // stochastic heating (adiabatic wall)
                    xflux(i,j,k,nvars+1) = 0.; // stochastic viscous heating (normal velocity zero at membrane)
                    xflux(i,j,k,nvars+2) = 0.; // stochastic viscous heating (slip BC)
                    xflux(i,j,k,nvars+3) = 0.; // stochastic dufour
                }
            });
        }

    }

    // Set transverse momentum at the membrane according to the full slip condition
    // XY
    for (MFIter mfi(edgeflux_y_in[0]); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();
        const Array4<Real>& edgey_u = edgeflux_y_in[0].array(mfi);
        const Array4<Real>& edgex_v = edgeflux_x_in[0].array(mfi);

        if (bx.smallEnd(0) == membrane_cell) {
              amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                  if (i == bx.smallEnd(0)) {
                      edgey_u(i,j,k) = 0.0;
                      edgex_v(i,j,k) = 0.0;
                  }
              });
        }
        if (bx.bigEnd(0) == membrane_cell) {
              amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                  if (i == bx.bigEnd(0)) {
                      edgey_u(i,j,k) = 0.0;
                      edgex_v(i,j,k) = 0.0;
                  }
              });
        }
    }

    for (MFIter mfi(edgeflux_z_in[0]); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.validbox();
        const Array4<Real>& edgez_u = edgeflux_z_in[0].array(mfi);
        const Array4<Real>& edgex_w = edgeflux_x_in[1].array(mfi);

        if (bx.smallEnd(0) == membrane_cell) {
              amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                  if (i == bx.smallEnd(0)) {
                      edgez_u(i,j,k) = 0.0;
                      edgex_w(i,j,k) = 0.0;
                  }
              });
        }
        if (bx.bigEnd(0) == membrane_cell) {
              amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
              {
                  if (i == bx.bigEnd(0)) {
                      edgez_u(i,j,k) = 0.0;
                      edgex_w(i,j,k) = 0.0;
                  }
              });
        }
    }
}
