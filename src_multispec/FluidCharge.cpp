#include "multispec_functions.H"

void DotWithZ(const MultiFab& mf,
              MultiFab& mfdotz,
              int abs_z)
{
    BL_PROFILE_VAR("DotWithZ()",DotWithZ);

    amrex::Vector<amrex::Real> z_temp;
    z_temp.resize(nspecies);

    for (int comp=0; comp<nspecies; ++comp) {
        z_temp[comp] = charge_per_mass[comp];
    }

    // Dot with abs(z) to estimate norms for relative errors (default=false)
    if (abs_z == 1) {
        for (int comp=0; comp<nspecies; ++comp) {
            z_temp[comp] = amrex::Math::abs(z_temp[comp]);
        }
    }

    mfdotz.setVal(0.);
    for (int comp=0; comp<nspecies; ++comp) {
        MultiFab::Saxpy(mfdotz,z_temp[comp],mf,comp,0,1,mfdotz.nGrow());
    }
}

void DotWithZFace(std::array< const MultiFab, AMREX_SPACEDIM >& mf,
                  std::array< MultiFab, AMREX_SPACEDIM >& mfdotz,
                  int abs_z)
{
    BL_PROFILE_VAR("DotWithZFace()",DotWithZFace);

    amrex::Vector<amrex::Real> z_temp;
    z_temp.resize(nspecies);

    for (int comp=0; comp<nspecies; ++comp) {
        z_temp[comp] = charge_per_mass[comp];
    }

    // Dot with abs(z) to estimate norms for relative errors (default=false)
    if (abs_z == 1) {
        for (int comp=0; comp<nspecies; ++comp) {
            z_temp[comp] = amrex::Math::abs(z_temp[comp]);
        }
    }

    for (int dir=0; dir<AMREX_SPACEDIM; ++dir) {
        mfdotz[dir].setVal(0.);
        for (int comp=0; comp<nspecies; ++comp) {
            MultiFab::Saxpy(mfdotz[dir],z_temp[comp],mf[dir],comp,0,1,mfdotz[dir].nGrow());
        }
    }
}

void ComputeChargeCoef(const MultiFab& rho_in,
                       const MultiFab& Temp_in,
                       MultiFab& charge_coef_in)
{
    BL_PROFILE_VAR("ComputeChargeCoef()",ComputeChargeCoef);

    for ( MFIter mfi(charge_coef_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.growntilebox(1);

        const Array4<const Real> & rho   = rho_in.array(mfi);
        const Array4<const Real> & Temp  = Temp_in.array(mfi);
        const Array4<Real> & charge_coef = charge_coef_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real n=0.;
            for (int comp=0; comp<nspecies; ++comp) {
                n += rho(i,j,k,comp) / molmass[comp];
            }

            for (int comp=0; comp<nspecies; ++comp) {
                charge_coef(i,j,k,comp) = rho(i,j,k,comp)*charge_per_mass[comp]
                    / (n*k_B*Temp(i,j,k));
            }
        });
    }
}

void EnforceChargeNeutrality()
{
    BL_PROFILE_VAR("EnforceChargeNeutrality",EnforceChargeNeutrality);

    Abort("EnforceChargeNeutrality() not written yet");
}

// compute face-centered A_\Phi (an nspecies vector)
// compute vector rho W z / (n k_B T) on cell centers and average to faces
// compute tensor rho W chi on cell centers and average to faces
// multiply them together and store the resulting vector in A_Phi
void ImplicitPotentialCoef()
{
    BL_PROFILE_VAR("ImplicitPotentialCoef",ImplicitPotentialCoef);

    Abort("ImplicitPotentialCoef() not written yet");
}

void ModifyS()
{
    BL_PROFILE_VAR("ModifyS",ModifyS);

    Abort("ModifyS() not written yet; needed for AdvanceTimestepIterative");
}

void ComputePermittivity()
{
    BL_PROFILE_VAR("ComputePermittivity",ComputePermittivity);

    if (dielectric_type == 0) {
        return;
    } else {
        Abort("ComputePermittivity() for dielectric_type !=0 not written yet");
    }
}

void ComputeLorentzForce(std::array< MultiFab, AMREX_SPACEDIM >& Lorentz_force,
                         std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot,
                         const MultiFab& permittivity,
                         const MultiFab& charge,
                         const Geometry& geom)
{
    BL_PROFILE_VAR("ComputeLorentzForce",ComputeLorentzForce);

    bool use_qE_Lorentz = false;

    if (use_qE_Lorentz) {

        AverageCCToFace(charge,Lorentz_force,0,1,SPEC_BC_COMP,geom);

        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            Lorentz_force[i].mult(-1.,0,1);
            MultiFab::Add(Lorentz_force[i],grad_Epot[i],0,0,1,0);
        }

        return;
    }

    // temporary multifab to hold div(-eps*E)
    MultiFab temp_cc(permittivity.boxArray(),permittivity.DistributionMap(),1,1);

    // Lorentz force = E div (eps*E) - (1/2) E^2 grad(eps)

    // discretize E div(eps*E)
    // start by averaging epsilon to faces
    // store this in Lorentz_force temporarily
    if (dielectric_type == 0) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            Lorentz_force[i].setVal(dielectric_const);
        }
    } else {
        AverageCCToFace(permittivity,Lorentz_force,0,1,SPEC_BC_COMP,geom);
        if (zero_eps_on_wall_type > 0) {
            // set beta to set to zero on certain boundary faces
            ZeroEpsOnWall(Lorentz_force);
        }
    }

    // multiply by grad_Epot = -E to get -eps*E on faces
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        MultiFab::Multiply(Lorentz_force[i],grad_Epot[i],0,0,1,0);
    }

    // take divergence of -eps*E and put it in temp_cc
    ComputeDiv(temp_cc,Lorentz_force,0,0,1,geom,0);

    // fill ghost cells for div(-eps*E)
    temp_cc.FillBoundary(geom.periodicity());
    MultiFabPhysBC(temp_cc,geom,0,1,SPEC_BC_COMP);

    // average div(-eps*E) to faces, store in Lorentz_force
    AverageCCToFace(temp_cc,Lorentz_force,0,1,SPEC_BC_COMP,geom);

    // multiply by -E to get E*div(eps*E) on faces
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        MultiFab::Multiply(Lorentz_force[i],grad_Epot[i],0,0,1,0);
    }

    if (dielectric_type != 0) {
        // spatially-varying permittivity
        // add -(1/2) E^2 grad(eps)
        Abort("ComputeLorentzForce() spatially-varying epsilon not implemented yet");
    }

}

void ComputeE_ext(std::array< MultiFab, AMREX_SPACEDIM >& E_ext) {

    BL_PROFILE_VAR("ComputeE_ext",ComputeE_ext);

    if (E_ext_type == 1) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            E_ext[i].setVal(E_ext_value[i]);
        }
    }

}

void ZeroEpsOnWall(std::array< MultiFab, AMREX_SPACEDIM >& beta) {

    BL_PROFILE_VAR("ZeroEpsOnWall",ZeroEpsOnWall);

    Abort("ZeroEpsOnWall not written yet");
}
