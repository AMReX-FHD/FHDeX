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

    int nspecies_gpu = nspecies;
    Real k_B_gpu = k_B;
    
    GpuArray<Real,MAX_SPECIES> molmass_gpu;
    GpuArray<Real,MAX_SPECIES> charge_per_mass_gpu;
    for (int comp=0; comp<nspecies_gpu; ++comp) {
        molmass_gpu[comp] = molmass[comp];
        charge_per_mass_gpu[comp] = charge_per_mass[comp];
    }

    for ( MFIter mfi(charge_coef_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.growntilebox(1);
        
        const Array4<const Real> & rho   = rho_in.array(mfi);
        const Array4<const Real> & Temp  = Temp_in.array(mfi);
        const Array4<Real> & charge_coef = charge_coef_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real n=0.;
            for (int comp=0; comp<nspecies_gpu; ++comp) {
                n += rho(i,j,k,comp) / molmass_gpu[comp];
            }

            for (int comp=0; comp<nspecies_gpu; ++comp) {
                charge_coef(i,j,k,comp) = rho(i,j,k,comp)*charge_per_mass_gpu[comp]
                    / (n*k_B_gpu*Temp(i,j,k));
            }
        });
    }
}

void EnforceChargeNeutrality()
{
    Abort("EnforceChargeNeutrality() not written yet");
}

// compute face-centered A_\Phi (an nspecies vector)
// compute vector rho W z / (n k_B T) on cell centers and average to faces
// compute tensor rho W chi on cell centers and average to faces
// multiply them together and store the resulting vector in A_Phi
void ImplicitPotentialCoef()
{
    Abort("ImplicitPotentialCoef() not written yet");
}

void ModifyS()
{
    Abort("ModifyS() not written yet; needed for AdvanceTimestepIterative");
}

void ComputePermittivity()
{
    if (dielectric_type == 0) {
        return;
    } else {
        Abort("ComputePermittivity() for dielectric_type !=0 not written yet");
    }
}

void ComputeLorentzForce(std::array< MultiFab, AMREX_SPACEDIM >& Lorentz_force,
                         std::array< const MultiFab, AMREX_SPACEDIM >& grad_Epot,
                         const MultiFab& permittivity,
                         const MultiFab& charge,
                         const Geometry& geom)
{

}

void ComputeE_ext(MultiFab& E_ext) {

    if (E_ext_type == 1) {
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            E_ext[i].setVal(E_ext_value[i]);
        }
    }
    
}

void ZeroEpsOnWall(std::array< MultiFab, AMREX_SPACEDIM >& beta) {
    Abort("ZeroEpsOnWall not written yet");
}
