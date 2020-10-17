#include "multispec_functions.H"

void DotWithZ(MultiFab& mf,
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

void DotWithZFace(MultiFab& mf,
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
            MultiFab::Saxpy(mfdotz[dir],z_temp[comp],mf,comp,0,1,mfdotz[dir].nGrow());
        }
    }
}

void ComputeChargeCoef(MultiFab& rho_in,
                       MultiFab& Temp_in,
                       MultiFab& charge_coef_in)
{
    BL_PROFILE_VAR("ComputeChargeCoef()",ComputeChargeCoef);

    int nspecies_gpu = nspecies;

    GpuArray<Real,MAX_SPECIES> molmass_gpu;
    for (int comp=0; comp<nspecies_gpu; ++comp) {
        molmass_gpu[comp] = molmass[comp];
    }


    for ( MFIter mfi(charge_coef_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.growntilebox(1);
        
        const Array4<Real> & rho         = rho_in.array(mfi);
        const Array4<Real> & Temp        = Temp_in.array(mfi);
        const Array4<Real> & charge_coef = charge_coef_in.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real n=0.;
            for (int comp=0; comp<nspecies_gpu; ++comp) {
                n += rho(i,j,k,comp) / molmass_gpu[comp];
            }

            for (int comp=0; comp<nspecies_gpu; ++comp) {

            }
            
        });
    }
}
