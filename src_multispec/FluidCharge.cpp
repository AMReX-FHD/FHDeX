#include "multispec_functions.H"

void DotWithZ(MultiFab& mf,
              MultiFab& mfdotz,
              int abs_z) {

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
                  int abs_z) {

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
