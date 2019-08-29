#include "compressible_functions.H"
#include "compressible_functions_F.H"
#include "compressible_namespace.H"

using namespace compressible;

void InitializeCompressibleNamespace() {

  bc_mass_lo.resize(AMREX_SPACEDIM);
  bc_mass_hi.resize(AMREX_SPACEDIM);
  bc_therm_lo.resize(AMREX_SPACEDIM);
  bc_therm_hi.resize(AMREX_SPACEDIM);
  bc_vel_lo.resize(AMREX_SPACEDIM);
  bc_vel_hi.resize(AMREX_SPACEDIM);

  bc_Yk.resize(MAX_SPECIES*LOHI*AMREX_SPACEDIM);
  bc_Xk.resize(MAX_SPECIES*LOHI*AMREX_SPACEDIM);

  initialize_compressible_namespace(bc_mass_lo.dataPtr(), bc_mass_hi.dataPtr(),
                                    bc_therm_lo.dataPtr(), bc_therm_hi.dataPtr(),
                                    bc_vel_lo.dataPtr(), bc_vel_hi.dataPtr(),
                                    bc_Yk.dataPtr(), bc_Xk.dataPtr());
}
