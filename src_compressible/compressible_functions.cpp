#include "compressible_functions.H"
#include "compressible_functions_F.H"
#include "compressible_namespace.H"

using namespace compressible;

void InitializeCompressibleNamespace() {

  mass_bc_lo.resize(AMREX_SPACEDIM);
  mass_bc_hi.resize(AMREX_SPACEDIM);
  therm_bc_lo.resize(AMREX_SPACEDIM);
  therm_bc_hi.resize(AMREX_SPACEDIM);
  vel_bc_lo.resize(AMREX_SPACEDIM);
  vel_bc_hi.resize(AMREX_SPACEDIM);
  
  Yk_bc.resize(MAX_SPECIES*LOHI*AMREX_SPACEDIM);
  Xk_bc.resize(MAX_SPECIES*LOHI*AMREX_SPACEDIM);

  initialize_compressible_namespace(mass_bc_lo.dataPtr(), mass_bc_hi.dataPtr(),
				    therm_bc_lo.dataPtr(), therm_bc_hi.dataPtr(),
				    vel_bc_lo.dataPtr(), vel_bc_hi.dataPtr(),
				    Yk_bc.dataPtr(), Xk_bc.dataPtr());
 
}
