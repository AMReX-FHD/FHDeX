#include "compressible_functions.H"
#include "compressible_functions_F.H"
#include "compressible_namespace.H"

using namespace compressible;

void InitializeCompressibleNamespace() {
  
  Yk_lo.resize(MAX_SPECIES*AMREX_SPACEDIM);
  Yk_hi.resize(MAX_SPECIES*AMREX_SPACEDIM);
  Xk_lo.resize(MAX_SPECIES*AMREX_SPACEDIM);
  Xk_hi.resize(MAX_SPECIES*AMREX_SPACEDIM);

  // Yk_lo_x.resize(MAX_SPECIES);
  // Yk_hi_x.resize(MAX_SPECIES);
  // Yk_lo_y.resize(MAX_SPECIES);
  // Yk_hi_y.resize(MAX_SPECIES);
  // Yk_lo_z.resize(MAX_SPECIES);
  // Yk_hi_z.resize(MAX_SPECIES);
  // Xk_lo_x.resize(MAX_SPECIES);
  // Xk_hi_x.resize(MAX_SPECIES);
  // Xk_lo_y.resize(MAX_SPECIES);
  // Xk_hi_y.resize(MAX_SPECIES);
  // Xk_lo_z.resize(MAX_SPECIES);
  // Xk_hi_z.resize(MAX_SPECIES);

  initialize_compressible_namespace(Yk_lo.dataPtr(), Yk_hi.dataPtr(), Xk_lo.dataPtr(), Xk_hi.dataPtr());
  
  // initialize_compressible_namespace(Yk_lo_x.dataPtr(), Yk_hi_x.dataPtr(), Yk_lo_y.dataPtr(), Yk_hi_y.dataPtr(), Yk_lo_z.dataPtr(), Yk_hi_z.dataPtr(),
  // 				    Xk_lo_x.dataPtr(), Xk_hi_x.dataPtr(), Xk_lo_y.dataPtr(), Xk_hi_y.dataPtr(), Xk_lo_z.dataPtr(), Xk_hi_z.dataPtr());
  
}
