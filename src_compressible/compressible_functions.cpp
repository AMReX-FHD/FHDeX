#include "compressible_functions.H"
#include "compressible_functions_F.H"
#include "compressible_namespace.H"

using namespace compressible;

void InitializeCompressibleNamespace() {
  
  Yk_bc.resize(MAX_SPECIES*LOHI*AMREX_SPACEDIM);
  Xk_bc.resize(MAX_SPECIES*LOHI*AMREX_SPACEDIM);

  initialize_compressible_namespace(Yk_bc.dataPtr(), Xk_bc.dataPtr());
 
}
