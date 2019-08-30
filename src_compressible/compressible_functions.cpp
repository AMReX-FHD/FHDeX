#include "compressible_functions.H"
#include "compressible_functions_F.H"
#include "compressible_namespace.H"

using namespace compressible;

void InitializeCompressibleNamespace() {

  bc_Yk.resize(MAX_SPECIES*LOHI*AMREX_SPACEDIM);
  bc_Xk.resize(MAX_SPECIES*LOHI*AMREX_SPACEDIM);

  initialize_compressible_namespace(bc_Yk.dataPtr(), bc_Xk.dataPtr());
}
