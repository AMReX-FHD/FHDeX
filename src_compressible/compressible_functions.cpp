#include "compressible_functions.H"


void InitializeCompressibleNamespace() {

    BL_PROFILE_VAR("InitializeCompressibleNamespace()",InitializeCompressibleNamespace);

    bc_Yk.resize(MAX_SPECIES*LOHI*AMREX_SPACEDIM);
    bc_Xk.resize(MAX_SPECIES*LOHI*AMREX_SPACEDIM);

    initialize_compressible_namespace(bc_Yk.dataPtr(), bc_Xk.dataPtr(), &plot_means, &plot_vars);
}
