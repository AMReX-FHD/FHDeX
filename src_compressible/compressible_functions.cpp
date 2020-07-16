#include "compressible_functions.H"


void InitializeCompressibleNamespace() {

    BL_PROFILE_VAR("InitializeCompressibleNamespace()",InitializeCompressibleNamespace);

    bc_Yk.resize(MAX_SPECIES*LOHI*AMREX_SPACEDIM);
    bc_Xk.resize(MAX_SPECIES*LOHI*AMREX_SPACEDIM);

    initialize_compressible_namespace(bc_Yk.dataPtr(), bc_Xk.dataPtr(), &plot_means, &plot_vars);
}


void InitConsVar(MultiFab& cons, const MultiFab& prim,
                 const amrex::Geometry geom) {

    const Real* dx = geom.CellSize();
    const RealBox& realDomain = geom.ProbDomain();

    // initialize conserved variables
    for ( MFIter mfi(cons); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();
        init_consvar(BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_ANYD(cons[mfi]),
                     BL_TO_FORTRAN_ANYD(prim[mfi]),
                     dx, ZFILL(realDomain.lo()), ZFILL(realDomain.hi()));
    }

}
