#include "particle_functions.H"

void ComputeDryMobility(MultiFab & dryMobility, species* particleInfo, const Geometry & Geom)
{
    BL_PROFILE_VAR("ComputeDryMobility()",ComputeDryMobility);

    const Real* dx = Geom.CellSize();
    const Real* plo = Geom.ProbLo();
    const Real* phi = Geom.ProbHi();
    const int ngc = dryMobility.nGrow();

    for ( MFIter mfi(dryMobility); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.validbox();
        compute_dry_mobility(AMREX_ARLIM_3D(bx.loVect()), AMREX_ARLIM_3D(bx.hiVect()), BL_TO_FORTRAN_3D(dryMobility[mfi]), AMREX_ZFILL(dx), AMREX_ZFILL(plo), AMREX_ZFILL(phi), &ngc, particleInfo);

    }
}


