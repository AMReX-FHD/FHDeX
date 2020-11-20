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
        compute_dry_mobility(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()), BL_TO_FORTRAN_3D(dryMobility[mfi]), ZFILL(dx), ZFILL(plo), ZFILL(phi), &ngc, particleInfo);
        
    }
}


