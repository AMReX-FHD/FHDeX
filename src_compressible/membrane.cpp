#include "compressible_functions.H"
#include "compressible_functions_F.H"
#include "common_functions.H"


void doMembrane(MultiFab& cons, MultiFab& prim, std::array<MultiFab, AMREX_SPACEDIM>& flux, const amrex::Geometry geom, const amrex::Real* dx, const amrex::Real dt)
{

    AMREX_D_TERM(flux[0].setVal(0.0);,
                 flux[1].setVal(0.0);,
                 flux[2].setVal(0.0););

    // Loop over boxes
    for ( MFIter mfi(cons); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

//        do_ssa(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
//                   cons[mfi].dataPtr(),  
//                   prim[mfi].dataPtr(), flux[0][mfi].dataPtr(), dx, &dt);

        do_langevin(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                   cons[mfi].dataPtr(),  
                   prim[mfi].dataPtr(), flux[0][mfi].dataPtr(), dx, &dt);          
    }

    flux[0].OverrideSync(geom.periodicity());

    for ( MFIter mfi(cons); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        apply_effusion(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                   cons[mfi].dataPtr(),  
                   flux[0][mfi].dataPtr(), dx, &dt);       
    }

    conservedToPrimitive(prim, cons);
    cons.FillBoundary(geom.periodicity());
    prim.FillBoundary(geom.periodicity());

}

