#include "common_functions.H"
#include "common_functions_F.H"


//Computes divergence at cell centres from velcocities at cell faces
void ComputeDiv(MultiFab& div, std::array<MultiFab, AMREX_SPACEDIM>& umac, const Geometry geom, int increment)
{
    int inc = increment; 

    for ( MFIter mfi(div); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

		compute_div(BL_TO_FORTRAN_BOX(bx),
					BL_TO_FORTRAN_ANYD(umac[0][mfi]), 
                    BL_TO_FORTRAN_ANYD(umac[1][mfi]),
#if AMREX_SPACEDIM == 3
                    BL_TO_FORTRAN_ANYD(umac[2][mfi]),
#endif
                    BL_TO_FORTRAN_ANYD(div[mfi]),					
                    geom.CellSize(), &inc);
    }
}

//Computes gradient at cell faces of cell centred scalar
void ComputeGrad(MultiFab& phi, std::array<MultiFab, AMREX_SPACEDIM>& umac, const Geometry geom)
{
    for ( MFIter mfi(phi); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

		compute_grad(BL_TO_FORTRAN_BOX(bx),
					BL_TO_FORTRAN_ANYD(umac[0][mfi]), 
                    BL_TO_FORTRAN_ANYD(umac[1][mfi]),
#if AMREX_SPACEDIM == 3
                    BL_TO_FORTRAN_ANYD(umac[2][mfi]),
#endif
                    BL_TO_FORTRAN_ANYD(phi[mfi]),					
                    geom.CellSize());

    }
}

