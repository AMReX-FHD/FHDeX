#include "common_functions.H"
#include "common_functions_F.H"


//Computes divergence at cell centres from velcocities at cell faces
void ComputeDiv(MultiFab& div, std::array<MultiFab, AMREX_SPACEDIM>& phi_fc, 
                int start_incomp, int start_outcomp, int ncomp, 
                const Geometry geom, int increment)
{
    for ( MFIter mfi(div); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        for (int incomp=start_incomp; incomp<start_incomp+ncomp; ++incomp) {

            int outcomp = incomp + start_outcomp - start_incomp;

            compute_div(BL_TO_FORTRAN_BOX(bx),
                        BL_TO_FORTRAN_N_ANYD(phi_fc[0][mfi],incomp), 
                        BL_TO_FORTRAN_N_ANYD(phi_fc[1][mfi],incomp),
#if (AMREX_SPACEDIM==3)
                        BL_TO_FORTRAN_N_ANYD(phi_fc[2][mfi],incomp),
#endif
                        BL_TO_FORTRAN_N_ANYD(div[mfi],outcomp),					
                        geom.CellSize(), &increment);
        }
    }
}

//Computes gradient at cell faces of cell centred scalar
void ComputeGrad(MultiFab& phi, std::array<MultiFab, AMREX_SPACEDIM>& phi_fc, const Geometry geom)
{
    for ( MFIter mfi(phi); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        compute_grad(BL_TO_FORTRAN_BOX(bx),
                     BL_TO_FORTRAN_ANYD(phi_fc[0][mfi]), 
                     BL_TO_FORTRAN_ANYD(phi_fc[1][mfi]),
#if (AMREX_SPACEDIM==3)
                     BL_TO_FORTRAN_ANYD(phi_fc[2][mfi]),
#endif
                     BL_TO_FORTRAN_ANYD(phi[mfi]),					
                     geom.CellSize());
    }
}

