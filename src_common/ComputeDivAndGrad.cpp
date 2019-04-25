#include "common_functions.H"
#include "common_functions_F.H"


//Computes divergence at cell centres from velcocities at cell faces
void ComputeDiv(MultiFab & div, const std::array<MultiFab, AMREX_SPACEDIM> & phi_fc,
                int start_incomp, int start_outcomp, int ncomp,
                const Geometry & geom, int increment)
{

    BL_PROFILE_VAR("ComputeDiv()",ComputeDiv);

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
                        geom.CellSize(), & increment);
        }
    }
}



//Computes gradient at cell faces of cell centred scalar
void ComputeGrad(const MultiFab& phi, std::array<MultiFab, AMREX_SPACEDIM>& gphi,
                 int start_incomp, int start_outcomp, int ncomp, const Geometry& geom)
{
    for ( MFIter mfi(phi); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        for (int incomp=start_incomp; incomp<start_incomp+ncomp; ++incomp) {

            int outcomp = incomp + start_outcomp - start_incomp;

            compute_grad(BL_TO_FORTRAN_BOX(bx),
                         BL_TO_FORTRAN_N_ANYD(gphi[0][mfi],outcomp),
                         BL_TO_FORTRAN_N_ANYD(gphi[1][mfi],outcomp),
#if (AMREX_SPACEDIM==3)
                         BL_TO_FORTRAN_N_ANYD(gphi[2][mfi],outcomp),
#endif
                         BL_TO_FORTRAN_N_ANYD(phi[mfi],incomp),
                         geom.CellSize());
        }
    }
}

//Computes gradient at cell centres from cell centred data - oututs to a three component mf.
void ComputeCentredGrad(const MultiFab& phi, std::array<MultiFab, AMREX_SPACEDIM>& gphi, const Geometry& geom)
{
    for ( MFIter mfi(phi); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();


            compute_grad_cc(BL_TO_FORTRAN_BOX(bx),
                         BL_TO_FORTRAN_3D(gphi[0][mfi]),
                         BL_TO_FORTRAN_3D(gphi[1][mfi]),
#if (AMREX_SPACEDIM==3)
                         BL_TO_FORTRAN_3D(gphi[2][mfi]),
#endif
                         BL_TO_FORTRAN_3D(phi[mfi]),
                         geom.CellSize());
        }    
}
