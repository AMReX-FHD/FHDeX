#include "common_functions.H"
#include "common_functions_F.H"


//Computes divergence at cell centres from velcocities at cell faces
void ComputeDiv(MultiFab& div,
                const std::array<MultiFab, AMREX_SPACEDIM>& phi_fc,
                int start_incomp, int start_outcomp, int ncomp,
                const Geometry& geom, int increment)
{

    BL_PROFILE_VAR("ComputeDiv()",ComputeDiv);

    const Real* dx = geom.CellSize(); 
    
    for ( MFIter mfi(div); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();

        const auto& div_fab = (&div)->array(mfi);
        const auto& phix_fab = (&phi_fc[0])->array(mfi);
        const auto& phiy_fab = (&phi_fc[1])->array(mfi);
#if (AMREX_SPACEDIM == 3)        
        const auto& phiz_fab = (&phi_fc[2])->array(mfi);
#endif

        if (increment == 0) {
            AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
            {
                div_fab(i,j,k,start_outcomp+n) =
                      (phix_fab(i+1,j,k,start_incomp+n) - phix_fab(i,j,k,start_incomp+n)) / dx[0]
                    + (phiy_fab(i,j+1,k,start_incomp+n) - phiy_fab(i,j,k,start_incomp+n)) / dx[1]
#if (AMREX_SPACEDIM == 3)
                    + (phiz_fab(i,j,k+1,start_incomp+n) - phiz_fab(i,j,k,start_incomp+n)) / dx[2]
#endif
                    ;
            });
        }
        else
        {
            AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
            {
                div_fab(i,j,k,start_outcomp+n) +=
                      (phix_fab(i+1,j,k,start_incomp+n) - phix_fab(i,j,k,start_incomp+n)) / dx[0]
                    + (phiy_fab(i,j+1,k,start_incomp+n) - phiy_fab(i,j,k,start_incomp+n)) / dx[1]
#if (AMREX_SPACEDIM == 3)
                    + (phiz_fab(i,j,k+1,start_incomp+n) - phiz_fab(i,j,k,start_incomp+n)) / dx[2]
#endif
                    ;
            });
        }
    }
}



//Computes gradient at cell faces of cell centred scalar
void ComputeGrad(const MultiFab& phi, std::array<MultiFab, AMREX_SPACEDIM>& gphi,
                 int start_incomp, int start_outcomp, int ncomp, const Geometry& geom)
{
    const Real* dx = geom.CellSize(); 
    
    for ( MFIter mfi(phi); mfi.isValid(); ++mfi ) {

        const Box& bx = mfi.validbox();
        
        const auto& phi_fab = (&phi)->array(mfi);
        const auto& gphix_fab = (&gphi[0]) -> array(mfi);
        const auto& gphiy_fab = (&gphi[1]) -> array(mfi);
        const Box& bx_x = amrex::growHi(bx,0,1);
        const Box& bx_y = amrex::growHi(bx,1,1);
#if (AMREX_SPACEDIM == 3)        
        const auto& gphiz_fab = (&gphi[2]) -> array(mfi);
        const Box& bx_z = amrex::growHi(bx,2,1);
#endif

        AMREX_HOST_DEVICE_FOR_4D(bx_x, ncomp, i, j, k, n,
        {
            gphix_fab(i,j,k,start_outcomp+n) =
                (phi_fab(i,j,k,start_incomp+n) - phi_fab(i-1,j,k,start_incomp+n) ) / dx[0];
        });

        AMREX_HOST_DEVICE_FOR_4D(bx_y, ncomp, i, j, k, n,
        {
            gphiy_fab(i,j,k,start_outcomp+n) =
                (phi_fab(i,j,k,start_incomp+n) - phi_fab(i,j-1,k,start_incomp+n) ) / dx[1];
        });

#if (AMREX_SPACEDIM == 3)
        AMREX_HOST_DEVICE_FOR_4D(bx_z, ncomp, i, j, k, n,
        {
            gphiz_fab(i,j,k,start_outcomp+n) =
                (phi_fab(i,j,k,start_incomp+n) - phi_fab(i,j,k-1,start_incomp+n) ) / dx[2];
        });
#endif

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
