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

        Array4<Real> const& div_fab = div.array(mfi);
        AMREX_D_TERM(Array4<Real const> const& phix_fab = phi_fc[0].array(mfi);,
                     Array4<Real const> const& phiy_fab = phi_fc[1].array(mfi);,
                     Array4<Real const> const& phiz_fab = phi_fc[2].array(mfi););

        if (increment == 0) {
            AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
            {                    
                div_fab(i,j,k,start_outcomp+n) =
                    AMREX_D_TERM(  (phix_fab(i+1,j,k,start_incomp+n) - phix_fab(i,j,k,start_incomp+n)) / dx[0],
                                 + (phiy_fab(i,j+1,k,start_incomp+n) - phiy_fab(i,j,k,start_incomp+n)) / dx[1],
                                 + (phiz_fab(i,j,k+1,start_incomp+n) - phiz_fab(i,j,k,start_incomp+n)) / dx[2]);;
            });
        }
        else
        {
            AMREX_HOST_DEVICE_FOR_4D(bx, ncomp, i, j, k, n,
            {                 
                div_fab(i,j,k,start_outcomp+n) +=
                    AMREX_D_TERM(  (phix_fab(i+1,j,k,start_incomp+n) - phix_fab(i,j,k,start_incomp+n)) / dx[0],
                                 + (phiy_fab(i,j+1,k,start_incomp+n) - phiy_fab(i,j,k,start_incomp+n)) / dx[1],
                                 + (phiz_fab(i,j,k+1,start_incomp+n) - phiz_fab(i,j,k,start_incomp+n)) / dx[2]);;
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
        
        Array4<Real const> const& phi_fab = phi.array(mfi);
        
        AMREX_D_TERM(Array4<Real> const& gphix_fab = gphi[0].array(mfi);,
                     Array4<Real> const& gphiy_fab = gphi[1].array(mfi);,
                     Array4<Real> const& gphiz_fab = gphi[2].array(mfi););
        
        AMREX_D_TERM(const Box& bx_x = mfi.nodaltilebox(0);,
                     const Box& bx_y = mfi.nodaltilebox(1);,
                     const Box& bx_z = mfi.nodaltilebox(2););

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

    const Real* dx = geom.CellSize(); 

    for ( MFIter mfi(phi); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.validbox();
        
        Array4<Real const> const& phi_fab = phi.array(mfi);

        AMREX_D_TERM(Array4<Real> const& gphix_fab = gphi[0].array(mfi);,
                     Array4<Real> const& gphiy_fab = gphi[1].array(mfi);,
                     Array4<Real> const& gphiz_fab = gphi[2].array(mfi););

        AMREX_HOST_DEVICE_FOR_3D(bx, i, j, k,
        {
            AMREX_D_TERM(gphix_fab(i,j,k) = (phi_fab(i+1,j,k) - phi_fab(i-1,j,k) ) / (2.*dx[0]);,
                         gphiy_fab(i,j,k) = (phi_fab(i,j+1,k) - phi_fab(i,j-1,k) ) / (2.*dx[1]);,
                         gphiz_fab(i,j,k) = (phi_fab(i,j,k+1) - phi_fab(i,j,k-1) ) / (2.*dx[2]););
        });
    }
}
