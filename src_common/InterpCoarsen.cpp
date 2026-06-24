#include "common_functions.H"


// Divergence free coarsening of vel, assumes two levels with ratio 2

void FaceFillCoarse(Vector<std::array< MultiFab, AMREX_SPACEDIM >>& mf, int map)
{
    std::array< MultiFab*, AMREX_SPACEDIM > mfCoarse;
    std::array< MultiFab*, AMREX_SPACEDIM > mfFine;

    mfCoarse = GetArrOfPtrs(mf[0]);
    mfFine = GetArrOfPtrs(mf[1]);

    FaceFillCoarse(mfCoarse,mfFine,map);
}

void FaceFillCoarse(std::array<MultiFab, AMREX_SPACEDIM>* & mf, int map)
{

    std::array< MultiFab*, AMREX_SPACEDIM > mfCoarse;
    std::array< MultiFab*, AMREX_SPACEDIM > mfFine;

    mfCoarse = GetArrOfPtrs(mf[0]);
    mfFine = GetArrOfPtrs(mf[1]);

    FaceFillCoarse(mfCoarse,mfFine,map);
}

void FaceFillCoarse(std::array<MultiFab*, AMREX_SPACEDIM> & mfCoarse,std::array<MultiFab*, AMREX_SPACEDIM> & mfFine, int map)
{
    BL_PROFILE_VAR("FaceFillCoarse()",FaceFillCoarse);

    const IntVect ratio(AMREX_D_DECL(2,2,2));

    if(map==0)
    {
        const std::array<MultiFab*, AMREX_SPACEDIM> coarse = {AMREX_D_DECL(mfCoarse[0],mfCoarse[1],mfCoarse[2])};
        const std::array<const MultiFab*, AMREX_SPACEDIM> fine = {AMREX_D_DECL(mfFine[0],mfFine[1],mfFine[2])};
        amrex::average_down_faces(fine,coarse,ratio);
    }
    else if(map == 1)
    {

        std::array< MultiFab, AMREX_SPACEDIM > coarseTemp;
        const BoxArray& fine_BA = mfFine[0]->boxArray();
        BoxArray coarsened_BA = fine_BA;
        coarsened_BA.coarsen(2);

        for(int d=0;d<AMREX_SPACEDIM;d++)
        {

            coarseTemp[d].define(convert(coarsened_BA,nodal_flag_dir[d]), mfFine[d]->DistributionMap(), 1, mfFine[d]->nGrow());

            coarseTemp[d].setVal(0.0);
        }

        for ( MFIter mfi(coarseTemp[0],TilingIfNotGPU()); mfi.isValid(); ++mfi ) {


        AMREX_D_TERM(Box bx_x = mfi.tilebox(nodal_flag_x);,
                     Box bx_y = mfi.tilebox(nodal_flag_y);,
                     Box bx_z = mfi.tilebox(nodal_flag_z););

        AMREX_D_TERM(Array4<Real> const& phix_c_fab = coarseTemp[0].array(mfi);,
                     Array4<Real> const& phiy_c_fab = coarseTemp[1].array(mfi);,
                     Array4<Real> const& phiz_c_fab = coarseTemp[2].array(mfi););

        AMREX_D_TERM(Array4<Real const> const& phix_f_fab = mfFine[0]->array(mfi);,
                     Array4<Real const> const& phiy_f_fab = mfFine[1]->array(mfi);,
                     Array4<Real const> const& phiz_f_fab = mfFine[2]->array(mfi););


#if (AMREX_SPACEDIM == 3)
        amrex::ParallelFor(bx_x, bx_y, bx_z, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            phix_c_fab(i,j,k) = 0.125* ( phix_f_fab(2*i  ,2*j,2*k  ) + phix_f_fab(2*i  ,2*j+1,2*k  )
                                        +phix_f_fab(2*i  ,2*j,2*k+1) + phix_f_fab(2*i  ,2*j+1,2*k+1) )
                             + 0.0625* ( phix_f_fab(2*i+1,2*j,2*k  ) + phix_f_fab(2*i+1,2*j+1,2*k  )
                                        +phix_f_fab(2*i+1,2*j,2*k+1) + phix_f_fab(2*i+1,2*j+1,2*k+1) )
                             + 0.0625* ( phix_f_fab(2*i-1,2*j,2*k  ) + phix_f_fab(2*i-1,2*j+1,2*k  )
                                        +phix_f_fab(2*i-1,2*j,2*k+1) + phix_f_fab(2*i-1,2*j+1,2*k+1) );
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            phiy_c_fab(i,j,k) = 0.125* ( phiy_f_fab(2*i,2*j  ,2*k  ) + phiy_f_fab(2*i+1,2*j  ,2*k  )
                                        +phiy_f_fab(2*i,2*j  ,2*k+1) + phiy_f_fab(2*i+1,2*j  ,2*k+1) )
                             + 0.0625* ( phiy_f_fab(2*i,2*j+1,2*k  ) + phiy_f_fab(2*i+1,2*j+1,2*k  )
                                        +phiy_f_fab(2*i,2*j+1,2*k+1) + phiy_f_fab(2*i+1,2*j+1,2*k+1) )
                             + 0.0625* ( phiy_f_fab(2*i,2*j-1,2*k  ) + phiy_f_fab(2*i+1,2*j-1,2*k  )
                                        +phiy_f_fab(2*i,2*j-1,2*k+1) + phiy_f_fab(2*i+1,2*j-1,2*k+1) );
        },
                           [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            phiz_c_fab(i,j,k) = 0.125* ( phiz_f_fab(2*i,2*j  ,2*k  ) + phiz_f_fab(2*i+1,2*j  ,2*k  )
                                        +phiz_f_fab(2*i,2*j+1,2*k  ) + phiz_f_fab(2*i+1,2*j+1,2*k  ) )
                             + 0.0625* ( phiz_f_fab(2*i,2*j  ,2*k+1) + phiz_f_fab(2*i+1,2*j  ,2*k+1)
                                        +phiz_f_fab(2*i,2*j+1,2*k+1) + phiz_f_fab(2*i+1,2*j+1,2*k+1) )
                             + 0.0625* ( phiz_f_fab(2*i,2*j  ,2*k-1) + phiz_f_fab(2*i+1,2*j  ,2*k-1)
                                        +phiz_f_fab(2*i,2*j+1,2*k-1) + phiz_f_fab(2*i+1,2*j+1,2*k-1) );
        });
#else
        amrex::Abort("FaceFillCoarse for 2D needs to be written");
#endif

        }
        for(int d=0;d<AMREX_SPACEDIM;d++)
        {
            mfCoarse[d]->ParallelCopy(coarseTemp[d], 0, 0, 1, 0, 0);
        }
    }
}

void CellFillCoarse(Vector<MultiFab>& mf, Vector<Geometry> geom)
{

    MultiFab* mfp = &mf[0];
    Geometry* geomp = &geom[0];
    CellFillCoarse(mfp,geomp);
}

void CellFillCoarse(Vector<MultiFab>& mf, Geometry* geom)
{
    MultiFab* mfp = &mf[0];

    CellFillCoarse(mfp,geom);
}

void CellFillCoarse(MultiFab* & mf, Geometry* geom)
{
    BL_PROFILE_VAR("CellFillCoarse()",CellFillCoarse);

    const IntVect ratio(AMREX_D_DECL(2,2,2));
    amrex::average_down(mf[1],mf[0], geom[1], geom[0],0,1,ratio);

}

void CellFillCoarse(const MultiFab & mfFine, Geometry geomFine, MultiFab & mfCoarse, Geometry geomCoarse)
{
    BL_PROFILE_VAR("CellFillCoarse()",CellFillCoarse);

    const IntVect ratio(AMREX_D_DECL(2,2,2));
    amrex::average_down(mfFine,mfCoarse, geomFine, geomCoarse,0,1,ratio);

}


