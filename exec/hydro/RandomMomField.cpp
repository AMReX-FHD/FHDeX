#include "common_functions.H"

void RandomMomField(std::array< MultiFab, AMREX_SPACEDIM >& rand_mom_add,
                    const std::array< MultiFab, AMREX_SPACEDIM >& temp_face,
                    const amrex::Geometry& geom,
                    const amrex::Real& dt) {


    // Fill random number from N(0,sqrt(dt))
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        amrex::FillRandomNormal (rand_mom_add[d],0,1,0.0,std::sqrt(dt));
    }

    const Real* dx = geom.CellSize();
    Real dVol = (AMREX_SPACEDIM==2) ? dx[0]*dx[1]*cell_depth : dx[0]*dx[1]*dx[2];

    // rand_mom_add equilibrium standard deviation
    std::array< MultiFab, AMREX_SPACEDIM >  eq_mom_std;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
            eq_mom_std[d].define( rand_mom_add[d].boxArray(),  rand_mom_add[d].DistributionMap(), 1, 0);
            eq_mom_std[d].setVal(std::sqrt(variance_coef_mom_scaling*k_B/dVol));
    }

    for ( MFIter mfi(rand_mom_add[0],false); mfi.isValid(); ++mfi ) {

        AMREX_D_TERM(const Array4<Real const> & Tx = temp_face[0].array(mfi);,
                     const Array4<Real const> & Ty = temp_face[1].array(mfi);,
                     const Array4<Real const> & Tz = temp_face[2].array(mfi););
        
        AMREX_D_TERM(const Array4<Real const> & eqx = eq_mom_std[0].array(mfi);,
                     const Array4<Real const> & eqy = eq_mom_std[1].array(mfi);,
                     const Array4<Real const> & eqz = eq_mom_std[2].array(mfi););
        
        AMREX_D_TERM(const Array4<Real>       & momx = rand_mom_add[0].array(mfi);,
                     const Array4<Real>       & momy = rand_mom_add[1].array(mfi);,
                     const Array4<Real>       & momz = rand_mom_add[2].array(mfi););

        AMREX_D_TERM(const Box & bx_x = mfi.tilebox(nodal_flag_x);,
                     const Box & bx_y = mfi.tilebox(nodal_flag_y);,
                     const Box & bx_z = mfi.tilebox(nodal_flag_z););

        amrex::ParallelFor(bx_x, bx_y,
#if (AMREX_SPACEDIM == 3)
                           bx_z,
#endif
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                momx(i,j,k) *= eqx(i,j,k) * std::sqrt(Tx(i,j,k));
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                momy(i,j,k) *= eqy(i,j,k) * std::sqrt(Ty(i,j,k));
            }
#if (AMREX_SPACEDIM == 3)
          , [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                momz(i,j,k) *= eqz(i,j,k) * std::sqrt(Tz(i,j,k));
            }
#endif
            );
    }

    // sync up random numbers at boundaries and ghost cells
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        rand_mom_add[d].OverrideSync(geom.periodicity());
        rand_mom_add[d].FillBoundary(geom.periodicity());
    }
    
    // Ensure zero total momentum
    Vector<Real> av_mom(AMREX_SPACEDIM, 0.0);
    // take staggered sum & divide by number of cells
    SumStag(rand_mom_add,av_mom,true);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        // subtract off average
        rand_mom_add[d].plus(-av_mom[d],0,1,0);
        rand_mom_add[d].OverrideSync(geom.periodicity());
        rand_mom_add[d].FillBoundary(geom.periodicity());
    }
}

