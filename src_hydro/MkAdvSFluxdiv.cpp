#include "hydro_functions.H"
#include "hydro_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"
#include "common_namespace.H"

using namespace amrex;
using namespace common;

void MkAdvSFluxdiv(const std::array<MultiFab, AMREX_SPACEDIM>& umac,
		   const MultiFab& m,
		   MultiFab& m_update,
		   const amrex::Real* dx,
		   const Geometry& geom,
		   const int& increment)
{

     BL_PROFILE_VAR("MkAdvSFluxdiv()",MkAdvSFluxdiv);

     std::array<MultiFab, AMREX_SPACEDIM> m_fc;
     DistributionMapping dmap = m.DistributionMap();
     AMREX_D_TERM(m_fc[0].define(convert(m.boxArray(),nodal_flag_x), dmap, 1, 1);,
      	          m_fc[1].define(convert(m.boxArray(),nodal_flag_y), dmap, 1, 1);,
	          m_fc[2].define(convert(m.boxArray(),nodal_flag_z), dmap, 1, 1););

     AverageCCToFace(m, 0, m_fc, 0, 1);

     // AMREX_D_TERM(m_fc[0].FillBoundary(geom.periodicity());,
     //    	  m_fc[1].FillBoundary(geom.periodicity());,
     //    	  m_fc[2].FillBoundary(geom.periodicity()););

     // //TODO: this will break in 2D
     // setBC(m_fc[0], m_fc[1], m_fc[2]);

     for (int d=0; d<AMREX_SPACEDIM; ++d) {
         m_fc[d].FillBoundary(geom.periodicity());
         MultiFABPhysBC(m_fc[d]);
     }

     // Loop over boxes
     for (MFIter mfi(umac[0]); mfi.isValid(); ++mfi) {

         // Create cell-centered box from semi-nodal box
         const Box& validBox_cc = enclosedCells(mfi.validbox());

         mk_advective_s_fluxdiv(ARLIM_3D(validBox_cc.loVect()), ARLIM_3D(validBox_cc.hiVect()),
                                BL_TO_FORTRAN_ANYD(umac[0][mfi]),
                                BL_TO_FORTRAN_ANYD(umac[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                                BL_TO_FORTRAN_ANYD(umac[2][mfi]),
#endif
                                BL_TO_FORTRAN_ANYD(m_fc[0][mfi]),
                                BL_TO_FORTRAN_ANYD(m_fc[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                                BL_TO_FORTRAN_ANYD(m_fc[2][mfi]),
#endif
                                BL_TO_FORTRAN_ANYD(m[mfi]),
                                BL_TO_FORTRAN_ANYD(m_update[mfi]),
                                dx, &increment);
     }

}
