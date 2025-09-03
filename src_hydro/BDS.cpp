#include "hydro_functions.H"

void BDS(MultiFab& s_update,
         const int ncomp,
         const int bccomp,
         MultiFab const& state,
         std::array< MultiFab, AMREX_SPACEDIM >& umac,
         MultiFab const& fq,
         Geometry const& geom,
         const Real dt,
         const int proj_type)
{

    BL_PROFILE("BDS_ComputeAofs()");

#if (AMREX_SPACEDIM==2)
    if ( geom.IsRZ() )
    {
        Abort("BDS does not currently support cylindrical coordinates");
    }
#endif

    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    GpuArray<Real, AMREX_SPACEDIM> dxinv;
    AMREX_D_TERM( dxinv[0] = 1./dx[0];,
                  dxinv[1] = 1./dx[1];,
                  dxinv[2] = 1./dx[2];);

    BoxArray ba = s_update.boxArray();
    DistributionMapping dmap = s_update.DistributionMap();

    std::array< MultiFab, AMREX_SPACEDIM > sedge;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
      sedge[d].define(convert(ba,nodal_flag_dir[d]), dmap, ncomp, 1);
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(s_update,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {

        const Box& bx   = mfi.tilebox();

        //
        // Get handlers to Array4
        //
        AMREX_D_TERM( const auto& sedgex = sedge[0].array(mfi);,
                      const auto& sedgey = sedge[1].array(mfi);,
                      const auto& sedgez = sedge[2].array(mfi););

        AMREX_D_TERM( const auto& u = umac[0].const_array(mfi);,
                      const auto& v = umac[1].const_array(mfi);,
                      const auto& w = umac[2].const_array(mfi););

        BDS_ComputeEdgeState(bx, ncomp, bccomp,
                             state.array(mfi),
                             AMREX_D_DECL(sedgex, sedgey, sedgez),
                             AMREX_D_DECL(u, v, w),
                             fq.array(mfi),
                             geom, dt, proj_type);

        // flip the sign to return div
        auto const& s_update_arr  = s_update.array(mfi);
        amrex::ParallelFor(bx, ncomp, [=]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
          s_update_arr( i, j, k, n ) -= (
                  (sedgex(i+1,j,k,n)*u(i+1,j,k)-sedgex(i,j,k,n)*u(i,j,k))*dxinv[0]
               +(sedgey(i,j+1,k,n)*v(i,j+1,k)-sedgey(i,j,k,n)*v(i,j,k))*dxinv[1]
#if (AMREX_SPACEDIM == 3)
               +(sedgez(i,j,k+1,n)*w(i,j,k+1)-sedgez(i,j,k,n)*w(i,j,k))*dxinv[2]
#endif
                                         );
        });

        //
        // NOTE this sync cannot protect temporaries in ComputeEdgeState, ComputeFluxes
        // or ComputeDivergence, since functions have their own scope. As soon as the
        // CPU hits the end of the function, it will call the destructor for all
        // temporaries created in that function.
        //
        Gpu::streamSynchronize();  // otherwise we might be using too much memory
    }

}

#if (AMREX_SPACEDIM == 2)

/**
 * Uses the Bell-Dawson-Shubin (BDS) algorithm, a higher order Godunov
 * method for scalar conservation laws in two dimensions, to compute
 * edge states.
 *
 * \param [in]     bx          Current grid patch
 * \param [in]     ncomp       Number of components to work on
 * \param [in]     q           Array4 of state, starting at component of interest
 * \param [in,out] xedge       Array4 containing x-edges, starting at component of interest
 * \param [in,out] yedge       Array4 containing y-edges, starting at component of interest
 * \param [in]     umac        x-Face velocities.
 * \param [in]     vmac        y-Face velocities.
 * \param [in]     fq          Array4 for forces, starting at component of interest
 * \param [in]     geom        Level geometry.
 * \param [in]     l_dt        Time step.
 *
 */

void BDS_ComputeEdgeState(Box const& bx,
                          int ncomp,
                          int bccomp,
                          Array4<Real const> const& q,
                          Array4<Real      > const& xedge,
                          Array4<Real      > const& yedge,
                          Array4<Real const> const& umac,
                          Array4<Real const> const& vmac,
                          Array4<Real const> const& fq,
                          Geometry geom,
                          Real l_dt,
                          const int proj_type)
{
    // For now, loop on components here
    for( int icomp = 0; icomp < ncomp; ++icomp)
    {
        // Temporary slope container per components
        Box const& bxg1 = amrex::grow(bx,1);
        FArrayBox slopefab(bxg1,3);
        Elixir slopeeli = slopefab.elixir();

        BDS_ComputeSlopes(bx, geom, icomp, bccomp,
                           q, slopefab.array());

        BDS_ComputeConc(bx, geom, icomp, bccomp,
                         q, xedge, yedge, slopefab.array(),
                         umac, vmac, fq, l_dt);
    }

    // project edge states to satisfy EOS
    BDS_Proj(bx,ncomp,xedge,yedge,proj_type);
}

/**
 * Compute bilinear slopes for BDS algorithm.
 *
 * \param [in]  bx      Current grid patch
 * \param [in]  geom    Level geometry.
 * \param [in]  icomp   Component of the state Array4.
 * \param [in]  s       Array4<const> of state vector.
 * \param [out] slopes  Array4 to store slope information.
 *
 */

void BDS_ComputeSlopes(Box const& bx,
                       const Geometry& geom,
                       int icomp, int bccomp,
                       Array4<Real const> const& s,
                       Array4<Real      > const& slopes)
{
    constexpr amrex::Real eps = 1.0e-8;
    bool limit_slopes = (advection_type == 2) ? true : false;

    // Define container for the nodal interpolated state
    Box const& ngbx = amrex::grow(amrex::convert(bx,IntVect(AMREX_D_DECL(1,1,1))),1);
    FArrayBox tmpnodefab(ngbx,1);
    Elixir tmpeli = tmpnodefab.elixir();
    auto const& sint = tmpnodefab.array();

    Box const& gbx = amrex::grow(bx,1);
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    Real hx = dx[0];
    Real hy = dx[1];

    Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    // compute mathematical boundary conditions
    Vector<int> bc_lo(AMREX_SPACEDIM);
    Vector<int> bc_hi(AMREX_SPACEDIM);
    BCPhysToMath(bccomp,bc_lo,bc_hi);

    bool lo_x_physbc = (bc_lo[0] == amrex::BCType::foextrap || bc_lo[0] == amrex::BCType::ext_dir) ? true : false;
    bool hi_x_physbc = (bc_hi[0] == amrex::BCType::foextrap || bc_hi[0] == amrex::BCType::ext_dir) ? true : false;
    bool lo_y_physbc = (bc_lo[1] == amrex::BCType::foextrap || bc_lo[1] == amrex::BCType::ext_dir) ? true : false;
    bool hi_y_physbc = (bc_hi[1] == amrex::BCType::foextrap || bc_hi[1] == amrex::BCType::ext_dir) ? true : false;

    // bicubic interpolation to corner points
    // (i,j,k) refers to lower corner of cell
    // Added k index -- placeholder for 2d
    ParallelFor(ngbx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

        // set node values equal to the average of the ghost cell values since they store the physical condition on the boundary
        if ( i<=dlo.x && lo_x_physbc ) {
            sint(i,j,k) = 0.5*(s(dlo.x-1,j,k,icomp) + s(dlo.x-1,j-1,k,icomp));
            return;
        }
        if ( i>=dhi.x+1 && hi_x_physbc ) {
            sint(i,j,k) = 0.5*(s(dhi.x+1,j,k,icomp) + s(dhi.x+1,j-1,k,icomp));
            return;
        }
        if ( j<=dlo.y && lo_y_physbc ) {
            sint(i,j,k) = 0.5*(s(i,dlo.y-1,k,icomp) + s(i-1,dlo.y-1,k,icomp));
            return;
        }
        if ( j>=dhi.y+1 && hi_y_physbc ) {
            sint(i,j,k) = 0.5*(s(i,dhi.y+1,k,icomp) + s(i-1,dhi.y+1,k,icomp));
            return;
        }

        // one cell inward from any physical boundary, revert to 4-point average
        if ( (i==dlo.x+1 && lo_x_physbc) ||
             (i==dhi.x   && hi_x_physbc) ||
             (j==dlo.y+1 && lo_y_physbc) ||
             (j==dhi.y   && hi_y_physbc) ) {

            sint(i,j,k) = 0.25* (s(i,j,k,icomp) + s(i-1,j,k,icomp) + s(i,j-1,k,icomp) + s(i-1,j-1,k,icomp));
            return;
        }

        sint(i,j,k) = (s(i-2,j-2,k,icomp) + s(i-2,j+1,k,icomp) + s(i+1,j-2,k,icomp) + s(i+1,j+1,k,icomp)
                - 7.0*(s(i-2,j-1,k,icomp) + s(i-2,j  ,k,icomp) + s(i-1,j-2,k,icomp) + s(i  ,j-2,k,icomp) +
                       s(i-1,j+1,k,icomp) + s(i  ,j+1,k,icomp) + s(i+1,j-1,k,icomp) + s(i+1,j  ,k,icomp))
               + 49.0*(s(i-1,j-1,k,icomp) + s(i  ,j-1,k,icomp) + s(i-1,j  ,k,icomp) + s(i  ,j  ,k,icomp)) ) / 144.0;
    });

    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
        // compute initial estimates of slopes from unlimited corner points

        // local variables
        Real sumloc, redfac, redmax, div, kdp, sumdif, sgndif;

        Array1D<Real, 1, 4> diff;
        Array1D<Real, 1, 4> smin;
        Array1D<Real, 1, 4> smax;
        Array1D<Real, 1, 4> sc;

        Array1D<bool, 1, 4> allow_change;
        for (int mm=1; mm<=4; ++mm) {
            allow_change(mm) = true;
        }

        if ( i<=dlo.x && lo_x_physbc ) {
            allow_change(1) = false;
            allow_change(2) = false;
        }
        if ( i>=dhi.x+1 && hi_x_physbc ) {
            allow_change(3) = false;
            allow_change(4) = false;
        }
        if ( j<=dlo.y && lo_y_physbc ) {
            allow_change(1) = false;
            allow_change(3) = false;
        }
        if ( j>=dhi.y+1 && hi_y_physbc ) {
            allow_change(2) = false;
            allow_change(4) = false;
        }

        // compute initial estimates of slopes from unlimited corner points
        // sx
        slopes(i,j,k,0) = 0.5*(sint(i+1,j+1,k) + sint(i+1,j,k) - sint(i,j+1,k) - sint(i,j,k)) / hx;
        // sy
        slopes(i,j,k,1) = 0.5*(sint(i+1,j+1,k) - sint(i+1,j,k) + sint(i,j+1,k) - sint(i,j,k)) / hy;
        // sxy
        slopes(i,j,k,2) =     (sint(i+1,j+1,k) - sint(i+1,j,k) - sint(i,j+1,k) + sint(i,j,k)) / (hx*hy);

        if (limit_slopes) {

            // ++ / sint(i+1,j+1)
            sc(4) = s(i,j,k,icomp) + 0.5*(hx*slopes(i,j,k,0) + hy*slopes(i,j,k,1)) + 0.25*hx*hy*slopes(i,j,k,2);

            // +- / sint(i+1,j  )
            sc(3) = s(i,j,k,icomp) + 0.5*(hx*slopes(i,j,k,0) - hy*slopes(i,j,k,1)) - 0.25*hx*hy*slopes(i,j,k,2);

            // -+ / sint(i  ,j+1)
            sc(2) = s(i,j,k,icomp) - 0.5*(hx*slopes(i,j,k,0) - hy*slopes(i,j,k,1)) - 0.25*hx*hy*slopes(i,j,k,2);

            // -- / sint(i  ,j  )
            sc(1) = s(i,j,k,icomp) - 0.5*(hx*slopes(i,j,k,0) + hy*slopes(i,j,k,1)) + 0.25*hx*hy*slopes(i,j,k,2);

            // enforce max/min bounds
            smin(4) = amrex::min(s(i,j,k,icomp), s(i+1,j,k,icomp), s(i,j+1,k,icomp), s(i+1,j+1,k,icomp));
            smax(4) = amrex::max(s(i,j,k,icomp), s(i+1,j,k,icomp), s(i,j+1,k,icomp), s(i+1,j+1,k,icomp));

            smin(3) = amrex::min(s(i,j,k,icomp), s(i+1,j,k,icomp), s(i,j-1,k,icomp), s(i+1,j-1,k,icomp));
            smax(3) = amrex::max(s(i,j,k,icomp), s(i+1,j,k,icomp), s(i,j-1,k,icomp), s(i+1,j-1,k,icomp));

            smin(2) = amrex::min(s(i,j,k,icomp), s(i-1,j,k,icomp), s(i,j+1,k,icomp), s(i-1,j+1,k,icomp));
            smax(2) = amrex::max(s(i,j,k,icomp), s(i-1,j,k,icomp), s(i,j+1,k,icomp), s(i-1,j+1,k,icomp));

            smin(1) = amrex::min(s(i,j,k,icomp), s(i-1,j,k,icomp), s(i,j-1,k,icomp), s(i-1,j-1,k,icomp));
            smax(1) = amrex::max(s(i,j,k,icomp), s(i-1,j,k,icomp), s(i,j-1,k,icomp), s(i-1,j-1,k,icomp));

            for(int mm=1; mm<=4; ++mm){
               if (allow_change(mm)) {
                   sc(mm) = amrex::max(amrex::min(sc(mm), smax(mm)), smin(mm));
               }
            }

            // iterative loop
            for(int ll=1; ll<=3; ++ll){

               // compute the amount by which the average of the nodal values differs from cell-center value
               sumloc = 0.25*(sc(4) + sc(3) + sc(2) + sc(1));
               sumdif = (sumloc - s(i,j,k,icomp))*4.0;

               // sgndif = +(-)1 if the node average is too large(small)
               sgndif = std::copysign(1.0,sumdif);

               // compute how much each node is larger(smaller) than the cell-centered value
               for(int mm=1; mm<=4; ++mm){
                  diff(mm) = (sc(mm) - s(i,j,k,icomp))*sgndif;
               }

               kdp = 0;

               // count how many nodes are larger(smaller) than the cell-centered value
               for(int mm=1; mm<=4; ++mm){
                  if (diff(mm) > eps && allow_change(mm)) {
                     kdp = kdp+1;
                  }
               }

               // adjust node values
               for(int mm=1; mm<=4; ++mm){

                  // don't allow boundary nodes to change value
                  if (!allow_change(mm)) continue;

                  // how many node values are left to potentially adjust
                  if (kdp<1) {
                     div = 1.0;
                  } else {
                     div = kdp;
                  }

                  // if the node needs adjusting, figure out by how much the remaining sum is divy'ed up
                  if (diff(mm)>eps) {
                     redfac = sumdif*sgndif/div;
                     kdp = kdp-1;
                  } else {
                     redfac = 0.0;
                  }

                  // don't let the adjustment introduce any new extrema
                  if (sgndif > 0.0) {
                     redmax = sc(mm) - smin(mm);
                  } else {
                     redmax = smax(mm) - sc(mm);
                  }
                  redfac = amrex::min(redfac,redmax);

                  // adjust nodal value and decrement the excess
                  sumdif = sumdif - redfac*sgndif;
                  sc(mm) = sc(mm) - redfac*sgndif;
               }
            }

            // final slopes
            // sx
            slopes(i,j,k,0) = 0.5*( sc(4) + sc(3) -sc(1) - sc(2) )/hx;
            // sy
            slopes(i,j,k,1) = 0.5*( sc(4) + sc(2) -sc(1) - sc(3) )/hy;
            // sxy
            slopes(i,j,k,2) =     ( sc(1) + sc(4) -sc(2) - sc(3) )/(hx*hy);
        }
    });
}

/**
 * Returns updated edge value.
 *
 * \param [in] s
 * \param [in] slope
 * \param [in] del
 *
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real eval(const Real s,
          Array1D<Real,1,3>& slope,
          Array1D<Real,1,2>& del)
{
    Real val = s + del(1)*slope(1) + del(2)*slope(2) + del(1)*del(2)*slope(3);

    return val;
}

/**
 * Compute Conc for BDS algorithm.
 *
 * \param [in]     bx          Current grid patch
 * \param [in]     geom        Level geometry.
 * \param [in]     icomp       Component of the Array4s.
 * \param [in]     s           Array4 of state.
 * \param [in,out] sedgex      Array4 containing x-edges.
 * \param [in,out] sedgey      Array4 containing y-edges.
 * \param [in]     slopes      Array4 containing slope information.
 * \param [in]     umac        Array4 for u-face velocity.
 * \param [in]     vmac        Array4 for v-face velocity.
 * \param [in]     force       Array4 for forces.
 * \param [in]     dt          Time step.
 *
 */

void BDS_ComputeConc(Box const& bx,
                     const Geometry& geom,
                     int icomp, int bccomp,
                     Array4<Real const> const& s,
                     Array4<Real      > const& sedgex,
                     Array4<Real      > const& sedgey,
                     Array4<Real const> const& slopes,
                     Array4<Real const> const& umac,
                     Array4<Real const> const& vmac,
                     Array4<Real const> const& force,
                     const Real dt)
{
    Box const& gbx = amrex::grow(bx,1);
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    FArrayBox ux_fab(gbx,1);
    FArrayBox vy_fab(gbx,1);
    Elixir uxeli = ux_fab.elixir();
    Elixir vyeli = vy_fab.elixir();
    auto const& ux    = ux_fab.array();
    auto const& vy    = vy_fab.array();

    Real hx = dx[0];
    Real hy = dx[1];

    Real dt2 = dt/2.0;
    Real dt3 = dt/3.0;

    Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    // compute mathematical boundary conditions
    Vector<int> bc_lo(AMREX_SPACEDIM);
    Vector<int> bc_hi(AMREX_SPACEDIM);
    BCPhysToMath(bccomp,bc_lo,bc_hi);

    bool lo_x_physbc = (bc_lo[0] == amrex::BCType::foextrap || bc_lo[0] == amrex::BCType::ext_dir) ? true : false;
    bool hi_x_physbc = (bc_hi[0] == amrex::BCType::foextrap || bc_hi[0] == amrex::BCType::ext_dir) ? true : false;
    bool lo_y_physbc = (bc_lo[1] == amrex::BCType::foextrap || bc_lo[1] == amrex::BCType::ext_dir) ? true : false;
    bool hi_y_physbc = (bc_hi[1] == amrex::BCType::foextrap || bc_hi[1] == amrex::BCType::ext_dir) ? true : false;

    // compute cell-centered ux, vy
    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
        ux(i,j,k) = (umac(i+1,j,k) - umac(i,j,k)) / hx;
        vy(i,j,k) = (vmac(i,j+1,k) - vmac(i,j,k)) / hy;
    });

    // compute sedgex on x-faces
    Box const& xbx = amrex::surroundingNodes(bx,0);
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

        // set edge values equal to the ghost cell value since they store the physical condition on the boundary
        if ( i==dlo.x && lo_x_physbc ) {
            sedgex(i,j,k,icomp) = s(i-1,j,k,icomp);
            // This will return from computing (i,j,k), not exit from ComputeConc.
            // Recall that ParallelFor is a macro that sets up a multidimensional loop
            // and then calls this lambda on the loop interior.
            return;
        }
        if ( i==dhi.x+1 && hi_x_physbc ) {
            sedgex(i,j,k,icomp) = s(i,j,k,icomp);
            return;
        }

        // local variables
        Array1D<Real, 1, 2> del;
        Array1D<Real, 1, 2> p1;
        Array1D<Real, 1, 2> p2;
        Array1D<Real, 1, 2> p3;

        Array1D<Real, 1, 3> slope_tmp;

        int ioff, joff;
        Real isign,jsign;
        Real val1,val2,val3;
        Real u;
        Real gamma;

        // Since faces of xbx can overlap between tileboxes, need a temporary
        // to hold intermediate edgestate value
        Real xedge_tmp;

        ///////////////////////////////////////////////
        // compute sedgex without transverse corrections
        ///////////////////////////////////////////////

        if (umac(i,j,k) > 0) {
            isign = 1.;
            ioff = -1;
        } else {
            isign = -1.;
            ioff = 0;
        }

        for(int n=1; n<=3; ++n){
            slope_tmp(n) = slopes(i+ioff,j,k,n-1);
        }

        // centroid of rectangular volume
        del(1) = isign*0.5*hx - 0.5*umac(i,j,k)*dt;
        del(2) = 0.;
        xedge_tmp = eval(s(i+ioff,j,k,icomp),slope_tmp,del);

        // source term
        xedge_tmp = xedge_tmp*(1. - dt2*ux(i+ioff,j,k));
        if (force) {
            xedge_tmp += dt2*force(i+ioff,j,k,icomp);
        }

        ///////////////////////////////////////////////
        // compute \Gamma^{y+}
        ///////////////////////////////////////////////

        if (vmac(i+ioff,j+1,k) > 0.) {
            jsign = 1.;
            joff = 0;
        } else {
            jsign = -1.;
            joff = 1;
        }

        u = 0.;
        if (umac(i,j,k)*umac(i,j+joff,k) > 0.) {
            u = umac(i,j+joff,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;

        p3(1) = isign*0.5*hx - u*dt;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j+1,k)*dt;

        for(int n=1; n<=3; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k,n-1);
        }

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p2(ll)+p3(ll))/2.;
        }
        val1 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p3(ll))/2.;
        }
        val2 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p2(ll))/2.;
        }
        val3 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.;

        // source term
        gamma = gamma*(1. - dt3*(ux(i+ioff,j+joff,k)+vy(i+ioff,j+joff,k)));

        ///////////////////////////////////////////////
        // correct sedgex with \Gamma^{y+}
        ///////////////////////////////////////////////

        gamma = gamma * vmac(i+ioff,j+1,k);

        xedge_tmp = xedge_tmp - dt*gamma/(2.*hy);

        ///////////////////////////////////////////////
        // compute \Gamma^{y-}
        ///////////////////////////////////////////////

        if (vmac(i+ioff,j,k) > 0.) {
            jsign = 1.;
            joff = -1;
        } else {
            jsign = -1.;
            joff = 0;
                }

        u = 0.;
        if (umac(i,j,k)*umac(i,j+joff,k) > 0.) {
            u = umac(i,j+joff,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;

        p3(1) = isign*0.5*hx - u*dt;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j,k)*dt;

        for(int n=1; n<=3; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k,n-1);
        }

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p2(ll)+p3(ll))/2.;
        }
        val1 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p3(ll))/2.;
        }
        val2 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p2(ll))/2.;
        }
        val3 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.;

        // source term
        gamma = gamma*(1. - dt3*(ux(i+ioff,j+joff,k)+vy(i+ioff,j+joff,k)));

        ///////////////////////////////////////////////
        // correct sedgex with \Gamma^{y-}
        ///////////////////////////////////////////////

        gamma = gamma * vmac(i+ioff,j,k);
        sedgex(i,j,k,icomp) = xedge_tmp + dt*gamma/(2.*hy);
    });

    // compute sedgey on y-faces
    Box const& ybx = amrex::surroundingNodes(bx,1);
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

        // set edge values equal to the ghost cell value since they store the physical condition on the boundary
        if ( j==dlo.y && lo_y_physbc ) {
            sedgey(i,j,k,icomp) = s(i,j-1,k,icomp);
            return;
        }
        if ( j==dhi.y+1 && hi_y_physbc ) {
            sedgey(i,j,k,icomp) = s(i,j,k,icomp);
            return;
        }

        // local variables
        Array1D<Real, 1, 2> del;
        Array1D<Real, 1, 2> p1;
        Array1D<Real, 1, 2> p2;
        Array1D<Real, 1, 2> p3;

        Array1D<Real, 1, 3> slope_tmp;

        int ioff, joff;
        Real isign,jsign;
        Real val1,val2,val3;
        Real v;
        Real gamma;

        // To hold intermediate edgestate value
        Real yedge_tmp;

        ///////////////////////////////////////////////
        // compute sedgey without transverse corrections
        ///////////////////////////////////////////////

        // centroid of rectangular volume
        if (vmac(i,j,k) > 0.) {
            jsign = 1.;
            joff = -1;
        } else {
            jsign = -1.;
            joff = 0;
        }

        for(int n=1; n<=3; ++n){
            slope_tmp(n) = slopes(i,j+joff,k,n-1);
        }

        del(1) = 0.;
        del(2) = jsign*0.5*hy - 0.5*vmac(i,j,k)*dt;
        yedge_tmp = eval(s(i,j+joff,k,icomp),slope_tmp,del);

        // source term
        yedge_tmp = yedge_tmp*(1. - dt2*vy(i,j+joff,k));
        if (force) {
            yedge_tmp += dt2*force(i,j+joff,k,icomp);
        }

        ///////////////////////////////////////////////
        // compute \Gamma^{x+}
        ///////////////////////////////////////////////

        if (umac(i+1,j+joff,k) > 0.) {
            isign = 1.;
            ioff = 0;
        } else {
            isign = -1.;
            ioff = 1;
        }

        v = 0.;
        if (vmac(i,j,k)*vmac(i+ioff,j,k) > 0.) {
            v = vmac(i+ioff,j,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx - umac(i+1,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy - v*dt;

        for(int n=1; n<=3; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k,n-1);
        }

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p2(ll)+p3(ll))/2.;
        }
        val1 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p3(ll))/2.;
        }
        val2 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p2(ll))/2.;
        }
        val3 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.;

        // source term
        gamma = gamma*(1. - dt3*(ux(i+ioff,j+joff,k)+vy(i+ioff,j+joff,k)));

        ///////////////////////////////////////////////
        // correct sedgey with \Gamma^{x+}
        ///////////////////////////////////////////////

        gamma = gamma * umac(i+1,j+joff,k);
        yedge_tmp = yedge_tmp - dt*gamma/(2.*hx);

        ///////////////////////////////////////////////
        // compute \Gamma^{x-}
        ///////////////////////////////////////////////

        if (umac(i,j+joff,k) > 0.) {
            isign = 1.;
            ioff = -1;
        } else {
            isign = -1.;
            ioff = 0;
        }

        v = 0.;
        if (vmac(i,j,k)*vmac(i+ioff,j,k) > 0.) {
            v = vmac(i+ioff,j,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx - umac(i,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy - v*dt;

        for(int n=1; n<=3; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k,n-1);
        }

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p2(ll)+p3(ll))/2.;
        }
        val1 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p3(ll))/2.;
        }
        val2 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for (int ll=1; ll<=2; ++ll) {
            del(ll) = (p1(ll)+p2(ll))/2.;
        }
        val3 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.;

        // source term
        gamma = gamma*(1. - dt3*(ux(i+ioff,j+joff,k)+vy(i+ioff,j+joff,k)));

        ///////////////////////////////////////////////
        // correct sedgey with \Gamma^{x-}
        ///////////////////////////////////////////////

        gamma = gamma * umac(i,j+joff,k);
        sedgey(i,j,k,icomp) = yedge_tmp + dt*gamma/(2.*hx);
    });
}

#elif (AMREX_SPACEDIM==3)

/**
 * Uses the Bell-Dawson-Shubin (BDS) algorithm, a higher order Godunov
 * method for scalar conservation laws in three dimensions, to compute
 * edge states.
 *
 * \param [in]     bx          Current grid patch
 * \param [in]     ncomp       Number of components to work on
 * \param [in]     q           Array4 of state, starting at component of interest
 * \param [in,out] xedge       Array4 containing x-edges, starting at component of interest
 * \param [in,out] yedge       Array4 containing y-edges, starting at component of interest
 * \param [in,out] zedge       Array4 containing z-edges, starting at component of interest
 * \param [in]     umac        x-Face velocities.
 * \param [in]     vmac        y-Face velocities.
 * \param [in]     wmac        z-Face velocities.
 * \param [in]     fq          Array4 for forces, starting at component of interest
 * \param [in]     geom        Level geometry.
 * \param [in]     l_dt        Time step.
 */

void BDS_ComputeEdgeState(Box const& bx,
                          int ncomp,
                          int bccomp,
                          Array4<Real const> const& q,
                          Array4<Real      > const& xedge,
                          Array4<Real      > const& yedge,
                          Array4<Real      > const& zedge,
                          Array4<Real const> const& umac,
                          Array4<Real const> const& vmac,
                          Array4<Real const> const& wmac,
                          Array4<Real const> const& fq,
                          Geometry geom,
                          Real l_dt,
                          const int proj_type)
{
    // For now, loop on components here
    for( int icomp = 0; icomp < ncomp; ++icomp)
    {
        // Temporary slope container per components
        Box const& bxg1 = amrex::grow(bx,1);
        FArrayBox slopefab(bxg1,7);
        Elixir slopeeli = slopefab.elixir();

        BDS_ComputeSlopes(bx, geom, icomp, bccomp,
                           q, slopefab.array());

        BDS_ComputeConc(bx, geom, icomp, bccomp,
                         q, xedge, yedge, zedge,
                         slopefab.array(),
                         umac, vmac, wmac, fq,
                         l_dt);
    }

    // project edge states to satisfy EOS
    BDS_Proj(bx,ncomp,xedge,yedge,zedge,proj_type);
}

/**
 * Compute bilinear slopes for BDS algorithm.
 *
 * \param [in]  bx      Current grid patch
 * \param [in]  geom    Level geometry.
 * \param [in]  icomp   Component of the state Array4.
 * \param [in]  s       Array4<const> of state vector.
 * \param [out] slopes  Array4 to store slope information.
 *
 */

void BDS_ComputeSlopes(Box const& bx,
                       const Geometry& geom,
                       int icomp, int bccomp,
                       Array4<Real const> const& s,
                       Array4<Real      > const& slopes)
{
    constexpr amrex::Real eps = 1.0e-8;
    bool limit_slopes = (advection_type == 2) ? true : false;

    // Define container for the nodal interpolated state
    Box const& ngbx = amrex::grow(amrex::convert(bx,IntVect(AMREX_D_DECL(1,1,1))),1);
    FArrayBox tmpnodefab(ngbx,1);
    Elixir tmpeli = tmpnodefab.elixir();
    auto const& sint = tmpnodefab.array();

    Box const& gbx = amrex::grow(bx,1);
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    Real hx = dx[0];
    Real hy = dx[1];
    Real hz = dx[2];

    Real c1 = (343.0/1728.0);
    Real c2 = (49.0 /1728.0);
    Real c3 = (7.0  /1728.0);
    Real c4 = (1.0  /1728.0);

    Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    // compute mathematical boundary conditions
    Vector<int> bc_lo(AMREX_SPACEDIM);
    Vector<int> bc_hi(AMREX_SPACEDIM);
    BCPhysToMath(bccomp,bc_lo,bc_hi);

    bool lo_x_physbc = (bc_lo[0] == amrex::BCType::foextrap || bc_lo[0] == amrex::BCType::ext_dir) ? true : false;
    bool hi_x_physbc = (bc_hi[0] == amrex::BCType::foextrap || bc_hi[0] == amrex::BCType::ext_dir) ? true : false;
    bool lo_y_physbc = (bc_lo[1] == amrex::BCType::foextrap || bc_lo[1] == amrex::BCType::ext_dir) ? true : false;
    bool hi_y_physbc = (bc_hi[1] == amrex::BCType::foextrap || bc_hi[1] == amrex::BCType::ext_dir) ? true : false;
    bool lo_z_physbc = (bc_lo[2] == amrex::BCType::foextrap || bc_lo[2] == amrex::BCType::ext_dir) ? true : false;
    bool hi_z_physbc = (bc_hi[2] == amrex::BCType::foextrap || bc_hi[2] == amrex::BCType::ext_dir) ? true : false;

    // tricubic interpolation to corner points
    // (i,j,k) refers to lower corner of cell
    ParallelFor(ngbx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

        // set node values equal to the average of the ghost cell values since they store the physical condition on the boundary
        if ( i<=dlo.x && lo_x_physbc ) {
            sint(i,j,k) = 0.25*(s(dlo.x-1,j,k,icomp) + s(dlo.x-1,j-1,k,icomp) + s(dlo.x-1,j,k-1,icomp) + s(dlo.x-1,j-1,k-1,icomp));
            return;
        }
        if ( i>=dhi.x+1 && hi_x_physbc ) {
            sint(i,j,k) = 0.25*(s(dhi.x+1,j,k,icomp) + s(dhi.x+1,j-1,k,icomp) + s(dhi.x+1,j,k-1,icomp) + s(dhi.x+1,j-1,k-1,icomp));
            return;
        }
        if ( j<=dlo.y && lo_y_physbc ) {
            sint(i,j,k) = 0.25*(s(i,dlo.y-1,k,icomp) + s(i-1,dlo.y-1,k,icomp) + s(i,dlo.y-1,k-1,icomp) + s(i-1,dlo.y-1,k-1,icomp));
            return;
        }
        if ( j>=dhi.y+1 && hi_y_physbc ) {
            sint(i,j,k) = 0.25*(s(i,dhi.y+1,k,icomp) + s(i-1,dhi.y+1,k,icomp) + s(i,dhi.y+1,k-1,icomp) + s(i-1,dhi.y+1,k-1,icomp));
            return;
        }
        if ( k<=dlo.z && lo_z_physbc ) {
            sint(i,j,k) = 0.25*(s(i,j,dlo.z-1,icomp) + s(i-1,j,dlo.z-1,icomp) + s(i,j-1,dlo.z-1,icomp) + s(i-1,j-1,dlo.z-1,icomp));
            return;
        }
        if ( k>=dhi.z+1 && hi_z_physbc ) {
            sint(i,j,k) = 0.25*(s(i,j,dhi.z+1,icomp) + s(i-1,j,dhi.z+1,icomp) + s(i,j-1,dhi.z+1,icomp) + s(i-1,j-1,dhi.z+1,icomp));
            return;
        }

        // one cell inward from any physical boundary, revert to 8-point average
        if ( (i==dlo.x+1 && lo_x_physbc) ||
             (i==dhi.x   && hi_x_physbc) ||
             (j==dlo.y+1 && lo_y_physbc) ||
             (j==dhi.y   && hi_y_physbc) ||
             (k==dlo.z+1 && lo_z_physbc) ||
             (k==dhi.z   && hi_z_physbc) ) {

            sint(i,j,k) = 0.125* (s(i,j,k  ,icomp) + s(i-1,j,k  ,icomp) + s(i,j-1,k  ,icomp) + s(i-1,j-1,k  ,icomp) +
                                  s(i,j,k-1,icomp) + s(i-1,j,k-1,icomp) + s(i,j-1,k-1,icomp) + s(i-1,j-1,k-1,icomp));
            return;
        }

        sint(i,j,k) = c1*( s(i  ,j  ,k  ,icomp) + s(i-1,j  ,k  ,icomp) + s(i  ,j-1,k  ,icomp)
                          +s(i  ,j  ,k-1,icomp) + s(i-1,j-1,k  ,icomp) + s(i-1,j  ,k-1,icomp)
                          +s(i  ,j-1,k-1,icomp) + s(i-1,j-1,k-1,icomp) )
                     -c2*( s(i-1,j  ,k+1,icomp) + s(i  ,j  ,k+1,icomp) + s(i-1,j-1,k+1,icomp)
                          +s(i  ,j-1,k+1,icomp) + s(i-1,j+1,k  ,icomp) + s(i  ,j+1,k  ,icomp)
                          +s(i-2,j  ,k  ,icomp) + s(i+1,j  ,k  ,icomp) + s(i-2,j-1,k  ,icomp)
                          +s(i+1,j-1,k  ,icomp) + s(i-1,j-2,k  ,icomp) + s(i  ,j-2,k  ,icomp)
                          +s(i-1,j+1,k-1,icomp) + s(i  ,j+1,k-1,icomp) + s(i-2,j  ,k-1,icomp)
                          +s(i+1,j  ,k-1,icomp) + s(i-2,j-1,k-1,icomp) + s(i+1,j-1,k-1,icomp)
                          +s(i-1,j-2,k-1,icomp) + s(i  ,j-2,k-1,icomp) + s(i-1,j  ,k-2,icomp)
                          +s(i  ,j  ,k-2,icomp) + s(i-1,j-1,k-2,icomp) + s(i  ,j-1,k-2,icomp) )
                     +c3*( s(i-1,j+1,k+1,icomp) + s(i  ,j+1,k+1,icomp) + s(i-2,j  ,k+1,icomp)
                          +s(i+1,j  ,k+1,icomp) + s(i-2,j-1,k+1,icomp) + s(i+1,j-1,k+1,icomp)
                          +s(i-1,j-2,k+1,icomp) + s(i  ,j-2,k+1,icomp) + s(i-2,j+1,k  ,icomp)
                          +s(i+1,j+1,k  ,icomp) + s(i-2,j-2,k  ,icomp) + s(i+1,j-2,k  ,icomp)
                          +s(i-2,j+1,k-1,icomp) + s(i+1,j+1,k-1,icomp) + s(i-2,j-2,k-1,icomp)
                          +s(i+1,j-2,k-1,icomp) + s(i-1,j+1,k-2,icomp) + s(i  ,j+1,k-2,icomp)
                          +s(i-2,j  ,k-2,icomp) + s(i+1,j  ,k-2,icomp) + s(i-2,j-1,k-2,icomp)
                          +s(i+1,j-1,k-2,icomp) + s(i-1,j-2,k-2,icomp) + s(i  ,j-2,k-2,icomp) )
                     -c4*( s(i-2,j+1,k+1,icomp) + s(i+1,j+1,k+1,icomp) + s(i-2,j-2,k+1,icomp)
                          +s(i+1,j-2,k+1,icomp) + s(i-2,j+1,k-2,icomp) + s(i+1,j+1,k-2,icomp)
                          +s(i-2,j-2,k-2,icomp) + s(i+1,j-2,k-2,icomp) );
    });

    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
        // compute initial estimates of slopes from unlimited corner points

        // local variables
        Real sumloc, redfac, redmax, div, kdp, sumdif, sgndif;

        // Variables local to this loop
        Array1D<Real, 1, 8> diff;
        Array1D<Real, 1, 8> smin;
        Array1D<Real, 1, 8> smax;
        Array1D<Real, 1, 8> sc;

        Array1D<bool, 1, 8> allow_change;
        for (int mm=1; mm<=8; ++mm) {
            allow_change(mm) = true;
        }

        if ( i==dlo.x && lo_x_physbc ) {
            allow_change(1) = false;
            allow_change(2) = false;
            allow_change(3) = false;
            allow_change(4) = false;
        }
        if ( i==dhi.x+1 && hi_x_physbc ) {
            allow_change(5) = false;
            allow_change(6) = false;
            allow_change(7) = false;
            allow_change(8) = false;
        }
        if ( j==dlo.y && lo_y_physbc ) {
            allow_change(1) = false;
            allow_change(2) = false;
            allow_change(5) = false;
            allow_change(6) = false;
        }
        if ( j==dhi.y+1 && hi_y_physbc ) {
            allow_change(3) = false;
            allow_change(4) = false;
            allow_change(7) = false;
            allow_change(8) = false;
        }
        if ( k==dlo.z && lo_z_physbc ) {
            allow_change(1) = false;
            allow_change(3) = false;
            allow_change(5) = false;
            allow_change(7) = false;
        }
        if ( k==dhi.z+1 && hi_z_physbc ) {
            allow_change(2) = false;
            allow_change(4) = false;
            allow_change(6) = false;
            allow_change(8) = false;
        }

         // compute initial estimates of slopes from unlimited corner points
         // sx
         slopes(i,j,k,0) = 0.25*(( sint(i+1,j  ,k  ) + sint(i+1,j+1,k  )
                                  +sint(i+1,j  ,k+1) + sint(i+1,j+1,k+1) )
                                -( sint(i  ,j  ,k  ) + sint(i  ,j+1,k  )
                                  +sint(i  ,j  ,k+1) + sint(i  ,j+1,k+1) )) / hx;
         // sy
         slopes(i,j,k,1) = 0.25*(( sint(i  ,j+1,k  ) + sint(i+1,j+1,k  )
                                  +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1) )
                                -( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  )
                                  +sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1) )) / hy;

         // sz
         slopes(i,j,k,2) = 0.25*(( sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1)
                                  +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1) )
                                -( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  )
                                  +sint(i  ,j+1,k  ) + sint(i+1,j+1,k  ) )) / hz;

         // sxy
         slopes(i,j,k,3) = 0.5*( ( sint(i  ,j  ,k  ) + sint(i  ,j  ,k+1)
                                  +sint(i+1,j+1,k  ) + sint(i+1,j+1,k+1) )
                                -( sint(i+1,j  ,k  ) + sint(i+1,j  ,k+1)
                                  +sint(i  ,j+1,k  ) + sint(i  ,j+1,k+1) )) / (hx*hy);

         // sxz
         slopes(i,j,k,4) = 0.5*( ( sint(i  ,j  ,k  ) + sint(i  ,j+1,k  )
                                  +sint(i+1,j  ,k+1) + sint(i+1,j+1,k+1) )
                                -( sint(i+1,j  ,k  ) + sint(i+1,j+1,k  )
                                  +sint(i  ,j  ,k+1) + sint(i  ,j+1,k+1) )) / (hx*hz);

         // syz
         slopes(i,j,k,5) = 0.5*( ( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  )
                                  +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1) )
                                -( sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1)
                                  +sint(i  ,j+1,k  ) + sint(i+1,j+1,k  ) )) / (hy*hz);

         // sxyz
         slopes(i,j,k,6) =       (-sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) + sint(i  ,j+1,k  )
                                  +sint(i  ,j  ,k+1) - sint(i+1,j+1,k  ) - sint(i+1,j  ,k+1)
                                  -sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1) ) / (hx*hy*hz);

         if (limit_slopes) {

             // +++ / sint(i+1,j+1,k+1)
             sc(8) = s(i,j,k,icomp)
                  +0.5  *(     hx*slopes(i,j,k,0)+   hy*slopes(i,j,k,1)+   hz*slopes(i,j,k,2))
                  +0.25 *(  hx*hy*slopes(i,j,k,3)+hx*hz*slopes(i,j,k,4)+hy*hz*slopes(i,j,k,5))
                  +0.125*hx*hy*hz*slopes(i,j,k,6);

             // ++- / sint(i+1,j+1,k  )
             sc(7) = s(i,j,k,icomp)
                  +0.5  *(     hx*slopes(i,j,k,0)+   hy*slopes(i,j,k,1)-   hz*slopes(i,j,k,2))
                  +0.25 *(  hx*hy*slopes(i,j,k,3)-hx*hz*slopes(i,j,k,4)-hy*hz*slopes(i,j,k,5))
                  -0.125*hx*hy*hz*slopes(i,j,k,6);

             // +-+ / sint(i+1,j  ,k+1)
             sc(6) = s(i,j,k,icomp)
                  +0.5  *(     hx*slopes(i,j,k,0)-   hy*slopes(i,j,k,1)+   hz*slopes(i,j,k,2))
                  +0.25 *( -hx*hy*slopes(i,j,k,3)+hx*hz*slopes(i,j,k,4)-hy*hz*slopes(i,j,k,5))
                  -0.125*hx*hy*hz*slopes(i,j,k,6);

             // +-- / sint(i+1,j  ,k  )
             sc(5) = s(i,j,k,icomp)
                  +0.5  *(     hx*slopes(i,j,k,0)-   hy*slopes(i,j,k,1)-   hz*slopes(i,j,k,2))
                  +0.25 *( -hx*hy*slopes(i,j,k,3)-hx*hz*slopes(i,j,k,4)+hy*hz*slopes(i,j,k,5))
                  +0.125*hx*hy*hz*slopes(i,j,k,6);

             // -++ / sint(i  ,j+1,k+1)
             sc(4) = s(i,j,k,icomp)
                  +0.5  *(    -hx*slopes(i,j,k,0)+   hy*slopes(i,j,k,1)+   hz*slopes(i,j,k,2))
                  +0.25 *( -hx*hy*slopes(i,j,k,3)-hx*hz*slopes(i,j,k,4)+hy*hz*slopes(i,j,k,5))
                  -0.125*hx*hy*hz*slopes(i,j,k,6);

             // -+- / sint(i  ,j+1,k  )
             sc(3) = s(i,j,k,icomp)
                  +0.5  *(    -hx*slopes(i,j,k,0)+   hy*slopes(i,j,k,1)-   hz*slopes(i,j,k,2))
                  +0.25 *( -hx*hy*slopes(i,j,k,3)+hx*hz*slopes(i,j,k,4)-hy*hz*slopes(i,j,k,5))
                  +0.125*hx*hy*hz*slopes(i,j,k,6);

             // --+ / sint(i  ,j  ,k+1)
             sc(2) = s(i,j,k,icomp)
                  +0.5  *(    -hx*slopes(i,j,k,0)-   hy*slopes(i,j,k,1)+   hz*slopes(i,j,k,2))
                  +0.25 *(  hx*hy*slopes(i,j,k,3)-hx*hz*slopes(i,j,k,4)-hy*hz*slopes(i,j,k,5))
                  +0.125*hx*hy*hz*slopes(i,j,k,6);

             // ---/ sint(i  ,j  ,k  )
             sc(1) = s(i,j,k,icomp)
                  +0.5  *(    -hx*slopes(i,j,k,0)-   hy*slopes(i,j,k,1)-   hz*slopes(i,j,k,2))
                  +0.25 *(  hx*hy*slopes(i,j,k,3)+hx*hz*slopes(i,j,k,4)+hy*hz*slopes(i,j,k,5))
                  -0.125*hx*hy*hz*slopes(i,j,k,6);

             // enforce max/min bounds
             smin(8) = min(s(i  ,j  ,k  ,icomp),s(i+1,j  ,k  ,icomp),s(i  ,j+1,k  ,icomp),s(i  ,j  ,k+1,icomp),
                           s(i+1,j+1,k  ,icomp),s(i+1,j  ,k+1,icomp),s(i  ,j+1,k+1,icomp),s(i+1,j+1,k+1,icomp));
             smax(8) = max(s(i  ,j  ,k  ,icomp),s(i+1,j  ,k  ,icomp),s(i  ,j+1,k  ,icomp),s(i  ,j  ,k+1,icomp),
                           s(i+1,j+1,k  ,icomp),s(i+1,j  ,k+1,icomp),s(i  ,j+1,k+1,icomp),s(i+1,j+1,k+1,icomp));

             smin(7) = min(s(i  ,j  ,k-1,icomp),s(i+1,j  ,k-1,icomp),s(i  ,j+1,k-1,icomp),s(i  ,j  ,k  ,icomp),
                           s(i+1,j+1,k-1,icomp),s(i+1,j  ,k  ,icomp),s(i  ,j+1,k  ,icomp),s(i+1,j+1,k  ,icomp));
             smax(7) = max(s(i  ,j  ,k-1,icomp),s(i+1,j  ,k-1,icomp),s(i  ,j+1,k-1,icomp),s(i  ,j  ,k  ,icomp),
                           s(i+1,j+1,k-1,icomp),s(i+1,j  ,k  ,icomp),s(i  ,j+1,k  ,icomp),s(i+1,j+1,k  ,icomp));

             smin(6) = min(s(i  ,j-1,k  ,icomp),s(i+1,j-1,k  ,icomp),s(i  ,j  ,k  ,icomp),s(i  ,j-1,k+1,icomp),
                           s(i+1,j  ,k  ,icomp),s(i+1,j-1,k+1,icomp),s(i  ,j  ,k+1,icomp),s(i+1,j  ,k+1,icomp));
             smax(6) = max(s(i  ,j-1,k  ,icomp),s(i+1,j-1,k  ,icomp),s(i  ,j  ,k  ,icomp),s(i  ,j-1,k+1,icomp),
                           s(i+1,j  ,k  ,icomp),s(i+1,j-1,k+1,icomp),s(i  ,j  ,k+1,icomp),s(i+1,j  ,k+1,icomp));

             smin(5) = min(s(i  ,j-1,k-1,icomp),s(i+1,j-1,k-1,icomp),s(i  ,j  ,k-1,icomp),s(i  ,j-1,k  ,icomp),
                           s(i+1,j  ,k-1,icomp),s(i+1,j-1,k  ,icomp),s(i  ,j  ,k  ,icomp),s(i+1,j  ,k  ,icomp));
             smax(5) = max(s(i  ,j-1,k-1,icomp),s(i+1,j-1,k-1,icomp),s(i  ,j  ,k-1,icomp),s(i  ,j-1,k  ,icomp),
                           s(i+1,j  ,k-1,icomp),s(i+1,j-1,k  ,icomp),s(i  ,j  ,k  ,icomp),s(i+1,j  ,k  ,icomp));

             smin(4) = min(s(i-1,j  ,k  ,icomp),s(i  ,j  ,k  ,icomp),s(i-1,j+1,k  ,icomp),s(i-1,j  ,k+1,icomp),
                           s(i  ,j+1,k  ,icomp),s(i  ,j  ,k+1,icomp),s(i-1,j+1,k+1,icomp),s(i  ,j+1,k+1,icomp));
             smax(4) = max(s(i-1,j  ,k  ,icomp),s(i  ,j  ,k  ,icomp),s(i-1,j+1,k  ,icomp),s(i-1,j  ,k+1,icomp),
                           s(i  ,j+1,k  ,icomp),s(i  ,j  ,k+1,icomp),s(i-1,j+1,k+1,icomp),s(i  ,j+1,k+1,icomp));

             smin(3) = min(s(i-1,j  ,k-1,icomp),s(i  ,j  ,k-1,icomp),s(i-1,j+1,k-1,icomp),s(i-1,j  ,k  ,icomp),
                           s(i  ,j+1,k-1,icomp),s(i  ,j  ,k  ,icomp),s(i-1,j+1,k  ,icomp),s(i  ,j+1,k  ,icomp));
             smax(3) = max(s(i-1,j  ,k-1,icomp),s(i  ,j  ,k-1,icomp),s(i-1,j+1,k-1,icomp),s(i-1,j  ,k  ,icomp),
                           s(i  ,j+1,k-1,icomp),s(i  ,j  ,k  ,icomp),s(i-1,j+1,k  ,icomp),s(i  ,j+1,k  ,icomp));

             smin(2) = min(s(i-1,j-1,k  ,icomp),s(i  ,j-1,k  ,icomp),s(i-1,j  ,k  ,icomp),s(i-1,j-1,k+1,icomp),
                           s(i  ,j  ,k  ,icomp),s(i  ,j-1,k+1,icomp),s(i-1,j  ,k+1,icomp),s(i  ,j  ,k+1,icomp));
             smax(2) = max(s(i-1,j-1,k  ,icomp),s(i  ,j-1,k  ,icomp),s(i-1,j  ,k  ,icomp),s(i-1,j-1,k+1,icomp),
                           s(i  ,j  ,k  ,icomp),s(i  ,j-1,k+1,icomp),s(i-1,j  ,k+1,icomp),s(i  ,j  ,k+1,icomp));

             smin(1) = min(s(i-1,j-1,k-1,icomp),s(i  ,j-1,k-1,icomp),s(i-1,j  ,k-1,icomp),s(i-1,j-1,k  ,icomp),
                           s(i  ,j  ,k-1,icomp),s(i  ,j-1,k  ,icomp),s(i-1,j  ,k  ,icomp),s(i  ,j  ,k  ,icomp));
             smax(1) = max(s(i-1,j-1,k-1,icomp),s(i  ,j-1,k-1,icomp),s(i-1,j  ,k-1,icomp),s(i-1,j-1,k  ,icomp),
                           s(i  ,j  ,k-1,icomp),s(i  ,j-1,k  ,icomp),s(i-1,j  ,k  ,icomp),s(i  ,j  ,k  ,icomp));

             for(int mm=1; mm<=8; ++mm){
                if (allow_change(mm)) {
                    sc(mm) = max(min(sc(mm), smax(mm)), smin(mm));
                }
             }

             // iterative loop
             for(int ll = 1; ll<=6; ++ll){

                // compute the amount by which the average of the nodal values differs from cell-center value
                sumloc = 0.125*(sc(1)+sc(2)+sc(3)+sc(4)+sc(5)+sc(6)+sc(7)+sc(8));
                sumdif = (sumloc - s(i,j,k,icomp))*8.0;

                // sgndif = +(-)1 if the node average is too large(small)
                sgndif = std::copysign(1.0_rt,sumdif);

                // compute how much each node is larger(smaller) than the cell-centered value
                for(int mm=1; mm<=8; ++mm){
                   diff(mm) = (sc(mm) - s(i,j,k,icomp))*sgndif;
                }

                kdp = 0;

               // count how many nodes are larger(smaller) than the cell-centered value
                for(int mm=1; mm<=8; ++mm){
                   if (diff(mm) > eps && allow_change(mm)) {
                      kdp = kdp+1;
                   }
                }

                // adjust node values
                for(int mm=1; mm<=8; ++mm){

                   // don't allow boundary nodes to change value
                   if (!allow_change(mm)) continue;

                   if (kdp<1) {
                      div = 1.0;
                   } else {
                      div = kdp;
                   }

                   if (diff(mm)>eps) {
                      redfac = sumdif*sgndif/div;
                      kdp = kdp-1;
                   } else {
                      redfac = 0.0;
                   }

                   if (sgndif > 0.0) {
                      redmax = sc(mm) - smin(mm);
                   } else {
                      redmax = smax(mm) - sc(mm);
                   }

                   redfac = min(redfac,redmax);
                   sumdif = sumdif - redfac*sgndif;
                   sc(mm) = sc(mm) - redfac*sgndif;
                }
             }

             // final slopes

             // sx
             slopes(i,j,k,0) = 0.25*( ( sc(5) + sc(7)
                                       +sc(6) + sc(8))
                                     -( sc(1) + sc(3)
                                       +sc(2) + sc(4)) ) / hx;

             // sy
             slopes(i,j,k,1) = 0.25*( ( sc(3) + sc(7)
                                       +sc(4) + sc(8))
                                     -( sc(1) + sc(5)
                                       +sc(2) + sc(6)) ) / hy;

             // sz
             slopes(i,j,k,2) = 0.25*( ( sc(2) + sc(6)
                                       +sc(4) + sc(8))
                                     -( sc(1) + sc(5)
                                       +sc(3) + sc(7)) ) / hz;

             // sxy
             slopes(i,j,k,3) = 0.5*( ( sc(1) + sc(2)
                                      +sc(7) + sc(8))
                                    -( sc(5) + sc(6)
                                      +sc(3) + sc(4)) ) / (hx*hy);

             // sxz
             slopes(i,j,k,4) = 0.5*( ( sc(1) + sc(3)
                                      +sc(6) + sc(8))
                                    -( sc(5) + sc(7)
                                      +sc(2) + sc(4)) ) / (hx*hz);

             // syz
             slopes(i,j,k,5) = 0.5*( ( sc(1) + sc(5)
                                      +sc(4) + sc(8))
                                    -( sc(2) + sc(6)
                                      +sc(3) + sc(7)) ) / (hy*hz);

             // sxyz
             slopes(i,j,k,6) = (-sc(1) + sc(5) + sc(3)
                                +sc(2) - sc(7) - sc(6)
                                -sc(4) + sc(8) ) / (hx*hy*hz);

         }
    });
}

/**
 * Returns updated edge value.
 *
 * \param [in] s
 * \param [in] slope
 * \param [in] del
 *
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real eval (const Real s,
           Array1D<Real,1,7>& slope,
           Array1D<Real,1,3>& del)
{
    Real val = s + del(1)*slope(1)        + del(2)*slope(2)        + del(3)*slope(3)
                 + del(1)*del(2)*slope(4) + del(1)*del(3)*slope(5) + del(2)*del(3)*slope(6)
                 + del(1)*del(2)*del(3)*slope(7);

    return val;
}

/**
 * Compute Conc for BDS algorithm.
 *
 * \param [in]     bx          Current grid patch
 * \param [in]     geom        Level geometry.
 * \param [in]     icomp       Component of the Array4s.
 * \param [in]     s           Array4 of state.
 * \param [in,out] sedgex      Array4 containing x-edges.
 * \param [in,out] sedgey      Array4 containing y-edges.
 * \param [in,out] sedgez      Array4 containing z-edges.
 * \param [in]     slopes      Array4 containing slope information.
 * \param [in]     umac        Array4 for u-face velocity.
 * \param [in]     vmac        Array4 for v-face velocity.
 * \param [in]     wmac        Array4 for z-face velocity.
 * \param [in]     force       Array4 for forces.
 * \param [in]     dt          Time step.
 *
 */

void BDS_ComputeConc(Box const& bx,
                     const Geometry& geom,
                     int icomp, int bccomp,
                     Array4<Real const> const& s,
                     Array4<Real      > const& sedgex,
                     Array4<Real      > const& sedgey,
                     Array4<Real      > const& sedgez,
                     Array4<Real const> const& slopes,
                     Array4<Real const> const& umac,
                     Array4<Real const> const& vmac,
                     Array4<Real const> const& wmac,
                     Array4<Real const> const& force,
                     const Real dt)
{
    Box const& gbx = amrex::grow(bx,1);
    GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

    FArrayBox ux_fab(gbx,1);
    FArrayBox vy_fab(gbx,1);
    FArrayBox wz_fab(gbx,1);
    Elixir uxeli = ux_fab.elixir();
    Elixir vyeli = vy_fab.elixir();
    Elixir wzeli = wz_fab.elixir();
    auto const& ux    = ux_fab.array();
    auto const& vy    = vy_fab.array();
    auto const& wz    = wz_fab.array();

    Real hx = dx[0];
    Real hy = dx[1];
    Real hz = dx[2];

    Real dt2 = dt/2.0;
    Real dt3 = dt/3.0;
    Real dt4 = dt/4.0;

    Real sixth = 1.0/6.0;

    Box const& domain = geom.Domain();
    const auto dlo = amrex::lbound(domain);
    const auto dhi = amrex::ubound(domain);

    // compute mathematical boundary conditions
    Vector<int> bc_lo(AMREX_SPACEDIM);
    Vector<int> bc_hi(AMREX_SPACEDIM);
    BCPhysToMath(bccomp,bc_lo,bc_hi);

    bool lo_x_physbc = (bc_lo[0] == amrex::BCType::foextrap || bc_lo[0] == amrex::BCType::ext_dir) ? true : false;
    bool hi_x_physbc = (bc_hi[0] == amrex::BCType::foextrap || bc_hi[0] == amrex::BCType::ext_dir) ? true : false;
    bool lo_y_physbc = (bc_lo[1] == amrex::BCType::foextrap || bc_lo[1] == amrex::BCType::ext_dir) ? true : false;
    bool hi_y_physbc = (bc_hi[1] == amrex::BCType::foextrap || bc_hi[1] == amrex::BCType::ext_dir) ? true : false;
    bool lo_z_physbc = (bc_lo[2] == amrex::BCType::foextrap || bc_lo[2] == amrex::BCType::ext_dir) ? true : false;
    bool hi_z_physbc = (bc_hi[2] == amrex::BCType::foextrap || bc_hi[2] == amrex::BCType::ext_dir) ? true : false;

    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
          ux(i,j,k) = (umac(i+1,j,k) - umac(i,j,k)) / hx;
          vy(i,j,k) = (vmac(i,j+1,k) - vmac(i,j,k)) / hy;
          wz(i,j,k) = (wmac(i,j,k+1) - wmac(i,j,k)) / hz;
    });

    // compute sedgex on x-faces
    Box const& xbx = amrex::surroundingNodes(bx,0);
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

        // set edge values equal to the ghost cell value since they store the physical condition on the boundary
        if ( i==dlo.x && lo_x_physbc ) {
            sedgex(i,j,k,icomp) = s(i-1,j,k,icomp);
            // This will return from computing (i,j,k), not exit from ComputeConc.
            // Recall that ParallelFor is a macro that sets up a multidimensional loop
            // and then calls this lambda on the loop interior.
            return;
        }
        if ( i==dhi.x+1 && hi_x_physbc ) {
            sedgex(i,j,k,icomp) = s(i,j,k,icomp);
            return;
        }

        //local variables
        Array1D<Real, 1, 3> del;
        Array1D<Real, 1, 3> p1;
        Array1D<Real, 1, 3> p2;
        Array1D<Real, 1, 3> p3;
        Array1D<Real, 1, 3> p4;

        Array1D<Real, 1, 7> slope_tmp;

        int ioff, joff, koff;

        Real isign,jsign,ksign;
        Real val1,val2,val3,val4,val5;
        Real u;
        Real uu,vv,ww;
        Real gamma,gamma2;

        // Since faces of xbx can overlap between tileboxes, need a temporary
        // to hold intermediate edgestate value
        Real xedge_tmp;

        ////////////////////////////////////////////////
        // compute sedgex without transverse corrections
        ////////////////////////////////////////////////

        if (umac(i,j,k) > 0.0) {
           isign = 1.0;
           ioff = -1;
        } else {
           isign = -1.0;
           ioff = 0;
        }

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j,k,n-1);
        }


        // centroid of rectangular volume
        del(1) = isign*0.5*hx - 0.5*umac(i,j,k)*dt;
        del(2) = 0.0;
        del(3) = 0.0;
        xedge_tmp = eval(s(i+ioff,j,k,icomp),slope_tmp,del);

        // source term
        xedge_tmp = xedge_tmp*(1. - dt2*ux(i+ioff,j,k));
        if (force) {
            xedge_tmp += dt2*force(i+ioff,j,k,icomp);
        }

        ////////////////////////////////////////////////
        // compute \Gamma^{y+} without corner corrections
        ////////////////////////////////////////////////

        if (vmac(i+ioff,j+1,k) > 0.0) {
           jsign = 1.0;
           joff = 0;
        } else {
           jsign = -1.0;
           joff = 1;
        }


        u = 0.0;
        if (umac(i,j,k)*umac(i,j+joff,k) > 0.0) {
           u = umac(i,j+joff,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = 0.0;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;
        p2(3) = 0.0;

        p3(1) = isign*0.5*hx - u*dt;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j+1,k)*dt;
        p3(3) = 0.0;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p2(ll)+p3(ll))/2.0;
        }
        val1 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p3(ll))/2.0;
        }
        val2 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll))/2.0;
        }
        val3 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.0;

        // source term
        gamma = gamma*(1. - dt3*(ux(i+ioff,j+joff,k)+vy(i+ioff,j+joff,k)));

        ////////////////////////////////////////////////
        // correct \Gamma^{y+} with \Gamma^{y+,z+}
        ////////////////////////////////////////////////

        if (wmac(i+ioff,j+joff,k+1) > 0.0) {
           ksign = 1.0;
           koff = 0;
        } else {
           ksign = -1.0;
           koff = 1;
        }

        uu = 0.0;
        if (umac(i,j,k)*umac(i,j+joff,k+koff) > 0.0) {
           uu = umac(i,j+joff,k+koff);
        }

        vv = 0.0;
        if (vmac(i+ioff,j+1,k)*vmac(i+ioff,j+1,k+koff) > 0.0) {
           vv = vmac(i+ioff,j+1,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j+1,k)*dt;
        p3(3) = ksign*0.5*hz;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - wmac(i+ioff,j+joff,k+1)*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * wmac(i+ioff,j+joff,k+1);

        gamma = gamma - dt*gamma2/(3.0*hz);

        ////////////////////////////////////////////////
        // correct \Gamma^{y+} with \Gamma^{y+,z-}
        ////////////////////////////////////////////////

        if (wmac(i+ioff,j+joff,k) > 0.0) {
           ksign = 1.0;
           koff = -1;
        } else {
           ksign = -1.0;
           koff = 0;
        }

        uu = 0.0;
        if (umac(i,j,k)*umac(i,j+joff,k+koff) > 0.0) {
           uu = umac(i,j+joff,k+koff);
        }

        vv = 0.0;
        if (vmac(i+ioff,j+1,k)*vmac(i+ioff,j+1,k+koff) > 0.0) {
           vv = vmac(i+ioff,j+1,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j+1,k)*dt;
        p3(3) = ksign*0.5*hz;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - wmac(i+ioff,j+joff,k)*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * wmac(i+ioff,j+joff,k);

        gamma = gamma + dt*gamma2/(3.0*hz);

        ////////////////////////////////////////////////
        // correct sedgex with \Gamma^{y+}
        ////////////////////////////////////////////////

        gamma = gamma * vmac(i+ioff,j+1,k);
        xedge_tmp = xedge_tmp - dt*gamma/(2.0*hy);

        ////////////////////////////////////////////////
        // compute \Gamma^{y-} without corner corrections
        ////////////////////////////////////////////////

        if (vmac(i+ioff,j,k) > 0.0) {
           jsign = 1.0;
           joff = -1;
        } else {
           jsign = -1.0;
           joff = 0;
        }

        u = 0.0;
        if (umac(i,j,k)*umac(i,j+joff,k) > 0.0) {
           u = umac(i,j+joff,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = 0.0;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;
        p2(3) = 0.0;

        p3(1) = isign*0.5*hx - u*dt;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j,k)*dt;
        p3(3) = 0.0;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p2(ll)+p3(ll))/2.0;
        }
        val1 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p3(ll))/2.0;
        }
        val2 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll))/2.0;
        }
        val3 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.0;

        // source term
        gamma = gamma*(1. - dt3*(ux(i+ioff,j+joff,k)+vy(i+ioff,j+joff,k)));

        ////////////////////////////////////////////////
        // correct \Gamma^{y-} with \Gamma^{y-,z+}
        ////////////////////////////////////////////////

        if (wmac(i+ioff,j+joff,k+1) > 0.0) {
           ksign = 1.0;
           koff = 0;
        } else {
           ksign = -1.0;
           koff = 1;
        }

        uu = 0.0;
        if (umac(i,j,k)*umac(i,j+joff,k+koff) > 0.0) {
           uu = umac(i,j+joff,k+koff);
        }

        vv = 0.0;
        if (vmac(i+ioff,j,k)*vmac(i+ioff,j,k+koff) > 0.0) {
           vv = vmac(i+ioff,j,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j,k)*dt;
        p3(3) = ksign*0.5*hz;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - wmac(i+ioff,j+joff,k+1)*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * wmac(i+ioff,j+joff,k+1);

        gamma = gamma - dt*gamma2/(3.0*hz);

        ////////////////////////////////////////////////
        // correct \Gamma^{y-} with \Gamma^{y-,z-}
        ////////////////////////////////////////////////

        if (wmac(i+ioff,j+joff,k) > 0.0) {
           ksign = 1.0;
           koff = -1;
        } else {
           ksign = -1.0;
           koff = 0;
        }

        uu = 0.0;
        if (umac(i,j,k)*umac(i,j+joff,k+koff) > 0.0) {
           uu = umac(i,j+joff,k+koff);
        }

        vv = 0.0;
        if (vmac(i+ioff,j,k)*vmac(i+ioff,j,k+koff) > 0.0) {
           vv = vmac(i+ioff,j,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j,k)*dt;
        p3(3) = ksign*0.5*hz;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - wmac(i+ioff,j+joff,k)*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * wmac(i+ioff,j+joff,k);

        gamma = gamma + dt*gamma2/(3.0*hz);

        ////////////////////////////////////////////////
        // correct sedgex with \Gamma^{y-}
        ////////////////////////////////////////////////

        gamma = gamma * vmac(i+ioff,j,k);
        xedge_tmp = xedge_tmp + dt*gamma/(2.0*hy);

        ////////////////////////////////////////////////
        // compute \Gamma^{z+} without corner corrections
        ////////////////////////////////////////////////

        if (wmac(i+ioff,j,k+1) > 0.0) {
           ksign = 1.0;
           koff = 0;
        } else {
           ksign = -1.0;
           koff = 1;
        }

        u = 0.0;
        if (umac(i,j,k)*umac(i,j,k+koff) > 0.0) {
           u = umac(i,j,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = 0.0;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = 0.0;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - u*dt;
        p3(2) = 0.0;
        p3(3) = ksign*0.5*hz - wmac(i+ioff,j,k+1)*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p2(ll)+p3(ll))/2.0;
        }
        val1 = eval(s(i+ioff,j,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p3(ll))/2.0;
        }
        val2 = eval(s(i+ioff,j,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll))/2.0;
        }
        val3 = eval(s(i+ioff,j,k+koff,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.0;

        // source term
        gamma = gamma*(1. - dt3*(ux(i+ioff,j,k+koff)+wz(i+ioff,j,k+koff)));

        ////////////////////////////////////////////////
        // correct \Gamma^{z+} with \Gamma^{z+,y+}
        ////////////////////////////////////////////////

        if (vmac(i+ioff,j+1,k+koff) > 0.0) {
           jsign = 1.0;
           joff = 0;
        } else {
           jsign = -1.0;
           joff = 1;
        }

        uu = 0.0;
        if (umac(i,j,k)*umac(i,j+joff,k+koff) > 0.0) {
           uu = umac(i,j+joff,k+koff);
        }

        ww = 0.0;
        if (wmac(i+ioff,j,k+1)*wmac(i+ioff,j+joff,k+1) > 0.0) {
           ww = wmac(i+ioff,j+joff,k+1);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p3(2) = jsign*0.5*hy;
        p3(3) = ksign*0.5*hz - wmac(i+ioff,j,k+1)*dt;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vmac(i+ioff,j+1,k+koff)*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * vmac(i+ioff,j+1,k+koff);

        gamma = gamma - dt*gamma2/(3.0*hy);

        ////////////////////////////////////////////////
        // correct \Gamma^{z+} with \Gamma^{z+,y-}
        ////////////////////////////////////////////////

        if (vmac(i+ioff,j,k+koff) > 0.0) {
           jsign = 1.0;
           joff = -1;
        } else {
           jsign = -1.0;
           joff = 0;
        }

        uu = 0.0;
        if (umac(i,j,k)*umac(i,j+joff,k+koff) > 0.0) {
           uu = umac(i,j+joff,k+koff);
        }

        ww = 0.0;
        if (wmac(i+ioff,j,k+1)*wmac(i+ioff,j+joff,k+1) > 0.0) {
           ww = wmac(i+ioff,j+joff,k+1);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p3(2) = jsign*0.5*hy;
        p3(3) = ksign*0.5*hz - wmac(i+ioff,j,k+1)*dt;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vmac(i+ioff,j,k+koff)*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * vmac(i+ioff,j,k+koff);

        gamma = gamma + dt*gamma2/(3.0*hy);

        ////////////////////////////////////////////////
        // correct sedgex with \Gamma^{z+}
        ////////////////////////////////////////////////

        gamma = gamma * wmac(i+ioff,j,k+1);
        xedge_tmp = xedge_tmp - dt*gamma/(2.0*hz);

        ////////////////////////////////////////////////
        // compute \Gamma^{z-} without corner corrections
        ////////////////////////////////////////////////

        if (wmac(i+ioff,j,k) > 0.0) {
           ksign = 1.0;
           koff = -1;
        } else {
           ksign = -1.0;
           koff = 0;
        }

        u = 0.0;
        if (umac(i,j,k)*umac(i,j,k+koff) > 0.0) {
           u = umac(i,j,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = 0.0;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = 0.0;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - u*dt;
        p3(2) = 0.0;
        p3(3) = ksign*0.5*hz - wmac(i+ioff,j,k)*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p2(ll)+p3(ll))/2.0;
        }
        val1 = eval(s(i+ioff,j,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p3(ll))/2.0;
        }
        val2 = eval(s(i+ioff,j,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll))/2.0;
        }
        val3 = eval(s(i+ioff,j,k+koff,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.0;

        // source term
        gamma = gamma*(1. - dt3*(ux(i+ioff,j,k+koff)+wz(i+ioff,j,k+koff)));

        ////////////////////////////////////////////////
        // correct \Gamma^{z-} with \Gamma^{z-,y+}
        ////////////////////////////////////////////////

        if (vmac(i+ioff,j+1,k+koff) > 0.0) {
           jsign = 1.0;
           joff = 0;
        } else {
           jsign = -1.0;
           joff = 1;
        }

        uu = 0.0;
        if (umac(i,j,k)*umac(i,j+joff,k+koff) > 0.0) {
           uu = umac(i,j+joff,k+koff);
        }

        ww = 0.0;
        if (wmac(i+ioff,j,k)*wmac(i+ioff,j+joff,k) > 0.0) {
           ww = wmac(i+ioff,j+joff,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p3(2) = jsign*0.5*hy;
        p3(3) = ksign*0.5*hz - wmac(i+ioff,j,k)*dt;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vmac(i+ioff,j+1,k+koff)*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * vmac(i+ioff,j+1,k+koff);

        gamma = gamma - dt*gamma2/(3.0*hy);

        ////////////////////////////////////////////////
        // correct \Gamma^{z-} with \Gamma^{z-,y-}
        ////////////////////////////////////////////////

        if (vmac(i+ioff,j,k+koff) > 0.0) {
           jsign = 1.0;
           joff = -1;
        } else {
           jsign = -1.0;
           joff = 0;
        }

        uu = 0.0;
        if (umac(i,j,k)*umac(i,j+joff,k+koff) > 0.0) {
           uu = umac(i,j+joff,k+koff);
        }

        ww = 0.0;
        if (wmac(i+ioff,j,k)*wmac(i+ioff,j+joff,k) > 0.0) {
           ww = wmac(i+ioff,j+joff,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - umac(i,j,k)*dt;
        p3(2) = jsign*0.5*hy;
        p3(3) = ksign*0.5*hz - wmac(i+ioff,j,k)*dt;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vmac(i+ioff,j,k+koff)*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * vmac(i+ioff,j,k+koff);

        gamma = gamma + dt*gamma2/(3.0*hy);

        ////////////////////////////////////////////////
        // correct sedgex with \Gamma^{z-}
        ////////////////////////////////////////////////

        gamma = gamma * wmac(i+ioff,j,k);
        sedgex(i,j,k,icomp) = xedge_tmp + dt*gamma/(2.0*hz);
    });

    // compute sedgey on y-faces
    Box const& ybx = amrex::surroundingNodes(bx,1);
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

        // set edge values equal to the ghost cell value since they store the physical condition on the boundary
        if ( j==dlo.y && lo_y_physbc ) {
            sedgey(i,j,k,icomp) = s(i,j-1,k,icomp);
            return;
        }
        if ( j==dhi.y+1 && hi_y_physbc ) {
            sedgey(i,j,k,icomp) = s(i,j,k,icomp);
            return;
        }

        //local variables
        Array1D<Real, 1, 3> del;
        Array1D<Real, 1, 3> p1;
        Array1D<Real, 1, 3> p2;
        Array1D<Real, 1, 3> p3;
        Array1D<Real, 1, 3> p4;

        Array1D<Real, 1, 7> slope_tmp;

        int ioff, joff, koff;

        Real isign,jsign,ksign;
        Real val1,val2,val3,val4,val5;
        Real v;
        Real uu,vv,ww;
        Real gamma,gamma2;

        // Hold intermediate edgestate value
        Real yedge_tmp;

        ////////////////////////////////////////////////
        // compute sedgey without transverse corrections
        ////////////////////////////////////////////////

        // centroid of rectangular volume
        if (vmac(i,j,k) > 0.0) {
           jsign = 1.0;
           joff = -1;
        } else {
           jsign = -1.0;
           joff = 0;
        }

        del(1) = 0.0;
        del(2) = jsign*0.5*hy - 0.5*vmac(i,j,k)*dt;
        del(3) = 0.0;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i,j+joff,k,n-1);
        }

        yedge_tmp = eval(s(i,j+joff,k,icomp),slope_tmp,del);

        // source term
        yedge_tmp = yedge_tmp*(1. - dt2*vy(i,j+joff,k));
        if (force) {
            yedge_tmp += dt2*force(i,j+joff,k,icomp);
        }

        ////////////////////////////////////////////////
        // compute \Gamma^{x+} without corner corrections
        ////////////////////////////////////////////////

        if (umac(i+1,j+joff,k) > 0.0) {
           isign = 1.0;
           ioff = 0;
        } else {
           isign = -1.0;
           ioff = 1;
        }

        v = 0.0;
        if (vmac(i,j,k)*vmac(i+ioff,j,k) > 0.0) {
           v = vmac(i+ioff,j,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = 0.0;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p2(3) = 0.0;

        p3(1) = isign*0.5*hx - umac(i+1,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy - v*dt;
        p3(3) = 0.0;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p2(ll)+p3(ll))/2.0;
        }
        val1 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p3(ll))/2.0;
        }
        val2 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll))/2.0;
        }
        val3 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.0;

        // source term
        gamma = gamma*(1. - dt3*(vy(i+ioff,j+joff,k)+ux(i+ioff,j+joff,k)));

        ////////////////////////////////////////////////
        // correct \Gamma^{x+} with \Gamma^{x+,z+}
        ////////////////////////////////////////////////

        if (wmac(i+ioff,j+joff,k+1) > 0.0) {
           ksign = 1.0;
           koff = 0;
        } else {
           ksign = -1.0;
           koff = 1;
        }

        vv = 0.0;
        if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) > 0.0) {
           vv = vmac(i+ioff,j,k+koff);
        }

        uu = 0.0;
        if (umac(i+1,j+joff,k)*umac(i+1,j+joff,k+koff) > 0.0) {
           uu = umac(i+1,j+joff,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - umac(i+1,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p3(3) = ksign*0.5*hz;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - wmac(i+ioff,j+joff,k+1)*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * wmac(i+ioff,j+joff,k+1);

        gamma = gamma - dt*gamma2/(3.0*hz);

        ////////////////////////////////////////////////
        // correct \Gamma^{x+} with \Gamma^{x+,z-}
        ////////////////////////////////////////////////

        if (wmac(i+ioff,j+joff,k) > 0.0) {
           ksign = 1.0;
           koff = -1;
        } else {
           ksign = -1.0;
           koff = 0;
        }

        vv = 0.0;
        if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) > 0.0) {
           vv = vmac(i+ioff,j,k+koff);
        }

        uu = 0.0;
        if (umac(i+1,j+joff,k)*umac(i+1,j+joff,k+koff) > 0.0) {
           uu = umac(i+1,j+joff,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - umac(i+1,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p3(3) = ksign*0.5*hz;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - wmac(i+ioff,j+joff,k)*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * wmac(i+ioff,j+joff,k);

        gamma = gamma + dt*gamma2/(3.0*hz);

        ////////////////////////////////////////////////
        // correct sedgey with \Gamma^{x+}
        ////////////////////////////////////////////////

        gamma = gamma * umac(i+1,j+joff,k);
        yedge_tmp = yedge_tmp - dt*gamma/(2.0*hx);

        ////////////////////////////////////////////////
        // compute \Gamma^{x-} without corner corrections
        ////////////////////////////////////////////////

        if (umac(i,j+joff,k) > 0.0) {
           isign = 1.0;
           ioff = -1;
        } else {
           isign = -1.0;
           ioff = 0;
        }

        v = 0.0;
        if (vmac(i,j,k)*vmac(i+ioff,j,k) > 0.0) {
           v = vmac(i+ioff,j,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = 0.0;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p2(3) = 0.0;

        p3(1) = isign*0.5*hx - umac(i,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy - v*dt;
        p3(3) = 0.0;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p2(ll)+p3(ll))/2.0;
        }
        val1 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p3(ll))/2.0;
        }
        val2 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll))/2.0;
        }
        val3 = eval(s(i+ioff,j+joff,k,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.0;

        // source term
        gamma = gamma*(1. - dt3*(vy(i+ioff,j+joff,k)+ux(i+ioff,j+joff,k)));

        ////////////////////////////////////////////////
        // correct \Gamma^{x-} with \Gamma^{x-,z+}
        ////////////////////////////////////////////////

        if (wmac(i+ioff,j+joff,k+1) > 0.0) {
           ksign = 1.0;
           koff = 0;
        } else {
           ksign = -1.0;
           koff = 1;
        }

        vv = 0.0;
        if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) > 0.0) {
           vv = vmac(i+ioff,j,k+koff);
        }

        uu = 0.0;
        if (umac(i,j+joff,k)*umac(i,j+joff,k+koff) > 0.0) {
           uu = umac(i,j+joff,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - umac(i,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p3(3) = ksign*0.5*hz;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - wmac(i+ioff,j+joff,k+1)*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * wmac(i+ioff,j+joff,k+1);

        gamma = gamma - dt*gamma2/(3.0*hz);

        ////////////////////////////////////////////////
        // correct \Gamma^{x-} with \Gamma^{x-,z-}
        ////////////////////////////////////////////////

        if (wmac(i+ioff,j+joff,k) > 0.0) {
           ksign = 1.0;
           koff = -1;
        } else {
           ksign = -1.0;
           koff = 0;
        }

        vv = 0.0;
        if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) > 0.0) {
           vv = vmac(i+ioff,j,k+koff);
        }

        uu = 0.0;
        if (umac(i,j+joff,k)*umac(i,j+joff,k+koff) > 0.0) {
           uu = umac(i,j+joff,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx - umac(i,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p3(3) = ksign*0.5*hz;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - wmac(i+ioff,j+joff,k)*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * wmac(i+ioff,j+joff,k);

        gamma = gamma + dt*gamma2/(3.0*hz);

        ////////////////////////////////////////////////
        // correct sedgey with \Gamma^{x-}
        ////////////////////////////////////////////////

        gamma = gamma * umac(i,j+joff,k);
        yedge_tmp = yedge_tmp + dt*gamma/(2.0*hx);

        ////////////////////////////////////////////////
        // compute \Gamma^{z+} without corner corrections
        ////////////////////////////////////////////////

        if (wmac(i,j+joff,k+1) > 0.0) {
           ksign = 1.0;
           koff = 0;
        } else {
           ksign = -1.0;
           koff = 1;
        }

        v = 0.0;
        if (vmac(i,j,k)*vmac(i,j,k+koff) > 0.0) {
           v = vmac(i,j,k+koff);
        }

        p1(1) = 0.0;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = 0.0;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p2(3) = ksign*0.5*hz;

        p3(1) = 0.0;
        p3(2) = jsign*0.5*hy - v*dt;
        p3(3) = ksign*0.5*hz - wmac(i,j+joff,k+1)*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p2(ll)+p3(ll))/2.0;
        }
        val1 = eval(s(i,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p3(ll))/2.0;
        }
        val2 = eval(s(i,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll))/2.0;
        }
        val3 = eval(s(i,j+joff,k+koff,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.0;

        // source term
        gamma = gamma*(1. - dt3*(vy(i,j+joff,k+koff)+wz(i,j+joff,k+koff)));

        ////////////////////////////////////////////////
        // correct \Gamma^{z+} with \Gamma^{z+,x+}
        ////////////////////////////////////////////////

        if (umac(i+1,j+joff,k+koff) > 0.0) {
           isign = 1.0;
           ioff = 0;
        } else {
           isign = -1.0;
           ioff = 1;
        }

        vv = 0.0;
        if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) > 0.0) {
           vv = vmac(i+ioff,j,k+koff);
        }

        ww = 0.0;
        if (wmac(i,j+joff,k+1)*wmac(i+ioff,j+joff,k+1) > 0.0) {
           ww = wmac(i+ioff,j+joff,k+1);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx;
        p3(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p3(3) = ksign*0.5*hz - wmac(i,j+joff,k+1)*dt;

        p4(1) = isign*0.5*hx - umac(i+1,j+joff,k+koff)*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * umac(i+1,j+joff,k+koff);

        gamma = gamma - dt*gamma2/(3.0*hx);

        ////////////////////////////////////////////////
        // correct \Gamma^{z+} with \Gamma^{z+,x-}
        ////////////////////////////////////////////////

        if (umac(i,j+joff,k+koff) > 0.0) {
           isign = 1.0;
           ioff = -1;
        } else {
           isign = -1.0;
           ioff = 0;
        }

        vv = 0.0;
        if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) > 0.0) {
           vv = vmac(i+ioff,j,k+koff);
        }

        ww = 0.0;
        if (wmac(i,j+joff,k+1)*wmac(i+ioff,j+joff,k+1) > 0.0) {
           ww = wmac(i+ioff,j+joff,k+1);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx;
        p3(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p3(3) = ksign*0.5*hz - wmac(i,j+joff,k+1)*dt;

        p4(1) = isign*0.5*hx - umac(i,j+joff,k+koff)*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * umac(i,j+joff,k+koff);

        gamma = gamma + dt*gamma2/(3.0*hx);

        ////////////////////////////////////////////////
        // correct sedgey with \Gamma^{z+}
        ////////////////////////////////////////////////

        gamma = gamma * wmac(i,j+joff,k+1);
        yedge_tmp = yedge_tmp - dt*gamma/(2.0*hz);

        ////////////////////////////////////////////////
        // compute \Gamma^{z-} without corner corrections
        ////////////////////////////////////////////////

        if (wmac(i,j+joff,k) > 0.0) {
           ksign = 1.0;
           koff = -1;
        } else {
           ksign = -1.0;
           koff = 0;
        }

        v = 0.0;
        if (vmac(i,j,k)*vmac(i,j,k+koff) > 0.0) {
           v = vmac(i,j,k+koff);
        }

        p1(1) = 0.0;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = 0.0;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p2(3) = ksign*0.5*hz;

        p3(1) = 0.0;
        p3(2) = jsign*0.5*hy - v*dt;
        p3(3) = ksign*0.5*hz - wmac(i,j+joff,k)*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p2(ll)+p3(ll))/2.0;
        }
        val1 = eval(s(i,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p3(ll))/2.0;
        }
        val2 = eval(s(i,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll))/2.0;
        }
        val3 = eval(s(i,j+joff,k+koff,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.0;

        // source term
        gamma = gamma*(1. - dt3*(vy(i,j+joff,k+koff)+wz(i,j+joff,k+koff)));

        ////////////////////////////////////////////////
        // correct \Gamma^{z-} with \Gamma^{z-,x+}
        ////////////////////////////////////////////////

        if (umac(i+1,j+joff,k+koff) > 0.0) {
           isign = 1.0;
           ioff = 0;
        } else {
           isign = -1.0;
           ioff = 1;
        }

        vv = 0.0;
        if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) > 0.0) {
           vv = vmac(i+ioff,j,k+koff);
        }

        ww = 0.0;
        if (wmac(i,j+joff,k)*wmac(i+ioff,j+joff,k) > 0.0) {
           ww = wmac(i+ioff,j+joff,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx;
        p3(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p3(3) = ksign*0.5*hz - wmac(i,j+joff,k)*dt;

        p4(1) = isign*0.5*hx - umac(i+1,j+joff,k+koff)*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * umac(i+1,j+joff,k+koff);

        gamma = gamma - dt*gamma2/(3.0*hx);

        ////////////////////////////////////////////////
        // correct \Gamma^{z-} with \Gamma^{z-,x-}
        ////////////////////////////////////////////////

        if (umac(i,j+joff,k+koff) > 0.0) {
           isign = 1.0;
           ioff = -1;
        } else {
           isign = -1.0;
           ioff = 0;
        }

        vv = 0.0;
        if (vmac(i,j,k)*vmac(i+ioff,j,k+koff) > 0.0) {
           vv = vmac(i+ioff,j,k+koff);
        }

        ww = 0.0;
        if (wmac(i,j+joff,k)*wmac(i+ioff,j+joff,k) > 0.0) {
           ww = wmac(i+ioff,j+joff,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p2(3) = ksign*0.5*hz;

        p3(1) = isign*0.5*hx;
        p3(2) = jsign*0.5*hy - vmac(i,j,k)*dt;
        p3(3) = ksign*0.5*hz - wmac(i,j+joff,k)*dt;

        p4(1) = isign*0.5*hx - umac(i,j+joff,k+koff)*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * umac(i,j+joff,k+koff);

        gamma = gamma + dt*gamma2/(3.0*hx);

        ////////////////////////////////////////////////
        // correct sedgey with \Gamma^{z-}
        ////////////////////////////////////////////////

        gamma = gamma * wmac(i,j+joff,k);
        sedgey(i,j,k,icomp) = yedge_tmp + dt*gamma/(2.0*hz);
    });

    // compute sedgez on z-faces
    Box const& zbx = amrex::surroundingNodes(bx,2);
    ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

        // set edge values equal to the ghost cell value since they store the physical condition on the boundary
        if ( k==dlo.z && lo_z_physbc ) {
            sedgez(i,j,k,icomp) = s(i,j,k-1,icomp);
            return;
        }
        if ( k==dhi.z+1 && hi_z_physbc ) {
            sedgez(i,j,k,icomp) = s(i,j,k,icomp);
            return;
        }

        //local variables
        Array1D<Real, 1, 3> del;
        Array1D<Real, 1, 3> p1;
        Array1D<Real, 1, 3> p2;
        Array1D<Real, 1, 3> p3;
        Array1D<Real, 1, 3> p4;

        Array1D<Real, 1, 7> slope_tmp;

        int ioff, joff, koff;

        Real isign,jsign,ksign;
        Real val1,val2,val3,val4,val5;
        Real w;
        Real uu,vv,ww;
        Real gamma,gamma2;

        // Hold intermediate edgestate value
        Real zedge_tmp;

        ////////////////////////////////////////////////
        // compute sedgez without transverse corrections
        ////////////////////////////////////////////////

        // centroid of rectangular volume
        if (wmac(i,j,k) > 0.0) {
           ksign = 1.0;
           koff = -1;
        } else {
           ksign = -1.0;
           koff = 0;
        }

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i,j,k+koff,n-1);
        }

        del(1) = 0.0;
        del(2) = 0.0;
        del(3) = ksign*0.5*hz - 0.5*wmac(i,j,k)*dt;
        zedge_tmp = eval(s(i,j,k+koff,icomp),slope_tmp,del);


        // source term
        zedge_tmp = zedge_tmp*(1. - dt2*wz(i,j,k+koff));
        if (force) {
            zedge_tmp += dt2*force(i,j,k+koff,icomp);
        }


        ////////////////////////////////////////////////
        // compute \Gamma^{x+} without corner corrections
        ////////////////////////////////////////////////

        if (umac(i+1,j,k+koff) > 0.0) {
           isign = 1.0;
           ioff = 0;
        } else {
           isign = -1.0;
           ioff = 1;
        }

        w = 0.0;
        if (wmac(i,j,k)*wmac(i+ioff,j,k) > 0.0) {
           w = wmac(i+ioff,j,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = 0.0;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = 0.0;
        p2(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx - umac(i+1,j,k+koff)*dt;
        p3(2) = 0.0;
        p3(3) = ksign*0.5*hz - w*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p2(ll)+p3(ll))/2.0;
        }
        val1 = eval(s(i+ioff,j,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p3(ll))/2.0;
        }
        val2 = eval(s(i+ioff,j,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll))/2.0;
        }
        val3 = eval(s(i+ioff,j,k+koff,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.0;

        // source term
        gamma = gamma*(1. - dt3*(wz(i+ioff,j,k+koff)+ux(i+ioff,j,k+koff)));

        ////////////////////////////////////////////////
        // correct \Gamma^{x+} with \Gamma^{x+,y+}
        ////////////////////////////////////////////////

        if (vmac(i+ioff,j+1,k+koff) > 0.0) {
           jsign = 1.0;
           joff = 0;
        } else {
           jsign = -1.0;
           joff = 1;
        }

        ww = 0.0;
        if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) > 0.0) {
           ww = wmac(i+ioff,j+joff,k);
        }

        uu = 0.0;
        if (umac(i+1,j,k+koff)*umac(i+1,j+joff,k+koff) > 0.0) {
           uu = umac(i+1,j+joff,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx - umac(i+1,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy;
        p3(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vmac(i+ioff,j+1,k+koff)*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * vmac(i+ioff,j+1,k+koff);

        gamma = gamma - dt*gamma2/(3.0*hy);

        ////////////////////////////////////////////////
        // correct \Gamma^{x+} with \Gamma^{x+,y-}
        ////////////////////////////////////////////////

        if (vmac(i+ioff,j,k+koff) > 0.0) {
           jsign = 1.0;
           joff = -1;
        } else {
           jsign = -1.0;
           joff = 0;
        }

        ww = 0.0;
        if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) > 0.0) {
           ww = wmac(i+ioff,j+joff,k);
        }

        uu = 0.0;
        if (umac(i+1,j,k+koff)*umac(i+1,j+joff,k+koff) > 0) {
           uu = umac(i+1,j+joff,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx - umac(i+1,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy;
        p3(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vmac(i+ioff,j,k+koff)*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * vmac(i+ioff,j,k+koff);

        gamma = gamma + dt*gamma2/(3.0*hy);

        ////////////////////////////////////////////////
        // correct sedgez with \Gamma^{x+}
        ////////////////////////////////////////////////

        gamma = gamma * umac(i+1,j,k+koff);
        zedge_tmp = zedge_tmp - dt*gamma/(2.0*hx);

        ////////////////////////////////////////////////
        // compute \Gamma^{x-} without corner corrections
        ////////////////////////////////////////////////

        if (umac(i,j,k+koff) > 0.0) {
           isign = 1.0;
           ioff = -1;
        } else {
           isign = -1.0;
           ioff = 0;
        }

        w = 0.0;
        if (wmac(i,j,k)*wmac(i+ioff,j,k) > 0.0) {
           w = wmac(i+ioff,j,k);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = 0.0;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = 0.0;
        p2(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx - umac(i,j,k+koff)*dt;
        p3(2) = 0.0;
        p3(3) = ksign*0.5*hz - w*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p2(ll)+p3(ll))/2.0;
        }
        val1 = eval(s(i+ioff,j,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p3(ll))/2.0;
        }
        val2 = eval(s(i+ioff,j,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll))/2.0;
        }
        val3 = eval(s(i+ioff,j,k+koff,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.0;

        // source term
        gamma = gamma*(1. - dt3*(wz(i+ioff,j,k+koff)+ux(i+ioff,j,k+koff)));

        ////////////////////////////////////////////////
        // correct \Gamma^{x-} with \Gamma^{x-,y+}
        ////////////////////////////////////////////////

        if (vmac(i+ioff,j+1,k+koff) > 0.0) {
           jsign = 1.0;
           joff = 0;
        } else {
           jsign = -1.0;
           joff = 1;
        }

        ww = 0.0;
        if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) > 0.0) {
           ww = wmac(i+ioff,j+joff,k);
        }

        uu = 0.0;
        if (umac(i,j,k+koff)*umac(i,j+joff,k+koff) > 0.0) {
           uu = umac(i,j+joff,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx - umac(i,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy;
        p3(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vmac(i+ioff,j+1,k+koff)*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * vmac(i+ioff,j+1,k+koff);

        gamma = gamma - dt*gamma2/(3.0*hy);

        ////////////////////////////////////////////////
        // correct \Gamma^{x-} with \Gamma^{x-,y-}
        ////////////////////////////////////////////////

        if (vmac(i+ioff,j,k+koff) > 0.0) {
           jsign = 1.0;
           joff = -1;
        } else {
           jsign = -1.0;
           joff = 0;
        }

        ww = 0.0;
        if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) > 0.0) {
           ww = wmac(i+ioff,j+joff,k);
        }

        uu = 0.0;
        if (umac(i,j,k+koff)*umac(i,j+joff,k+koff) > 0.0) {
           uu = umac(i,j+joff,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx - umac(i,j+joff,k)*dt;
        p3(2) = jsign*0.5*hy;
        p3(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p4(1) = isign*0.5*hx - uu*dt;
        p4(2) = jsign*0.5*hy - vmac(i+ioff,j,k+koff)*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * vmac(i+ioff,j,k+koff);

        gamma = gamma + dt*gamma2/(3.0*hy);

        ////////////////////////////////////////////////
        // correct sedgez with \Gamma^{x-}
        ////////////////////////////////////////////////

        gamma = gamma * umac(i,j,k+koff);
        zedge_tmp = zedge_tmp + dt*gamma/(2.0*hx);

        ////////////////////////////////////////////////
        // compute \Gamma^{y+} without corner corrections
        ////////////////////////////////////////////////

        if (vmac(i,j+1,k+koff) > 0.0) {
           jsign = 1.0;
           joff = 0;
        } else {
           jsign = -1.0;
           joff = 1;
        }

        w = 0.0;
        if (wmac(i,j,k)*wmac(i,j+joff,k) > 0.0) {
           w = wmac(i,j+joff,k);
        }

        p1(1) = 0.0;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = 0.0;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p3(1) = 0.0;
        p3(2) = jsign*0.5*hy - vmac(i,j+1,k+koff)*dt;
        p3(3) = ksign*0.5*hz - w*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p2(ll)+p3(ll))/2.0;
        }
        val1 = eval(s(i,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p3(ll))/2.0;
        }
        val2 = eval(s(i,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll))/2.0;
        }
        val3 = eval(s(i,j+joff,k+koff,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.0;

        // source term
        gamma = gamma*(1. - dt3*(wz(i,j+joff,k+koff)+vy(i,j+joff,k+koff)));

        ////////////////////////////////////////////////
        // correct \Gamma^{y+} with \Gamma^{y+,x+}
        ////////////////////////////////////////////////

        if (umac(i+1,j+joff,k+koff) > 0.0) {
           isign = 1.0;
           ioff = 0;
        } else {
           isign = -1.0;
           ioff = 1;
        }

        ww = 0.0;
        if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) > 0.0) {
           ww = wmac(i+ioff,j+joff,k);
        }

        vv = 0.0;
        if (vmac(i,j+1,k+koff)*vmac(i+ioff,j+1,k+koff) > 0.0) {
           vv = vmac(i+ioff,j+1,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j+1,k)*dt;
        p3(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p4(1) = isign*0.5*hx - umac(i+1,j+joff,k+koff)*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * umac(i+1,j+joff,k+koff);

        gamma = gamma - dt*gamma2/(3.0*hx);

        ////////////////////////////////////////////////
        // correct \Gamma^{y+} with \Gamma^{y+,x-}
        ////////////////////////////////////////////////

        if (umac(i,j+joff,k+koff) > 0.0) {
           isign = 1.0;
           ioff = -1;
        } else {
           isign = -1.0;
           ioff = 0;
        }

        ww = 0.0;
        if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) > 0.0) {
           ww = wmac(i+ioff,j+joff,k);
        }

        vv = 0.0;
        if (vmac(i,j+1,k+koff)*vmac(i+ioff,j+1,k+koff) > 0.0) {
           vv = vmac(i+ioff,j+1,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j+1,k)*dt;
        p3(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p4(1) = isign*0.5*hx - umac(i,j+joff,k+koff)*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * umac(i,j+joff,k+koff);

        gamma = gamma + dt*gamma2/(3.0*hx);

        ////////////////////////////////////////////////
        // correct sedgez with \Gamma^{y+}
        ////////////////////////////////////////////////

        gamma = gamma * vmac(i,j+1,k+koff);
        zedge_tmp = zedge_tmp - dt*gamma/(2.0*hy);

        ////////////////////////////////////////////////
        // compute \Gamma^{y-} without corner corrections
        ////////////////////////////////////////////////

        if (vmac(i,j,k+koff) > 0.0) {
           jsign = 1.0;
           joff = -1;
        } else {
           jsign = -1.0;
           joff = 0;
        }

        w = 0.0;
        if (wmac(i,j,k)*wmac(i,j+joff,k) > 0.0) {
           w = wmac(i,j+joff,k);
        }

        p1(1) = 0.0;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = 0.0;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p3(1) = 0.0;
        p3(2) = jsign*0.5*hy - vmac(i,j,k+koff)*dt;
        p3(3) = ksign*0.5*hz - w*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p2(ll)+p3(ll))/2.0;
        }
        val1 = eval(s(i,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p3(ll))/2.0;
        }
        val2 = eval(s(i,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll))/2.0;
        }
        val3 = eval(s(i,j+joff,k+koff,icomp),slope_tmp,del);

        // average these centroid values to get the average value
        gamma = (val1+val2+val3)/3.0;

        // source term
        gamma = gamma*(1. - dt3*(wz(i,j+joff,k+koff)+vy(i,j+joff,k+koff)));

        ////////////////////////////////////////////////
        // correct \Gamma^{y-} with \Gamma^{y-,x+};
        ////////////////////////////////////////////////

        if (umac(i+1,j+joff,k+koff) > 0.0) {
           isign = 1.0;
           ioff = 0;
        } else {
           isign = -1.0;
           ioff = 1;
        }

        ww = 0.0;
        if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) > 0.0) {
           ww = wmac(i+ioff,j+joff,k);
        }

        vv = 0.0;
        if (vmac(i,j,k+koff)*vmac(i+ioff,j,k+koff) > 0.0) {
           vv = vmac(i+ioff,j,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j,k)*dt;
        p3(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p4(1) = isign*0.5*hx - umac(i+1,j+joff,k+koff)*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * umac(i+1,j+joff,k+koff);

        gamma = gamma - dt*gamma2/(3.0*hx);

        ////////////////////////////////////////////////
        // correct \Gamma^{y-} with \Gamma^{y-,x-}
        ////////////////////////////////////////////////

        if (umac(i,j+joff,k+koff) > 0.0) {
           isign = 1.0;
           ioff = -1;
        } else {
           isign = -1.0;
           ioff = 0;
        }

        ww = 0.0;
        if (wmac(i,j,k)*wmac(i+ioff,j+joff,k) > 0.0) {
           ww = wmac(i+ioff,j+joff,k);
        }

        vv = 0.0;
        if (vmac(i,j,k+koff)*vmac(i+ioff,j,k+koff) > 0.0) {
           vv = vmac(i+ioff,j,k+koff);
        }

        p1(1) = isign*0.5*hx;
        p1(2) = jsign*0.5*hy;
        p1(3) = ksign*0.5*hz;

        p2(1) = isign*0.5*hx;
        p2(2) = jsign*0.5*hy;
        p2(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p3(1) = isign*0.5*hx;
        p3(2) = jsign*0.5*hy - vmac(i+ioff,j,k)*dt;
        p3(3) = ksign*0.5*hz - wmac(i,j,k)*dt;

        p4(1) = isign*0.5*hx - umac(i,j+joff,k+koff)*dt;
        p4(2) = jsign*0.5*hy - vv*dt;
        p4(3) = ksign*0.5*hz - ww*dt;

        for(int n=1; n<=7; ++n){
            slope_tmp(n) = slopes(i+ioff,j+joff,k+koff,n-1);
        }

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.0;
        }
        val1 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll));
        }
        val2 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll));
        }
        val3 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll));
        }
        val4 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        for(int ll=1; ll<=3; ++ll ){
           del(ll) = 0.5*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll));
        }
        val5 = eval(s(i+ioff,j+joff,k+koff,icomp),slope_tmp,del);

        gamma2 = -0.8*val1 + 0.45*(val2+val3+val4+val5);

        // divu source term
        gamma2 = gamma2*(1. - dt4*(ux(i+ioff,j+joff,k+koff)+vy(i+ioff,j+joff,k+koff)+wz(i+ioff,j+joff,k+koff)));

        gamma2 = gamma2 * umac(i,j+joff,k+koff);

        gamma = gamma + dt*gamma2/(3.0*hx);

        ////////////////////////////////////////////////
        // correct sedgez with \Gamma^{y-}
        ////////////////////////////////////////////////

        gamma = gamma * vmac(i,j,k+koff);
        sedgez(i,j,k,icomp) = zedge_tmp + dt*gamma/(2.0*hy);
    });
}
#endif

// project edge states to satisfy EOS
void BDS_Proj(const Box& bx,
              int ncomp,
              Array4<Real> const& sedgex,
              Array4<Real> const& sedgey,
#if (AMREX_SPACEDIM == 3)
              Array4<Real> const& sedgez,
#endif
              const int proj_type)
{

   Box const& xbx = amrex::surroundingNodes(bx,0);
   ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

     GpuArray<Real, MAX_SPECIES> sedge;
     for (int n=0; n<ncomp; ++n) {
       sedge[n] = sedgex(i,j,k,n);
     }

     BDS_Proj_local(sedge,ncomp,proj_type);

     for (int n=0; n<ncomp; ++n) {
       sedgex(i,j,k,n) = sedge[n];
     }

   });

   Box const& ybx = amrex::surroundingNodes(bx,1);
   ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

     GpuArray<Real, MAX_SPECIES> sedge;
     for (int n=0; n<ncomp; ++n) {
       sedge[n] = sedgey(i,j,k,n);
     }

     BDS_Proj_local(sedge,ncomp,proj_type);

     for (int n=0; n<ncomp; ++n) {
       sedgey(i,j,k,n) = sedge[n];
     }


   });

#if (AMREX_SPACEDIM == 3)
   Box const& zbx = amrex::surroundingNodes(bx,2);
   ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

     GpuArray<Real, MAX_SPECIES> sedge;
     for (int n=0; n<ncomp; ++n) {
       sedge[n] = sedgez(i,j,k,n);
     }

     BDS_Proj_local(sedge,ncomp,proj_type);

     for (int n=0; n<ncomp; ++n) {
       sedgez(i,j,k,n) = sedge[n];
     }


   });
#endif

}
