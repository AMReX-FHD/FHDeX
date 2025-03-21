#include "myfunc.H"
#include "mykernel.H"
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>
#include "common_functions.H"
#include "rng_functions.H"

void advance (MultiFab& phi_old,
              MultiFab& phi_new,
              Array<MultiFab, AMREX_SPACEDIM>& flux,
              Array<MultiFab, AMREX_SPACEDIM>& stochFlux,
              Real dt,
	      Real npts_scale,
              Geometry const& geom,
              Vector<BCRec> const& BoundaryCondition,
              int Ncomp, amrex::Real phileft, amrex::Real phiright, int avg_type)
{

    // Fill the ghost cells of each grid from the other grids
    // includes periodic domain boundaries
    phi_old.FillBoundary(geom.periodicity());
    // There are physical boundaries to fill.
    FillDomainBoundary(phi_old, geom, BoundaryCondition);

    const BCRec& bc = BoundaryCondition[0];

    //
    // Note that this simple example is not optimized.
    // The following two MFIter loops could be merged
    // and we do not have to use flux MultiFab.
    //
    // =======================================================

    // This example supports both 2D and 3D.  Otherwise,
    // we would not need to use AMREX_D_TERM.
    AMREX_D_TERM(const Real dxinv = geom.InvCellSize(0);,
                 const Real dyinv = geom.InvCellSize(1);,
                 const Real dzinv = geom.InvCellSize(2););

    const Box& domain_bx = geom.Domain();
    const auto dom_lo = lbound(domain_bx);
    const auto dom_hi = ubound(domain_bx);

    //Real variance = dxinv*dyinv/(npts_scale*dt);
    Real variance = dxinv*dyinv/dt;

#if(AMREX_SPACEDIM > 2)
    variance *=dzinv;
#endif
//    amrex::Print() << dxinv << " " << dyinv <<  std::endl;
//    amrex::Print() << npts_scale  << std::endl;
//    amrex::Print() << " variance " << variance << std::endl;

    // fill random numbers (can skip density component 0)
//    for(int d=0;d<AMREX_SPACEDIM;d++) {
//            MultiFabFillRandom(stochFlux[d], 0, variance, geom);
//    }
            MultiFabFillRandom(stochFlux[0], 0, variance, geom);


    // Compute fluxes one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& xbx = mfi.nodaltilebox(0);
        auto const& fluxx = flux[0].array(mfi);

//        const Box& ybx = mfi.nodaltilebox(1);
//        auto const& fluxy = flux[1].array(mfi);
#if (AMREX_SPACEDIM > 2)
        const Box& zbx = mfi.nodaltilebox(2);
        auto const& fluxz = flux[2].array(mfi);
#endif

        auto const& stochfluxx = stochFlux[0].array(mfi);
//        auto const& stochfluxy = stochFlux[1].array(mfi);
#if (AMREX_SPACEDIM > 2)
        auto const& stochfluxz = stochFlux[2].array(mfi);
#endif
        const Box& bx = mfi.validbox();
        const auto lo = lbound(bx);
        const auto hi = ubound(bx);

        auto const& phi = phi_old.array(mfi);




        amrex::ParallelFor(xbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                compute_flux_x(i,j,k,fluxx,stochfluxx,phi,dxinv,
                               lo.x, hi.x, dom_lo.x, dom_hi.x, bc.lo(0), bc.hi(0),Ncomp,
                               phileft,phiright,avg_type);
            });

//        amrex::ParallelFor(ybx,
//            [=] AMREX_GPU_DEVICE (int i, int j, int k)
//            {
//                compute_flux_y(i,j,k,fluxy,stochfluxy,phi,dyinv,
//                               lo.y, hi.y, dom_lo.y, dom_hi.y, bc.lo(1), bc.hi(1),Ncomp);
//            });
#if (AMREX_SPACEDIM > 2)
        amrex::ParallelFor(zbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                compute_flux_z(i,j,k,fluxz,stochfluxz,phi,dzinv,
                               lo.z, hi.z, dom_lo.z, dom_hi.z, bc.lo(2), bc.hi(2),Ncomp);
            });
#endif
    }

    // Advance the solution one grid at a time
    for ( MFIter mfi(phi_old); mfi.isValid(); ++mfi )
    {
        const Box& vbx = mfi.validbox();
        auto const& fluxx = flux[0].array(mfi);
        auto const& fluxy = flux[1].array(mfi);
#if (AMREX_SPACEDIM > 2)
        auto const& fluxz = flux[2].array(mfi);
#endif
        auto const& phiOld = phi_old.array(mfi);
        auto const& phiNew = phi_new.array(mfi);

        amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            update_phi(i,j,k,phiOld,phiNew,
                       AMREX_D_DECL(fluxx,fluxy,fluxz),
                       dt,
                       AMREX_D_DECL(dxinv,dyinv,dzinv),
                       Ncomp);
        });
    }
}

void init_phi(amrex::MultiFab& phi_new, amrex::Geometry const& geom, amrex::Real npts_scale,int Ncomp,
              amrex::Real phileft, amrex::Real phiright){

    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
    GpuArray<Real,AMREX_SPACEDIM> prob_lo = geom.ProbLoArray();
    // =======================================
    // Initialize phi_new by calling a Fortran routine.
    // MFIter = MultiFab Iterator
    for (MFIter mfi(phi_new); mfi.isValid(); ++mfi)
    {
        const Box& vbx = mfi.validbox();
        auto const& phiNew = phi_new.array(mfi);
        amrex::ParallelFor(vbx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            init_phi(i,j,k,phiNew,dx,prob_lo, npts_scale,Ncomp,phileft, phiright);
        });
    }
}
