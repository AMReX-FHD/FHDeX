#include <AmrCoreAdv.H>
#include <Kernels.H>
#include <myfunc.H>

using namespace amrex;

// Advance a single level for a single time step, updates flux registers
void
AmrCoreAdv::AdvancePhiAtLevel (int lev, Real /*time*/, Real dt_lev, int /*iteration*/, int /*ncycle*/)
{
    Array<MultiFab, AMREX_SPACEDIM> fluxes;
    Array<MultiFab, AMREX_SPACEDIM> stochFluxes;
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
    {
        BoxArray ba = grids[lev];
        ba.surroundingNodes(i);
        fluxes[i].define(ba, dmap[lev], phi_new[lev].nComp(), 0);
        fluxes[i].setVal(0.);
        stochFluxes[i].define(ba, dmap[lev], phi_new[lev].nComp(), 0);
        stochFluxes[i].setVal(0.);
    }

    phi_old[lev].FillBoundary(Geom(lev).periodicity());

    // We do this here so we can print the FABs for debugging
    phi_new[lev].setVal(0.0);

    advance_phi(phi_old[lev], phi_new[lev], fluxes, stochFluxes, dt_lev, npts_scale, geom[lev], bcs, pot);

    // Increment or decrement the flux registers by area and time-weighted fluxes
    // Note that the fluxes have already been scaled by dt and area
    // The update in update_phi is phi_new = phi_old + dt/dx (F_{i+1/2} - F_{i-1/2}),
    // i.e. phi_t = +div(F) with F = (1/2) grad(phi) + phi grad(V)/gamma + noise.
    // CrseInit with scale +1 and Reflux with scale +1 then remove exactly this
    // coarse flux across every coarse/fine face from the adjacent uncovered
    // coarse cells (AMReX refluxes with -F on a low face and +F on a high face).
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    if (do_reflux) {
        if (flux_reg[lev+1]) {
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                // update the lev+1/lev flux register (index lev+1)
                flux_reg[lev+1]->CrseInit(fluxes[i],i,0,0,fluxes[i].nComp(),1.0);
            }
        }
        if (flux_reg[lev]) {
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                // update the lev/lev-1 flux register (index lev)
                flux_reg[lev]->FineAdd(fluxes[i],i,0,0,fluxes[i].nComp(),-1.0);
            }
        }
    }
}
