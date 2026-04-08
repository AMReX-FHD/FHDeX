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

    advance_phi(phi_old[lev], phi_new[lev], fluxes, stochFluxes, dt_lev, npts_scale, geom[lev], bcs,
            m_ensemble_dir, m_ext_pot, m_ext_pot_alpha, m_ext_pot_beta, m_ext_pot_gamma);

    // Increment or decrement the flux registers by area and time-weighted fluxes
    // Note that the fluxes have already been scaled by dt and area
    // In this example we are solving phi_t = -div(+F)
    // The fluxes contain, e.g., F_{i+1/2,j} = (phi*u)_{i+1/2,j}
    // Keep this in mind when considering the different sign convention for updating
    // the flux registers from the coarse or fine grid perspective
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
