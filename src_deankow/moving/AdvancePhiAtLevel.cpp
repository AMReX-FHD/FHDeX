#include <AmrCoreAdv.H>
#include <Kernels.H>
#include <myfunc.H>
#include <update_metric.H>
#include <AMReX_Print.H>

using namespace amrex;

// Advance a single level for a single time step, updates flux registers
void
AmrCoreAdv::AdvancePhiAtLevel (int lev, Real time, Real dt_lev, int /*iteration*/, int /*ncycle*/)
{
    Array<MultiFab, AMREX_SPACEDIM> fluxes;
    Array<MultiFab, AMREX_SPACEDIM> stochFluxes;
    for (int i = 0; i < AMREX_SPACEDIM; ++i)
    {
        BoxArray ba = grids[lev];
        stochFluxes[i].define(ba, dmap[lev], num_flux, 1);
        stochFluxes[i].setVal(0.);
        ba.surroundingNodes(i);
        fluxes[i].define(ba, dmap[lev], num_flux*phi_new[lev].nComp(), 0);
    }

    phi_old[lev].FillBoundary(Geom(lev).periodicity());

    // We do this here so we can print the FABs for debugging
    phi_new[lev].setVal(0.0);

    if( lev == 1 ) {
        amrex::Print() << "NOT SURE HOW I GOT HERE" << std::endl;
    }

        amrex::Print() << "entering advance_phi" << std::endl;
       amrex::Print() << " " << dt_lev << std::endl;
       amrex::Print() << " time " << time << " " << dt_lev << std::endl;

       const auto dx     = geom[lev].CellSizeArray();
       const auto problo = geom[lev].ProbLoArray();

       for (MFIter mfi(gmetric); mfi.isValid(); ++mfi)
        {
            const Box& gbx = mfi.validbox();
            auto const& gmet_arr = gmetric.array(mfi);
            auto const& gsqr_arr = sqrgmetric.array(mfi);
            auto const& detg_arr = detg.array(mfi);
            auto const& newgmet_arr = newgmetric.array(mfi);
            auto const& newgsqr_arr = newsqrgmetric.array(mfi);
            auto const& newdetg_arr = newdetg.array(mfi);
            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                new_dk_metric(i,j,k,gmet_arr,gsqr_arr,detg_arr,newgmet_arr,newgsqr_arr,newdetg_arr,dx,problo,time);
            });
        }


//    advance_phi(phi_old[lev], phi_new[lev], fluxes, stochFluxes, dt_lev, npts_scale, geom[lev], bcs);

        amrex::Print() << "calling advance_phi" << std::endl;

    advance_phi(phi_old[lev], phi_new[lev], fluxes, stochFluxes, gmetric, sqrgmetric, detg, newgmetric, newsqrgmetric, newdetg,
                dt_lev, num_part, dorand, num_flux, ext_pot, geom[lev], bcs, time);

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
