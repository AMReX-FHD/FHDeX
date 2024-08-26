#include "reactDiff_functions.H"
#include "chemistry_functions.H"

// this solves dn/dt = f(n) - g (note the minus sign for g)
// where f(n) are the chemical production rates (deterministic or stochastic)
// and g=ext_src is a constant (in time) *deterministic* source term.
// to model stochastic particle production (sources) include g in the definition of f instead.
// or add it as a reaction 0->products
void AdvanceReaction(MultiFab& n_old,
                     MultiFab& n_new,
                     const MultiFab& ext_src,
                     const Real& dt,
                     const Real& time,
                     const Geometry& geom) {

    BoxArray ba = n_old.boxArray();
    DistributionMapping dmap = n_old.DistributionMap();

    // if there are no reactions to process, copy n_old to n_new,
    // account for ext_src and return
    if (nreaction < 1) {
        MultiFab::LinComb(n_new,1,n_old,0,-dt,ext_src,0,0,nspecies,0);
        n_new.FillBoundary(geom.periodicity());
        MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);
        return;
    }

    MultiFab rate(ba,dmap,nspecies,0);

    if (reactDiff_reaction_type == 0) { // first-order det/tau-leaping/CLE, or SSA

        // ChemicalRates();
        
        MultiFab::LinComb(n_new,1,n_old,0,-dt,rate,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,-dt,ext_src,0,0,nspecies,0);
/*
      ! calculate rates
      ! rates could be deterministic or stochastic depending on use_Poisson_rng
      call chemical_rates(mla,n_old,rate,dx,dt,vol_fac_in=volume_factor)

      ! update
      do n=1,nlevs
        call multifab_copy_c(n_new(n),1,n_old(n),1,nspecies,0)
        call multifab_saxpy_3(n_new(n),dt,rate(n))
        call multifab_saxpy_3(n_new(n),-dt,ext_src(n))  ! note the negative sign

        call multifab_fill_boundary(n_new(n))
        call multifab_physbc(n_new(n),1,scal_bc_comp,nspecies, &
                             the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
      end do
*/
    } else if (reactDiff_reaction_type == 1) { // second-order det/tau-leaping/CLE

    } else {
        Abort("AdvanceReaction() - invalid reactDiff_reaction_type");
    }

    n_new.FillBoundary(geom.periodicity());
    MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);
}
