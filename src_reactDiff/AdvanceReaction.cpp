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

    // only used for reactDiff_reaction_type = 1
    Vector<Real> mattingly_lin_comb_coef(2);
    mattingly_lin_comb_coef[0] = 1.;
    mattingly_lin_comb_coef[1] = 0.;

    if (reactDiff_reaction_type == 0) { // first-order det/tau-leaping/CLE, or SSA

        // calculate rates
        // rates could be deterministic or stochastic depending on reaction_type
        ChemicalRates(n_old,rate,geom,dt,n_old,mattingly_lin_comb_coef,volume_factor);

        MultiFab::LinComb(n_new,1,n_old,0,dt,rate,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,-dt,ext_src,0,0,nspecies,0); //note the negative sign

        n_new.FillBoundary(geom.periodicity());
        MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);

    } else if (reactDiff_reaction_type == 1) { // second-order det/tau-leaping/CLE

        // here we use Mattingly's predictor-corrector with theta=0.5d0 (for rection_type=1).
        // with these parameters this is actually equivalent to a traditional midpoint scheme.
        Real theta = 0.5;
        Real alpha1 = 2.;
        Real alpha2 = 1.;

        //!!!!!!!!!!!!!!
        // predictor   !
        //!!!!!!!!!!!!!!

        // calculate rates from a(n_old)
        ChemicalRates(n_old,rate,geom,theta*dt,n_old,mattingly_lin_comb_coef,volume_factor);

        // update
        MultiFab::LinComb(n_new,1,n_old,0,theta*dt,rate,0,0,nspecies,0);
        MultiFab::Saxpy(n_new,-theta*dt,ext_src,0,0,nspecies,0); //note the negative sign
        n_new.FillBoundary(geom.periodicity());
        MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);

        //!!!!!!!!!!!!!!
        // corrector   !
        //!!!!!!!!!!!!!!

        // Here we write this in the form that Mattingly et al do
        //  where we just continue the second half of the time step from where we left

        mattingly_lin_comb_coef[0] = -alpha2;
        mattingly_lin_comb_coef[1] = alpha1;

        // calculate rates from 2*a(n_pred)-a(n_old)
        ChemicalRates(n_old,rate,geom,(1.-theta)*dt,n_new,mattingly_lin_comb_coef,volume_factor);

        // update
        MultiFab::Saxpy(n_new,(1.-theta)*dt,rate,0,0,nspecies,0);
        // note the negative sign
        // also note that ext_src does not change in the time interval (t,t+dt)
        MultiFab::Saxpy(n_new,-(1.-theta)*dt*(alpha1-alpha2),ext_src,0,0,nspecies,0);
        n_new.FillBoundary(geom.periodicity());
        MultiFabPhysBC(n_new, geom, 0, nspecies, SPEC_BC_COMP, time);

    } else {
        Abort("AdvanceReaction() - invalid reactDiff_reaction_type");
    }

}
