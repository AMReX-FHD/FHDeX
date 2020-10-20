#include "multispec_functions.H"
#include "common_functions.H"

void ElectroDiffusiveMassFluxdiv(const MultiFab& rho,
                                 const MultiFab& Temp,
                                 const MultiFab& rhoWchi,
                                 std::array< MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
                                 MultiFab& diff_mass_fluxdiv,
                                 std::array< const MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                                 MultiFab& charge,
                                 std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot,
                                 MultiFab& Epot,
                                 const MultiFab& permittivity,
                                 Real dt,
                                 int zero_initial_Epot,
                                 const Geometry& geom)
{
    BoxArray ba = rho.boxArray();
    DistributionMapping dmap = rho.DistributionMap();
    
    // electro_mass_flux(d) is face-centered, has nspecies component, zero ghost 
    // cells & nodal in direction d
    std::array< MultiFab, AMREX_SPACEDIM > electro_mass_flux;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        electro_mass_flux[d].define(convert(ba,nodal_flag_dir[d]), dmap, nspecies, 0);
    }
  
    // compute the face-centered electro_mass_flux (each direction: cells+1 faces while 
    // cells contain interior+2 ghost cells)
    ElectroDiffusiveMassFlux(rho,Temp,rhoWchi,electro_mass_flux,diff_mass_flux,
                             stoch_mass_flux,charge,grad_Epot,Epot,permittivity,
                             dt,zero_initial_Epot,geom);


    
    // add fluxes to diff_mass_flux
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        MultiFab::Add(diff_mass_flux[i],electro_mass_flux[i],0,0,nspecies,0);
    }

    // add flux divergence to diff_mass_fluxdiv
    ComputeDiv(diff_mass_fluxdiv,electro_mass_flux,0,0,nspecies,geom,1);    
}


void ElectroDiffusiveMassFlux(const MultiFab& rho,
                              const MultiFab& Temp,
                              const MultiFab& rhoWchi,
                              std::array< MultiFab, AMREX_SPACEDIM >& electro_mass_flux,
                              std::array< MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
                              std::array< const MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                              MultiFab& charge,
                              std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot,
                              MultiFab& Epot,
                              const MultiFab& permittivity,
                              Real dt,
                              int zero_initial_Epot,
                              const Geometry& geom)
{
    BoxArray ba = rho.boxArray();
    DistributionMapping dmap = rho.DistributionMap();

    MultiFab alpha      (ba,dmap,       1,0);
    MultiFab charge_coef(ba,dmap,nspecies,1);
    MultiFab rhs        (ba,dmap,       1,0);

    std::array< MultiFab, AMREX_SPACEDIM > beta;
    std::array< MultiFab, AMREX_SPACEDIM > rhoWchi_fc;
    std::array< MultiFab, AMREX_SPACEDIM > permittivity_fc;
    std::array< MultiFab, AMREX_SPACEDIM > E_ext;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        beta           [d].define(convert(ba,nodal_flag_dir[d]), dmap,                 1, 0);
        rhoWchi_fc     [d].define(convert(ba,nodal_flag_dir[d]), dmap, nspecies*nspecies, 0);
        permittivity_fc[d].define(convert(ba,nodal_flag_dir[d]), dmap,                 1, 0);
        E_ext          [d].define(convert(ba,nodal_flag_dir[d]), dmap,                 1, 0);
    }

    // if periodic, ensure charge sums to zero by subtracting off the average
    if (geom.isAllPeriodic()) {
        Real sum = charge.sum() / ba.numPts();
        charge.plus(-sum,0,1);
        if (amrex::Math::abs(sum) > 1.e-12) {
            Print() << "average charge = " << sum << std::endl;
            Warning("Warning: electrodiffusive_mass_flux - average charge is not zero");
        }
    }

    bool any_shift = false;
    for (int i=0; i<AMREX_SPACEDIM*LOHI; ++i) {
        if (shift_cc_to_boundary[i] == 1) any_shift=true;
    }

    // compute face-centered rhoWchi from cell-centered values 
    if (any_shift) {
        Abort("ShiftCCToBoundaryFace not written yet");        
    } else {
        AverageCCToFace(rhoWchi,rhoWchi_fc,0,nspecies*nspecies,1,geom);
    }
    
    // solve poisson equation for phi (the electric potential)
    // -del dot epsilon grad Phi = charge
    if (zero_initial_Epot) {
        Epot.setVal(0.);
    }

    // fill ghost cells for Epot at walls using Dirichlet value
    MultiFabPhysBC(Epot,geom,0,1,3);

    // set alpha=0
    alpha.setVal(0.);

    if (electroneutral==0 || E_ext_type != 0) {
        // permittivity on faces
        if (any_shift) {
            Abort("ShiftCCToBoundaryFace not written yet");
        } else {
            AverageCCToFace(permittivity,permittivity_fc,0,1,1,geom);            
        }
    }



    if (electroneutral == 1) {
        // For electroneutral we only support homogeneous Neumann BCs for potential
        // This is the correct Poisson BC for impermeable walls
        // For reservoirs, the BCs are actually inhomogeneous but computed on-the-fly by the code later on
        // Here we setup just the homogeneous Poisson problem -- this is all that the multigrid solver can handle     
        // Reactive walls are not yet supported
        Abort("ElectroDiffusiveMassFluxdiv.cpp: electroneutral not written yet");
    } else {

        // non-electroneutral

        // set beta=permittivity (epsilon)
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            MultiFab::Copy(beta[i],permittivity_fc[i],0,0,1,0);
        }

        if (zero_eps_on_wall_type > 0) {
            Abort("ElectroDiffusiveMassFluxdiv.cpp: zero_eps_on_wall_type > 0 not written yet");
        }

        // set rhs equal to charge
        MultiFab::Copy(rhs,charge,0,0,1,0);

        

    }


}
