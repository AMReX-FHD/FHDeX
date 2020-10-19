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
//    ElectroDiffusiveMassFlux();


    
    // add fluxes to diff_mass_flux
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        MultiFab::Add(diff_mass_flux[i],electro_mass_flux[i],0,0,nspecies,0);
    }

    // add flux divergence to diff_mass_fluxdiv
    ComputeDiv(diff_mass_fluxdiv,electro_mass_flux,0,0,nspecies,geom,1);    
}
