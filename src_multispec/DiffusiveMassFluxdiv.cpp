#include "multispec_functions.H"
#include "multispec_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "multispec_namespace.H"
#include "common_namespace.H"

using namespace multispec;
using namespace common;
using namespace amrex;

void DiffusiveMassFluxdiv(const MultiFab& rho,
			  const MultiFab& rhotot,
			  const MultiFab& molarconc,
			  const MultiFab& rhoWchi,
			  const MultiFab& Gamma,
			  MultiFab& diff_mass_fluxdiv,
			  std::array< MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
			  const Geometry& geom)
{

    BL_PROFILE_VAR("DiffusiveMassFluxdiv()",DiffusiveMassFluxdiv);

    // compute the face-centered flux (each direction: cells+1 faces while 
    // cells contain interior+2 ghost cells) 
    DiffusiveMassFlux(rho,rhotot,molarconc,rhoWchi,Gamma,diff_mass_flux,geom);

    // compute divergence of determinstic flux 
    ComputeDiv(diff_mass_flux,diff_mass_fluxdiv,geom,nspecies);

}

void DiffusiveMassFlux(const MultiFab& rho,
		       const MultiFab& rhotot,
		       const MultiFab& molarconc,
		       const MultiFab& rhoWchi,
		       const MultiFab& Gamma,
		       std::array< MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
		       const Geometry& geom)
{

    BL_PROFILE_VAR("DiffusiveMassFlux()",DiffusiveMassFluxdiv);

    int i;
    BoxArray ba = rho.boxArray();
    DistributionMapping dmap = rho.DistributionMap();
    int nspecies = rho.nComp();
    int nspecies2 = nspecies*nspecies;

    const Real* dx = geom.CellSize();

    // build local face-centered multifab with nspecies^2 component, zero ghost cells 
    // and nodal in direction i

    // rho*W*chi*Gamma
    std::array< MultiFab, AMREX_SPACEDIM > rhoWchi_face;
    AMREX_D_TERM(rhoWchi_face[0].define(convert(ba,nodal_flag_x), dmap, nspecies2, 0);,
                 rhoWchi_face[1].define(convert(ba,nodal_flag_y), dmap, nspecies2, 0);,
                 rhoWchi_face[2].define(convert(ba,nodal_flag_z), dmap, nspecies2, 0););

    // Gamma-matrix
    std::array< MultiFab, AMREX_SPACEDIM > Gamma_face;
    AMREX_D_TERM(Gamma_face[0].define(convert(ba,nodal_flag_x), dmap, nspecies2, 0);,
                 Gamma_face[1].define(convert(ba,nodal_flag_y), dmap, nspecies2, 0);,
                 Gamma_face[2].define(convert(ba,nodal_flag_z), dmap, nspecies2, 0););

    AverageCCToFace(rhoWchi, 0, rhoWchi_face, 0, nspecies2);
    AverageCCToFace(Gamma, 0, Gamma_face, 0, nspecies2);
    // Note: Add shifting option?

    //Computes gradient at cell faces of cell centred scalar
    ComputeGrad(molarconc, diff_mass_flux, 0, 0, nspecies, geom);

    // MatvecMul needs to add A*x result to x
    for(i=0;i<AMREX_SPACEDIM;i++) {
      MatvecMul(diff_mass_flux[i], Gama_face[i], nspecies);
    }

    for(i=0;i<AMREX_SPACEDIM;i++) {
      MatvecMul(diff_mass_flux[i], rhoWchi_face[i], nspecies);
    }

    // //correct fluxes to ensure mass conservation to roundoff
    // if (correct_flux==1 && (nspecies > 1)) {
    //   // Print() << "Checking conservation of deterministic fluxes \n";
    //   correction_flux(mla, rho, rhotot, diff_mass_flux, the_bc_tower%bc_tower_array);
    // }

}
