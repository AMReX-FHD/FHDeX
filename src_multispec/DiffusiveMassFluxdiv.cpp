#include "multispec_functions.H"
#include "common_functions.H"

// FIXME: Fill ghost cells

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
    ComputeDiv(diff_mass_fluxdiv,diff_mass_flux,0,0,nspecies,geom,0);  // increment = 0

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
    std::array< MultiFab, AMREX_SPACEDIM > rhoWchi_face;    // rho*W*chi*Gamma
    std::array< MultiFab, AMREX_SPACEDIM > Gamma_face;      // Gamma-matrix

    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        rhoWchi_face[d].define(convert(ba,nodal_flag_dir[d]), dmap, nspecies2, 0);
        Gamma_face[d]  .define(convert(ba,nodal_flag_dir[d]), dmap, nspecies2, 0);
    }

    // compute face-centered rhoWchi from cell-centered values 
    AverageCCToFace(rhoWchi, rhoWchi_face, 0, nspecies2, 1, geom);

    // calculate face-centrered grad(molarconc) 
    ComputeGrad(molarconc, diff_mass_flux, 0, 0, nspecies, geom);

    // compute face-centered Gama from cell-centered values 
    AverageCCToFace(Gamma, Gamma_face, 0, nspecies2, 1, geom);

    // compute Gama*grad(molarconc): Gama is nspecies^2 matrix; grad(x) is
    // nspecies component vector 
    for(i=0; i<AMREX_SPACEDIM; i++) {
      MatvecMul(diff_mass_flux[i], Gamma_face[i]);
    }

    if (is_nonisothermal) {
        //
        //
        //
    }

    if (barodiffusion_type > 0) {
        //
        //
        //
    }

    // compute -rhoWchi * (Gamma*grad(x) + ... ) on faces
    for(i=0;i<AMREX_SPACEDIM;i++) {
      MatvecMul(diff_mass_flux[i], rhoWchi_face[i]);
    }

    // If there are walls with zero-flux boundary conditions
    if (is_nonisothermal) {
        ZeroEdgevalWalls(diff_mass_flux, geom, 0, nspecies);
    }

    //correct fluxes to ensure mass conservation to roundoff
    if (correct_flux==1 && (nspecies > 1)) {
        CorrectionFlux(rho,rhotot,diff_mass_flux);
    }

}
