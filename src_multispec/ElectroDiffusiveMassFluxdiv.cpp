#include "multispec_functions.H"
#include "common_functions.H"
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

void ElectroDiffusiveMassFluxdiv(const MultiFab& rho,
                                 const MultiFab& Temp,
                                 const MultiFab& rhoWchi,
                                 std::array< MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
                                 MultiFab& diff_mass_fluxdiv,
                                 std::array< MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
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
                              std::array< MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
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
    for (int i=0; i<AMREX_SPACEDIM*2; ++i) {
        if (shift_cc_to_boundary[i] == 1) any_shift=true;
    }

    // compute face-centered rhoWchi from cell-centered values
    if (any_shift) {
        Abort("ShiftCCToBoundaryFace not written yet");
    } else {
        AverageCCToFace(rhoWchi,rhoWchi_fc,0,nspecies*nspecies,SPEC_BC_COMP,geom);
    }

    // solve poisson equation for phi (the electric potential)
    // -del dot epsilon grad Phi = charge
    if (zero_initial_Epot) {
        Epot.setVal(0.);
    }

    // fill ghost cells for Epot at walls using Dirichlet value
    MultiFabPhysBC(Epot,geom,0,1,EPOT_BC_COMP);

    // set alpha=0
    alpha.setVal(0.);

    if (electroneutral==0 || E_ext_type != 0) {
        // permittivity on faces
        if (any_shift) {
            Abort("ShiftCCToBoundaryFace not written yet");
        } else {
            AverageCCToFace(permittivity,permittivity_fc,0,1,SPEC_BC_COMP,geom);
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

        // for inhomogeneous Neumann bc's, there is no need to put in homogeneous form
        // the AMReX solver supports inhomogeneous Neumann provided the value in the
        // ghost cell contains the Neumann value.

    }

    if (E_ext_type != 0) {

        // compute external electric field on edges
        ComputeE_ext(E_ext);

        // compute epsilon*E_ext
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            MultiFab::Multiply(permittivity_fc[i],E_ext[i],0,0,1,0);
        }

        // compute div (epsilon*E_ext) and SUBTRACT it from solver rhs
        // this needs to be tested with spatially-varying E_ext OR epsilon
        rhs.mult(-1.,0,1,0);
        ComputeDiv(rhs,permittivity_fc,0,0,1,geom,1);
        rhs.mult(-1.,0,1,0);
    }

    // solve (alpha - del dot beta grad) Epot = charge (for electro-explicit)
    //   Inhomogeneous Dirichlet or homogeneous Neumann is OK
    // solve (alpha - del dot beta grad) Epot = z^T F (for electro-neutral)
    //   Only homogeneous Neumann BCs supported

    LinOpBCType lo_linop_bc[3];
    LinOpBCType hi_linop_bc[3];

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        if (bc_es_lo[i] == -1 && bc_es_hi[i] == -1) {
            lo_linop_bc[i] = LinOpBCType::Periodic;
            hi_linop_bc[i] = LinOpBCType::Periodic;
        }
        if(bc_es_lo[i] == 2)
        {
            lo_linop_bc[i] = LinOpBCType::inhomogNeumann;
        }
        if(bc_es_hi[i] == 2)
        {
            hi_linop_bc[i] = LinOpBCType::inhomogNeumann;
        }
        if(bc_es_lo[i] == 1)
        {
            lo_linop_bc[i] = LinOpBCType::Dirichlet;
        }
        if(bc_es_hi[i] == 1)
        {
            hi_linop_bc[i] = LinOpBCType::Dirichlet;
        }
    }

    // create solver opject
    // (alpha - del dot beta grad) Epot = charge (for electro-explicit)
    MLABecLaplacian linop({geom}, {ba}, {dmap});

    //set BCs
    linop.setDomainBC({AMREX_D_DECL(lo_linop_bc[0],
                                    lo_linop_bc[1],
                                    lo_linop_bc[2])},
                      {AMREX_D_DECL(hi_linop_bc[0],
                                    hi_linop_bc[1],
                                    hi_linop_bc[2])});

    // fill in ghost cells with Dirichlet/Neumann values
    // the ghost cells will hold the value ON the boundary
    MultiFabPotentialBC_solver(Epot,geom);

    // tell MLPoisson about these potentially inhomogeneous BC values
    linop.setLevelBC(0, &Epot);

    // this forces the solver to NOT enforce solvability
    // thus if there are Neumann conditions on phi they must
    // be correct or the Poisson solver won't converge
    linop.setEnforceSingularSolvable(false);

    // set alpha=0, beta=1 (will overwrite beta with epsilon next)
    linop.setScalars(0.0, 1.0);

    // set beta=epsilon
    linop.setBCoeffs(0, amrex::GetArrOfConstPtrs(beta));

    // Multi Level Multi Grid
    MLMG mlmg(linop);

    // Solver parameters
    mlmg.setMaxIter(poisson_max_iter);
    mlmg.setVerbose(poisson_verbose);
    mlmg.setBottomVerbose(poisson_bottom_verbose);

    // Do solve
    mlmg.solve({&Epot}, {&rhs}, poisson_rel_tol, 0.0);

    // restore original solver tolerance
    if (electroneutral == 1) {
        Abort("ElectroDiffusiveMassFluxdiv.cpp: electroneutral not written yet");
    }

    // for periodic problems subtract off the average of Epot
    // we can generalize this later for walls
    if (geom.isAllPeriodic()) {
        Real sum = Epot.sum() / ba.numPts();
        Epot.plus(-sum,0,1);
    }

    // fill ghost cells for electric potential
    if (electroneutral == 1) {
        // for electroneutral problems with reservoirs,
        // the inhomogeneous BC for phi must be computed here for each face:
        // grad(phi) = -z^T*F_diffstoch/(z^T*A_Phi)
        //  We only need this to fill in BCs for phi
        Abort("ElectroDiffusiveMassFluxdiv.cpp: electroneutral not written yet");
    }


    MultiFabPhysBC(Epot,geom,0,1,EPOT_BC_COMP);
    Epot.FillBoundary(geom.periodicity());

    // compute the gradient of the electric potential
    ComputeGrad(Epot,grad_Epot,0,0,1,EPOT_BC_COMP,geom,0);

    if (E_ext_type != 0) {
        // add external electric field
        // since E = -grad(Epot), we have to subtract the external field from grad_Epot
        for (int i=0; i<AMREX_SPACEDIM; ++i) {
            MultiFab::Subtract(grad_Epot[i],E_ext[i],0,0,1,0);
        }
    }

    if (zero_eps_on_wall_type > 0) {
        // Set E-field ie grad_Epot to be zero on certain boundary faces.
        // This enforces dphi/dn = 0 on the parts of the (Dirichlet) wall we want.
        Abort("ElectroDiffusiveMassFluxdiv.cpp: zero_eps_on_wall_type > 0 not written yet");
    }

    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        grad_Epot[i].FillBoundary(geom.periodicity());
    }

    // compute the charge flux coefficient
    ComputeChargeCoef(rho,Temp,charge_coef);

    // average charge flux coefficient to faces, store in flux
    for (int i=0; i<shift_cc_to_boundary.size(); ++i) {
        if (shift_cc_to_boundary[i] == 1) {
            Abort("ElectroDiffusiveMassFluxdiv.cpp: shift_cc_to_boundary not implemented");
        }
    }
    AverageCCToFace(charge_coef, electro_mass_flux, 0, nspecies, SPEC_BC_COMP, geom);

    // multiply flux coefficient by gradient of electric potential
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        for (int comp=0; comp<nspecies; ++comp) {
            MultiFab::Multiply(electro_mass_flux[i], grad_Epot[i], 0, comp, 1, 0);
        }
    }

    if (use_multiphase == 1) {
        LimitEMF(rho,electro_mass_flux);
    }

    // compute -rhoWchi * (... ) on faces
    for (int i=0; i<AMREX_SPACEDIM; ++i) {
        MatvecMul(electro_mass_flux[i],rhoWchi_fc[i]);
    }

    // for walls we need to zero the electro_mass_flux since we have already zero'd the diff and stoch
    // mass fluxes.  For inhomogeneous Neumann conditions on Epot, the physically correct thing would
    // have been to compute species gradients to exactly counterbalance the Neumann conditions on Epot
    // so the total mass flux (diff + stoch + Epot) is zero, but the numerical remedy here is to simply
    // zero them individually.
    ZeroEdgevalWalls(electro_mass_flux, geom, 0, nspecies);
}


void LimitEMF(const MultiFab& rho_in,
              std::array< MultiFab, AMREX_SPACEDIM >& electro_mass_flux) {


    for ( MFIter mfi(rho_in,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        const Box& bx = mfi.tilebox();

        const Array4<Real const>& rho = rho_in.array(mfi);

        AMREX_D_TERM(const Array4<Real>& emfx = electro_mass_flux[0].array(mfi);,
                     const Array4<Real>& emfy = electro_mass_flux[1].array(mfi);,
                     const Array4<Real>& emfz = electro_mass_flux[2].array(mfi););

        amrex::ParallelFor(bx, nspecies,  [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            if (rho(i,j,k,n) <= 0.) {

                if (emfx(i,j,k,n) > 0.) {
                    for (int m=0; m<nspecies; ++m) {
                        emfx(i,j,k,m) = 0.;
                    }
                }
                if (emfx(i+1,j,k,n) < 0.) {
                    for (int m=0; m<nspecies; ++m) {
                        emfx(i+1,j,k,m) = 0.;
                    }
                }

                if (emfy(i,j,k,n) > 0.) {
                    for (int m=0; m<nspecies; ++m) {
                        emfy(i,j,k,m) = 0.;
                    }
                }
                if (emfy(i,j+1,k,n) < 0.) {
                    for (int m=0; m<nspecies; ++m) {
                        emfy(i,j+1,k,m) = 0.;
                    }
                }
#if (AMREX_SPACEDIM == 3)
                if (emfz(i,j,k,n) > 0.) {
                    for (int m=0; m<nspecies; ++m) {
                        emfz(i,j,k,m) = 0.;
                    }
                }
                if (emfz(i,j,k+1,n) < 0.) {
                    for (int m=0; m<nspecies; ++m) {
                        emfz(i,j,k+1,m) = 0.;
                    }
                }
#endif
            }
        });
    } // end MFIter
}
