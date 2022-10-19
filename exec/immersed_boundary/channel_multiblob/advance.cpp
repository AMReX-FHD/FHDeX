#include "main_driver.H"

#include "hydro_functions.H"

#include "common_functions.H"


#include "gmres_functions.H"


#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

#include <IBMultiBlobContainer.H>
#include <IBMarkerMD.H>


#include <immbdy_namespace.H>

using namespace amrex;
using namespace immbdy_md;

using namespace ib_colloid;


using ParticleVector = typename IBMarkerContainer::ParticleVector;


// argv contains the name of the inputs file entered at the command line
void advance(std::array< MultiFab, AMREX_SPACEDIM >& umac,
             std::array< MultiFab, AMREX_SPACEDIM >& umacNew,
             MultiFab & pres, MultiFab & tracer,
             IBMultiBlobContainer & ib_mbc,
             const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_predict,
             const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_correct,
                   std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
             const MultiFab & beta, const MultiFab & gamma,
             const std::array< MultiFab, NUM_EDGE >& beta_ed,
             const Geometry geom, const Real & dt)
{

    BL_PROFILE_VAR("advance()",advance);

    const Real * dx  = geom.CellSize();
    const Real dtinv = 1.0/dt;
    Real theta_alpha = 1.;
    Real norm_pre_rhs;

    const BoxArray &              ba = beta.boxArray();
    const DistributionMapping & dmap = beta.DistributionMap();


    /****************************************************************************
     *                                                                          *
     * Define temporary MultiFabs                                               *
     *                                                                          *
     ***************************************************************************/

    // RHS pressure in GMRES
    MultiFab gmres_rhs_p(ba, dmap, 1, 1);
    gmres_rhs_p.setVal(0.);

    // RHS velocities in GMRES
    std::array< MultiFab, AMREX_SPACEDIM > gmres_rhs_u;
    // Velocity components updated by diffusion operator
    std::array< MultiFab, AMREX_SPACEDIM > Lumac;
    // Advective terms
    std::array< MultiFab, AMREX_SPACEDIM > advFluxdiv;
    // Advective terms (for predictor)
    std::array< MultiFab, AMREX_SPACEDIM > advFluxdivPred;
    // Staggered momentum
    std::array< MultiFab, AMREX_SPACEDIM > uMom;
    // Pressure gradient at inflow/outflow
    std::array< MultiFab, AMREX_SPACEDIM > pg;

    for (int i=0; i<AMREX_SPACEDIM; i++) {
           gmres_rhs_u[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                 Lumac[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 0);
            advFluxdiv[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
        advFluxdivPred[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                  uMom[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                    pg[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);

        // put in to fix fpe traps
        advFluxdivPred[i].setVal(0);
        advFluxdiv[i].setVal(0);
    }

    // Tracer concentration field for predictor
    MultiFab tracerPred(ba, dmap, 1, 1);
    tracerPred.setVal(0.);

    // Tracer advective terms
    MultiFab advFluxdivS(ba, dmap, 1, 1);
    advFluxdivS.setVal(0.);


    //___________________________________________________________________________
    // Scaled alpha, beta, gamma:
    //  * scaled by 1/2 => corrector/ crank-nicolson terms
    //  * scaled by -1  => move implicit term (LHS) to explicit term (RHS)

    // alpha_fc_0 arrays
    std::array< MultiFab, AMREX_SPACEDIM > alpha_fc_0;

    for (int i=0; i<AMREX_SPACEDIM; i++){
        alpha_fc_0[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
        alpha_fc_0[i].setVal(0.);
    }

    // Scaled by 1/2:
    // beta_wtd cell centered
    MultiFab beta_wtd(ba, dmap, 1, 1);
    MultiFab::Copy(beta_wtd, beta, 0, 0, 1, 1);
    beta_wtd.mult(0.5, 1);

    // beta_wtd on nodes in 2d, on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed_wtd;
#if (AMREX_SPACEDIM == 2)
    beta_ed_wtd[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    MultiFab::Copy(beta_ed_wtd[0], beta_ed[0], 0, 0, 1, 1);
    beta_ed_wtd[0].mult(0.5, 1);
#elif (AMREX_SPACEDIM == 3)
    for(int d=0; d<AMREX_SPACEDIM; d++) {
        beta_ed_wtd[d].define(convert(ba, nodal_flag_edge[d]), dmap, 1, 1);

        MultiFab::Copy(beta_ed_wtd[d], beta_ed[d], 0, 0, 1, 1);
        beta_ed_wtd[d].mult(0.5, 1);
    }
#endif

    // Scaled by 1/2:
    // gamma_wtd cell centered
    MultiFab gamma_wtd(ba, dmap, 1, 1);
    MultiFab::Copy(gamma_wtd, gamma, 0, 0, 1, 1);
    gamma_wtd.mult(-0.5, 1);

    // Scaled by -1/2:
    // beta_negwtd cell centered
    MultiFab beta_negwtd(ba, dmap, 1, 1);
    MultiFab::Copy(beta_negwtd, beta, 0, 0, 1, 1);
    beta_negwtd.mult(-0.5, 1);

    // beta_negwtd on nodes in 2d, on edges in 3d
    std::array< MultiFab, NUM_EDGE > beta_ed_negwtd;
#if (AMREX_SPACEDIM == 2)
    beta_ed_negwtd[0].define(convert(ba,nodal_flag), dmap, 1, 1);
    MultiFab::Copy(beta_ed_negwtd[0], beta_ed[0], 0, 0, 1, 1);
    beta_ed_negwtd[0].mult(-0.5, 1);
#elif (AMREX_SPACEDIM == 3)
    for(int d=0; d<AMREX_SPACEDIM; d++) {
        beta_ed_negwtd[d].define(convert(ba,nodal_flag_edge[d]), dmap, 1, 1);

        MultiFab::Copy(beta_ed_negwtd[d], beta_ed[d], 0, 0, 1, 1);
        beta_ed_negwtd[d].mult(-0.5, 1);
    }
#endif

    // Scaled by -1/2:
    // gamma cell centered
    MultiFab gamma_negwtd(ba, dmap, 1, 1);
    MultiFab::Copy(gamma_negwtd, gamma, 0, 0, 1, 1);
    gamma_negwtd.mult(-0.5, 1);



    /****************************************************************************
     *                                                                          *
     * Immersed-Marker parameters                                               *
     *                                                                          *
     ***************************************************************************/

    int ib_lev = 0;


    /****************************************************************************
     *                                                                          *
     * Apply deterministic boundary conditions                                  *
     *                                                                          *
     ***************************************************************************/

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        umac[i].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umac[i], geom, i);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umac[i], geom, i, is_inhomogeneous);
    }


    /****************************************************************************
     *                                                                          *
     * Advance tracer                                                           *
     * NOTE: this is a centered scheme => can be unstable of |Pe| > 2           *
     *                                                                          *
     ***************************************************************************/

    // Compute tracer:
    tracer.FillBoundary(geom.periodicity());
    MultiFabPhysBC(tracer, geom, 0, 1, SPEC_BC_COMP);

    MkAdvSFluxdiv_cc(umac, tracer, advFluxdivS, geom, 0, 1, 0);
    advFluxdivS.mult(dt, 1);

    // compute predictor
    MultiFab::Copy(tracerPred, tracer, 0, 0, 1, 0);
    MultiFab::Add(tracerPred, advFluxdivS, 0, 0, 1, 0);

    tracerPred.FillBoundary(geom.periodicity());
    MultiFabPhysBC(tracerPred, geom, 0, 1, SPEC_BC_COMP);

    MkAdvSFluxdiv_cc(umac, tracerPred, advFluxdivS, geom, 0, 1, 0);
    advFluxdivS.mult(dt, 1);

    // advance in time
    MultiFab::Add(tracer, tracerPred,  0, 0, 1, 0);
    MultiFab::Add(tracer, advFluxdivS, 0, 0, 1, 0);
    tracer.mult(0.5, 1);



    /****************************************************************************
     *                                                                          *
     * PREDICTOR (Heun's method) step:                                          *
     *                                                                          *
     ***************************************************************************/


    //___________________________________________________________________________
    // Interpolate immersed boundary predictor J(u^n)

    std::array<MultiFab, AMREX_SPACEDIM> umac_buffer;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        umac_buffer[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        umac_buffer[d].setVal(0.);
        MultiFab::Copy(umac_buffer[d], umac[d], 0, 0, 1, umac[d].nGrow());
        umac_buffer[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umac[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umac[d], geom, d, is_inhomogeneous);
    }

    ib_mbc.ResetPredictor(ib_lev);
    ib_mbc.InterpolatePredictor(ib_lev, umac_buffer);


    //___________________________________________________________________________
    // Move markers according to predictor velocity

    // MovePredictor needs the neighbor list to be filled
    ib_mbc.clearNeighbors();
    ib_mbc.fillNeighbors();

    ib_mbc.MovePredictor(ib_lev, dt);


    //___________________________________________________________________________
    // Calculate internal force model

    ib_mbc.PredictorForces(ib_lev);


    // //___________________________________________________________________________
    // // Update forces between markers
    // ib_mbc.clearNeighbors();
    // ib_mbc.fillNeighbors(); // Does ghost cells

    // ib_mbc.buildNeighborList(ib_mbc.CheckPair);


    // for (IBMBIter pti(ib_mbc, ib_lev); pti.isValid(); ++pti) {

    //     // Get marker data (local to current thread)
    //     TileIndex index(pti.index(), pti.LocalTileIndex());
    //     AoS & markers = ib_mbc.GetParticles(ib_lev).at(index).GetArrayOfStructs();
    //     long np = ib_mbc.GetParticles(ib_lev).at(index).numParticles();

    //     // m_index.second is used to keep track of the neighbor list
    //     for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

    //         ParticleType & mark = markers[m_index.first];
    //     }
    // }


    //___________________________________________________________________________
    // Spread forces to predictor
    // Remember: Spread, Fold, Fill, Sum

    std::array<MultiFab, AMREX_SPACEDIM> fc_force_pred;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        fc_force_pred[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        fc_force_pred[d].setVal(0.);
    }


    // Spread predictor forces
    ib_mbc.SpreadPredictor(ib_lev, fc_force_pred);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        fc_force_pred[d].SumBoundary(geom.periodicity());


    //___________________________________________________________________________
    // Compute predictor velocity field

    // Compute momentum fluxes: uMom = rho*u^n
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        MultiFab::Copy(uMom[d], umac[d], 0, 0, 1, 1);

        // let rho = 1
        uMom[d].mult(1.0, 1);

        uMom[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(uMom[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(uMom[d], geom, d, is_inhomogeneous);
    }

    // Compute advective fluxes: advFluxdiv = - D(\rho uu^n) = - D(u^n uMom)
    MkAdvMFluxdiv(umac, uMom, advFluxdiv, dx, 0);

    // (crank-nicolson terms) Explicit part of the diffusive operator Lu^n/2.
    // Note that we are using the weighted coefficients (to deal with the 1/2
    // part)
    StagApplyOp(geom, beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac,
                alpha_fc_0, dx, theta_alpha);


    //___________________________________________________________________________
    // Compute pressure gradient due to the BC: gp = Gp
    pres.setVal(0.); // Initial guess for pressure
    MultiFabPhysBC(pres, geom, 0, 1, PRES_BC_COMP, 0); // Apply pressure boundary conditions
    for (int d=0; d<AMREX_SPACEDIM; ++d) pg[d].setVal(0);
    ComputeGrad(pres, pg, 0, 0, 1, PRES_BC_COMP, geom);

    // Construct RHS of Navier Stokes Equation
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        Lumac[d].FillBoundary(geom.periodicity());

        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);

        gmres_rhs_u[d].mult(dtinv, 1);
        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_predict[d], 0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], fc_force_pred[d],    0, 0, 1, 0);

        // FillBoundary before adding boundary conditions to prevent
        // overwriting (which are defined on ghost cells)
        pg[d].FillBoundary(geom.periodicity());
        gmres_rhs_u[d].FillBoundary(geom.periodicity());

        MultiFab::Subtract(gmres_rhs_u[d], pg[d], 0, 0, 1, 1);
    }

    // initial guess for new solution
    for (int i=0; i<AMREX_SPACEDIM; i++)
        MultiFab::Copy(umacNew[i], umac[i], 0, 0, 1, 1);

    // Call GMRES to compute u^(n+1/2). Lu^(n+1/2) is computed implicitly. Note
    // that we are using the un-weighted coefficients.
    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, pres,
          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
          geom, norm_pre_rhs);

    // Apply boundary conditions to the solution
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umacNew[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umacNew[d], geom, d, is_inhomogeneous);
    }



    /****************************************************************************
     *                                                                          *
     * ADVANCE (Crank-Nicolson + Heun's method) step:                           *
     *                                                                          *
     ***************************************************************************/


    //___________________________________________________________________________
    // Interpolate immersed boundary J(u^*)

    std::array<MultiFab, AMREX_SPACEDIM> umacNew_buffer;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        umacNew_buffer[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        umacNew_buffer[d].setVal(0.);
        MultiFab::Copy(umacNew_buffer[d], umacNew[d], 0, 0, 1, umacNew[d].nGrow());
        umacNew_buffer[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(umacNew[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(umacNew[d], geom, d, is_inhomogeneous);
    }

    ib_mbc.ResetMarkers(ib_lev);
    ib_mbc.InterpolateMarkers(ib_lev, umacNew_buffer);


    //___________________________________________________________________________
    // Move markers according to velocity

    // MoveMarkers needs the neighbor list to be filled (TODO: this might be
    // redundant work, as we're already filling neighbors in the predictor
    // step)
    ib_mbc.clearNeighbors();
    ib_mbc.fillNeighbors();

    ib_mbc.MoveMarkers(ib_lev, dt);
    ib_mbc.RedistributeMarkers(); // Don't forget to send particles to the right CPU


    //___________________________________________________________________________
    // Calculate internal force model

    ib_mbc.MarkerForces(ib_lev);

    //___________________________________________________________________________
    // Move Multi-Blobs

    // Accumulate hydrodynamic drag forces
    ib_mbc.ResetDrag(ib_lev);
    ib_mbc.AccumulateDrag(ib_lev);
    ib_mbc.sumNeighbors(IBMBReal::dragx, AMREX_SPACEDIM, 0, 0);

    for (IBMBIter pti(ib_mbc, ib_lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        IBMultiBlobContainer::TileIndex index(pti.index(), pti.LocalTileIndex());
        IBMultiBlobContainer::AoS & markers =
            ib_mbc.GetParticles(ib_lev).at(index).GetArrayOfStructs();
        long np = ib_mbc.GetParticles(ib_lev).at(index).numParticles();

        for (int i=0; i<np; ++i) {

            ParticleType & mark = markers[i];
            for (int d=0; d<AMREX_SPACEDIM; ++d)
                std::cout << "drag_" << d << " = "
                          << mark.rdata(IBMBReal::dragx + d)
                          << std::endl;
        }
    }


    // Move Multi-Blobs
    ib_mbc.MoveBlob(ib_lev, dt);
    ib_mbc.Redistribute(); // Don't forget to send particles to the right CPU

    // Diagnostics: print COM positions
    for (IBMBIter pti(ib_mbc, ib_lev); pti.isValid(); ++pti) {

        // Get marker data (local to current thread)
        IBMultiBlobContainer::TileIndex index(pti.index(), pti.LocalTileIndex());
        IBMultiBlobContainer::AoS & markers =
            ib_mbc.GetParticles(ib_lev).at(index).GetArrayOfStructs();
        long np = ib_mbc.GetParticles(ib_lev).at(index).numParticles();

        for (int i=0; i<np; ++i) {

            ParticleType & mark = markers[i];
            for (int d=0; d<AMREX_SPACEDIM; ++d)
                std::cout << "pos_" << d << " = "
                           << mark.pos(d)
                           << std::endl;
        }
    }


    // //___________________________________________________________________________
    // // Update forces between markers (these repeated clear/fill/build neighbor
    // // calls might be redundant)
    // ib_mbc.clearNeighbors();
    // ib_mbc.fillNeighbors(); // Does ghost cells

    // ib_mbc.buildNeighborList(ib_mbc.CheckPair);


    // for (IBMBIter pti(ib_mbc, ib_lev); pti.isValid(); ++pti) {

    //     // Get marker data (local to current thread)
    //     TileIndex index(pti.index(), pti.LocalTileIndex());
    //     AoS & markers = ib_mbc.GetParticles(ib_lev).at(index).GetArrayOfStructs();
    //     long np = ib_mbc.GetParticles(ib_lev).at(index).numParticles();

    //     // m_index.second is used to keep track of the neighbor list
    //     for (MarkerListIndex m_index(0, 0); m_index.first<np; ++m_index.first) {

    //         ParticleType & mark = markers[m_index.first];
    //     }
    // }


    //___________________________________________________________________________
    // Spread forces to corrector
    // Remember: Spread, Fold Fill, Sum

    std::array<MultiFab, AMREX_SPACEDIM> fc_force_corr;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        fc_force_corr[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 6);
        fc_force_corr[d].setVal(0.);
    }

    // Spread to the `fc_force` multifab
    ib_mbc.SpreadMarkers(ib_lev, fc_force_corr);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        fc_force_corr[d].SumBoundary(geom.periodicity());


    //___________________________________________________________________________
    // Compute corrector velocity field

    // Compute momentum fluxes for corrector: uMom = rho u^*
    // let rho = 1
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        MultiFab::Copy(uMom[d], umacNew[d], 0, 0, 1, 0);
        //
        // let rho = 1
        uMom[d].mult(1.0, 1);

        uMom[d].FillBoundary(geom.periodicity());
        MultiFabPhysBCDomainVel(uMom[d], geom, d);
        int is_inhomogeneous = 1;
        MultiFabPhysBCMacVel(uMom[d], geom, d, is_inhomogeneous);
    }

    // Compute advective fluxes for the corrector:
    // advFluxdivPred = - D(\rho uu^*) = - D(u^* uMom)
    MkAdvMFluxdiv(umacNew,uMom,advFluxdivPred,dx,0);

    // trapezoidal advective terms
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        advFluxdiv[d].mult(0.5, 1);
        advFluxdivPred[d].mult(0.5, 1);
    }


    //___________________________________________________________________________
    // Compute pressure, and pressure gradient due to the BC: gp = Gp
    // Note that the pressure gradient due to the BC is left unchanged
    pres.setVal(0.); // Initial guess for pressure
    MultiFabPhysBC(pres, geom, 0, 1, PRES_BC_COMP, 0); // Apply pressure boundary conditions

    // Construct RHS of Navier Stokes Equation
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        Lumac[d].FillBoundary(geom.periodicity());

        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);
        gmres_rhs_u[d].mult(dtinv, 1);  // advance by dt

        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_correct[d], 0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], advFluxdivPred[d],   0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], fc_force_corr[d],    0, 0, 1, 0);

        // FillBoundary before adding boundary conditions to prevent
        // overwriting (which are defined on ghost cells)
        gmres_rhs_u[d].FillBoundary(geom.periodicity());

        MultiFab::Subtract(gmres_rhs_u[d], pg[d], 0, 0, 1, 1);
    }

    // Initial guess for new solution
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 1);

    // Call GMRES to compute u^(n+1). Lu^(n+1)/2 is computed implicitly. Note
    // that we are using the weighted coefficients (to deal witht he 1/2 part)
    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, pres,
          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
          geom, norm_pre_rhs);

    for (int i=0; i<AMREX_SPACEDIM; i++)
        MultiFab::Copy(umac[i], umacNew[i], 0, 0, 1, 0);
}
