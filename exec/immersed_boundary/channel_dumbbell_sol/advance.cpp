
#include "main_driver.H"

#include "hydro_functions.H"
#include "hydro_functions_F.H"

#include "common_functions.H"
#include "common_functions_F.H"

#include "common_namespace.H"

#include "gmres_functions.H"
#include "gmres_functions_F.H"

#include "gmres_namespace.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MultiFabUtil.H>
#include <IBMarkerContainer.H>
#include <IBMarkerMD.H>


using namespace amrex;
using namespace common;
using namespace gmres;
using namespace immbdy_md;


// argv contains the name of the inputs file entered at the command line
void advance(std::array< MultiFab, AMREX_SPACEDIM >& umac,
             std::array< MultiFab, AMREX_SPACEDIM >& umacNew,
             MultiFab& pres, MultiFab& tracer,
             IBMarkerContainer & ib_mc,
             const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_predict,
             const std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv_correct,
             const std::array< MultiFab, AMREX_SPACEDIM >& alpha_fc,
             const MultiFab& beta, const MultiFab& gamma,
             const std::array< MultiFab, NUM_EDGE >& beta_ed,
             const Geometry geom, const Real& dt)
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

    for (int i=0; i<AMREX_SPACEDIM; i++) {
           gmres_rhs_u[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                 Lumac[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
            advFluxdiv[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
        advFluxdivPred[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
                  uMom[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);
    }

    // Tracer concentration field for predictor
    MultiFab tracerPred(ba, dmap, 1, 1);

    // Tracer advective terms
    MultiFab advFluxdivS(ba, dmap, 1, 1);


    //___________________________________________________________________________
    // Scaled alpha, beta, gamma:

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
     * Apply non-stochastic boundary conditions                                 *
     *                                                                          *
     ***************************************************************************/

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        umac[i].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(umac[i], i, geom, i);
        MultiFABPhysBCMacVel(umac[i], i, geom, i);
    }


    /****************************************************************************
     *                                                                          *
     * Advance tracer                                                           *
     *                                                                          *
     ***************************************************************************/

    // Compute tracer:
    tracer.FillBoundary(geom.periodicity());
    MultiFABPhysBC(tracer, geom);

    MkAdvSFluxdiv(umac, tracer, advFluxdivS, dx, geom, 0);
    advFluxdivS.mult(dt, 1);

    // compute predictor
    MultiFab::Copy(tracerPred, tracer, 0, 0, 1, 0);
    MultiFab::Add(tracerPred, advFluxdivS, 0, 0, 1, 0);

    tracerPred.FillBoundary(geom.periodicity());
    MultiFABPhysBC(tracerPred, geom);

    MkAdvSFluxdiv(umac, tracerPred, advFluxdivS, dx, geom, 0);
    advFluxdivS.mult(dt, 1);

    // advance in time
    MultiFab::Add(tracer, tracerPred,  0, 0, 1, 0);
    MultiFab::Add(tracer, advFluxdivS, 0, 0, 1, 0);
    tracer.mult(0.5, 1);

    // amrex::Print() << "tracer L0 norm = " << tracer.norm0() << "\n";
    //////////////////////////

    //////////////////////////////////////////////////
    // ADVANCE velocity field
    //////////////////////////////////////////////////

    // PREDICTOR STEP (heun's method: part 1)
    // compute advective term

    for (int i=0; i<AMREX_SPACEDIM; i++)
        MultiFab::Copy(uMom[i], umac[i], 0, 0, 1, 1);

    // let rho = 1
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        uMom[d].mult(1.0, 1);
    }

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        uMom[i].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(uMom[i], i, geom, i);
        MultiFABPhysBCMacVel(uMom[i], i, geom, i);
    }

    MkAdvMFluxdiv(umac, uMom, advFluxdiv, dx, 0);

    // crank-nicolson terms
    StagApplyOp(beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac, alpha_fc_0, dx, theta_alpha);


    //___________________________________________________________________________
    // Interpolate immersed boundary predictor
    ib_mc.ResetPredictor(0);
    ib_mc.InterpolatePredictor(0, umac);

    //___________________________________________________________________________
    // Move markers according to predictor velocity
    ib_mc.MovePredictor(0, dt);


    //
    //
    // TODO: Update forces between markers
    //
    //
    ib_mc.clearNeighbors();
    ib_mc.fillNeighbors(); // Does ghost cells

    ib_mc.buildNeighborList(ib_mc.CheckPair);

    int ib_lev = 0;

    // Parameters for spring force calculation
    Real spr_k = 10.0 ; // spring constant
    Real init_dx = 0.01, init_dy = 0., init_dz = 0.; //initial distance btw markers. Need to update depending on initial coordinates.

    // Parameters for calling bending force calculation
    Real bend_k = 10.0; //bending stiffness
    Real cos_theta0 = 1.0; //initial cos_theta value

    RealVect f, f_p, f_m;
    RealVect r, r_p, r_m;


    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = ib_mc.GetParticles(ib_lev).at(index);
        long np = particle_data.size();

        AoS & particles = particle_data.GetArrayOfStructs();
        ParticleType * nbhd_data = ib_mc.get_neighbors(ib_lev, index);
        const Vector<int> & nbhd = ib_mc.get_neighbor_list(ib_lev, index);
        int nbhd_index = 0;


        //// Set all forces to zero. Get ready for updating ////////
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];
            part.rdata(IBM_realData::pred_forcex) = 0.;
            part.rdata(IBM_realData::pred_forcey) = 0.;
            part.rdata(IBM_realData::pred_forcez) = 0.;
        }

        // Set bending forces to zero
        f_p = RealVect{0., 0., 0.};
        f   = RealVect{0., 0., 0.};
        f_m = RealVect{0., 0., 0.};
 

        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];

            // position vector for the current particle
            r = RealVect{part.rdata(IBM_realData::pred_posx), part.rdata(IBM_realData::pred_posy), part.rdata(IBM_realData::pred_posz)};

            // std::cout << "my id = " << part.id() << std::endl;

            int nn = nbhd[nbhd_index];
            nbhd_index ++; // pointing at first neighbor

            // assuming there is no previous/minus neighbor before searching
            bool found_npart_m = false;

            // Loops over neighbor list
            for (int j=0; j < nn; ++j){
                int ni = nbhd[nbhd_index] - 1; // -1 <= neighbor list uses Fortran indexing
                // std::cout << "neighbor list index " << ni << std::endl;

                ParticleType * npart; 
                ParticleType * npart_m; // for pointing to the minus particle if found

                if (ni >= np) {
                    ni = ni - np;
                    npart = & nbhd_data[ni];
                    // std::cout << "neighbor = " << npart->pos(1) << std::endl;
                    // std::cout << "neighbor id = " << npart->id() << std::endl;
                } else {
                    npart = & particles[ni];
                    // std::cout << "particle = " << npart->ipos(1) << std::endl;
                    // std::cout << "particle id = " << npart->id() << std::endl;
                }

                //check if the neighbor candidate is the previous/minus neighbor created.
                //If so, compute and update spring forces for both current particle and its previous neighbor
                if ((npart->id()==part.idata(IBM_intData::id_0)) && (npart->cpu()==part.idata(IBM_intData::cpu_0))) {

                    // add on differential changes in spring forces 
                    part.rdata(IBM_realData::pred_forcex) -= spr_k * ((part.rdata(IBM_realData::pred_posx)-npart->rdata(IBM_realData::pred_posx) - init_dx));
                    part.rdata(IBM_realData::pred_forcey) -= spr_k * ((part.rdata(IBM_realData::pred_posy)-npart->rdata(IBM_realData::pred_posy) - init_dy));
                    part.rdata(IBM_realData::pred_forcez) -= spr_k * ((part.rdata(IBM_realData::pred_posz)-npart->rdata(IBM_realData::pred_posz) - init_dz));

                    // action and reaction. add same forces back to the neighbor particle.
                    npart->rdata(IBM_realData::pred_forcex) += spr_k * ((part.rdata(IBM_realData::pred_posx)-npart->rdata(IBM_realData::pred_posx) - init_dx));
                    npart->rdata(IBM_realData::pred_forcey) += spr_k * ((part.rdata(IBM_realData::pred_posy)-npart->rdata(IBM_realData::pred_posy) - init_dy));
                    npart->rdata(IBM_realData::pred_forcez) += spr_k * ((part.rdata(IBM_realData::pred_posz)-npart->rdata(IBM_realData::pred_posz) - init_dz));

                    // position vector for the previous (minus) neirghbor particle
                    npart_m = npart;
                    r_m = RealVect{npart->rdata(IBM_realData::pred_posx),npart->rdata(IBM_realData::pred_posy), npart->rdata(IBM_realData::pred_posz)};                     
                    found_npart_m = true;
                }

                // check if the neighbor candidate is the next (plus) neighbor created.
                // If so, call bending force calculation and update forces on the previous/minus, current, and next/plus particles. 
                if (npart->idata(IBM_intData::id_0) == part.id() && npart->idata(IBM_intData::cpu_0) == part.cpu()) {
                    // position vector for the next/plus neighbor particle
                    //ParticleType * npart_p = npart; //npart is the plus particle
                    r_p = RealVect{npart->rdata(IBM_realData::pred_posx),npart->rdata(IBM_realData::pred_posy), npart->rdata(IBM_realData::pred_posz)};
                    
                    // make sure the current particle has a minus neighbor as well
                    if (found_npart_m) {
                    //calling the bending force calculation
                    bending_f(f, f_p, f_m, r, r_p, r_m, bend_k, cos_theta0);
                    
                    Print() << "f= " << f_p << std::endl;
                    Print() << "f= " << f << std::endl;
                    Print() << "f= " << f_m << std::endl;
                    
                    // updating the force on the minus, current, and plus particles. 
                    part.rdata(IBM_realData::pred_forcex) += f[0]; 
                    part.rdata(IBM_realData::pred_forcey) += f[1];   
                    part.rdata(IBM_realData::pred_forcez) += f[2];              

                    npart_m->rdata(IBM_realData::pred_forcex) += f_m[0];
                    npart_m->rdata(IBM_realData::pred_forcey) += f_m[1];
                    npart_m->rdata(IBM_realData::pred_forcez) += f_m[2];

                    npart->rdata(IBM_realData::pred_forcex) += f_p[0];
                    npart->rdata(IBM_realData::pred_forcey) += f_p[1];
                    npart->rdata(IBM_realData::pred_forcez) += f_p[2];
                    }
                }
                nbhd_index ++;
            }
            
            
         }
    }


Abort();




    //___________________________________________________________________________
    // Spread forces to predictor
    std::array<MultiFab, AMREX_SPACEDIM> fc_force_pred;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        fc_force_pred[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        fc_force_pred[d].setVal(0.);
    }

    // Spread to the `fc_force` multifab
    ib_mc.fillNeighbors(); // Don't forget to fill neighbor particles
    ib_mc.SpreadMarkers(0, fc_force_pred);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        fc_force_pred[d].FillBoundary(geom.periodicity());


    for (int d=0; d<AMREX_SPACEDIM; d++) {
        Lumac[d].FillBoundary(geom.periodicity());

        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);

        gmres_rhs_u[d].mult(dtinv, 1);
        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_predict[d], 0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], fc_force_pred[d],    0, 0, 1, 0);
    }

    std::array< MultiFab, AMREX_SPACEDIM > pg;
    for (int i=0; i<AMREX_SPACEDIM; i++)
        pg[i].define(convert(ba, nodal_flag_dir[i]), dmap, 1, 1);

    pres.setVal(0.);  // initial guess
    SetPressureBC(pres, geom);

    ComputeGrad(pres, pg, 0, 0, 1, geom);

    for (int i=0; i<AMREX_SPACEDIM; i++) {
        pg[i].FillBoundary(geom.periodicity());
        gmres_rhs_u[i].FillBoundary(geom.periodicity());

        MultiFab::Subtract(gmres_rhs_u[i], pg[i], 0, 0, 1, 1);
    }

    // initial guess for new solution
    for (int i=0; i<AMREX_SPACEDIM; i++)
        MultiFab::Copy(umacNew[i], umac[i], 0, 0, 1, 1);

    // call GMRES to compute predictor
    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, pres,
          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
          geom, norm_pre_rhs);

    // Compute predictor advective term
    // let rho = 1
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        umacNew[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(umacNew[d], d, geom, d);
        MultiFABPhysBCMacVel(umacNew[d], d, geom, d);

        MultiFab::Copy(uMom[d], umacNew[d], 0, 0, 1, 0);
        uMom[d].mult(1.0, 1);

        uMom[d].FillBoundary(geom.periodicity());
        MultiFABPhysBCDomainVel(uMom[d], d, geom, d);
        MultiFABPhysBCMacVel(uMom[d], d, geom, d);
    }

    MkAdvMFluxdiv(umacNew,uMom,advFluxdivPred,dx,0);

    // ADVANCE STEP (crank-nicolson + heun's method)

    // Compute gmres_rhs

    // trapezoidal advective terms
    for (int d=0; d<AMREX_SPACEDIM; d++) {
        advFluxdiv[d].mult(0.5, 1);
        advFluxdivPred[d].mult(0.5, 1);
    }

    // crank-nicolson terms
    StagApplyOp(beta_negwtd, gamma_negwtd, beta_ed_negwtd, umac, Lumac, alpha_fc_0, dx, theta_alpha);



    //___________________________________________________________________________
    // Interpolate immersed boundary
    ib_mc.ResetMarkers(0);
    ib_mc.InterpolateMarkers(0, umacNew);

    //___________________________________________________________________________
    // Move markers according to velocity
    ib_mc.MoveMarkers(0, dt);
    ib_mc.Redistribute(); // Don't forget to send particles to the right CPU


    //
    //
    // TODO: Update forces between markers
    //
    //
    ib_mc.clearNeighbors();
    ib_mc.fillNeighbors(); // Does ghost cells

    ib_mc.buildNeighborList(ib_mc.CheckPair);

    for (IBMarIter pti(ib_mc, ib_lev); pti.isValid(); ++pti) {

        PairIndex index(pti.index(), pti.LocalTileIndex());
        auto & particle_data = ib_mc.GetParticles(ib_lev).at(index);
        long np = particle_data.size();

        AoS & particles = particle_data.GetArrayOfStructs();
        ParticleType * nbhd_data = ib_mc.get_neighbors(ib_lev, index);

        const Vector<int> & nbhd = ib_mc.get_neighbor_list(ib_lev, index);
        int nbhd_index = 0;


        //// Set all forces to zero. Get ready for updating ////////
        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];
            part.rdata(IBM_realData::forcex) = 0.;
            part.rdata(IBM_realData::forcey) = 0.;
            part.rdata(IBM_realData::forcez) = 0.;
        }

        for (int i = 0; i < np; ++i) {
            ParticleType & part = particles[i];

            // std::cout << "my id = " << part.id() << std::endl;

            int nn = nbhd[nbhd_index];
            nbhd_index ++; // pointing at first neighbor

            // Loops over neighbor list
            for (int j=0; j < nn; ++j){
                int ni = nbhd[nbhd_index] - 1; // -1 <= neighbor list uses Fortran indexing
                // std::cout << "neighbor list index " << ni << std::endl;

                ParticleType * npart;

                if (ni >= np) {
                    ni = ni - np;
                    npart = & nbhd_data[ni];
                    // std::cout << "neighbor = " << npart->pos(1) << std::endl;
                    // std::cout << "neighbor id = " << npart->id() << std::endl;
                } else {
                    npart = & particles[ni];
                    // std::cout << "particle = " << npart->pos(1) << std::endl;
                    // std::cout << "particle id = " << npart->id() << std::endl;
                }

                //check if the neighbor candidate(s) is the actual neighbor to apply force.
                //If so, compute and update forces for both current particle and its interacting neighbor
                if ((npart->id()==part.idata(IBM_intData::id_0)) && (npart->cpu()==part.idata(IBM_intData::cpu_0))) {

                    // add on differential changes in forces using velocities
                    // part.rdata(IBM_realData::pred_forcex) += spr_k * dt * (part.rdata(IBM_realData::pred_velx)-npart->rdata(IBM_realData::pred_velx));
                    // part.rdata(IBM_realData::pred_forcey) += spr_k * dt * (part.rdata(IBM_realData::pred_vely)-npart->rdata(IBM_realData::pred_vely));
                    // part.rdata(IBM_realData::pred_forcez) += spr_k * dt * (part.rdata(IBM_realData::pred_velz)-npart->rdata(IBM_realData::pred_velz));

                    // npart->rdata(IBM_realData::pred_forcex) += spr_k * dt * (part.rdata(IBM_realData::pred_velx)-npart->rdata(IBM_realData::pred_velx));
                    // npart->rdata(IBM_realData::pred_forcey) += spr_k * dt * (part.rdata(IBM_realData::pred_vely)-npart->rdata(IBM_realData::pred_vely));
                    // npart->rdata(IBM_realData::pred_forcez) += spr_k * dt * (part.rdata(IBM_realData::pred_velz)-npart->rdata(IBM_realData::pred_velz));

                    // Alternatively one can use position changes for calculating forces, but need to set all forces to zero earlier before updating 
                    part.rdata(IBM_realData::forcex) -= spr_k * ((part.pos(0)-npart->pos(0) - init_dx));
                    part.rdata(IBM_realData::forcey) -= spr_k * ((part.pos(1)-npart->pos(1) - init_dy));
                    part.rdata(IBM_realData::forcez) -= spr_k * ((part.pos(2)-npart->pos(2) - init_dz));

                    // action and reaction. add same forces back to the neighbor particle.
                    npart->rdata(IBM_realData::forcex) += spr_k * ((part.pos(0)-npart->pos(0) - init_dx));
                    npart->rdata(IBM_realData::forcey) += spr_k * ((part.pos(1)-npart->pos(1) - init_dy));
                    npart->rdata(IBM_realData::forcez) += spr_k * ((part.pos(2)-npart->pos(2) - init_dz));
                }

                nbhd_index ++;
            }

         }
    }



    //___________________________________________________________________________
    // Spread forces to corrector
    std::array<MultiFab, AMREX_SPACEDIM> fc_force_corr;
    for (int d=0; d<AMREX_SPACEDIM; ++d){
        fc_force_corr[d].define(convert(ba, nodal_flag_dir[d]), dmap, 1, 1);
        fc_force_corr[d].setVal(0.);
    }

    // Spread to the `fc_force` multifab
    ib_mc.fillNeighbors(); // Don't forget to fill neighbor particles
    ib_mc.SpreadMarkers(0, fc_force_corr);
    for (int d=0; d<AMREX_SPACEDIM; ++d)
        fc_force_corr[d].FillBoundary(geom.periodicity());



    for (int d=0; d<AMREX_SPACEDIM; d++) {
        MultiFab::Copy(gmres_rhs_u[d], umac[d], 0, 0, 1, 1);

        gmres_rhs_u[d].mult(dtinv, 1);
        MultiFab::Add(gmres_rhs_u[d], mfluxdiv_correct[d], 0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], Lumac[d],            0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], advFluxdiv[d],       0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], advFluxdivPred[d],   0, 0, 1, 0);
        MultiFab::Add(gmres_rhs_u[d], fc_force_corr[d],    0, 0, 1, 0);

        gmres_rhs_u[d].FillBoundary(geom.periodicity());

        MultiFab::Subtract(gmres_rhs_u[d], pg[d], 0, 0, 1, 1);

        MultiFab::Copy(umacNew[d], umac[d], 0, 0, 1, 0);
    }

    pres.setVal(0.);  // initial guess
    SetPressureBC(pres, geom);

    // call GMRES here
    GMRES(gmres_rhs_u, gmres_rhs_p, umacNew, pres,
          alpha_fc, beta_wtd, beta_ed_wtd, gamma_wtd, theta_alpha,
          geom, norm_pre_rhs);

    for (int i=0; i<AMREX_SPACEDIM; i++)
        MultiFab::Copy(umac[i], umacNew[i], 0, 0, 1, 0);

    //////////////////////////////////////////////////

}
