#include "rng_functions.H"
#include "gmres_functions.H"
#include "common_functions.H"
#include "hydro_functions_F.H"
//#include "analysis_functions_F.H"
//#include "StructFact_F.H"
#include "StochMFlux.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>


StochMFlux::StochMFlux(BoxArray ba_in, DistributionMapping dmap_in, Geometry geom_in,
		       int n_rngs_in) {

  BL_PROFILE_VAR("StochMFlux::StochMFlux()",StochMFlux);

  n_rngs = n_rngs_in;
  geom = geom_in;

  mflux_cc.resize(n_rngs);
  mflux_ed.resize(n_rngs);

  // Define mflux multifab vectors
  // TEMPORARY ASSUMPTION: filtering_width = 0
  for (int i=0; i<n_rngs; ++i) {
    mflux_cc[i].define(ba_in, dmap_in, AMREX_SPACEDIM, 1);
  }
#if (AMREX_SPACEDIM == 2)
  for (int i=0; i<n_rngs; ++i) {
    mflux_ed[i][0].define(convert(ba_in,nodal_flag), dmap_in, ncomp_ed, 0);
  }
#elif (AMREX_SPACEDIM == 3)
  for (int i=0; i<n_rngs; ++i) {
    mflux_ed[i][0].define(convert(ba_in,nodal_flag_xy), dmap_in, ncomp_ed, 0);
    mflux_ed[i][1].define(convert(ba_in,nodal_flag_xz), dmap_in, ncomp_ed, 0);
    mflux_ed[i][2].define(convert(ba_in,nodal_flag_yz), dmap_in, ncomp_ed, 0);
  }
#endif

  // Define weighted mflux multifab vectors
  mflux_cc_weighted.define(ba_in, dmap_in, AMREX_SPACEDIM, 1);
#if (AMREX_SPACEDIM == 2)
  mflux_ed_weighted[0].define(convert(ba_in,nodal_flag), dmap_in, ncomp_ed, 0);
#elif (AMREX_SPACEDIM == 3)
  mflux_ed_weighted[0].define(convert(ba_in,nodal_flag_xy), dmap_in, ncomp_ed, 0);
  mflux_ed_weighted[1].define(convert(ba_in,nodal_flag_xz), dmap_in, ncomp_ed, 0);
  mflux_ed_weighted[2].define(convert(ba_in,nodal_flag_yz), dmap_in, ncomp_ed, 0);
#endif
}


void StochMFlux::weightMflux(Vector< amrex::Real > weights) {
  mflux_cc_weighted.setVal(0.0);

  for (int i=0; i<n_rngs; ++i) {
    MultiFab::Saxpy(mflux_cc_weighted, weights[i], mflux_cc[i], 0, 0, AMREX_SPACEDIM, 1);
  }
  mflux_cc_weighted.FillBoundary(geom.periodicity());

  MultiFABPhysBC(mflux_cc_weighted, geom);

  for (int d=0; d<NUM_EDGE; ++d) {
    mflux_ed_weighted[d].setVal(0.0);
    for (int i=0; i<n_rngs; ++i) {
      MultiFab::Saxpy(mflux_ed_weighted[d], weights[i], mflux_ed[i][d], 0, 0, ncomp_ed, 0);
    }
    mflux_ed_weighted[d].FillBoundary(geom.periodicity());

    // TODO: is this the correct BC?
    MultiFABPhysBC(mflux_ed_weighted[d], d, geom);
  }
}

void StochMFlux::fillMStochastic() {

    BL_PROFILE_VAR("StochMFlux::StochMFlux()",StochMFlux);

    for (int i=0; i<n_rngs; ++i) {

        switch(stoch_stress_form) {

        case 0: // Non-symmetric
            // Print() << "Non-symmetric \n";
            for (int n=0; n<AMREX_SPACEDIM; ++n) {
                MultiFABFillRandom(mflux_cc[i],n,1.0,geom);
            }

            for (int d=0; d<NUM_EDGE; ++d) {
                for (int n=0; n<ncomp_ed; ++n) {
                    MultiFABFillRandom(mflux_ed[i][d],n,1.0,geom);
                }
            }
            break;

        default: // Symmetric
            // Print() << "Symmetric \n";
            for (int n=0; n<AMREX_SPACEDIM; ++n) {
                MultiFABFillRandom(mflux_cc[i],n,2.0,geom);
            }

            for (int d=0; d<NUM_EDGE; ++d) {
                MultiFABFillRandom(mflux_ed[i][d],0,1.0,geom);
                MultiFab::Copy(mflux_ed[i][d], mflux_ed[i][d], 0, 1, ncomp_ed-1, 0);
            }
            break;
        }
    }

    // TODO: Put stochastic BCs here ?
}

void StochMFlux::multbyVarSqrtEtaTemp(const MultiFab& eta_cc,
				      const std::array< MultiFab, NUM_EDGE >& eta_ed,
				      const MultiFab& temp_cc,
				      const std::array< MultiFab, NUM_EDGE >& temp_ed,
				      const amrex::Real& dt) {

  const Real* dx = geom.CellSize();

  Real dVol = dx[0]*dx[1];
  if (AMREX_SPACEDIM == 2) {
    dVol *= cell_depth;
  } else {
    if (AMREX_SPACEDIM == 3) {
      dVol *= dx[2];
    }
  }

  // Compute variance using computed differential volume
  Real variance = variance_coef_mom*sqrt(variance_coef_mom*2.0*k_B/(dVol*dt));
  //Real variance = variance_coef_mom*sqrt(variance_coef_mom*2.0*k_B/(dVol));

  // Scale mflux_weighted by variance
  mflux_cc_weighted.mult(variance, 1);
  for (int d=0; d<NUM_EDGE; d++) {
    mflux_ed_weighted[d].mult(variance, 0);
  }

  // Multiply mflux_weighted by sqrt(eta*temperature)
  // Loop over boxes
  for (MFIter mfi(mflux_cc_weighted); mfi.isValid(); ++mfi) {
    // Note: Make sure that multifab is cell-centered
    const Box& validBox = mfi.validbox();

    mult_by_sqrt_eta_temp(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
			  BL_TO_FORTRAN_FAB(mflux_cc_weighted[mfi]),
			  BL_TO_FORTRAN_FAB(mflux_ed_weighted[0][mfi]),
#if (AMREX_SPACEDIM == 3)
			  BL_TO_FORTRAN_FAB(mflux_ed_weighted[1][mfi]),
			  BL_TO_FORTRAN_FAB(mflux_ed_weighted[2][mfi]),
#endif
			  BL_TO_FORTRAN_ANYD(eta_cc[mfi]),
			  BL_TO_FORTRAN_ANYD(eta_ed[0][mfi]),
#if (AMREX_SPACEDIM == 3)
			  BL_TO_FORTRAN_ANYD(eta_ed[1][mfi]),
			  BL_TO_FORTRAN_ANYD(eta_ed[2][mfi]),
#endif
			  BL_TO_FORTRAN_ANYD(temp_cc[mfi]),
			  BL_TO_FORTRAN_ANYD(temp_ed[0][mfi])
#if (AMREX_SPACEDIM == 3)
			  , BL_TO_FORTRAN_ANYD(temp_ed[1][mfi]),
			  BL_TO_FORTRAN_ANYD(temp_ed[2][mfi])
#endif
			  );
  }

  mflux_cc_weighted.FillBoundary(geom.periodicity());
  for (int d=0; d<NUM_EDGE; ++d) {
    mflux_ed_weighted[d].FillBoundary(geom.periodicity());
  }
}

void StochMFlux::stochMforce(std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv,
			     const MultiFab& eta_cc,
			     const std::array< MultiFab, NUM_EDGE >& eta_ed,
			     const MultiFab& temp_cc,
			     const std::array< MultiFab, NUM_EDGE >& temp_ed,
			     const Vector< amrex::Real >& weights,
			     const amrex::Real& dt) {

  BL_PROFILE_VAR("StochMFlux::stochMforce()",stochMforce);

  // Take linear combination of mflux multifabs at each stage
  StochMFlux::weightMflux(weights);

  // Multiply weighted mflux (cc & edge) by sqrt(eta*temperature)
  StochMFlux::multbyVarSqrtEtaTemp(eta_cc,eta_ed,temp_cc,temp_ed,dt);

  const Real* dx = geom.CellSize();

  // Loop over boxes
  int increment = 0;
  for (MFIter mfi(mflux_cc_weighted); mfi.isValid(); ++mfi) {
    // Note: Make sure that multifab is cell-centered
    const Box& validBox = mfi.validbox();

    stoch_m_force(ARLIM_3D(validBox.loVect()), ARLIM_3D(validBox.hiVect()),
		  BL_TO_FORTRAN_FAB(mflux_cc_weighted[mfi]),
		  BL_TO_FORTRAN_FAB(mflux_ed_weighted[0][mfi]),
#if (AMREX_SPACEDIM == 3)
		  BL_TO_FORTRAN_FAB(mflux_ed_weighted[1][mfi]),
		  BL_TO_FORTRAN_FAB(mflux_ed_weighted[2][mfi]),
#endif
		  BL_TO_FORTRAN_ANYD(mfluxdiv[0][mfi]),
		  BL_TO_FORTRAN_ANYD(mfluxdiv[1][mfi]),
#if (AMREX_SPACEDIM == 3)
		  BL_TO_FORTRAN_ANYD(mfluxdiv[2][mfi]),
#endif
		  dx, &increment);
  }

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
      mfluxdiv[d].FillBoundary(geom.periodicity());
      // TODO: is this the right BC?
      MultiFABPhysBC(mfluxdiv[d], d, geom);
      // MultiFABPhysBCDomainVel(mfluxdiv[d], d);
      // MultiFABPhysBCMacVel(mfluxdiv[d], d);
  }
}

void StochMFlux::addMfluctuations(std::array< MultiFab, AMREX_SPACEDIM >& umac,
				  const MultiFab& rhotot, const MultiFab& Temp,
				  const amrex::Real& variance, Geometry geom) {

  std::array< MultiFab, AMREX_SPACEDIM > m_old;
  std::array< MultiFab, AMREX_SPACEDIM > rhotot_fc;
  std::array< MultiFab, AMREX_SPACEDIM > Temp_fc;

  for (int d=0; d<AMREX_SPACEDIM; d++) {
    m_old[d].define(     umac[d].boxArray(), umac[d].DistributionMap(), 1, 1);
    rhotot_fc[d].define( umac[d].boxArray(), umac[d].DistributionMap(), 1, 1);
    Temp_fc[d].define(   umac[d].boxArray(), umac[d].DistributionMap(), 1, 1);
  }

  // NOTE: these only operate on valid cells
  AverageCCToFace(rhotot, 0, rhotot_fc, 0, 1);
  AverageCCToFace(Temp,   0, Temp_fc,   0, 1);

  for (int d=0; d<AMREX_SPACEDIM; d++) {
      rhotot_fc[d].FillBoundary(geom.periodicity());
      Temp_fc[d].FillBoundary(geom.periodicity());

      MultiFABPhysBC(rhotot_fc[d], d, geom);
      MultiFABPhysBC(Temp_fc[d], d, geom);
  }

  // Convert umac to momenta, rho*umac
  for (int d=0; d<AMREX_SPACEDIM; d++) {
    MultiFab::Copy(     m_old[d], umac[d],      0, 0, 1, 1);
    MultiFab::Multiply( m_old[d], rhotot_fc[d], 0, 0, 1, 1);
  }

  addMfluctuations_stag(m_old, rhotot_fc, Temp_fc, variance,geom);


  // Convert momenta to umac, (1/rho)*momentum
  for (int d=0; d<AMREX_SPACEDIM; d++) {
      MultiFab::Copy(   umac[d], m_old[d],     0, 0, 1, 1);
      MultiFab::Divide( umac[d], rhotot_fc[d], 0, 0, 1, 1);
  }
}

void StochMFlux::addMfluctuations_stag(std::array< MultiFab, AMREX_SPACEDIM >& m_old,
				       const std::array< MultiFab, AMREX_SPACEDIM >& rhotot_fc,
				       const std::array< MultiFab, AMREX_SPACEDIM >& Temp_fc,
				       const amrex::Real& variance, Geometry geom) {

  const Real* dx = geom.CellSize();
  Real dVol = dx[0]*dx[1];
  if (AMREX_SPACEDIM == 2) {
    dVol *= cell_depth;
  } else {
    if (AMREX_SPACEDIM == 3) {
      dVol *= dx[2];
    }
  }

  // Initialize variances
  Real variance_mom = abs(variance)*k_B/dVol;
  std::array<MultiFab, AMREX_SPACEDIM> variance_mfab;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    variance_mfab[d].define(m_old[d].boxArray(), m_old[d].DistributionMap(),1,1);
  }

  std::array< MultiFab, AMREX_SPACEDIM > mac_temp;
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    mac_temp[d].define(m_old[d].boxArray(), m_old[d].DistributionMap(),1,1);
  }

  // Fill momentum multifab with random numbers, scaled by equilibrium variances
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    // Set variance multifab to sqrt(rho*temp)
    MultiFab::Copy(     variance_mfab[d],rhotot_fc[d],0,0,1,1);
    MultiFab::Multiply( variance_mfab[d],Temp_fc[d],  0,0,1,1);
    SqrtMF(variance_mfab[d]);

    // Fill momentum with random numbers, scaled by sqrt(var*k_B/dV)
    MultiFABFillRandom(mac_temp[d],0,variance_mom,geom);

    // TODO: add stochastic BCs here?


    // Scale random momenta further by factor of sqrt(rho*temp)
    MultiFab::Multiply(mac_temp[d],variance_mfab[d],0,0,1,1);

    MultiFab::Saxpy(m_old[d], 1.0, mac_temp[d],0,0,1,1);

    // For safety, although called by MultiFABFillRandom()
    m_old[d].OverrideSync(geom.periodicity());
    m_old[d].FillBoundary(geom.periodicity());
  }

  if (variance < 0.0) {
    // Ensure zero total momentum
    Vector<Real> av_mom;
    // take staggered sum & divide by number of cells
    SumStag(geom,m_old,0,av_mom,true);
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
      // subtract off average
      m_old[d].plus(-av_mom[d],1);
      m_old[d].OverrideSync(geom.periodicity());
      m_old[d].FillBoundary(geom.periodicity());
    }
  }

  for (int i=0; i<AMREX_SPACEDIM; i++) {
      m_old[i].FillBoundary(geom.periodicity());
      MultiFABPhysBCDomainVel(m_old[i], i, geom,i);
      MultiFABPhysBCMacVel(m_old[i], i, geom,i);
  }
}

void StochMFlux::writeMFs(std::array< MultiFab, AMREX_SPACEDIM >& mfluxdiv) {
  std::string plotfilename;
  std::string dimStr = "xyz";

  // Write out original fluxes
  for (int i=0; i<n_rngs; ++i){
    plotfilename = "a_mfluxcc_stage"+std::to_string(i);
    VisMF::Write(mflux_cc[i],plotfilename);

    for (int d=0; d<NUM_EDGE; ++d) {
      plotfilename = "a_mfluxnd_stage"+std::to_string(i)+"_";
      plotfilename += dimStr[d];
      VisMF::Write(mflux_ed[i][d],plotfilename);
    }
  }

  // Write out weighted fluxes
  plotfilename = "a_mfluxcc_weighted";
  VisMF::Write(mflux_cc_weighted,plotfilename);

  for (int d=0; d<NUM_EDGE; ++d) {
    plotfilename = "a_mfluxnd_weighted_";
    plotfilename += dimStr[d];
    VisMF::Write(mflux_ed_weighted[d],plotfilename);
  }

  // Write out fluxdiv
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    plotfilename = "a_mfluxdiv_";
    plotfilename += dimStr[d];
    VisMF::Write(mfluxdiv[d],plotfilename);
  }
}
