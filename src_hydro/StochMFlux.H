#ifndef _StochMFlux_H
#define _StochMFlux_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_VisMF.H>

#include <common_namespace.H>

using namespace amrex;
using namespace common;

class StochMFlux {

#if (AMREX_SPACEDIM == 2)
  const int ncomp_cc = 2;
#elif (AMREX_SPACEDIM == 3)
  const int ncomp_cc = 3;
#endif
  const int ncomp_ed = 2;
  
  int n_rngs = 1;
  int stoch_stress_form = 1;
  
  Vector< amrex::Real > weights;
  // weights = {std::sqrt(0.5), std::sqrt(0.5)};

  Geometry geom;
  BoxArray ba;
  DistributionMapping dmap;

  // mflux cell centered
  Vector< MultiFab > mflux_cc;

  // mflux on nodes in 2d
  // mflux on edges in 3d
  Vector<std::array< MultiFab, NUM_EDGE > >  mflux_ed;

  // weighted mflux cell centered & nodal
  MultiFab mflux_cc_weighted;
  std::array< MultiFab, NUM_EDGE >  mflux_ed_weighted;

  // eta & temperature cell centered
  MultiFab  eta_cc;
  MultiFab temp_cc;

  // eta & temperature nodal
  std::array< MultiFab, NUM_EDGE >   eta_ed;
  std::array< MultiFab, NUM_EDGE >  temp_ed;

  // mflux cell centred
  std::array< MultiFab, AMREX_SPACEDIM >  mfluxdiv;

public:

  StochMFlux(BoxArray,DistributionMapping,Geometry,int,int,Vector<amrex::Real>);
  void fillMStochastic();
  void weightMflux();
  void mult_by_sqrt_eta();
  void stochMforce();
  void writeMFs();

};

StochMFlux::StochMFlux(BoxArray ba_in, DistributionMapping dmap_in, Geometry geom_in, int stoch_stress_form_in=1, int n_rngs_in=1, Vector<amrex::Real> weights_in={1.0}) {
  n_rngs = n_rngs_in;
  stoch_stress_form = stoch_stress_form_in;
  geom = geom_in;
  ba = ba_in;
  dmap = dmap_in;
  
  weights.resize(n_rngs);
  mflux_cc.resize(n_rngs);
  mflux_ed.resize(n_rngs);

  for (int i=0; i<n_rngs; ++i) {
    weights[i] = weights_in[i];
  }

  // Define mfluxdiv multifabs
  mfluxdiv[0].define(convert(ba,nodal_flag_x), dmap, 1, 1);
#if (AMREX_SPACEDIM == 3)
  mfluxdiv[1].define(convert(ba,nodal_flag_y), dmap, 1, 1);
  mfluxdiv[2].define(convert(ba,nodal_flag_z), dmap, 1, 1);
#endif
  
  // Define mflux multifab vectors
  // TEMPORARY ASSUMPTION: filtering_width = 0
  for (int i=0; i<n_rngs; ++i) {
    mflux_cc[i].define(ba, dmap, ncomp_cc, 1);
  }
#if (AMREX_SPACEDIM == 2)
  for (int i=0; i<n_rngs; ++i) {
    mflux_ed[i][0].define(convert(ba,nodal_flag), dmap, ncomp_ed, 0);
  }
#elif (AMREX_SPACEDIM == 3)
  for (int i=0; i<n_rngs; ++i) {
    mflux_ed[i][0].define(convert(ba,nodal_flag_xy), dmap, ncomp_ed, 0);
    mflux_ed[i][1].define(convert(ba,nodal_flag_xz), dmap, ncomp_ed, 0);
    mflux_ed[i][2].define(convert(ba,nodal_flag_yz), dmap, ncomp_ed, 0);
  }
#endif

  // Define weighted mflux multifab vectors
  mflux_cc_weighted.define(ba, dmap, ncomp_cc, 1);
#if (AMREX_SPACEDIM == 2)
  mflux_ed_weighted[0].define(convert(ba,nodal_flag), dmap, ncomp_ed, 0);
#elif (AMREX_SPACEDIM == 3)
  mflux_ed_weighted[0].define(convert(ba,nodal_flag_xy), dmap, ncomp_ed, 0);
  mflux_ed_weighted[1].define(convert(ba,nodal_flag_xz), dmap, ncomp_ed, 0);
  mflux_ed_weighted[2].define(convert(ba,nodal_flag_yz), dmap, ncomp_ed, 0);
#endif
  
  // Define eta & temperature multifabs
  // eta cell-centered
  eta_cc.define(ba, dmap, 1, 1);
  eta_cc.setVal(1.0);
  // temperature cell-centered
  temp_cc.define(ba, dmap, 1, 1);
  temp_cc.setVal(1.0);
#if (AMREX_SPACEDIM == 2)
  // eta nodal
  eta_ed[0].define(convert(ba,nodal_flag), dmap, 1, 0);
  eta_ed[0].setVal(1.0);
  // temperature nodal
  temp_ed[0].define(convert(ba,nodal_flag), dmap, 1, 0);
  temp_ed[0].setVal(1.0);
#elif (AMREX_SPACEDIM == 3)
  // eta nodal
  eta_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
  eta_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
  eta_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
  eta_ed[0].setVal(1.0);
  eta_ed[1].setVal(1.0);
  eta_ed[2].setVal(1.0);
  // temperature nodal
  temp_ed[0].define(convert(ba,nodal_flag_xy), dmap, 1, 0);
  temp_ed[1].define(convert(ba,nodal_flag_xz), dmap, 1, 0);
  temp_ed[2].define(convert(ba,nodal_flag_yz), dmap, 1, 0);
  temp_ed[0].setVal(1.0);
  temp_ed[1].setVal(1.0);
  temp_ed[2].setVal(1.0);
#endif
}

void StochMFlux::weightMflux() {
  mflux_cc_weighted.setVal(0.0);
  for (int i=0; i<n_rngs; ++i) {
    MultiFab::Saxpy(mflux_cc_weighted, weights[i], mflux_cc[i], 0, 0, ncomp_cc, 1);
  }
  for (int d=0; d<NUM_EDGE; ++d) {
    mflux_ed_weighted[d].setVal(0.0);
    for (int i=0; i<n_rngs; ++i) {
      MultiFab::Saxpy(mflux_ed_weighted[d], weights[i], mflux_ed[i][d], 0, 0, ncomp_ed, 0);
    }
  }

}

void StochMFlux::fillMStochastic() {  
  for (int i=0; i<n_rngs; ++i){ 
    switch(stoch_stress_form) {
    case 0: // Non-symmetric
      // Print() << "Non-symmetric \n";
      for (int n=0; n<ncomp_cc; ++n) {
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
      for (int n=0; n<ncomp_cc; ++n) {
	MultiFABFillRandom(mflux_cc[i],n,2.0,geom);
      }

      for (int d=0; d<NUM_EDGE; ++d) {
	MultiFABFillRandom(mflux_ed[i][d],0,1.0,geom);
        MultiFab::Copy(mflux_ed[i][d], mflux_ed[i][d], 0, 1, ncomp_ed-1, 0);
      }
      break;
    }
  }

  // Take linear combination of mflux multifabs at each stage
  StochMFlux::weightMflux();
}

void StochMFlux::mult_by_sqrt_eta() { 
  // Multiply mflux_weighted by sqrt(eta*temperature)
}

void StochMFlux::stochMforce() { 
  // Multiply mflux (cc & edge) by sqrt(eta*temperature)
  StochMFlux::mult_by_sqrt_eta();

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

}

void StochMFlux::writeMFs() {
  std::string plotfilename;
  std::string dimStr = "xyz";
  for (int i=0; i<n_rngs; ++i){ 
    plotfilename = "a_mfluxcc_stage"+std::to_string(i);
    VisMF::Write(mflux_cc[i],plotfilename);
    for (int d=0; d<NUM_EDGE; ++d) {
      plotfilename = "a_mfluxnd_stage"+std::to_string(i)+"_";
      plotfilename += dimStr[d];
      VisMF::Write(mflux_ed[i][d],plotfilename);
    }
  }
} 

#endif
