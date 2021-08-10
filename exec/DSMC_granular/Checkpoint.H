#ifndef _Checkpoint_H_
#define _Checkpoint_H_

#include "AMReX_PlotFileUtil.H"
#include "AMReX_PlotFileDataImpl.H"
#include "DsmcParticleContainer.H"
#include <sys/stat.h>

#include "common_functions.H"
#include "common_namespace.H"

void WriteCheckPoint(
	int step,
	const Real time,
	const Real dt,
	int statsCount,
	int timeCount,
	const amrex::MultiFab& cuInst,
	const amrex::MultiFab& cuMeans,
	const amrex::MultiFab& cuVars,
	const amrex::MultiFab& primInst,
	const amrex::MultiFab& primMeans,
	const amrex::MultiFab& primVars,
	const amrex::MultiFab& coVars,
	const FhdParticleContainer& particles,
	const amrex::MultiFab& spatialCross1D,
	const amrex::MultiFab& rhotimeCross,
	const amrex::MultiFab& utimeCross,
	const amrex::MultiFab& KtimeCross,
	const amrex::MultiFab& rho_time,
	const amrex::MultiFab& u_time,
	const amrex::MultiFab& K_time,
	const amrex::MultiFab& vmom);

void ReadCheckPoint(
	int& step,
	Real& time,
	Real& dt,
	int& statsCount,
	int& timeCount,
	amrex::MultiFab& cuInst,
	amrex::MultiFab& cuMeans,
	amrex::MultiFab& cuVars,
	amrex::MultiFab& primInst,
	amrex::MultiFab& primMeans,
	amrex::MultiFab& primVars,
	amrex::MultiFab& coVars,
	amrex::MultiFab& spatialCross1D,
	amrex::MultiFab& rhotimeCross,
	amrex::MultiFab& utimeCross,
	amrex::MultiFab& KtimeCross,
	amrex::MultiFab& rho_time,
	amrex::MultiFab& u_time,
	amrex::MultiFab& K_time,
	amrex::MultiFab& vmom,
	const int ncuvars,
	const int nprimvars,
	const int ncovars,
	const int ncross,
	const int ntime,
	const int npart);

void ReadCheckPointParticles(FhdParticleContainer& particles);

void ReadFile(const std::string& filename, Vector<char>& charBuf, 
              bool bExitOnError);
              
void globalVars(Real time);

#endif
