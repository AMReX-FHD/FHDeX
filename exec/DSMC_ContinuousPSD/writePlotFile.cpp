#include "INS_functions.H"
#include "common_functions.H"
#include "species.H"
#include "DsmcParticleContainer.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;
using namespace std;

void FhdParticleContainer::writePlotFile(const MultiFab& mfcuInst,
						const MultiFab& mfcuMeans,
						const MultiFab& mfcuVars,
						const MultiFab& mfprimInst,
						const MultiFab& mfprimMeans,
						const MultiFab& mfprimVars,
						const MultiFab& mfcoVars,
						const Geometry& geom,
						Real time,
						int step) {
	BL_PROFILE_VAR("writePlotFile()",writePlotFile);

	int ncon    = 5;
	int nprim   = 14;
	int nvars   = 21 + ncon + nprim; // covariances + prim. vars + cons. vars

	amrex::BoxArray ba = mfcuInst.boxArray();
	amrex::DistributionMapping dmap = mfcuInst.DistributionMap();
	
	std::string pltcu    = amrex::Concatenate("pltcu",step,12);
	std::string pltprim  = amrex::Concatenate("pltprim",step,12);
	std::string pltvar   = amrex::Concatenate("pltvar",step,12);

	amrex::MultiFab mfcuplt(ba, dmap, 2*ncon, 0);
	amrex::MultiFab mfprimplt(ba, dmap, 2*nprim, 0);
	amrex::MultiFab mfvarplt(ba, dmap, nvars, 0);

	amrex::Vector<std::string> cuNames(2*ncon);
	amrex::Vector<std::string> primNames(2*nprim);
	amrex::Vector<std::string> varNames(nvars);

	//////////////////////////////////////
	// Conserved Means and Instants
	//////////////////////////////////////

	/*
		Conserved Vars:
		0  - rho = (1/V) += m
		1  - Jx  = (1/V) += mu
		2  - Jy  = (1/V) += mv
		3  - Jz  = (1/V) += mw
		4  - K   = (1/V) += m|v|^2
		... (repeat for each species)
	*/

	int cnt = 0;
	// Instant values
	cuNames[cnt++] = "rhoInstant";
	cuNames[cnt++] = "JxInstant";
	cuNames[cnt++] = "JyInstant";
	cuNames[cnt++] = "JzInstant";
	cuNames[cnt++] = "KInstant";
	MultiFab::Copy(mfcuplt, mfcuInst, 0, ncon*0, ncon, 0);

	// Mean values
	cuNames[cnt++] = "rhoMean";
	cuNames[cnt++] = "JxMean";
	cuNames[cnt++] = "JyMean";
	cuNames[cnt++] = "JzMean";
	cuNames[cnt++] = "KMean";
	MultiFab::Copy(mfcuplt, mfcuMeans, 0, ncon*1, ncon, 0);
	WriteSingleLevelPlotfile(pltcu, mfcuplt, cuNames, geom, time, step);

	//////////////////////////////////////
	// Primitive Means and Instants
	//////////////////////////////////////

	/*
	   Primitive Vars:
		0	- n  (X_ns)
		1  - rho(Y_ns)
		2  - u  (u_ns)
		3  - v  (v_ns)
		4  - w  (w_ns)
		5  - uu (uu_ns)
		6  - uv (uv_ns)
		7  - uw (uw_ns)
		8  - vv (vv_ns)
		9  - vw (vw_ns)
		10 - ww (ww_ns)
		11 - T  (P_ns)  = (1/3V) += m|v|^2
		12 - P  (T_ns)  = (1/3Nk_B) += m|v|^2
		13 - E  (E_ns)  = (1/2) += |v|^2 + c_v*T
		... (repeat for each species)
	*/

	cnt = 0;
	// Instant values
	primNames[cnt++] = "nInstant";
	primNames[cnt++] = "rhoInstant";
	primNames[cnt++] = "uInstant";
	primNames[cnt++] = "vInstant";
	primNames[cnt++] = "wInstant";
	primNames[cnt++] = "uuInstant";
	primNames[cnt++] = "uvInstant";
	primNames[cnt++] = "uwInstant";
	primNames[cnt++] = "vvInstant";
	primNames[cnt++] = "vwInstant";
	primNames[cnt++] = "wwInstant";
	primNames[cnt++] = "TInstant";
	primNames[cnt++] = "PInstant";
	primNames[cnt++] = "EInstant";
	MultiFab::Copy(mfprimplt, mfprimInst, 0, nprim*0, nprim, 0);

	// Mean values
	primNames[cnt++] = "nMean";
	primNames[cnt++] = "rhoMean";
	primNames[cnt++] = "uMean";
	primNames[cnt++] = "vMean";
	primNames[cnt++] = "wMean";
	primNames[cnt++] = "uuMean";
	primNames[cnt++] = "uvMean";
	primNames[cnt++] = "uwMean";
	primNames[cnt++] = "vvMean";
	primNames[cnt++] = "vwMean";
	primNames[cnt++] = "wwMean";
	primNames[cnt++] = "TMean";
	primNames[cnt++] = "PMean";
	primNames[cnt++] = "EMean";
	MultiFab::Copy(mfprimplt, mfprimMeans, 0, nprim*1, nprim, 0);
	WriteSingleLevelPlotfile(pltprim, mfprimplt, primNames, geom, time, step);

	//////////////////////////////////////
	// Variances
	//////////////////////////////////////

	// Covariances
	/*
		0  - drho.dJx
		1  - drho.dJy
		2  - drho.dJz
		3  - drho.dT
		4  - drho.d(rho*E)
		5  - dJx.dJy
		6  - dJx.dJz
		7  - dJy.dJz
		8  - dJx.d(rho*E)
		9  - dJy.d(rho*E)
		10 - dJz.d(rho*E)
		11 - drho.du
		12 - drho.dv
		13 - drho.dw
		14 - du.dv
		15 - du.dw
		16 - dv.dw
		17 - drho.dT
		18 - du.dT
		19 - dv.dT
		20 - dw.dT
	*/
	cnt = 0;
	int istart = 0;

	// Conserved Variances
	varNames[cnt++] = "rhoVar";
	varNames[cnt++] = "JxVar";
	varNames[cnt++] = "JyVar";
	varNames[cnt++] = "JzVar";
	varNames[cnt++] = "KVar";
	MultiFab::Copy(mfvarplt, mfcuVars, 0, istart, ncon, 0);
	istart += ncon;

	// Primitive Variances
	varNames[cnt++] = "nVar";
	varNames[cnt++] = "rhoVar";
	varNames[cnt++] = "uVar";
	varNames[cnt++] = "vVar";
	varNames[cnt++] = "wVar";
	varNames[cnt++] = "uuVar";
	varNames[cnt++] = "uvVar";
	varNames[cnt++] = "uwVar";
	varNames[cnt++] = "vvVar";
	varNames[cnt++] = "vwVar";
	varNames[cnt++] = "wwVar";
	varNames[cnt++] = "TVar";
	varNames[cnt++] = "PVar";
	varNames[cnt++] = "EVar";
	MultiFab::Copy(mfvarplt, mfprimVars, 0, istart, nprim, 0);
	istart += nprim;

	varNames[cnt++] = "rho.Jx";
	varNames[cnt++] = "rho.Jy";
	varNames[cnt++] = "rho.Jz";
	varNames[cnt++] = "rho.T";
	varNames[cnt++] = "rho.E";
	varNames[cnt++] = "Jx.Jy";
	varNames[cnt++] = "Jx.Jz";
	varNames[cnt++] = "Jy.Jz";
	varNames[cnt++] = "Jx.E";
	varNames[cnt++] = "Jy.E";
	varNames[cnt++] = "Jz.E";
	varNames[cnt++] = "rho.u";
	varNames[cnt++] = "rho.v";
	varNames[cnt++] = "rho.w";
	varNames[cnt++] = "u.v";
	varNames[cnt++] = "u.w";
	varNames[cnt++] = "v.w";
	varNames[cnt++] = "rho.T";
	varNames[cnt++] = "u.T";
	varNames[cnt++] = "v.T";
	varNames[cnt++] = "w.T";
	MultiFab::Copy(mfvarplt, mfcoVars, 0, istart, 21, 0);
	WriteSingleLevelPlotfile(pltvar, mfvarplt, varNames, geom, time, step);
}
