#include "INS_functions.H"
#include "common_functions.H"
#include "species.H"
#include "DsmcParticleContainer.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;
using namespace std;

void FhdParticleContainer::writePlotFile(const MultiFab& mfcuConInst,
   						 						  const MultiFab& mfcuConMeans,
   						 						  const MultiFab& mfcuPrimInst,
   						 						  const MultiFab& mfcuPrimMeans,
   						 						  const MultiFab& mfcuVars,
   						 						  const MultiFab& mfcoVars,
		                                   const Geometry& geom, 
		                                   Real time,
		                                   int step) {
	BL_PROFILE("writePlotFile()");
	
	int npltvar = 0;
	int ncovar  = 21;
	int nprim   = (nspecies+1)*13;
	int ncon    = (nspecies+1)*5;
	npltvar    += ncovar;	// covar
	npltvar    += 3*nprim;	// primitive (3x - Inst, Mean, Var)
	npltvar    += 3*ncon;	// conserved (3x)

	amrex::BoxArray ba = mfcuConInst.boxArray();
	amrex::DistributionMapping dmap = mfcuConInst.DistributionMap();

	std::string pltfname = amrex::Concatenate("plt",step,9);
	amrex::Print() << " Writing plotfile " << pltfname << "\n";
	
	amrex::MultiFab mfplot(ba, dmap, npltvar, 0);
	amrex::Vector<std::string> varNames(npltvar);
	
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
	
	int cnt = 0;
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
	
	int pltbegin=0;
	MultiFab::Copy(mfplot, mfcoVars, 0, pltbegin, ncovar, 0);
	pltbegin += ncovar;
	
	/*
		Conserved Vars:
		0  - rho = (1/V) += m
		1  - Jx  = (1/V) += mu
		2  - Jy  = (1/V) += mv
		3  - Jz  = (1/V) += mw
		4  - K   = (1/V) += m|v|^2
		... (repeat for each species)
	*/
	
   // Instant values
	varNames[cnt++] = "rhoInstant";
	varNames[cnt++] = "JxInstant";
	varNames[cnt++] = "JyInstant";
	varNames[cnt++] = "JzInstant";
	varNames[cnt++] = "KInstant";

	for(int ispec=0;ispec<nspecies;ispec++) {
		varNames[cnt++] = amrex::Concatenate("rhoInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("JxInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("JyInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("JzInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("KInstant_",ispec,2);
	}
	MultiFab::Copy(mfplot, mfcuConInst, 0, pltbegin, ncon, 0);
	pltbegin += ncon;

   // Mean values
	varNames[cnt++] = "rhoMean";
	varNames[cnt++] = "JxMean";
	varNames[cnt++] = "JyMean";
	varNames[cnt++] = "JzMean";
	varNames[cnt++] = "KMean";

	for(int ispec=0;ispec<nspecies;ispec++) {
		varNames[cnt++] = amrex::Concatenate("rhoMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("JxMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("JyMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("JzMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("KMean_",ispec,2);
	}
	MultiFab::Copy(mfplot, mfcuConMeans, 0, pltbegin, ncon, 0);
	pltbegin += ncon;
	
	/*
	   Primitive Vars:
		0	- n  (X_ns)
		1  - u  (u_ns)
		2  - v  (v_ns)
		3  - w  (w_ns)
		4  - uu (uu_ns)
		5  - uv (uv_ns)
		6  - uw (uw_ns)
		7  - vv (vv_ns)
		8  - vw (vw_ns)
		9  - ww (ww_ns)
		10 - P  (P_ns)  = (1/3V) += m|v|^2
		11 - T  (T_ns)  = (1/3Nk_B) += m|v|^2
		12 - E  (E_ns)  = (1/2) += |v|^2 + c_v*T
		... (repeat for each species)
	*/	
	
   // Instant values
	varNames[cnt++] = "nInstant";
	varNames[cnt++] = "uInstant";
	varNames[cnt++] = "vInstant";
	varNames[cnt++] = "wInstant";
	varNames[cnt++] = "uuInstant";
	varNames[cnt++] = "uvInstant";
	varNames[cnt++] = "uwInstant";
	varNames[cnt++] = "vvInstant";
	varNames[cnt++] = "vwInstant";
	varNames[cnt++] = "wwInstant";
	varNames[cnt++] = "PInstant";
	varNames[cnt++] = "TInstant";
	varNames[cnt++] = "EInstant";

	for(int ispec=0;ispec<nspecies;ispec++) {
		varNames[cnt++] = amrex::Concatenate("nInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("wInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uuInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uvInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uwInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vvInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vwInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("wwInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("PInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("TInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("EInstant_",ispec,2);
	}
	
	MultiFab::Copy(mfplot, mfcuPrimInst, 0, pltbegin, nprim, 0);
	pltbegin += nprim;
	
   // Mean values
	varNames[cnt++] = "nMean";
	varNames[cnt++] = "uMean";
	varNames[cnt++] = "vMean";
	varNames[cnt++] = "wMean";
	varNames[cnt++] = "uuMean";
	varNames[cnt++] = "uvMean";
	varNames[cnt++] = "uwMean";
	varNames[cnt++] = "vvMean";
	varNames[cnt++] = "vwMean";
	varNames[cnt++] = "wwMean";
	varNames[cnt++] = "PMean";
	varNames[cnt++] = "TMean";
	varNames[cnt++] = "EMean";

	for(int ispec=0;ispec<nspecies;ispec++) {
		varNames[cnt++] = amrex::Concatenate("nMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("wMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uuMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uvMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uwMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vvMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vwMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("wwMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("PMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("TMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("EMean_",ispec,2);
	}
	
	MultiFab::Copy(mfplot, mfcuPrimMeans, 0, pltbegin, nprim, 0);
	pltbegin += nprim;
	
   // Conserved Variances
	varNames[cnt++] = "rhoVar";
	varNames[cnt++] = "JxVar";
	varNames[cnt++] = "JyVar";
	varNames[cnt++] = "JzVar";
	varNames[cnt++] = "KVar";

	for(int ispec=0;ispec<nspecies;ispec++) {
		varNames[cnt++] = amrex::Concatenate("rhoVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("JxVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("JyVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("JzVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("KVar_",ispec,2);
	}

	// Primitive Variances
	varNames[cnt++] = "nVar";
	varNames[cnt++] = "uVar";
	varNames[cnt++] = "vVar";
	varNames[cnt++] = "wVar";
	varNames[cnt++] = "uuVar";
	varNames[cnt++] = "uvVar";
	varNames[cnt++] = "uwVar";
	varNames[cnt++] = "vvVar";
	varNames[cnt++] = "vwVar";
	varNames[cnt++] = "wwVar";
	varNames[cnt++] = "PVar";
	varNames[cnt++] = "TVar";
	varNames[cnt++] = "EVar";

	for(int ispec=0;ispec<nspecies;ispec++) {
		varNames[cnt++] = amrex::Concatenate("nVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("wVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uuVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uvVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uwVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vvVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vwVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("wwVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("PVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("TVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("EVar_",ispec,2);
	}
	MultiFab::Copy(mfplot, mfcuVars, 0, pltbegin, ncon+nprim, 0);
	pltbegin += (ncon+nprim);
	
	WriteSingleLevelPlotfile(pltfname, mfplot, varNames, geom, time, step);
}
