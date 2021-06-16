#include "common_functions.H"
#include "species.H"
#include "DsmcParticleContainer.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;
using namespace std;

void FhdParticleContainer::writePlotFile(const MultiFab& Cu,
                                         const MultiFab& CuMeans, 
                                         const MultiFab& CuVars,
                                         const MultiFab& CoVars,
                                         const Geometry& geom, 
                                         Real time, 
                                         int plotCount) {
	BL_PROFILE("writePlotFile()");
	
	int ncovar = 10;     			     // covar 
	int nstats  = (nspecies+1)*13; // variables
	int npltvar = ncovar + 3*nstats; // covars + instant + mean + variance

	amrex::BoxArray ba = Cu.boxArray();
	amrex::DistributionMapping dmap = Cu.DistributionMap();

	std::string pltfname = amrex::Concatenate("plt",plotCount,9);
	amrex::Print() << " Writing plotfile " << pltfname << "\n";
	
	amrex::MultiFab mfplot(ba, dmap, npltvar, 0);
	amrex::Vector<std::string> varNames(npltvar);
	// Co/variances
	/*
		0  - drho.drho
		1  - drho.dJx
		2  - drho.dJy
		3  - drho.dJz
		4  - drho.dE
		5  - dJx.dJx
		6  - dJx.dJy
		7  - dJx.dJz
		8  - dJy.dJy
		9  - dJy.dJz
		10 - dJz.dJz
		11 - dJx.dE
		12 - dJy.dE
		13 - dJz.dE
		14 - dE.dE
	*/
	
	int cnt = 0;
	varNames[cnt++] = "rhoJx";
	varNames[cnt++] = "rhoJy";
	varNames[cnt++] = "rhoJz";
	varNames[cnt++] = "rhoE";
	varNames[cnt++] = "JxJy";
	varNames[cnt++] = "JxJz";
	varNames[cnt++] = "JyJz";
	varNames[cnt++] = "JxE";
	varNames[cnt++] = "JyJE";
	varNames[cnt++] = "JzE";
	int pltbegin=0;
	MultiFab::Copy(mfplot, CoVars, 0, pltbegin, ncovar, 0);
	pltbegin += ncovar;
	/*
		0 	- n - number density
		1	- rho - mass density
		2 	- Jx - x-mom density
		3 	- Jy - y-mom density
		4 	- Jz - z-mom density
		5 	- tau_xx - xx shear stress
		6 	- tau_xy - xy shear stress
		7 	- tau_xz - xz shear stress
		8 	- tau_yy - yy shear stress
		9  - tau_yz - yz shear stress
		10 - tau_zz - zz shear stress
		11 - E - energy
		12 - T_g - granular temperature
		13 - X_i - mole fraction for spec. i
		14 - Y_i - mass fraction for spec. i
		15 - u_i - x-vel for spec. i
		16 - v_i - y-vel for spec. i
		17 - w_i - z-vel for spec. i
		18 - uu_i 
		19 - uv_i
		20 - uw_i
		21 - vv_i
		22 - vw_i
		23 - ww_i
		24 - E_i
		25 - T_g_i
		... (repeat for each add. species)
	*/
	
  // Instant values
	varNames[cnt++] = "nInstant";
	varNames[cnt++] = "rhoInstant";
	varNames[cnt++] = "JxInstant";
	varNames[cnt++] = "JyInstant";
	varNames[cnt++] = "JzInstant";
	varNames[cnt++] = "txxInstant";
	varNames[cnt++] = "txyInstant";
	varNames[cnt++] = "txzInstant";
	varNames[cnt++] = "tyyInstant";
	varNames[cnt++] = "tyzInstant";
	varNames[cnt++] = "tzzInstant";
	varNames[cnt++] = "EInstant";
	varNames[cnt++] = "TgInstant";

	for(int ispec=0;ispec<nspecies;ispec++) {
		varNames[cnt++] = amrex::Concatenate("XInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("YInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("wInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uuInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uvInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uwInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vvInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vwInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("wwInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("EInstant_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("TgInstant_",ispec,2);
	}
	MultiFab::Copy(mfplot, Cu, 0, pltbegin, nstats, 0);
	pltbegin += nstats;
	
  // Mean Values
	varNames[cnt++] = "nMean";
	varNames[cnt++] = "rhoMean";
	varNames[cnt++] = "JxMean";
	varNames[cnt++] = "JyMean";
	varNames[cnt++] = "JzMean";
	varNames[cnt++] = "txxMean";
	varNames[cnt++] = "txyMean";
	varNames[cnt++] = "txzMean";
	varNames[cnt++] = "tyyMean";
	varNames[cnt++] = "tyzMean";
	varNames[cnt++] = "tzzMean";
	varNames[cnt++] = "EMean";
	varNames[cnt++] = "TgMean";

	for(int ispec=0;ispec<nspecies;ispec++) {
		varNames[cnt++] = amrex::Concatenate("XMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("YMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("wMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uuMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uvMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uwMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vvMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vwMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("wwMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("EMean_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("TgMean_",ispec,2);
	}
	MultiFab::Copy(mfplot, CuMeans, 0, pltbegin, nstats, 0);
	pltbegin += nstats;
	
  // Variances
	varNames[cnt++] = "nVar";
	varNames[cnt++] = "rhoVar";
	varNames[cnt++] = "JxVar";
	varNames[cnt++] = "JyVar";
	varNames[cnt++] = "JzVar";
	varNames[cnt++] = "txxVar";
	varNames[cnt++] = "txyVar";
	varNames[cnt++] = "txzVar";
	varNames[cnt++] = "tyyVar";
	varNames[cnt++] = "tyzVar";
	varNames[cnt++] = "tzzVar";
	varNames[cnt++] = "EVar";
	varNames[cnt++] = "TgVar";

	for(int ispec=0;ispec<nspecies;ispec++) {
		varNames[cnt++] = amrex::Concatenate("XVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("YVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("wVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uuVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uvVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uwVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vvVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vwVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("wwVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("EVar_",ispec,2);
		varNames[cnt++] = amrex::Concatenate("TgVar_",ispec,2);
	}
	MultiFab::Copy(mfplot, CuVars, 0, pltbegin, nstats, 0);
	pltbegin += nstats;

	WriteSingleLevelPlotfile(pltfname, mfplot, varNames, geom, time, plotCount);
}
