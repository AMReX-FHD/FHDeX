#include "common_functions.H"
#include "species.H"
#include "DsmcParticleContainer.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;
using namespace std;

void FhdParticleContainer::writePlotFile(const MultiFab& covar, const MultiFab& cuMean,
               const Geometry& geom, Real time, int statCount){
	BL_PROFILE("WritePlotFile()");
	
	int ncovar = 15;     			  // covar 
	int nmeanvar = (nspecies+1)*13; // cuMean
	int npltvar = ncovar + nmeanvar;

	amrex::BoxArray ba = cuMean.boxArray();
	amrex::DistributionMapping dmap = cuMean.DistributionMap();

	std::string pltfname = amrex::Concatenate("plt",statCount,9);
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
	varNames[cnt++] = "rhorho";
	varNames[cnt++] = "rhoJx";
	varNames[cnt++] = "rhoJy";
	varNames[cnt++] = "rhoJz";
	varNames[cnt++] = "rhoE";
	varNames[cnt++] = "JxJx";
	varNames[cnt++] = "JxJy";
	varNames[cnt++] = "JxJz";
	varNames[cnt++] = "JxE";
	varNames[cnt++] = "JyJy";
	varNames[cnt++] = "JyJz";
	varNames[cnt++] = "JyJE";
	varNames[cnt++] = "JzJz";
	varNames[cnt++] = "JzE";
	varNames[cnt++] = "EE";
	int pltbegin=0;
	MultiFab::Copy(mfplot, covar, 0, pltbegin, ncovar, 0);
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
	
	varNames[cnt++] = "nm";
	varNames[cnt++] = "rhom";
	varNames[cnt++] = "Jxm";
	varNames[cnt++] = "Jym";
	varNames[cnt++] = "Jzm";
	varNames[cnt++] = "txxm";
	varNames[cnt++] = "txym";
	varNames[cnt++] = "txzm";
	varNames[cnt++] = "tyym";
	varNames[cnt++] = "tyzm";
	varNames[cnt++] = "tzzm";
	varNames[cnt++] = "Em";
	varNames[cnt++] = "Tgm";

	for(int ispec=0;ispec<nspecies;ispec++) {
		varNames[cnt++] = amrex::Concatenate("X",ispec,2);
		varNames[cnt++] = amrex::Concatenate("Y",ispec,2);
		varNames[cnt++] = amrex::Concatenate("u",ispec,2);
		varNames[cnt++] = amrex::Concatenate("v",ispec,2);
		varNames[cnt++] = amrex::Concatenate("w",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uu",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uv",ispec,2);
		varNames[cnt++] = amrex::Concatenate("uw",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vv",ispec,2);
		varNames[cnt++] = amrex::Concatenate("vw",ispec,2);
		varNames[cnt++] = amrex::Concatenate("ww",ispec,2);
		varNames[cnt++] = amrex::Concatenate("E",ispec,2);
		varNames[cnt++] = amrex::Concatenate("Tg",ispec,2);
	}
	MultiFab::Copy(mfplot, cuMean, 0, pltbegin, nmeanvar, 0);
	pltbegin += nmeanvar;
	WriteSingleLevelPlotfile(pltfname, mfplot, varNames, geom, time, statCount);
}
