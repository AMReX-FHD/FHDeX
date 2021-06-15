#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include "DsmcParticleContainer.H"
#include "INS_functions.H"
using namespace amrex;
using namespace std;

void FhdParticleContainer::WritePlotFile(const MultiFab& covar, const MultiFab& cuMean,
               const Geometry& geom, Real time, int step){
	BL_PROFILE("WritePlotFile()");

	BoxArray ba = covar.boxArray();

	const DistributionMapping& dmap = covar.DistributionMap();
	const std::string& plotfilename = amrex::Concatenate("plt", step,9);

	amrex::Print() << "  Writing plotfile " << plotfilename << "\n";

	Vector<std::string> varnames;

	int nPlotvar = 15;     // covar 
	nPlotvar += (nspecies+1)*13; // cuMean
	MultiFab mfplot(ba, dmap, nPlotvar, 0);
	Vector<std::string> varNames(nPlotvar);
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
	varNames[cnt++] = "n_m";
	varNames[cnt++] = "rho_m";
	varNames[cnt++] = "Jx_m";
	varNames[cnt++] = "Jy_m";
	varNames[cnt++] = "Jz_m";
	varNames[cnt++] = "tau_xx_m";
	varNames[cnt++] = "tau_xy_m";
	varNames[cnt++] = "tau_xz_m";
	varNames[cnt++] = "tau_yy_m";
	varNames[cnt++] = "tau_yz_m";
	varNames[cnt++] = "tau_zz_m";
	varNames[cnt++] = "E_m";
	varNames[cnt++] = "Tg_m";

	for(int ispec=0;ispec<nspecies;ispec++) {
		varNames[cnt++] = "X_" + std::to_string(ispec);
		varNames[cnt++] = "Y_" + std::to_string(ispec);
		varNames[cnt++] = "u_" + std::to_string(ispec);
		varNames[cnt++] = "v_" + std::to_string(ispec);
		varNames[cnt++] = "w_" + std::to_string(ispec);
		varNames[cnt++] = "uu_" + std::to_string(ispec);
		varNames[cnt++] = "uv_" + std::to_string(ispec);
		varNames[cnt++] = "uw_" + std::to_string(ispec);
		varNames[cnt++] = "vv_" + std::to_string(ispec);
		varNames[cnt++] = "vw_" + std::to_string(ispec);
		varNames[cnt++] = "ww_" + std::to_string(ispec);
		varNames[cnt++] = "E_" + std::to_string(ispec);
		varNames[cnt++] = "Tg_" + std::to_string(ispec);
	}

	int plotstart = 0;
	int plotend = 15;
	MultiFab::Copy(mfplot, covar,  0, plotstart, plotend, 0);
	plotstart = plotend; // can add more Multifabs here
	plotend += (nspecies+1)*13;
	MultiFab::Copy(mfplot, cuMean, 0, plotstart, plotend, 0);
	amrex::WriteSingleLevelPlotfile(plotfilename, mfplot, varnames, geom, time, step);
}
