#ifndef _multispec_functions_H_
#define _multispec_functions_H_

#include <AMReX_MultiFab.H>

#include "multispec_functions_F.H"
#include "multispec_namespace.H"

#include "common_functions.H"

#include "StochMassFlux.H"
#include "StochMomFlux.H"

using namespace multispec;
using namespace amrex;

/////////////////////////////////////////////////////////////////////////////////
// in multispec_functions.cpp

void InitializeMultispecNamespace();

/////////////////////////////////////////////////////////////////////////////////
// in ComputeDivReversibleStress.cpp
void ComputeDivReversibleStress(std::array<MultiFab,AMREX_SPACEDIM>& div_reversible_stress,
                                const MultiFab& rhotot_in,
                                MultiFab& rho_in,
                                const Geometry& geom);

/////////////////////////////////////////////////////////////////////////////////
// in ComputeMassFluxdiv.cpp

void ComputeMassFluxdiv(MultiFab& rho,
                        MultiFab& rhotot,
                        const MultiFab& Temp,
			MultiFab& diff_mass_fluxdiv,
                        MultiFab& stoch_mass_fluxdiv,
			std::array<MultiFab,AMREX_SPACEDIM>& diff_mass_flux,
                        std::array<MultiFab,AMREX_SPACEDIM>& stoch_mass_flux,
                        StochMassFlux& sMassFlux,
			const Real& dt, const Real& stage_time, const Geometry& geom,
                        Vector<Real>& weights,
                        MultiFab& charge,
                        std::array<MultiFab,AMREX_SPACEDIM>& grad_Epot,
                        MultiFab& Epot,
                        MultiFab& permittivity,
                        const int& zero_initial_Epot=1);

void ComputeHigherOrderTerm(const MultiFab& molarconc,
                            std::array<MultiFab,AMREX_SPACEDIM>& diff_mass_flux,
                            const Geometry& geom);

/////////////////////////////////////////////////////////////////////////////////
// in DiffusiveMassFlux.cpp

void DiffusiveMassFluxdiv(const MultiFab& rho,
			  const MultiFab& rhotot,
			  const MultiFab& molarconc,
			  const MultiFab& rhoWchi,
			  const MultiFab& Gamma,
			  MultiFab& diff_mass_fluxdiv,
			  std::array< MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
			  const Geometry& geom);

void DiffusiveMassFlux(const MultiFab& rho,
		       const MultiFab& rhotot,
		       const MultiFab& molarconc,
		       const MultiFab& rhoWchi,
		       const MultiFab& Gamma,
		       std::array< MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
		       const Geometry& geom);

/////////////////////////////////////////////////////////////////////////////////
// in ComputeMixtureProperties.cpp
void ComputeMixtureProperties(const MultiFab& rho,
			      const MultiFab& rhotot,
			      MultiFab& D_bar,
			      MultiFab& D_therm,
			      MultiFab& Hessian);

/////////////////////////////////////////////////////////////////////////////////
// in CorrectionFlux.cpp

void CorrectionFlux(const MultiFab& rho, const MultiFab& rhotot,
		    std::array< MultiFab, AMREX_SPACEDIM >& flux);

/////////////////////////////////////////////////////////////////////////////////
// in CorrectionFlux.cpp

void CorrectionFlux(const MultiFab& rho, const MultiFab& rhotot,
		    std::array< MultiFab, AMREX_SPACEDIM >& flux);

/////////////////////////////////////////////////////////////////////////////////
// in ElectroDiffusiveMassFluxdiv.cpp

void ElectroDiffusiveMassFluxdiv(const MultiFab& rho,
                                 const MultiFab& Temp,
                                 const MultiFab& rhoWchi,
                                 std::array< MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
                                 MultiFab& diff_mass_fluxdiv,
                                 std::array< MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                                 MultiFab& charge,
                                 std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot,
                                 MultiFab& Epot,
                                 const MultiFab& permittivity,
                                 Real dt,
                                 int zero_initial_Epot,
                                 const Geometry& geom);

void ElectroDiffusiveMassFlux(const MultiFab& rho,
                              const MultiFab& Temp,
                              const MultiFab& rhoWchi,
                              std::array< MultiFab, AMREX_SPACEDIM >& electro_mass_flux,
                              std::array< MultiFab, AMREX_SPACEDIM >& diff_mass_flux,
                              std::array< MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                              MultiFab& charge,
                              std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot,
                              MultiFab& Epot,
                              const MultiFab& permittivity,
                              Real dt,
                              int zero_initial_Epot,
                              const Geometry& geom);

void LimitEMF(const MultiFab& rho_in,
              std::array< MultiFab, AMREX_SPACEDIM >& electro_mass_flux);

/////////////////////////////////////////////////////////////////////////////////
// in FluicCharge.cpp

void DotWithZ(const MultiFab& mf,
              MultiFab& mfdotz,
              int abs_z=0);

void DotWithZFace(std::array< const MultiFab, AMREX_SPACEDIM >& mf,
                  std::array< MultiFab, AMREX_SPACEDIM >& mfdotz,
                  int abs_0);

void ComputeChargeCoef(const MultiFab& rho_in,
                       const MultiFab& Temp_in,
                       MultiFab& charge_coef_in);

void EnforceChargeNeutrality();

void ImplicitPotentialCoef();

void ModifyS();

void ComputePermittivity();

void ComputeLorentzForce(std::array< MultiFab, AMREX_SPACEDIM >& Lorentz_force,
                         std::array< MultiFab, AMREX_SPACEDIM >& grad_Epot,
                         const MultiFab& permittivity,
                         const MultiFab& charge,
                         const Geometry& geom);

void ComputeE_ext(std::array< MultiFab, AMREX_SPACEDIM >& E_ext);

void ZeroEpsOnWall(std::array< MultiFab, AMREX_SPACEDIM >& beta);

/////////////////////////////////////////////////////////////////////////////////
// in InitialProjection.cpp

void InitialProjection(std::array< MultiFab, AMREX_SPACEDIM >& umac,
                       MultiFab& rho, MultiFab& rhotot,
                       MultiFab& diff_mass_fluxdiv,
                       MultiFab& stoch_mass_fluxdiv,
                       std::array< MultiFab, AMREX_SPACEDIM >& stoch_mass_flux,
                       StochMassFlux& sMassFlux,
                       const MultiFab& Temp, const MultiFab& eta,
                       const std::array< MultiFab, NUM_EDGE >& eta_ed,
                       const Real& dt, const Real& time, const Geometry& geom,
                       MultiFab& charge_old,
                       std::array<MultiFab,AMREX_SPACEDIM>& grad_Epot_old,
                       MultiFab& Epot,
                       MultiFab& permittivity);

/////////////////////////////////////////////////////////////////////////////////
// in MassFluxUtil.cpp

void ComputeMolconcMolmtot(const MultiFab& rho,
			   const MultiFab& rhotot,
			   MultiFab& molarconc,
			   MultiFab& molmtot);

void ComputeGamma(const MultiFab& molarconc,
		  const MultiFab& Hessian,
		  MultiFab& Gamma);

void ComputeRhoWChi(const MultiFab& rho,
		    const MultiFab& rhotot,
		    const MultiFab& molarconc,
		    MultiFab& rhoWchi,
		    const MultiFab& D_bar);

void ComputeZetaByTemp(const MultiFab& molarconc,
 		       const MultiFab& D_bar,
 		       const MultiFab& Temp,
 		       MultiFab& zeta_by_Temp,
 		       const MultiFab& D_therm);

void ComputeSqrtLonsagerFC(const MultiFab& rho,
                           const MultiFab& rhotot,
                           std::array< MultiFab, AMREX_SPACEDIM >& sqrtLonsager_fc,
                           const Geometry& geom);

/////////////////////////////////////////////////////////////////////////////////
// in MatvecMul.cpp

void MatvecMul(MultiFab& x,
	       const MultiFab& A);

/////////////////////////////////////////////////////////////////////////////////
// in MkDiffusiveMFluxdiv.cpp

void MkDiffusiveMFluxdiv(std::array<MultiFab, AMREX_SPACEDIM> & m_update,
                         const std::array<MultiFab, AMREX_SPACEDIM> & umac,
                         const MultiFab& eta,
                         const std::array<MultiFab, NUM_EDGE> & eta_ed,
                         const MultiFab& kappa,
                         const Geometry& geom,
                         const Real* dx,
                         const int& increment);

/////////////////////////////////////////////////////////////////////////////////
// in ProjectOntoEOS.cpp

void ProjectOntoEOS(MultiFab& rho_in);

/////////////////////////////////////////////////////////////////////////////////
// in RhoUtil.cpp

void RhototBCInit();

void ComputeRhotot(const MultiFab& rho, MultiFab& rhotot, int include_ghost=0);

void ConvertRhoCToC(MultiFab& rho, const MultiFab& rhotot, MultiFab& conc, int rho_to_c);

void FillRhoRhototGhost(MultiFab& rho, MultiFab& rhotot, const Geometry& geom);

void FillRhototGhost(MultiFab& rhotot_in, const MultiFab& conc_in, const Geometry& geom);

/////////////////////////////////////////////////////////////////////////////////
// Device calls from  MassFluxUtil.cpp
void ComputeChi(int nspecies_in, 
               GpuArray<Real, MAX_SPECIES>& molmass_in,
               GpuArray<Real, MAX_SPECIES>& rhoN,
               Real rhotot,
               GpuArray<Real, MAX_SPECIES>& MolarConcN,
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& chi,
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& D_barN,
               int chi_iterations_in); 

/** 
 * \brief Compute total molar mass
 *
 * \param[in] W Mass fractions.
 * \param[out] molmtot Total molar mass
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetMolTot( GpuArray<Real, MAX_SPECIES>& W,
                GpuArray<Real, MAX_SPECIES>& molmass_in,
                Real& molmtot,
                int nspecies_in)
{
    Real sumOverN = 0;

    for (int n=0; n<nspecies_in; ++n){
        sumOverN += W[n] / molmass_in[n];                
    }  

    molmtot = 1.0 / sumOverN; 
}  

/**
 * \param[in] nspecies_in Number of species
 * \param[in] molmass_in 
 * \param[in] RhoN Density, last dim for number of species
 * \param[in] rhotot_in Total density in each cell
 * \param[out] MolarConcN Molar concentration
 * \param[out] molmtot_in Total molar mass 
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeMolconcMolmtotLocal( int nspecies_in,
                GpuArray<Real, MAX_SPECIES>& molmass_in,
                GpuArray<Real, MAX_SPECIES>& RhoN,
                Real rhotot_in,
                GpuArray<Real, MAX_SPECIES>& MolarConcN,
                Real& molmtot_in,
                int& func_call_counter){
        
    func_call_counter++;
//Print().SetPrecision(17) << " inner B molmtot_in: " << molmtot_in << std::endl;

    GpuArray<Real, MAX_SPECIES> W; /**< Mass fraction w_i = rho_i/rho */

    // calculate mass fraction and total molar mass (1/m=Sum(w_i/m_i))

    for (int n=0; n<nspecies_in; ++n){
        W[n] = RhoN[n] / rhotot_in;
    }


//Print().SetPrecision(17) << "inner B molmtot_in: " << molmtot_in << std::endl;
    GetMolTot(W, molmass_in, molmtot_in, nspecies_in); 
//Print().SetPrecision(17) << "inner A molmtot_in: " << molmtot_in << std::endl;

    // calculate molar concentrations in each cell (x_i=m*w_i/m_i) 

    for (int n=0; n<nspecies_in; ++n){
        MolarConcN[n] = molmtot_in * W[n] / molmass_in[n]; //HACK error found
    }
    
//Print().SetPrecision(17) << "inner A molmtot_in: " << molmtot_in;
Print().SetPrecision(17) << "molarconc(1): " << MolarConcN[0] << std::endl;
Print().SetPrecision(17) << "molarconc(2): " << MolarConcN[1] << std::endl;
Print().SetPrecision(17) << "molarconc(3): " << MolarConcN[2] << std::endl;

}


/**
 * \brief Set all values in 2D Array to zero.
 *
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void MatrixToZeros(int nspecies_in,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& matrixIn) {

    for (int n=1; n<=nspecies_in; ++n ){
        for (int m=1; m<=nspecies_in; ++m ){
            matrixIn(n,m) = 0.0;
        }
    }
}

/**
 * \brief Populate X_xxT
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PopulateX_xxT(int nspecies,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& X_xxT,
                GpuArray<Real, MAX_SPECIES>& MolarConcN) {

    for (int row=1; row<=nspecies;++row ){ 
        // diagonal entries
        X_xxT(row,row) = MolarConcN[row-1] - std::pow(MolarConcN[row-1],2);
        for (int column=1; column<=row-1; ++column ){
            // off-diagonal entries
            // form x*traspose(x) off diagonals 
            X_xxT(row,column) = -MolarConcN[row-1]*MolarConcN[column-1];
            //symmetric
            X_xxT(column,row) = X_xxT(row,column);
        }
    }
}

/**
 * \brief Compute Gamma = I + matmul(X_xxT, Hessian). j
 *
 * Adds the identity matrix during initialization. Then adds product of X_xxT and Hessian matrices.
 *
 * \param[in] nspecies_in Number of species.
 * \param[out] GammaN 
 * \param[in] X_xxT 
 * \param[in] HessianN 
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GammaIPlusMatmul(int nspecies_in,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& GammaN,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& X_xxT,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& HessianN) {
        
    //Compute Gamma 
    for (int row=1; row<=nspecies_in; ++row){
        for (int column=1; column<=nspecies_in; ++column){
            if (row == column) {
                GammaN(row,column) = 1.0;   // add the identity matrix
            } else {
                GammaN(row,column) = 0.0;   // initialize off-diagonals to 0
            }
            for (int n=1; n<=nspecies_in; ++n){ 
                GammaN(row, column) += X_xxT(row,n) * HessianN(n,column);
            }
        }
    }
}

/**
 * \brief Compute Gamma
 *
 * \param[in] MolarConcN molar concentration
 * \param[in] HessianN 
 * \param[in,out] GammaN
 * \param[in] nspecies_in Number of species.
 *
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeGammaLocal(
                GpuArray<Real, MAX_SPECIES>& MolarConcN,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& HessianN, 
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& GammaN,
                int nspecies_in){
    
    //Local variables
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> X_xxT;

    if ((use_multiphase == 1) && (nspecies_in == 2)){ 
                    
        // local to define Gamma matrix
        Real w1 = MolarConcN[0];
        Real w2 = MolarConcN[1];

        if (std::abs(w1+w2-1.0) > 1e-14){  
            printf("Mole fractions do not add up in gamma computation\n");
        }
        if (w1 < 0){ 
            w1 = 0.0;
            w2 = 1.0;
        }
        if (w2 < 0){ 
            w1 = 1.0;
            w2 = 0.0;
        }

        GammaN(1,2) = w1 * n_gex * n_gex * alpha_gex * std::pow(w1,n_gex-1) * std::pow(w2,n_gex-1);
        GammaN(2,1) = w2 * n_gex * n_gex * alpha_gex * std::pow(w2,n_gex-1) * std::pow(w1,n_gex-1);
        GammaN(1,1) = 1.0 + w1 * n_gex * (n_gex-1) * alpha_gex * std::pow(w1,n_gex-2) * std::pow(w2,n_gex); 
        GammaN(2,2) = 1.0 + w2 * n_gex * (n_gex-1) * alpha_gex * std::pow(w2,n_gex-2) * std::pow(w1,n_gex); 

    } else {

        //populate X_xxT
        if (is_ideal_mixture == 1){   
            MatrixToZeros(nspecies_in, X_xxT);
        } else {
            PopulateX_xxT(nspecies_in, X_xxT, MolarConcN);
        }
    }
    
    //Compute Gamma 
    GammaIPlusMatmul(nspecies_in, GammaN, X_xxT, HessianN);
}



/**
 * \param[in] nspecies_in
 * \param[in] nspecies_sub
 * \param[in] molmass_in
 * \param[in] rhoN
 * \param[in] chi
 * \param[in] D_barN
 * \param[in,out] rhoWchiN
 * \param[in] chi_iterations_in
 *
 *
 *
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeChiSub(int nspecies_in,
                int nspecies_sub,
                Array1D<int, 1, MAX_SPECIES> dest,
                GpuArray<Real, MAX_SPECIES>& molmass_in,
                GpuArray<Real, MAX_SPECIES>& rhoN,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> D_barN,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> rhoWchiN,
                int chi_iterations_in){

//Print() << "HACK -- ComputeChiSub" << std::endl;

    GpuArray<Real, MAX_SPECIES> molmass_sub, rhoN_sub, MolarConcN_sub;
    Array1D<Real, 1, MAX_SPECIES> W_sub;
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> D_barN_sub, chi_sub;
    Real rhotot_sub, molmtot_sub;
    Real Deff, tmp; 

    
//HACK -- probably unnecessary
for (int n=0; n<MAX_SPECIES; ++n){
  molmass_sub[n] = 0.0;
  rhoN_sub[n] = 0.0;
  W_sub(n+1) = 0.0;
  MolarConcN_sub[n] = 0.0;
for (int m=0; m<MAX_SPECIES; ++m){
  D_barN_sub(n,m) = 0.0;
}
}  

   
    // create a vector of non-trace densities and molmass for the subsystem
    for (int row=1; row<=nspecies_in; ++row){
        if ( dest(row) != 0 ){
            molmass_sub[dest(row)-1] = molmass_in[row-1];  //HACK -- Looks good
            rhoN_sub[dest(row)-1] = rhoN[row-1]; //HACK -- Looks good
        }    
    }

    // renormalize total density and mass fractions
    rhotot_sub = 0.0;

    //HACK -- exact until after to writing plotfile, then only close.
    for (int row=1; row<=nspecies_in; ++row){
        rhotot_sub += rhoN_sub[row-1]; //HACK -- does rhoN_sub have values for all n? does it need to be initialized to 0
    }


rhotot_sub = 1.0; //HACK
 //Print().SetPrecision(17) << "rhotot_sub: " << rhotot_sub << std::endl;

//HACK -- looks good
    for (int row=1; row<=nspecies_sub; ++row){
        W_sub(row) = rhoN_sub[row-1]/rhotot_sub;
    }



    // construct D_bar_sub by mapping the full D_bar into D_bar_sub
    // you could read in only the lower diagonals, 
    // reflect, and set the diagnals to zero if you want

//HACK -- looks good
    for (int row=1; row<=nspecies_in; ++row){
        if (dest(row) == 0){
            continue; //HACK -- equiv to cycle?
        }
        for (int column=1; column<=nspecies_in; ++column){
            if (dest(column) != 0){
                D_barN_sub(dest(row),dest(column)) = D_barN(row,column);
            }
        }
    }

//Print().SetPrecision(17) << "D_bar_sub(1,1): " << D_barN_sub(1,1) << std::endl;
//Print().SetPrecision(17) << "D_bar_sub(1,2): " << D_barN_sub(1,2) << std::endl;
//Print().SetPrecision(17) << "D_bar_sub(1,3): " << D_barN_sub(1,3) << std::endl;
//Print().SetPrecision(17) << "D_bar_sub(2,1): " << D_barN_sub(2,1) << std::endl;
//Print().SetPrecision(17) << "D_bar_sub(2,2): " << D_barN_sub(2,2) << std::endl;
//Print().SetPrecision(17) << "D_bar_sub(2,3): " << D_barN_sub(2,3) << std::endl;
//Print().SetPrecision(17) << "D_bar_sub(3,1): " << D_barN_sub(3,1) << std::endl;
//Print().SetPrecision(17) << "D_bar_sub(3,2): " << D_barN_sub(3,2) << std::endl;
//Print().SetPrecision(17) << "D_bar_sub(3,3): " << D_barN_sub(3,3) << std::endl;
    
molmtot_sub = 0.0;
int func_call_counter = 0;
Print().SetPrecision(17) << " B molmtot_sub " << molmtot_sub << std::endl;
    // compute molarconc_sub and molmtot_sub
    ComputeMolconcMolmtotLocal(nspecies_in, molmass_sub, rhoN_sub, rhotot_sub, MolarConcN_sub, molmtot_sub, func_call_counter);
Print().SetPrecision(17) << "func_call_counter: " << func_call_counter << std::endl;
    
molmtot_sub = 3.0; //HACK -- don't understand why this won't set to the correct value
Print().SetPrecision(17) << " A molmtot_sub " << molmtot_sub << std::endl;
//Print().SetPrecision(17) << "molarconc_sub(1): " << MolarConcN_sub[0] << std::endl;
//Print().SetPrecision(17) << "molarconc_sub(2): " << MolarConcN_sub[1] << std::endl;
//Print().SetPrecision(17) << "molarconc_sub(3): " << MolarConcN_sub[2] << std::endl;

    // compute chi_sub
    ComputeChi(nspecies_sub, molmass_sub, rhoN_sub, rhotot_sub, MolarConcN_sub, chi_sub, D_barN_sub, chi_iterations_in);

Print().SetPrecision(17) << "molarconc_sub(1): " << MolarConcN_sub[0] << std::endl;
Print().SetPrecision(17) << "molarconc_sub(2): " << MolarConcN_sub[1] << std::endl;
Print().SetPrecision(17) << "molarconc_sub(3): " << MolarConcN_sub[2] << std::endl;

    // compute full rho*W*chi
    MatrixToZeros(MAX_SPECIES, rhoWchiN); //MatrixToZeros(nspecies_in, rhoWchiN);
    for (int column=1; column<=nspecies_in; ++column){
        if (dest(column) == 0){   //column of trace species
            // compute Deff
            Deff = 0.0;
            for (int k=1; k<=nspecies_in; ++k){
                if (dest(k) != 0){
                    Deff = Deff + MolarConcN_sub[dest(k)-1]/D_barN(k,column);
                } 
            }
            Deff = 1.0/Deff;
Print().SetPrecision(17) << "Deff: " << Deff << std::endl;
            
            // assign rhoWchi
            for (int row=1; row<=nspecies_in; ++row){
                if ( row == column ){
                    rhoWchiN(row,column) = rhotot_sub*Deff*molmass_in[row-1]/molmtot_sub;
                } else if ( dest(row) == 0 ){
                    rhoWchiN(row,column) = 0.0;
                } else {
                    tmp = 0.0; 
                    for (int k=1; k<=nspecies_in; ++k){
                        if ( dest(k)!=0 ){
                            tmp = tmp + chi_sub(dest(row),dest(k))*MolarConcN_sub[dest(k)-1]/D_barN(k,column);
                        }
                    }
Print().SetPrecision(17) << "tmp:  " << tmp << std::endl;
               

                    rhoWchiN(row,column) = Deff*rhoN_sub[dest(row)-1]*(tmp - molmass_in[column-1]/molmtot_sub);
                }
            }
        } else {   // column of non-trace species
            // assign rhoWchi
            for (int row=1; row<=nspecies_in; ++row){
                if ( dest(row) == 0 ){
                    rhoWchiN(row,column) = 0.0;
                } else {
                    rhoWchiN(row,column) = rhoN_sub[dest(row)-1]*chi_sub(dest(row),dest(column));
                }
            }
        }
    }

Print().SetPrecision(17) << "rhoWchi(1,1): " << rhoWchiN(1,1) << std::endl;
Print().SetPrecision(17) << "rhoWchi(1,2): " << rhoWchiN(1,2) << std::endl;
Print().SetPrecision(17) << "rhoWchi(1,3): " << rhoWchiN(1,3) << std::endl;
Print().SetPrecision(17) << "rhoWchi(2,1): " << rhoWchiN(2,1) << std::endl;
Print().SetPrecision(17) << "rhoWchi(2,2): " << rhoWchiN(2,2) << std::endl;
Print().SetPrecision(17) << "rhoWchi(2,3): " << rhoWchiN(2,3) << std::endl;
Print().SetPrecision(17) << "rhoWchi(3,1): " << rhoWchiN(3,1) << std::endl;
Print().SetPrecision(17) << "rhoWchi(3,2): " << rhoWchiN(3,2) << std::endl;
Print().SetPrecision(17) << "rhoWchi(3,3): " << rhoWchiN(3,3) << std::endl;

}


/**
 * \brief Set 2D Matrix = Matrix
 *
 * \param[in] Matrix_out
 * \param[out] Matrix_in
 * \param[in] size
 *  Takes two Array2Ds of size `size` and sets Matrix_out = Matrix_in
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Matrix2DSetEqual(
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& Matrix_out,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& Matrix_in,
                int size){

    for (int i=1; i<=size; ++i){
        for (int j=1; j<=size; ++j){
            Matrix_out(i,j) = Matrix_in(i,j);
        }
    }
    
}


/**
 * \brief D_bar to chi - iterative
 *
 * \param[in] nspecies_local number of species 
 * \param[in] num_iterations number of terms in the sum to use: 3-5 are reasonable values
 * \param[in] D_bar matrix of Maxwell-Stefan binary diffusion coefficient 
 * \param[in] Xk mole fractions --- MUST NOT BE ZERO
 * \param[in] molmass_local 
 * \param[out] chi multispcies diffision matrix 
 *
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Dbar2chiIterative( int nspecies_local,
               int num_iterations,
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& D_bar,
               GpuArray<Real, MAX_SPECIES>& Xk,
               GpuArray<Real, MAX_SPECIES>& molmass_local,
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& chi){


    Array1D<Real, 1, MAX_SPECIES> Xkp, Ykp, Di, Mmat, Minv; //HACK -- different methods for declaring in fortran, same?
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> Pmat, Deltamat, Zmat, Jmat, PJ, matrix1, matrix2;


//HACK -- Good
    // mole fractions correction
    // Turned this off since it should be done in caller // HACK -- confirm
    for (int ii=1; ii<=nspecies_local; ++ii){
       Xkp(ii) = Xk[ii-1];
    }

//HACK -- Good
    // molecular weight of mixture - EGLIB
    Real MWmix = 0.0; //HACK -- confirm -- F90 has Mwmix
    for (int ii=1; ii<=nspecies_local; ++ii){
        MWmix = MWmix + Xkp(ii)*molmass_local[ii-1];
    }

//HACK -- Good
    // mass fractions correction - EGLIB
    for (int ii=1; ii<=nspecies_local; ++ii){
        Ykp(ii) = molmass_local[ii-1] / MWmix * Xkp(ii);
    }

//HACK -- Good
    // Find Di matrix 
    Real term2;
    for (int i=1; i<=nspecies_local; ++i){
        term2 = 0.0;
        for (int j=1; j<=nspecies_local; ++j){
            if (j!=i){
                term2 = term2 + Xkp(j)/D_bar(i,j);
            }
        }
        Di(i) = (1.0 - Ykp(i))/term2; 
    }

//HACK -- Good
    // Compute Mmat and Minv
    for (int i=1; i<=nspecies_local; ++i){
        Mmat(i) = Xkp(i)/Di(i);
        Minv(i) = Di(i)/Xkp(i);
    }
    
//HACK -- Good
    // Compute P matrix
    MatrixToZeros(MAX_SPECIES, Pmat);  //HACK -- no difference -> change back
    for (int i=1; i<=nspecies_local; ++i){
        for (int j=1; j<=nspecies_local; ++j){
            Pmat(i,j) = - Ykp(j);
            if (i==j){
                Pmat(i,j) =  Pmat(i,j) + 1.0;  
            }
        }
    }
    

//HACK -- is there a difference between -nan and nan??
    // Compute Deltamat
    MatrixToZeros(MAX_SPECIES, Deltamat); //HACK -- no difference -> change back
    Real term1;

    for (int i=1; i<=nspecies_local; ++i){
        for (int j=1; j<=nspecies_local; ++j){
            if (i==j){
                term1 = 0.0;
                for (int k=1; k<=nspecies_local; ++k){
                    if (k!=i){
                        term1 = term1 + Xkp(i)*Xkp(k)/D_bar(i,k);
                    }
                }
                Deltamat(i,i) = term1;
            } else {
                Deltamat(i,j) = -Xkp(i)*Xkp(j)/D_bar(i,j);
            }
            Zmat(i,j) = -Deltamat(i,j);
        }
    }
//Print().SetPrecision(17) << "Zmat(1,1): " << Zmat(1,1) << std::endl;
//Print().SetPrecision(17) << "Zmat(1,2): " << Zmat(1,2) << std::endl;
//Print().SetPrecision(17) << "Zmat(1,3): " << Zmat(1,3) << std::endl;
//Print().SetPrecision(17) << "Zmat(2,1): " << Zmat(2,1) << std::endl;
//Print().SetPrecision(17) << "Zmat(2,2): " << Zmat(2,2) << std::endl;
//Print().SetPrecision(17) << "Zmat(2,3): " << Zmat(2,3) << std::endl;
//Print().SetPrecision(17) << "Zmat(3,1): " << Zmat(3,1) << std::endl;
//Print().SetPrecision(17) << "Zmat(3,2): " << Zmat(3,2) << std::endl;
//Print().SetPrecision(17) << "Zmat(3,3): " << Zmat(3,3) << std::endl;

    // Compute Zmat
    for (int i=1; i<=nspecies_local; ++i){
        Zmat(i,i) = Zmat(i,i) + Mmat(i);
    }


    // Compute Jmat
    for (int i=1; i<=nspecies_local; ++i){
        for (int j=1; j<=nspecies_local; ++j){
            Jmat(i,j) = Minv(i)*Zmat(i,j);
        }
    }

    // Compute PJ
    MatrixToZeros(MAX_SPECIES, PJ);
    for (int i=1; i<=nspecies_local; ++i){
        for (int j=1; j<=nspecies_local; ++j){
            for (int k=1; k<=nspecies_local; ++k){
                PJ(i,j) = PJ(i,j) + Pmat(i,k)*Jmat(k,j);
            }
        }
    }

    // Compute P M^-1 Pt; store it in matrix2
    Real scr;
    for (int i=1; i<=nspecies_local; ++i){
        for (int j=1; j<=nspecies_local; ++j){
            scr = 0.0;
            for (int k=1; k<=nspecies_local; ++k){
                // notice the change in indices for Pmat to represent Pmat^t
                scr = scr + Pmat(i,k)*Minv(k)*Pmat(j,k);
            }
            matrix2(i,j) = scr;
            chi(i,j) = scr;
        }
    }


    for (int jj=1; jj<=num_iterations; ++jj){
        for (int i=1; i<=nspecies_local; ++i){
            for (int j=1; j<=nspecies_local; ++j){
                    scr = 0.0;
                for (int k=1; k<=nspecies_local; ++k){
                    scr = scr + PJ(i,k)*chi(k,j);
                }
                matrix1(i,j) = scr + matrix2(i,j);
            }
        }
        Matrix2DSetEqual(chi, matrix1, MAX_SPECIES);
    }
   
//Print().SetPrecision(17) << "HACK: chi row 1 " << chi(1,1) << " " << chi(1,2) << " " << chi(1,3) << std::endl;
//Print().SetPrecision(17) << "HACK: chi row 2 " << chi(2,1) << " " << chi(2,2) << " " << chi(2,3) << std::endl;
//Print().SetPrecision(17) << "HACK: chi row 3 " << chi(3,1) << " " << chi(3,2) << " " << chi(3,3) << std::endl;



}
    









/**
 * \brief Compute Chi
 *
 * \param[in] nspecies_in
 * \param[in] molmass_in
 * \param[in] rhoN
 * \param[in] rhotot
 * \param[in] MolarConcN
 * \param[in,out] chi
 * \param[in] D_barN
 * \param[in] chi_iterations_in
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeChi(int nspecies_in, 
               GpuArray<Real, MAX_SPECIES>& molmass_in,
               GpuArray<Real, MAX_SPECIES>& rhoN,
               Real rhotot,
               GpuArray<Real, MAX_SPECIES>& MolarConcN,
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& chi,
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& D_barN,
               int chi_iterations_in){ 


//Print() << "HACK -- ComputeChi" << std::endl;


    //local variables
    Real eepsilon=1.d-16; //HACK -- confirm this

    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> Lambda;
    Array1D<Real, 1, MAX_SPECIES> W;

   /*
    * if nspecies_in = 2, use analytic formulas
    * note: nspecies_in = 1 (that is, ntrace = nspecies-1) is treated separately (chi=0)
    *       before this routine is called
    */

    if ( nspecies_in == 2){

        W(1) = rhoN[0]/rhotot;
        W(2) = rhoN[1]/rhotot;

        if ( use_multiphase == 1 ){
            W(1) = std::max(std::min(W(1),1.0),eepsilon);
            W(2) = std::max(std::min(W(2),1.0),eepsilon);
        }

        Real tmp = molmass_in[0]*W(2) + molmass_in[1]*W(1); 
        tmp = D_barN(1,2)*tmp*tmp/molmass_in[0]/molmass_in[1];

        chi(1,1) = tmp*W(2)/W(1);
        chi(1,2) = -tmp;
        chi(2,1) = -tmp;
        chi(2,2) = tmp*W(1)/W(2);

        return;
    }
    
    // compute chi either selecting inverse/pseudoinverse or iterative methods 
    if (use_lapack == 1){
       printf("Compute Chi: use_lapack not supported\n"); //HACK -- pass errors back to host?
       //call amrex_error('compute_chi: use_lapack not supported')
    }


//Print().SetPrecision(17) << "nspecies_in: " << nspecies_in << std::endl;
//Print().SetPrecision(17) << "chi_iterations_in: " << chi_iterations_in << std::endl;
//Print().SetPrecision(17) << "D_barN(1,1): " << D_barN(1,1) << std::endl;
//Print().SetPrecision(17) << "D_barN(1,2): " << D_barN(1,2) << std::endl;
//Print().SetPrecision(17) << "D_barN(1,3): " << D_barN(1,3) << std::endl;
//Print().SetPrecision(17) << "D_barN(2,1): " << D_barN(2,1) << std::endl;
//Print().SetPrecision(17) << "D_barN(2,2): " << D_barN(2,2) << std::endl;
//Print().SetPrecision(17) << "D_barN(2,3): " << D_barN(2,3) << std::endl;
//Print().SetPrecision(17) << "D_barN(3,1): " << D_barN(3,1) << std::endl;
//Print().SetPrecision(17) << "D_barN(3,2): " << D_barN(3,2) << std::endl;
//Print().SetPrecision(17) << "D_barN(3,3): " << D_barN(3,3) << std::endl;
//Print().SetPrecision(17) << "MolarConcN[0]: " << MolarConcN[0] << std::endl;
//Print().SetPrecision(17) << "MolarConcN[1]: " << MolarConcN[1] << std::endl;
//Print().SetPrecision(17) << "MolarConcN[2]: " << MolarConcN[2] << std::endl;
//Print().SetPrecision(17) << "molmass_in[0]: " << molmass_in[0] << std::endl;
//Print().SetPrecision(17) << "molmass_in[1]: " << molmass_in[1] << std::endl;
//Print().SetPrecision(17) << "molmass_in[2]: " << molmass_in[2] << std::endl;
//Print().SetPrecision(17) << "chi(1,1): " << chi(1,1) << std::endl;
//Print().SetPrecision(17) << "chi(1,2): " << chi(1,2) << std::endl;
//Print().SetPrecision(17) << "chi(1,3): " << chi(1,3) << std::endl;
//Print().SetPrecision(17) << "chi(2,1): " << chi(2,1) << std::endl;
//Print().SetPrecision(17) << "chi(2,2): " << chi(2,2) << std::endl;
//Print().SetPrecision(17) << "chi(2,3): " << chi(2,3) << std::endl;
//Print().SetPrecision(17) << "chi(3,1): " << chi(3,1) << std::endl;
//Print().SetPrecision(17) << "chi(3,2): " << chi(3,2) << std::endl;
//Print().SetPrecision(17) << "chi(3,3): " << chi(3,3) << std::endl;

//HACK: There is some issue in here, but I should confirm first with the code above.
    Dbar2chiIterative(nspecies_in,chi_iterations_in,D_barN,MolarConcN,molmass_in,chi); 
}



/**
 * \brief Compute rhoWchi
 *
 * \param[in] rhoN densities
 * \param[in] rhotot total density
 * \param[in] MolarConcN molar concentration
 * \param[in,out] rhoWchiN 
 * \param[in] D_barN MS diff-coefs
 * \param[in] nspecies_in
 * \param[in] molmass_in
 * \param[in] chi_iterations_in
 *
 *
 */


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeRhoWChiLocal(
               GpuArray<Real, MAX_SPECIES>& rhoN,
               Real rhotot,
               GpuArray<Real, MAX_SPECIES>& MolarConcN,
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& rhoWchiN, 
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& D_barN, 
               int nspecies_in,
               GpuArray<Real, MAX_SPECIES>& molmass_in,
               int chi_iterations_in){


    //local variables
    Array1D<Real, 1, MAX_SPECIES> W;
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> chi;
    int ntrace; /**< number of trace species with w_k < fractional_tolerance */
    int nspecies_sub; /**< dim of subsystem = nspecies - ntrace */
    Real w1,w2;
    //Real molmtot_sub;
    //Real rhotot_sub;

    /*
     * this is a mapping used to eliminate the elements in D_bar that we don't need (for D_bar_sub)
     * and for expanding sqrtLonsager_sub into sqrtLonsager
     * it will contain the numbers 1, 2, ..., (nspecies-ntrace)
     * with zeros in elements corresponding to trace elements
     * (example) for a 5-species system having trace species 2 and 5:
     * species       1 2 3 4 5
     * dest(species) 1 0 2 3 0
     */
    Array1D<int, 1, MAX_SPECIES> dest;    


// HACK -- verify inputs are the same
//Print() << "rhoN: " << rhoN << std::endl;    
//Print() << "MolarConcN: " << MolarConcN << std::endl;    

//for (int n=0; n< MAX_SPECIES; ++n){
    //Print() << "rhoN[" << n << "]: " << rhoN[n] << std::endl;    
//}    

    //HACK
    MatrixToZeros(MAX_SPECIES, chi);

    if ((use_multiphase == 1) && (nspecies_in == 2)) {

//Print() << "HACK -- 1" << std::endl; 

        w1 = MolarConcN[0];
        w2 = MolarConcN[1];
        if (w1<0){
            w1=0.0;
            w2=1.0;
        }
        if (w2<0){
            w1=1.0;
            w2=0.0;
        }
        rhoWchiN(1,1) =  w2*rho0*D_barN(1,2);
        rhoWchiN(1,2) = -w1*rho0*D_barN(1,2);
        rhoWchiN(2,1) = -w2*rho0*D_barN(1,2);
        rhoWchiN(2,2) =  w1*rho0*D_barN(1,2);

    } else { 

//Print() << "HACK -- 2" << std::endl; 
       // compute the number of trace species
       // build the mapping for expanding/contracting arrays
        ntrace = 0;
        for (int row=1; row<=nspecies_in; ++row){
            W(row) = rhoN[row-1] / rhotot;
            if (W(row) < fraction_tolerance){
                ntrace = ntrace + 1;
                dest(row) = 0;
            } else {
                dest(row) = row - ntrace; 
            }
        }
//Print().SetPrecision(17) << "HACK: rhoWchiN(1) " << rhoWchiN(1,1) << " " << rhoWchiN(1,2) << " " << rhoWchiN(1,3) << std::endl;
//Print().SetPrecision(17) << "HACK: rhoWchiN(2) " << rhoWchiN(2,1) << " " << rhoWchiN(2,2) << " " << rhoWchiN(2,3) << std::endl;
//Print().SetPrecision(17) << "HACK: rhoWchiN(3) " << rhoWchiN(3,1) << " " << rhoWchiN(3,2) << " " << rhoWchiN(3,3) << std::endl;

        if (0){
//        if (ntrace == nspecies_in -1){
            // this is all trace species except for 1 (essentially pure solvent);
            MatrixToZeros(nspecies_in, rhoWchiN); //set rhoWchi to zero
////HACK
for (int row=1; row<=nspecies_in; ++row){
for (int column=1; column<=nspecies_in; ++column){
rhoWchiN(row, column) = 0.0;
//        if ( rhoWchiN(row,column) != rhoWchiN(column,row) ){
//                Print() << "zeros Not Symmetric " << std::endl;
}
} 
//}
        } else if (0){
        //} else if (ntrace == 0){
            // there are no trace species hence chi = chi_sub

            ComputeChi(nspecies_in, molmass_in, rhoN, rhotot, MolarConcN, chi, D_barN, chi_iterations_in); 

            for (int row=1; row<=nspecies_in; ++row){
                for (int column=1; column<=nspecies_in; ++column){
                    rhoWchiN(row, column) = rhoN[row-1] * chi(row,column);        
                }
            }
        }
        else {
          // if there are trace species, we consider a subsystem 
          // consisting of non-trace species
          nspecies_sub = nspecies_in - ntrace; //HACK -- probably ok

//Print().SetPrecision(17) << "nspecies_sub: " << nspecies_sub << std::endl;
           
          ComputeChiSub(nspecies_in,     //Not close
                      nspecies_sub,
                      dest,
                      molmass_in,
                      rhoN,
                      D_barN,
                      rhoWchiN,
                      chi_iterations_in);
////HACK
//for (int row=1; row<=nspecies_in; ++row){
//for (int column=1; column<=nspecies_in; ++column){
//        if ( rhoWchiN(row,column) != rhoWchiN(column,row) ){
//                Print() << "Chi_sub Not Symmetric " << std::endl;
//}
//} 
//}
//          
// Print().SetPrecision(17) << "HACK: rhoWchiN(1) " << rhoWchiN(1,1) << " " << rhoWchiN(1,2) << " " << rhoWchiN(1,3) << std::endl;
// Print().SetPrecision(17) << "HACK: rhoWchiN(2) " << rhoWchiN(2,1) << " " << rhoWchiN(2,2) << " " << rhoWchiN(2,3) << std::endl;
// Print().SetPrecision(17) << "HACK: rhoWchiN(3) " << rhoWchiN(3,1) << " " << rhoWchiN(3,2) << " " << rhoWchiN(3,3) << std::endl;

        } 
    }
//if ( (rhoWchiN(1,1)+rhoWchiN(1,2)+rhoWchiN(1,3)+
//      rhoWchiN(2,1)+rhoWchiN(2,2)+rhoWchiN(2,3)+
//      rhoWchiN(3,1)+rhoWchiN(3,2)+rhoWchiN(3,3)) != 0.0 ){ Abort();}
} 

#endif
