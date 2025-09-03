#ifndef _multispec_functions_H_
#define _multispec_functions_H_

#include <AMReX_MultiFab.H>

#include "multispec_namespace.H"

#include "common_functions.H"

#include "StochMassFlux.H"
#include "StochMomFlux.H"

#include <limits>

using namespace multispec;
using namespace amrex;

/////////////////////////////////////////////////////////////////////////////////
// in multispec_functions.cpp

void InitializeMultispecNamespace();

/////////////////////////////////////////////////////////////////////////////////
// in ComputeDivReversibleStress.cpp
void ComputeDivFHReversibleStress(std::array<MultiFab,AMREX_SPACEDIM>& div_reversible_stress,
                                const MultiFab& rhotot_in,
                                MultiFab& rho_in,
                                const Geometry& geom);


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

void ComputeFHHigherOrderTerm(const MultiFab& molarconc,
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

void ComputeEta(const MultiFab& rho_in,
                const MultiFab& rhotot_in,
                MultiFab& eta_in);

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

void ComputeMassfrac(const MultiFab& rho,
               const MultiFab& rhotot,
               MultiFab& massfrac);

void ComputeGamma(const MultiFab& molarconc,
          const MultiFab& Hessian,
          MultiFab& Gamma);

void ComputeFHGamma(const MultiFab& massfrac,
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


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeChi(int nspecies_in,
                GpuArray<Real, MAX_SPECIES>& molmass_in,
                GpuArray<Real, MAX_SPECIES>& rhoN,
                Real rhotot,
                GpuArray<Real, MAX_SPECIES>& MolarConcN,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& chi,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& D_barN);

/**
 * \brief Compute total molar mass
 *
 * \param[in] W Mass fractions.
 * \param[in] molmass_in
 * \param[in] nspecies_in Number of species
 * \param[out] molmtot Total molar mass
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetMolTot ( GpuArray<Real, MAX_SPECIES>& W,
                GpuArray<Real, MAX_SPECIES>& molmass_in,
                Real& molmtot_in,
                int nspecies_in)
{
    Real sumWOverN = 0.0;

    for (int n=0; n<nspecies_in; ++n){
        sumWOverN += W[n] / molmass_in[n];
    }

    molmtot_in = 1.0 / sumWOverN;
}

/**
 * \brief Compute total molar mass
 *
 * \param[in] W Mass fractions.
 * \param[in] molmass_in
 * \param[in] nspecies_in Number of species
 * \param[out] molmtot Total molar mass
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GetMonomerTot ( GpuArray<Real, MAX_SPECIES>& W,
                GpuArray<Real, MAX_SPECIES>& monomer_in,
                Real& monomertot,
                int nspecies_in)
{
    Real sumWOverN = 0.0;

    for (int n=0; n<nspecies_in; ++n){
        sumWOverN += W[n] / monomer_in[n];
    }

    monomertot = 1.0 / sumWOverN;
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
void ComputeMolconcMolmtotLocal ( int nspecies_in,
                GpuArray<Real, MAX_SPECIES>& molmass_in,
                GpuArray<Real, MAX_SPECIES>& RhoN,
                Real rhotot_in,
                GpuArray<Real, MAX_SPECIES>& MolarConcN,
                Real& molmtot_in){

    GpuArray<Real, MAX_SPECIES> W; /**< Mass fraction w_i = rho_i/rho */

    // calculate mass fraction and total molar mass (1/m=Sum(w_i/m_i))
    for (int n=0; n<nspecies_in; ++n){   //HACK: Ask Andy about nspecies -> nspecies_in
        W[n] = RhoN[n] / rhotot_in;
    }

    GetMolTot(W, molmass_in, molmtot_in, nspecies_in);

    // calculate molar concentrations in each cell (x_i=m*w_i/m_i)
    for (int n=0; n<nspecies_in; ++n){
        MolarConcN[n] = molmtot_in * W[n] / molmass_in[n];
    }

}


/**
 * \brief Set all values in 2D Array to zero.
 *
 * \param[in] size upper limit of Array2D index (1 to size) to set to zeros.
 * \param[in,out] matrixIn Array2D to be set to zeros.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void MatrixToZeros (int size,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& matrixIn) {
//HACK:required

    for (int n=1; n<=size; ++n ){
        for (int m=1; m<=size; ++m ){
            matrixIn(n,m) = 0.0;
        }
    }
}

/**
 * \brief Populate X_xxT
 *
 * \param[out] X_xxT
 * \param[in] MolarConcN
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void PopulateX_xxT (
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& X_xxT,
                GpuArray<Real, MAX_SPECIES>& MolarConcN) {

    for (int row=1; row<=nspecies;++row ){
        // diagonal entries
        X_xxT(row,row) = MolarConcN[row-1] - MolarConcN[row-1]*MolarConcN[row-1];
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
 * \brief Compute Gamma = I + matmul(X_xxT, Hessian).
 *
 * Adds the identity matrix during initialization. Then adds product of X_xxT and Hessian matrices.
 *
 * \param[out] GammaN
 * \param[in] X_xxT
 * \param[in] HessianN
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void GammaIPlusMatmul (
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& GammaN,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& X_xxT,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& HessianN) {

    //Compute Gamma
    for (int row=1; row<=nspecies; ++row){
        for (int column=1; column<=nspecies; ++column){
            if (row == column) {
                GammaN(row,column) = 1.0;   // add the identity matrix
            } else {
                GammaN(row,column) = 0.0;   // initialize off-diagonals to 0
            }
            for (int n=1; n<=nspecies; ++n){
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
 *
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeGammaLocal (
                GpuArray<Real, MAX_SPECIES>& MolarConcN,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& HessianN,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& GammaN){

    //Local variables
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> X_xxT;

    if ((use_multiphase == 1) && (nspecies == 2)){

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

        if( n_gex == 1 ){

           GammaN(1,2) = w1 * n_gex * n_gex * alpha_gex * std::pow(w1,n_gex-1) * std::pow(w2,n_gex-1);
           GammaN(2,1) = w2 * n_gex * n_gex * alpha_gex * std::pow(w2,n_gex-1) * std::pow(w1,n_gex-1);
           GammaN(1,1) = 1.0 ;
           GammaN(2,2) = 1.0 ;

        } else {

           GammaN(1,2) = w1 * n_gex * n_gex * alpha_gex * std::pow(w1,n_gex-1) * std::pow(w2,n_gex-1);
           GammaN(2,1) = w2 * n_gex * n_gex * alpha_gex * std::pow(w2,n_gex-1) * std::pow(w1,n_gex-1);
           GammaN(1,1) = 1.0 + w1 * n_gex * (n_gex-1) * alpha_gex * std::pow(w1,n_gex-2) * std::pow(w2,n_gex);
           GammaN(2,2) = 1.0 + w2 * n_gex * (n_gex-1) * alpha_gex * std::pow(w2,n_gex-2) * std::pow(w1,n_gex);

        }

    } else {

        //populate X_xxT
        if (is_ideal_mixture == 1){
            MatrixToZeros(nspecies, X_xxT);
        } else {
            PopulateX_xxT(X_xxT, MolarConcN);
        }

        //Compute Gamma
        GammaIPlusMatmul(GammaN, X_xxT, HessianN);

    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeFHGammaLocal (
                GpuArray<Real, MAX_SPECIES>& massFrac,
                Array2D<Real, 0, MAX_SPECIES-1, 0, MAX_SPECIES-1>& GammaN){

        Real Nbar;

        GetMonomerTot ( massFrac, fh_monomers, Nbar, nspecies);

        for (int n=0; n<nspecies; n++){
           for (int m=0; m<nspecies; m++){
              Real wn = std::max(std::min( massFrac[n],1.),0.);
              GammaN(n,m) = Nbar*wn*fh_chi(n,m);
           }
           GammaN(n,n) += Nbar/ fh_monomers[n];
        }

}



/**
 * \brief Compute rhoWchi using a subsystem of species above the fractional weight tolerance.
 *
 * \param[in] nspecies_sub Number of species in the subsystem (= nspecies - ntrace).
 * \param[in] dest Mapping between system and subsystem
 * \param[in] molmass_in
 * \param[in] rhoN
 * \param[in] D_barN
 * \param[in,out] rhoWchiN
 *
 *
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeChiSub(int nspecies_sub,
                   Array1D<int, 1, MAX_SPECIES> dest,
                   GpuArray<Real, MAX_SPECIES>& molmass_in,
                   GpuArray<Real, MAX_SPECIES>& rhoN,
                   Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& D_barN,
                   Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& rhoWchiN){

    GpuArray<Real, MAX_SPECIES> molmass_sub, rhoN_sub, MolarConcN_sub;
    Array1D<Real, 1, MAX_SPECIES> W_sub;
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> D_barN_sub, chi_sub;
    Real rhotot_sub, molmtot_sub;
    Real Deff, tmp;

    // create a vector of non-trace densities and molmass for the subsystem
    for (int row=1; row<=nspecies; ++row){
        if ( dest(row) != 0 ){
            molmass_sub[dest(row)-1] = molmass_in[row-1];
            rhoN_sub[dest(row)-1] = rhoN[row-1];
        }
    }

    // renormalize total density and mass fractions
    rhotot_sub = 0.0;
    for (int row=1; row<=nspecies_sub; ++row){
        rhotot_sub += rhoN_sub[row-1];
    }

    for (int row=1; row<=nspecies_sub; ++row){
        W_sub(row) = rhoN_sub[row-1]/rhotot_sub;
    }

    // construct D_bar_sub by mapping the full D_bar into D_bar_sub
    // you could read in only the lower diagonals,
    // reflect, and set the diagnals to zero if you want

    for (int row=1; row<=nspecies; ++row){
        if (dest(row) == 0){
            continue;
        }
        for (int column=1; column<=nspecies; ++column){
            if (dest(column) != 0){
                D_barN_sub(dest(row),dest(column)) = D_barN(row,column);
            }
        }
    }

    // compute molarconc_sub and molmtot_sub
    ComputeMolconcMolmtotLocal(nspecies_sub, molmass_sub, rhoN_sub, rhotot_sub, MolarConcN_sub, molmtot_sub);

    // compute chi_sub
//      amrex::Print() << " entering reduced system "  << nspecies_sub << std::endl;
//      for (int n=1; n<=nspecies_sub; ++n){
//          amrex::Print() << "n , dest,molmass rho massf " << n << " " << dest(n) << " " << molmass_sub[n-1] << " " <<
//              rhoN_sub[n-1] << " " << MolarConcN_sub[n-1] << std::endl;
//              for (int m=1; m<=nspecies_sub; ++m){
//              amrex::Print() << "m d chi " << m << " " << D_barN_sub(n,m) << " " <<  chi_sub(n,m) << std::endl;
//          }
//      }
    ComputeChi(nspecies_sub, molmass_sub, rhoN_sub, rhotot_sub, MolarConcN_sub, chi_sub, D_barN_sub);


    // compute full rho*W*chi
    MatrixToZeros(nspecies, rhoWchiN);
    for (int column=1; column<=nspecies; ++column){
        if (dest(column) == 0){   //column of trace species
            // compute Deff
            Deff = 0.0;
            for (int k=1; k<=nspecies; ++k){
                if (dest(k) != 0){
                    Deff = Deff + MolarConcN_sub[dest(k)-1]/D_barN(k,column);
                }
            }
            Deff = 1.0/Deff;

            // assign rhoWchi
            for (int row=1; row<=nspecies; ++row){
                if ( row == column ){
                    rhoWchiN(row,column) = rhotot_sub*Deff*molmass_in[row-1]/molmtot_sub;
                } else if ( dest(row) == 0 ){
                    rhoWchiN(row,column) = 0.0;
                } else {
                    tmp = 0.0;
                    for (int k=1; k<=nspecies; ++k){
                        if ( dest(k)!=0 ){
                            tmp = tmp + chi_sub(dest(row),dest(k))*MolarConcN_sub[dest(k)-1]/D_barN(k,column);
                        }
                    }


                    rhoWchiN(row,column) = Deff*rhoN_sub[dest(row)-1]*(tmp - molmass_in[column-1]/molmtot_sub);
                }
            }
        } else {   // column of non-trace species
            // assign rhoWchi
            for (int row=1; row<=nspecies; ++row){
                if ( dest(row) == 0 ){
                    rhoWchiN(row,column) = 0.0;
                } else {
                    rhoWchiN(row,column) = rhoN_sub[dest(row)-1]*chi_sub(dest(row),dest(column));
                }
            }
        }
    }
}


/**
 * \brief Set 2-D square Matrix = 2-D squareMatrix
 *
 * \param[out] Matrix_out
 * \param[in] Matrix_in
 * \param[in] size
 *  Takes two Array2Ds of size `size` and sets Matrix_out = Matrix_in
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Matrix2DSetEqual (
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
 * \param[in] D_bar matrix of Maxwell-Stefan binary diffusion coefficient
 * \param[in] Xk mole fractions --- MUST NOT BE ZERO
 * \param[in] molmass_local
 * \param[out] chi multispecies diffusion matrix
 *
 * The number of terms in the sum is determined by the global variable `chi_iterations`.
 * 3-5 is a reasonable range for the value.
 *
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Dbar2chiIterative ( int nspecies_in,
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& D_bar,
               GpuArray<Real, MAX_SPECIES>& Xk,
               GpuArray<Real, MAX_SPECIES>& molmass_local,
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& chi){

    Array1D<Real, 1, MAX_SPECIES> Xkp, Ykp, Di, Mmat, Minv;
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> Pmat, Deltamat, Zmat, Jmat, PJ, matrix1, matrix2;

    // mole fractions correction
    Real sum = 0.;
    for (int ii=1; ii<=nspecies_in; ++ii){
       Xkp(ii) = Xk[ii-1];
       sum += Xkp(ii);
    }
//    if(amrex::abs(sum-1.) > 1.e-12){
//        amrex::Print() << " bad x in dbarier " << sum  << std::endl;
//    }

    // molecular weight of mixture - EGLIB
    Real MWmix = 0.0;
    for (int ii=1; ii<=nspecies_in; ++ii){
        MWmix = MWmix + Xkp(ii)*molmass_local[ii-1];

    }

    // mass fractions correction - EGLIB
    for (int ii=1; ii<=nspecies_in; ++ii){
        Ykp(ii) = molmass_local[ii-1] / MWmix * Xkp(ii);
    }

    // Find Di matrix
    Real term2;
    for (int i=1; i<=nspecies_in; ++i){
        term2 = 0.0;
        for (int j=1; j<=nspecies_in; ++j){
            if (j!=i){
                term2 = term2 + Xkp(j)/D_bar(i,j);
            }
        }
        Di(i) = (1.0 - Ykp(i))/term2;
    }

    // Compute Mmat and Minv
    for (int i=1; i<=nspecies_in; ++i){
        Mmat(i) = Xkp(i)/Di(i);
        Minv(i) = Di(i)/Xkp(i);
    }

    // Compute P matrix
    MatrixToZeros(nspecies_in, Pmat);
    for (int i=1; i<=nspecies_in; ++i){
        for (int j=1; j<=nspecies_in; ++j){
            Pmat(i,j) = - Ykp(j);
            if (i==j){
                Pmat(i,j) =  Pmat(i,j) + 1.0;
            }
        }
    }

    // Compute Deltamat
    MatrixToZeros(nspecies_in, Deltamat);
    Real term1;

    for (int i=1; i<=nspecies_in; ++i){
        for (int j=1; j<=nspecies_in; ++j){
            if (i==j){
                term1 = 0.0;
                for (int k=1; k<=nspecies_in; ++k){
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

    // Compute Zmat
    for (int i=1; i<=nspecies_in; ++i){
        Zmat(i,i) = Zmat(i,i) + Mmat(i);
    }

    // Compute Jmat
    for (int i=1; i<=nspecies_in; ++i){
        for (int j=1; j<=nspecies_in; ++j){
            Jmat(i,j) = Minv(i)*Zmat(i,j);
        }
    }

    // Compute PJ
    MatrixToZeros(nspecies_in, PJ);
    for (int i=1; i<=nspecies_in; ++i){
        for (int j=1; j<=nspecies_in; ++j){
            for (int k=1; k<=nspecies_in; ++k){
                PJ(i,j) = PJ(i,j) + Pmat(i,k)*Jmat(k,j);
            }
        }
    }

    // Compute P M^-1 Pt; store it in matrix2
    Real scr;
    for (int i=1; i<=nspecies_in; ++i){
        for (int j=1; j<=nspecies_in; ++j){
            scr = 0.0;
            for (int k=1; k<=nspecies_in; ++k){
                // notice the change in indices for Pmat to represent Pmat^t
                scr = scr + Pmat(i,k)*Minv(k)*Pmat(j,k);
            }
            matrix2(i,j) = scr;
            chi(i,j) = scr;
        }
    }


    for (int jj=1; jj<=chi_iterations; ++jj){
        for (int i=1; i<=nspecies_in; ++i){
            for (int j=1; j<=nspecies_in; ++j){
                    scr = 0.0;
                for (int k=1; k<=nspecies_in; ++k){
                    scr = scr + PJ(i,k)*chi(k,j);
                }
                matrix1(i,j) = scr + matrix2(i,j);
            }
        }
        Matrix2DSetEqual(chi, matrix1, MAX_SPECIES);
    }
}






/**
 * \brief Compute Chi
 *
 * \param[in] nspecies_in Number of species
 * \param[in] molmass_in
 * \param[in] rhoN
 * \param[in] rhotot
 * \param[in] MolarConcN
 * \param[in,out] chi multispcies diffusion matrix.
 * \param[in] D_barN matrix of Maxwell-Stefan binary diffusion coefficient
 *
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeChi (int nspecies_in,
               GpuArray<Real, MAX_SPECIES>& molmass_in,
               GpuArray<Real, MAX_SPECIES>& rhoN,
               Real rhotot,
               GpuArray<Real, MAX_SPECIES>& MolarConcN,
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& chi,
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& D_barN){

    //local variables
    Real eepsilon=1.e-16;

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
       printf("Compute Chi: use_lapack not supported\n");
       //call amrex_error('compute_chi: use_lapack not supported')
    }

    Dbar2chiIterative(nspecies_in, D_barN,MolarConcN,molmass_in,chi);

}



/**
 * \brief Compute rhoWchi
 *
 * \param[in] rhoN densities
 * \param[in] rhotot total density
 * \param[in] MolarConcN molar concentration
 * \param[in,out] rhoWchiN
 * \param[in] D_barN MS diff-coefs
 * \param[in] molmass_in
 *
 */


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeRhoWChiLocal (
               GpuArray<Real, MAX_SPECIES>& rhoN,
               Real rhotot,
               GpuArray<Real, MAX_SPECIES>& MolarConcN,
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& rhoWchiN,
               Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& D_barN,
               GpuArray<Real, MAX_SPECIES>& molmass_in){

    //local variables
    Array1D<Real, 1, MAX_SPECIES> W;
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> chi;
    int ntrace; /**< number of trace species with w_k < fractional_tolerance */
    int nspecies_sub; /**< dim of subsystem = nspecies - ntrace */
    Real w1,w2;

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

    if ((use_multiphase == 1) && (nspecies == 2)) {

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

       // compute the number of trace species
       // build the mapping for expanding/contracting arrays
        ntrace = 0;
        for (int row=1; row<=nspecies; ++row){
            W(row) = rhoN[row-1] / rhotot;
            if (W(row) < fraction_tolerance){
                ntrace = ntrace + 1;
                dest(row) = 0;
            } else {
                dest(row) = row - ntrace;
            }
        }

        if (ntrace == nspecies -1){
            // this is all trace species except for 1 (essentially pure solvent);
            MatrixToZeros(nspecies, rhoWchiN); //set rhoWchi to zero

        } else if (ntrace == 0){
            // there are no trace species, hence chi = chi_sub

            ComputeChi(nspecies, molmass_in, rhoN, rhotot, MolarConcN, chi, D_barN);

            for (int row=1; row<=nspecies; ++row){
                for (int column=1; column<=nspecies; ++column){
                    rhoWchiN(row, column) = rhoN[row-1] * chi(row,column);
                }
            }
        }
        else {
          // if there are trace species, we consider a subsystem
          // consisting of non-trace species
          nspecies_sub = nspecies - ntrace;

//          amrex::Print() << " entering subspecies " << ntrace << " " << nspecies_sub << std::endl;
//          for (int n=1; n<=nspecies; ++n){
//              amrex::Print() << "n , dest,molmass rho " << n << " " << dest(n) << " " << molmass_in[n-1] << " " << rhoN[n-1] << std::endl;
//              for (int m=1; m<=nspecies; ++m){
//                  amrex::Print() << "m d rhowchi " << m << " " << D_barN(n,m) << " " <<  rhoWchiN(n,m) << std::endl;
//              }
//          }
          ComputeChiSub(nspecies_sub,
                        dest,
                        molmass_in,
                        rhoN,
                        D_barN,
                        rhoWchiN);

        }
    }
}

/**
 * \brief Compute nonnegative rho face-centered average
 *
 * \param[in] Rho1 density from a neighboring cell
 * \param[in] Rho2 density from a neighboring cell
 * \param[in] dx
 * \param[out] RhoAv face-centered average
 *
 */


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeNonnegativeRhoAv (
                GpuArray<Real, MAX_SPECIES>& Rho1,
                GpuArray<Real, MAX_SPECIES>& Rho2,
                const GpuArray<Real, AMREX_SPACEDIM>& dx,
                GpuArray<Real, MAX_SPECIES>& molmass_in,
                GpuArray<Real, MAX_SPECIES>& RhoAv){

    Real dv, value1, value2, tmp1, tmp2;

#if (AMREX_SPACEDIM == 2)
    dv = 1.0;
    for (int dim=0; dim<2; ++dim){
        dv *= dx[dim];
    }
    dv *= cell_depth;
#elif (AMREX_SPACEDIM == 3)
    dv = 1.0;
    for (int dim=0; dim<3; ++dim){
        dv *= dx[dim];
    }
#endif

    if ((use_multiphase == 1) && (nspecies == 2)){
       value1 = Rho1[0];
       value2 = Rho2[0];
       if ((value1 <= 0.0) || (value2 <= 0.0)){
          RhoAv[0] = 0.0;
       } else {
          RhoAv[0] = std::min( 0.5*(value1 + value2), rho0);
       }
       RhoAv[1] = rho0 - RhoAv[0];

    } else {

        for (int comp=1; comp<=nspecies; ++comp){
            value1 = Rho1[comp-1]/molmass_in[comp-1]; //Convert to number density
            value2 = Rho2[comp-1]/molmass_in[comp-1];

            switch (avg_type) {

            case 1 : {    // Arithmetic with a C0-smoothed Heaviside

                if ( (value1 <= 0.0) || (value2 <= 0.0) ) {
                   RhoAv[comp-1] = 0.0;
                } else {
                   tmp1 = std::min(dv*value1,1.0);
                   tmp2 = std::min(dv*value2,1.0);
                   RhoAv[comp-1] = molmass_in[comp-1]*(value1+value2)/2.0*tmp1*tmp2;
                }

            } break;

            case 2 :   // Geometric
                RhoAv[comp-1] = molmass_in[comp-1]*std::sqrt(std::max(value1,0.0)*std::max(value2,0.0));
            break;

            case 3 : {   // Harmonic
                // What we want here is the harmonic mean of max(value1,0) and max(value2,0)
                // Where we define the result to be zero if either one is zero
                // But numerically we want to avoid here division by zero
                Real tiny = std::numeric_limits<Real>::min(); //replacement for fortran tiny(1.d0) function call
                if ( (value1 <= 10.0*tiny) || (value2 <= 10.0*tiny) ) {
                   RhoAv[comp-1] = 0.0;
                } else {
                   RhoAv[comp-1] = molmass_in[comp-1]*2.0/(1.0/value1+1.0/value2);
                }
            } break;

            case 10 : {  // Arithmetic with (discontinuous) Heaviside
                if ( (value1 <= 0.0) || (value2 <= 0.0) ) {
                   RhoAv[comp-1] = 0.0;
                } else {
                   RhoAv[comp-1] = molmass_in[comp-1]*(value1+value2)/2.0;
                }
            } break;

            case 11 : {  // Arithmetic with C1-smoothed Heaviside
                if ( (value1 <= 0.0) || (value2 <= 0.0) ) {
                   RhoAv[comp-1] = 0.0;
                } else {
                    tmp1 = dv*value1;
                    if (tmp1 < 1.0) {
                        tmp1 = (3.0-2.0*tmp1)*tmp1*tmp1;
                    } else {
                        tmp1 = 1.0;
                    }
                    tmp2 = dv*value2;
                    if (tmp2 < 1.0) {
                        tmp2 = (3.0 - 2.0*tmp2)*tmp2*tmp2;
                    } else {
                        tmp2 = 1.0;
                    }
                    RhoAv[comp-1] = molmass_in[comp-1]*(value1 + value2)/2.0*tmp1*tmp2;
                }
            } break;

            case 12 : {  // Arithmetic with C2-smoothed Heaviside
                if ( (value1 <= 0.0) || (value2 <= 0.0) ) {
                    RhoAv[comp-1] = 0.0;
                } else {
                    tmp1 = dv*value1;
                    if (tmp1 < 1.0) {
                        tmp1=(10.0 - 15.0*tmp1 + 6.0*tmp1*tmp1)*std::pow(tmp1,3);
                    } else {
                        tmp1=1.0;
                    }
                    tmp2 = dv*value2;
                    if (tmp2 < 1.0) {
                        tmp2 = (10.0-15.0*tmp2+6.0*tmp2*tmp2)*std::pow(tmp2,3);
                    } else {
                        tmp2 = 1.0;
                    }
                    RhoAv[comp-1] = molmass_in[comp-1]*(value1+value2)/2.0*tmp1*tmp2;
                }
            } break;

            default :
                printf("compute_nonnegative_rho_av: invalid avg_type\n");
                Abort();
            }
        }
    }
}



/**
 * \brief Water-glycerol mixtures near room temperature
 *
 * \param[in,out] chi
 * \param[in] RhoN
 * \param[in] rhotot
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ChiWaterGlycerol (
                Real& chi,
                GpuArray<Real, MAX_SPECIES>& RhoN,
                Real rhotot_in){

    // local
    Real c_loc;

    // mass fraction of glycerol
    c_loc = RhoN[0]/rhotot_in;

    // chi = chi0 * rational function
    chi = Dbar[0]*(1.024-1.001692692*c_loc)/(1.0+0.6632641981*c_loc);

}


/**
 * \brief compute D_bar
 *
 * \param[in] RhoN
 * \param[in] rhotot
 * \param[in,out] D_bar
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeDBarLocal (
          GpuArray<Real, MAX_SPECIES>& RhoN,
          Real rhotot_in,
          Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& D_bar){

    // off-diagonal components of symmetric matrices
    Array1D<Real, 1, MAX_SPECIES*(MAX_SPECIES-1)/2> DBar_loc;

    switch (std::abs(mixture_type)) {

    case 1 : { // water-glycerol

        //Note: This case could not be verified.
        printf("mixture_type=1 is not currently supported.\n");
        Abort();

       // we require nspecies=2
       // Dbar(1) = chi0 in the binary notation
        if (nspecies != 2){
            printf("mixture_properties_mass_local assumes nspecies=2 if mixture_type=3 (water-glycerol)\n");
            Abort();
        }

        ChiWaterGlycerol(DBar_loc(1), RhoN, rhotot_in);

    } break;

    case 2 : { // Electrolyte mixture

        //Note: The case could not be verified.
        printf("mixture_type=2 is not currently supported.\n");
        Abort();

        if (nspecies != 3){
            printf("mixture_properties_mass_local assumes nspecies=3 if mixture_type=2 (water-glycerol)\n");
            Abort();
        }

        // This is the leading-order correction for dilute solutions
        // In particular, the cross-diffusion coefficient ~sqrt(concentration) as per renormalization theory
        // The ordering of the values is D_12; D_13, D_23

        DBar_loc(1) = Dbar[0]*std::sqrt(RhoN[0]/rhotot_in); // counter-ion cross coefficient D_12~sqrt(w)
        DBar_loc(2) = Dbar[1]; // D_13 = self diffusion of first ion
        DBar_loc(3) = Dbar[2]; // D_23 = self diffusion of second ion

    } break;

    default : //Keep it constant

        for ( int n=1; n <= nspecies*(nspecies-1)/2; ++n){
            DBar_loc(n) = Dbar[n-1];
        }

    }

    // Complete the process by filling the matrices using generic formulae -- this part should not change
    // populate D_bar and Hessian matrix
    int n=0;
    for (int row=1; row<=nspecies; ++row){
        for (int column=1; column<=row-1; ++column){
            n = n+1;
            D_bar(row, column) = DBar_loc(n);         // SM-diffcoeff's read from input
            D_bar(column, row) = D_bar(row, column);     // symmetric
        }
        // populate diagonals
        D_bar(row, row) = 0.0;           // as self-diffusion is zero
    }
}

/**
 * \brief lower triangle and diagonal are overwritten by cholesky factor
 *
 * \param[in] A input matrix
 * \param[in] np
 *
 * A is ths input matrix. Upon return the lower triangle and diagonal
 * are overwritten by the cholesky factor.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void CholeskyDC (Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& A,
                int np){

    int ising;
    Real sum1;

    Array1D<Real, 1, MAX_SPECIES> P;

    Real small_number = 0.0; // Some tolerance

    for (int i=1; i<=np; ++i){
        ising = 0;

        for (int j=i; j<=np; ++j){
            sum1 = A(i,j);

            for (int k=i-1; k>=1; --k){
                sum1 = sum1 - A(i,k)*A(j,k);
            }
            if(i==j){

                if(sum1<=small_number){
                    P(i) = 0.0;
                    ising = 1;

                } else {
                    P(i) = std::sqrt(sum1);
                }
            } else {

                if( ising==0 ){
                    A(j,i) = sum1/P(i);
                } else {
                    A(j,i) = 0.0;
                }
            }
        }
    }

    for (int i = 1; i<=np; ++i){
        for (int j = i+1; j<=np; ++j){
            A(i,j) = 0.0; // Zero upper triangle
        }
    A(i,i) = P(i);
    }

}

/**
 * \brief Compute Sqrt L-Onsager on subset
 *
 * \param[in] dest map from full system to subset
 * \param[in] Molmass_sub
 * \param[in] Rho_sub
 * \param[in] rhotot_sub
 * \param[in] nspecies_sub
 * \param[in,out] SqrtLOnsager_sub
 *
 * dest is a mapping used to eliminate elements in D_bar we don't need (for D_bar_sub)
 * and for expanding sqrtLonsager_sub into sqrtLonsager
 * it will contain the numbers 1, 2, ..., (nspecies-ntrace)
 * with zeros in elements corresponding to trace elements
 * (example) for a 5-species system having trace species 2 and 5:
 * species       1 2 3 4 5
 * dest(species) 1 0 2 3 0
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeSqrtLOnsagerSub (
                Array1D<int, 1, MAX_SPECIES>& dest,
                GpuArray<Real, MAX_SPECIES>& Molmass_in,
                GpuArray<Real, MAX_SPECIES>& Rho_in,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& DBar_in,
                Real rhotot_in,
                int nspecies_sub,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& SqrtLOnsager){

    // local variables
    GpuArray<Real, MAX_SPECIES> Molmass_sub;
    GpuArray<Real, MAX_SPECIES> Rho_sub;
    GpuArray<Real, MAX_SPECIES> Molarconc_sub;

    Real molmtot_sub;
    Array1D<Real, 1, MAX_SPECIES> W_sub;
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> Chi_sub;
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> DBar_sub;
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> SqrtLOnsager_sub;

    // create a vector of non-trace densities and molmass for the subsystem
    for (int row=1; row<=nspecies; ++row){
        if (dest(row) != 0){
            Molmass_sub[dest(row)-1] = Molmass_in[row-1];
            Rho_sub[dest(row)-1] = Rho_in[row-1];
        }
    }

    // renormalize total density and mass fractions
    Real rhotot_sub = 0.0;
    for (int n=0; n<nspecies_sub; ++n){
        rhotot_sub += Rho_sub[n];
    }

    for (int row=1; row<= nspecies_sub; ++row){
        W_sub(row) = Rho_sub[row-1]/rhotot_sub;
    }

    // first, compute the full D_bar
    // then, construct D_bar_sub by mapping the full D_bar into D_bar_sub
    // you could read in only the lower diagonals,
    // reflect, and set the diagnals to zero if you want

    ComputeDBarLocal(Rho_in, rhotot_in, DBar_in);

    for (int row=1; row<= nspecies; ++row){
        if (dest(row) == 0){
            continue;
        }
        for (int column=1; column<=nspecies; ++column){
            if (dest(column) != 0){
                DBar_sub(dest(row),dest(column)) = DBar_in(row,column);
            }
        }
    }

    // compute molarconc_sub and molmtot_sub
    ComputeMolconcMolmtotLocal(nspecies_sub, Molmass_sub, Rho_sub, rhotot_sub, Molarconc_sub, molmtot_sub);



    // compute chi_sub
    ComputeChi(nspecies_sub, Molmass_sub, Rho_sub, rhotot_sub, Molarconc_sub, Chi_sub, DBar_sub);

    // compute Onsager matrix L_sub (store in sqrtLonsager_sub)
    for (int column=1; column<= nspecies_sub; ++column){
        for (int row=1; row<= nspecies_sub; ++row){
            SqrtLOnsager_sub(row, column) = molmtot_sub*rhotot_sub*W_sub(row)*Chi_sub(row,column)*W_sub(column)/k_B;
        }
    }

    // compute cell-centered Cholesky factor, sqrtLonsager_sub
    CholeskyDC(SqrtLOnsager_sub, nspecies_sub);

    // expand sqrtLonsager_sub into sqrtLonsager
    MatrixToZeros(nspecies, SqrtLOnsager);
    for (int row=1; row<= nspecies; ++row){
        if (dest(row) == 0){
           continue;
        }
        for (int column=1; column<= nspecies; ++column){
            if (dest(column) != 0){
                SqrtLOnsager(row,column) = SqrtLOnsager_sub(dest(row),dest(column));
            }
        }
    }

}

/**
 * \brief Compute Sqrt L-Onsager
 *
 * \param[in] molmass_in
 * \param[in] RhoN
 * \param[in] rhotot_in
 * \param[in,out] sqrtLonsager
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeSqrtLOnsagerLocal (
                GpuArray<Real, MAX_SPECIES>& molmass_in,
                GpuArray<Real, MAX_SPECIES>& RhoN,
                Real rhotot_in,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& SqrtLOnsager){

    // local variables
    Array1D<Real, 1, MAX_SPECIES> W;

    GpuArray<Real, MAX_SPECIES> Molarconc_loc;
    Real molmtot_loc;
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> Chi_loc;
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> DBar_loc;

    int ntrace;                /**<   number of trace species with w_k < fractional_tolerance */
    int nspecies_sub;          /**<   dim of subsystem = nspecies-ntrace */

   /*
    * this is a mapping used to eliminate elements in D_bar we don't need (for D_bar_sub)
    * and for expanding sqrtLonsager_sub into sqrtLonsager
    * it will contain the numbers 1, 2, ..., (nspecies-ntrace)
    * with zeros in elements corresponding to trace elements
    * (example) for a 5-species system having trace species 2 and 5:
    *  species       1 2 3 4 5
    *  dest(species) 1 0 2 3 0
    */
    Array1D<int, 1, MAX_SPECIES> dest;

    // compute the number of trace species
    // build the mapping for expanding/contracting arrays
    ntrace = 0;
    for (int row=1; row <= nspecies; ++row){
        W(row) = RhoN[row-1]/rhotot_in;
        if (W(row) < fraction_tolerance){
           ntrace = ntrace + 1;
           dest(row) = 0;
        } else {
           dest(row) = row - ntrace;
        }
    }

    if (ntrace == nspecies-1){

        // this is all trace species except for 1 (essentially pure solvent);
        // set sqrtLonsager to zero
        MatrixToZeros(MAX_SPECIES, SqrtLOnsager);

    } else if (ntrace == 0) {    // there are no trace species

        // compute molarconc and molmtot
        ComputeMolconcMolmtotLocal(nspecies, molmass_in, RhoN, rhotot_in, Molarconc_loc, molmtot_loc);

        // compute D_bar
        ComputeDBarLocal(RhoN, rhotot_in, DBar_loc);

        // compute chi
        ComputeChi(nspecies, molmass_in, RhoN, rhotot_in, Molarconc_loc, Chi_loc, DBar_loc);

        // compute Onsager matrix L (store in sqrtLonsager)
        for (int column=1; column<= nspecies; ++column){
            for (int row=1; row<=nspecies; ++row){
                SqrtLOnsager(row, column) = molmtot_loc*rhotot_in*W(row)*Chi_loc(row,column)*W(column)/k_B;
            }
        }

        // compute cell-centered Cholesky factor, sqrtLonsager
        CholeskyDC(SqrtLOnsager,nspecies);

    } else {

        // if there are trace species, we consider a subsystem
        // consisting of non-trace species

        nspecies_sub = nspecies - ntrace;
        ComputeSqrtLOnsagerSub( dest,
                        molmass_in,
                        RhoN,
                        DBar_loc,
                        rhotot_in,
                        nspecies_sub,
                        SqrtLOnsager);

    }

}

/**
 * \brief Compute Zeta by temperature
 *
 * \param[in] MolarConc molar concentration
 * \param[in] DBar MS diff-coeffs
 * \param[in] Temp Temperature
 * \param[in,out] ZetaByTemp zeta/T
 * \param[in] DTerm Thermo diff-coeffs
 *
 *
 */

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ComputeZetaByTempLocal (
                GpuArray<Real, MAX_SPECIES>& MolarConc,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& DBar,
                Real Temp,
                GpuArray<Real, MAX_SPECIES>& ZetaByTemp,
                GpuArray<Real, MAX_SPECIES>& DTherm){

    // local variables
    Real Sum_knoti;
    Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES> Lambda;

    // compute zeta_by_Temp for thermodiffusion
    if (is_nonisothermal == 1){

       // compute Lambda_ij matrix; molarconc is
       // expressed in terms of molmtot,mi,rhotot etc.
        for (int row=1; row<=nspecies; ++row){
            for (int column=1; column<=row-1; ++column){
                Lambda(row, column) = -MolarConc[row-1]*MolarConc[column-1]/DBar(row,column);
                Lambda(column, row) = Lambda(row, column);
            }
        }
        for (int row=1; row<=nspecies; ++row){
            Sum_knoti = 0.0;
            for (int column=1; column<=nspecies; ++column){
                if (column != row){
                    Sum_knoti = Sum_knoti + Lambda(row,column)*(DTherm[row-1] - DTherm[column-1]);
                }
                ZetaByTemp[row-1] = Sum_knoti/Temp;
            }
        }
    }

}

/**
 * \brief Mixture Mass Properties
 *
 * \params[in] Rho
 * \params[in] rhotot
 * \params[out] DBar
 * \params[out] DTherm
 * \params[out] Hessian
 *
 * The default case should be to simply set Dbar, Dtherm and Hessian to constants, and read from the input file.
 *
 *
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void MixturePropsMassLocal (GpuArray<Real, MAX_SPECIES>& Rho,
                Real rhotot,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& DBar,
                GpuArray<Real, MAX_SPECIES>& DThermLoc,
                Array2D<Real, 1, MAX_SPECIES, 1, MAX_SPECIES>& Hessian){

    // Local values of transport and thermodynamic coefficients (functions of composition!):
    GpuArray<Real, MAX_SPECIES*(MAX_SPECIES-1)/2> HOffdiagLoc; //off-diagonal components of symmetric matrices
    GpuArray<Real, MAX_SPECIES> HDiagLoc; // Diagonal component

    // populate D_bar
    ComputeDBarLocal(Rho, rhotot, DBar);

    // For now we only encode constant Hessian matrices since we do not have any thermodynamic
    // models coded up (Wilson, NTLR, UNIQUAC, etc.)

    for (int n=0; n<nspecies; ++n){
        HDiagLoc[n] = H_diag[n];
        DThermLoc[n]  = Dtherm[n];    // thermal diffcoeff's read from input
    }

    for (int n=0; n < nspecies*(nspecies-1)/2; ++n){
        HOffdiagLoc[n] = H_offdiag[n];
    }

    // Complete the process by filling the matrices using generic formulae -- this part should not change
    // populate Hessian matrix
    int n=0;

    if (is_ideal_mixture == 0){
        for (int row=1; row<=nspecies; ++row){
            for (int column=1; column<=row-1; ++column){
                ++n;
                Hessian(row, column) = HOffdiagLoc[n-1];  // positive semidefinite matrix read from input
                Hessian(column, row) = Hessian(row,column);  // Hessian is symmetric
            }
        Hessian(row, row) = HDiagLoc[row-1];
        }
    } else {
        MatrixToZeros(nspecies,Hessian); // set matrix to zero for ideal-mixture
    }

}


#endif
